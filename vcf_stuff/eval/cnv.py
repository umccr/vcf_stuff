import glob
from os.path import join, basename, splitext, dirname

from ngs_utils.file_utils import safe_mkdir
from ngs_utils.logger import critical, info
import numpy as np
from pybedtools import BedTool
import csv
from cyvcf2 import VCF, Writer

"""
==> HCC2218_cnvkit-call.cns <==
chromosome  start    end       gene            log2       baf   cn  cn1  cn2  depth    probes  weight
1           30334    622042    MIR1302-10      0.397727         3             186.342  15      3.32545
1           721404   1277866   RP11-206L10.9   -1.28897   0     0   0    0    27.6628  264     66.6316

==> HCC2218_facets_cncf.tsv <==
chrom  seg  num.mark  nhet  cnlr.median          mafR                  segclust  cnlr.median.clust   mafR.clust            start      end        cf.em              tcn.em  lcn.em
1      1    3611      318   -1.07727964252529    7.12566902807843      1         -1.06980799758558   9.84653017083575      13483      31905860   0.958261957063053  1       0
1      2    6288      539   -0.0951260579848125  -0.00186540158968537  13        -0.116727310455426  0.0012470526861517    31960791   142826805  1                  2       1

==> HCC2218_manta.vcf <==
1    44195682    MantaDEL:369:0:1:0:0:0  T   <DEL>  .   PASS   BPI_AF=0.787,0.402;BPI_END=44201846;BPI_START=44195682;END=44201846;SOMATIC;SOMATICSCORE=205;SVLEN=-6164;SVTYPE=DEL     PR:SR   41,0:45,0       14,25:31,50
1    211021468   MantaDEL:1421:0:1:0:1:0 A   <DEL>  .   PASS   BPI_AF=0.694,0.073;BPI_END=211033724;BPI_START=211021468;END=211033724;SOMATIC;SOMATICSCORE=77;SVLEN=-12256;SVTYPE=DEL  PR:SR   46,0:155,0      35,6:193,14

==> HCC2218_purple.cnv.tsv <==
#chromosome  start      end        copyNumber          bafCount  observedBAF         actualBAF           segmentStartSupport  segmentEndSupport  method
1            1          36946000   0.8035740452186546  2         1.0                 1.0                 TELOMERE             NONE               BAF_WEIGHTED
1            36946001   121485000  1.5790169587407195  3         0.5291133717401169  0.5291133717401169  NONE                 CENTROMERE         BAF_WEIGHTED

==> HCC2218_truthset_cnv_bcbio.tsv <==
chrom  start      end        tot_cn
1      10001      31920394   1
1      144810725  145092947  5
"""


# Check with peter for what he looks in CNV
class CnvCall:
    def __init__(self, chrom=None, start=None, end=None, gene=None, cn=None, event=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.gene = gene
        self.cn = cn
        self.event = event

    def get_bed_raw(self):
        return [self.chrom, str(self.start), str(self.end), self.gene, self.cn, self.event]


header_by_caller = {
    'cnvkit': 'chromosome	start	end	gene	log2	baf	cn	cn1	cn2	depth	probes	weight'.split(),
    'facets': 'chrom	seg	num.mark	nhet	cnlr.median	mafR	segclust	cnlr.median.clust	mafR.clust	start	end	cf.em	tcn.em	lcn.em'.split(),
    'purple': '#chromosome	start	end	copyNumber	bafCount	observedBAF	actualBAF	segmentStartSupport	segmentEndSupport	method'.split(),
    'hcc2218truthset': 'chrom	start	end	tot_cn'.split(),
}

parse_row_by_caller = {
    'cnvkit': lambda r: CnvCall(
        chrom=r.get('chromosome'),
        start=r.get('start'),
        end  =r.get('end'),
        gene =r.get('gene'),
        cn   =r.get('cn'),
    ),
    'facets': lambda r: CnvCall(
        chrom=r.get('chrom').replace('23', 'X').replace('24', 'Y'),
        start=r.get('start'),
        end  =r.get('end'),
        cn   =r.get('tcn.em'),
    ),
    'purple': lambda r: CnvCall(
        chrom=r.get('#chromosome'),
        start=r.get('start'),
        end  =r.get('end'),
        cn   =round(float(r.get('copyNumber'))),
    ),
    'hcc2218truthset': lambda r: CnvCall(
        chrom=r.get('chrom'),
        start=r.get('start'),
        end  =r.get('end'),
        cn   =r.get('tot_cn'),
    ),
}



def iter_manta(fh):
    vcf = VCF(fh)
    for r in vcf:
        if r.INFO['SVTYPE'] == 'DEL' and r.INFO.get('BPI_AF') is not None:
            cn = round(float(np.mean(r.INFO['BPI_AF'])) * 2)
            yield CnvCall(
                chrom=r.CHROM,
                start=r.INFO['BPI_START'],
                end  =r.INFO['BPI_END'],
                cn   =cn
            )


def get_iter_cnv(header, parse_row_fn):
    def parse_row(cnv_path):
        with open(cnv_path) as fh:
            next(fh)
            reader = csv.DictReader(fh, fieldnames=header, delimiter='\t')
            for row in reader:
                yield parse_row_fn(row)
    return parse_row


def cnv_to_bed(cnv_path, out_bed_path):
    with open(cnv_path) as fh:
        parse_fn = None
        header = next(fh).strip().split('\t')

        if header[0].startswith('##fileformat=VCF'):
            # Manta
            info(f'Detected {cnv_path} as caller "manta"')
            parse_fn = iter_manta

        else:
            for caller, hdr in header_by_caller.items():
                if header == hdr:
                    print(f'Parsing {cnv_path} as caller "{caller}" with header {hdr}')
                    parse_fn = get_iter_cnv(header, parse_row_by_caller[caller])

        if not parse_fn:
            critical(f'Cannot detect CNV file format in {cnv_path}')

    with open(out_bed_path, 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        for i, call in enumerate(parse_fn(cnv_path)):
            bed_row = call.get_bed_raw()
            writer.writerow(bed_row)
            if i == 0:
                print(bed_row)
                print()


if __name__ == '__main__':
    for fp in glob.glob('/Users/vsaveliev/git/umccr/vcf_stuff/tests/data/cnv/*'):
        if 'titan' not in fp:
            res = join('/Users/vsaveliev/git/umccr/vcf_stuff/tests/results/eval_cnv/',
                       basename(splitext(fp)[0] + '.bed'))
            safe_mkdir(dirname(res))
            cnv_to_bed(fp, res)




