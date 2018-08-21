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

==> hartwig_truthset_clean.tsv <==
Truth_type  Truth_orientation  Truth_sv_len  Truth_chr1  Truth_bp1  Truth_chr2  Truth_bp2  Manta_filter     Manta_chr1  Manta_bp1  Manta_chr2  Manta_bp2  BPI_adj_Manta_chr1  BPI_adj_Manta_bp1  BPI_adj_Manta_chr2  BPI_adj_Manta_bp2  Manual_check_1  Manual_check_2  Gridss_filter1   Gridss_diff_BPI1  Gridss_filter2   Gridss_diff_BPI2  Gridss_chr1  Gridss_bp1  Gridss_chr2  Gridss_bp2  Purple_type1  Purple_type2  Purple_chr1  Purple_bp1  Purple_chr2  Purple_bp2  Conserting_chr1  Conserting_bp1  Conserting_chr2  Conserting_bp2
DEL         INNIE              33588         1           207981233  1           208014821  NA               1           207981229  1           208014817  1                   207981233          1                   208014821          NA              NA              NA               0                 NA               1                 1            207981233   1            208014822   DEL           DEL           1            207981001   1            208015001   1                207981233       1                208014818
DUP         OUTIE              153518        1           224646602  1           224800120  NA               1           224646602  1           224800120  1                   224646602          1                   224800120          NA              NA              NA               1                 NA               0                 1            224646603   1            224800120   DUP           DUP           1            224647001   1            224800001   1                224646603       1                224800120
DEL         INNIE              3238          1           224782795  1           224786033  NA               1           224782795  1           224786033  1                   224782795          1                   224786033          NA              NA              NA               0                 NA               1                 1            224782795   1            224786034   DEL           DEL           1            224782795   1            224786033   1                224782795       1                224786034
BND         TANDEM_RIGHT       NA            1           87337011   10          36119127   MinSomaticScore  1           87336976   10          36119258   1                   87337011           10                  36119127           NA              NA              SINGLE_ASSEMBLY  0                 SINGLE_ASSEMBLY  0                 1            87337011    10           36119127    BND           NA            1            87337001    NA           NA          NA               NA              NA               NA
DEL         INNIE              2732596       10          33386464   10          36119060   MinSomaticScore  10          33386392   10          36118937   10                  33386464           10                  36119060           TRUE            TRUE            NA               NA                NA               NA                NA           NA          NA           NA          BND           NA            10           33386001    NA           NA          NA               NA              NA               NA
BND         TANDEM_RIGHT       NA            10          60477422   12          72667074   NA               10          60477422   12          72667074   10                  60477422           12                  72667074           NA              NA              SINGLE_ASSEMBLY  1                 SINGLE_ASSEMBLY  0                 10           60477423    12           72667074    NA            NA            NA           NA          NA           NA          NA               NA              NA               NA
BND         OUTIE              NA            10          7059510    19          17396810   NA               10          7059511    19          17396810   10                  7059510            19                  17396810           NA              NA              NA               1                 NA               0                 10           7059511     19           17396810    NA            NA            NA           NA          NA           NA          NA               NA              19               17396811
BND         INNIE              NA            10          7132876    19          17397643   NA               10          7132872    19          17397640   10                  7132876            19                  17397643           NA              NA              NA               0                 NA               1                 10           7132876     19           17397644    BND           NA            10           7133001     NA           NA          NA               NA              19               17397637
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
    'hcc2218_truth': 'chrom	start	end	tot_cn'.split(),
    'colo829_hartwig_truth': 'Truth_type	Truth_orientation	Truth_sv_len	Truth_chr1	Truth_bp1	Truth_chr2	Truth_bp2	Manta_filter	Manta_chr1	Manta_bp1	Manta_chr2	Manta_bp2	BPI_adj_Manta_chr1	BPI_adj_Manta_bp1	BPI_adj_Manta_chr2	BPI_adj_Manta_bp2	Manual_check_1	Manual_check_2	Gridss_filter1	Gridss_diff_BPI1	Gridss_filter2	Gridss_diff_BPI2	Gridss_chr1	Gridss_bp1	Gridss_chr2	Gridss_bp2	Purple_type1	Purple_type2	Purple_chr1	Purple_bp1	Purple_chr2	Purple_bp2	Conserting_chr1	Conserting_bp1	Conserting_chr2	Conserting_bp2'.split(),
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
    'hcc2218_truth': lambda r: CnvCall(
        chrom=r.get('chrom'),
        start=r.get('start'),
        end  =r.get('end'),
        cn   =r.get('tot_cn'),
    ),
    'colo829_hartwig_truth': lambda r: CnvCall(
        chrom=r.get('Truth_chr1'),
        start=r.get('Truth_bp1'),
        end  =r.get('Truth_bp2'),
        event={'DEL': 'Del', 'DUP': 'Amp'}[r.get('Truth_type')],
    ) if r['Truth_chr1'] == r['Truth_chr2'] and r['Truth_type'] in ['DEL', 'DUP'] else None,
}


def iter_manta(fh):
    vcf = VCF(fh)
    for r in vcf:
        if r.INFO['SVTYPE'] in ['DEL', 'DUP']:
            yield CnvCall(
                chrom=r.CHROM,
                start=r.INFO['BPI_START'],
                end  =r.INFO['BPI_END'],
                event={'DEL': 'Del', 'DUP': 'Amp'}[r.INFO['SVTYPE']],
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
            if call:
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




