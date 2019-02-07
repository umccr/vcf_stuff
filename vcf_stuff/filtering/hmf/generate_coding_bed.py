#!/usr/bin/env python
from bed_annotation import get_canonical_transcripts_ids
from hpc_utils.hpc import get_ref_file
# from pybedtools import BedTool
from ngs_utils import gtf
import sys

''' Generates coding_regions BED file for SAGE
    https://github.com/hartwigmedical/hmftools/tree/master/sage

    Usage: 
    python {__file__} | sort -k1,1V -k2,2n | grep -v ^MT | grep -v ^GL | bedtools merge -i - > coding_regions.canonical.sort.merged.bed
'''

GENOME = 'GRCh37'
if len(sys.argv) > 1:
    out = open(sys.argv[1], 'w')
else:
    out = sys.stdout

gtf_fpath = get_ref_file(GENOME, key='gtf')
db = gtf.get_gtf_db(gtf_fpath)

canon_by_gname = get_canonical_transcripts_ids(GENOME)
# def get_only_canonical_filter():
#     return lambda x: x[BedCols.ENSEMBL_ID] == canon_tx_by_gname.get(x[BedCols.HUGO]) or \
#                      x[BedCols.ENSEMBL_ID] == canon_tx_by_gname.get(x[BedCols.GENE])

def _get_attr(_rec, _key):
    val = _rec.attributes.get(_key)
    if val is None:
        return None
    assert len(val) == 1, (_key, str(val))
    return val[0]

features = []
for rec in db.all_features(order_by=('seqid', 'start', 'end')):
    if rec.end - rec.start < 0: continue
    if rec.featuretype == 'CDS' and _get_attr(rec, 'transcript_id') in canon_by_gname.values():
        features.append([rec.chrom,
                         str(rec.start - 1),
                         str(rec.end)])

for f in features:
    out.write('\t'.join(f) + '\n')

try:
    out.close()
except:
    pass

