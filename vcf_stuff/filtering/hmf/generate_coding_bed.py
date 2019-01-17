#!/usr/bin/env python
from bed_annotation import get_canonical_transcripts_ids
from hpc_utils.hpc import get_ref_file
# from pybedtools import BedTool
from ngs_utils import gtf

''' Generates coding_regions BED file for SAGE
    https://github.com/hartwigmedical/hmftools/tree/master/sage
'''

GENOME = 'GRCh37'
OUTPUT= 'coding_regions.canonical.bed'

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

with open(OUTPUT, 'w') as out:
    for f in features:
        out.write('\t'.join(f) + '\n')
#BedTool(features).saveas(OUTPUT)


# Then sort and clean:
# sort -k1,1V -k2,2n coding_regions.canonical.bed | grep -v ^MT | grep -v ^GL > coding_regions.sort.bed