#!/usr/bin/env python

from hpc_utils.hpc import get_ref_file
# from pybedtools import BedTool
from ngs_utils import gtf

''' Generates coding_regions BED file for SAGE
    https://github.com/hartwigmedical/hmftools/tree/master/sage
'''

GENOME = 'GRCh37'
OUTPUT= 'coding_regions.bed'

gtf_fpath = get_ref_file(GENOME, key='gtf')
db = gtf.get_gtf_db(gtf_fpath)

features = []
for rec in db.all_features(order_by=('seqid', 'start', 'end')):
    if rec.end - rec.start < 0: continue
    if rec.featuretype == 'CDS':
        features.append([rec.chrom,
                         str(rec.start - 1),
                         str(rec.end)])

with open(OUTPUT, 'w') as out:
    for f in features:
        out.write('\t'.join(f) + '\n')
#BedTool(features).saveas(OUTPUT)
