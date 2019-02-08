#!/usr/bin/env python
import os
import sys
from os.path import dirname
from bed_annotation import get_canonical_transcripts_ids
from ngs_utils.file_utils import open_gzipsafe
from ngs_utils.logger import warn
from hpc_utils.hpc import get_ref_file

''' Generates coding_regions BED file for SAGE
    https://github.com/hartwigmedical/hmftools/tree/master/sage

    Usage: 
    python {__file__} -g GRCh37 | sort -k1,1V -k2,2n | grep -v ^MT | grep -v ^GL | bedtools merge -i - > coding_regions.canonical.sort.merged.bed
'''

GENOME = 'GRCh37'
out = sys.stdout


canon_by_gname = get_canonical_transcripts_ids(GENOME)


genes_set = set()
genes_without_canon = set()
gtf_path = os.path.join(get_ref_file(key='pyensembl_data'), 'GRCh37/ensembl75/Homo_sapiens.GRCh37.75.gtf.gz')
warn(f'Parsing {gtf_path}')
with open_gzipsafe(gtf_path) as f:
    lines_cnt = 0
    exons_cnt = 0
    for l in f:
        if not l.startswith('#') and l.strip():
            lines_cnt += 1
            fields = l.strip().split('\t')
            try:
                chrom, biotype, feature, start, end, _, _, _, annotations = fields
            except:
                warn(f'Cannot read fields{str(fields)}')
                raise
            if (biotype == 'protein_coding' or 'decay' in biotype) and feature in ['CDS']:
                annotations = {kv.split()[0].strip().strip('"'): kv.split()[1].strip().strip('"') for kv in annotations.split('; ')}
                gene_name = annotations['gene_name']
                # import pdb; pdb.set_trace()
                # TODO: check canonical transcripts against cancer genes
                transcript_id = annotations['transcript_id']
                canon_transcript_id = canon_by_gname.get(gene_name)
                if not canon_transcript_id:
                    genes_without_canon.add(gene_name)
                elif transcript_id == canon_transcript_id:
                    start = int(start) - 1
                    end = int(end)
                    if end - start >= 3:
                        print('\t'.join([chrom, str(start), str(end), gene_name]))
                        genes_set.add(gene_name)
                        exons_cnt += 1
    if exons_cnt % 10000 == 0:
        warn(f'Processed {len(genes_set)} genes, written {exons_cnt} CDS and stop_codon regions...')

warn(f'Done. Processed {len(genes_set)} genes, written {exons_cnt} exons')
if genes_without_canon:
    warn(f'No canonical transcript for {len(genes_without_canon)} genes: {genes_without_canon}')

sys.exit(1)



# os.environ['PYENSEMBL_CACHE_DIR'] = dirname(get_ref_file(key='pyensembl_data'))
# from pyensembl import EnsemblRelease
# data = EnsemblRelease(95 if '38' in GENOME else 75)
#
# canon_by_gname = get_canonical_transcripts_ids(GENOME)
# canon_ids = canon_by_gname.values()
#
#
# # gene_cnt = 0
# # exon_cnt = 0
# # for transcript in data.transcripts():
# #     gene_cnt += 1
# #     if transcript.id in canon_ids:
# #         for exon in transcript.exons:
# #             exon_cnt += 1
# #             print('\t'.join([exon.contig, str(exon.start - 1), str(exon.end), exon.gene_name]))
# #     if gene_cnt % 1000 == 0:
# #         warn(f'Processed {gene_cnt} genes, {exon_cnt} exons...')
# # warn(f'Done. Processed {gene_cnt} genes, {exon_cnt} exons')
#
#
# gene_cnt = 0
# exon_cnt = 0
# trx_cnt = 0
# for gene in data.genes():
#     gene_cnt += 1
#     import pdb; pdb.set_trace()
#     # if gene.transcript_biotype == 'protein_coding'
#     gene_canon_ids   = [t for t in data.transcript_ids_of_gene_id(gene.id) if t in canon_ids]
#     if len(gene_canon_ids) >= 2:
#         warn(f'Warning: multiple canonical transcripts for gene {gene}: {gene_canon_ids}')
#     if len(gene_canon_ids) == 0:
#         warn(f'Warning: no canonical transcripts for gene {gene}')
#     else:
#         trx_cnt += 1
#         # Check if coding. Use CDS and not exons.
#         for exon_id in data.exon_ids_of_transcript_id(gene_canon_ids[0]):
#             exon_cnt += 1
#             exon = data.exon_by_id(exon_id)
#             print('\t'.join([exon.contig, str(exon.start - 1), str(exon.end), exon.gene_name]))
#     if gene_cnt % 1000 == 0:
#         warn(f'Processed {gene_cnt} genes, written {trx_cnt} tanscripts, {exon_cnt} exons...')
# warn(f'Done. Processed {gene_cnt} genes, written {trx_cnt} tanscripts, {exon_cnt} exons')
#
#
#
#
# # gtf_fpath = get_ref_file(GENOME, key='gtf')
# # db = gtf.get_gtf_db(gtf_fpath)
#
# # from BCBio import GFF
# # with open(gtf_fpath) as f:
# #     for rec in GFF.parse(f):
# #         print(rec)
#
# # def get_only_canonical_filter():
# #     return lambda x: x[BedCols.ENSEMBL_ID] == canon_tx_by_gname.get(x[BedCols.HUGO]) or \
# #                      x[BedCols.ENSEMBL_ID] == canon_tx_by_gname.get(x[BedCols.GENE])
#
# # def _get_attr(_rec, _key):
# #     val = _rec.attributes.get(_key)
# #     if val is None:
# #         return None
# #     assert len(val) == 1, (_key, str(val))
# #     return val[0]
# #
# # features = []
# # for rec in db.all_features(order_by=('seqid', 'start', 'end')):
# #     if rec.end - rec.start < 0: continue
# #     if rec.featuretype == 'CDS':
# #         gene = _get_attr(rec, 'transcript_id')
# #         trx = _get_attr(rec, 'transcript_id')
# #         if canon_by_gname.get(gene) == trx:
# #             features.append([rec.chrom,
# #                              str(rec.start - 1),
# #                              str(rec.end)])
# #
# # for f in features:
# #     out.write('\t'.join(f) + '\n')

