import os
import gzip
import csv
import pandas as pd

from ngs_utils.file_utils import add_suffix
from ngs_utils.vcf_utils import get_tumor_sample_name
from ngs_utils.bed_utils import get_chrom_order
from hpc_utils.hpc import get_loc
from vcf_stuff.eval import dislay_stats_df
from vcf_stuff.eval.cnv import cnv_to_bed
from pybedtools import BedTool
from collections import defaultdict



rule all:
    input:
        report = 'eval/report.tsv',
        table = 'eval/table.tsv'
    output:
        report = 'report.tsv',
        table = 'table.tsv'
    shell:
        'ln -s {input.report} {output.report} && '
        'ln -s {input.table} {output.table}'


def merge_regions():
    samples_regions = config.get('sample_regions')
    truth_regions = config['truth_regions'] if 'truth_regions' in config else None

    output = 'narrow/regions.bed'
    if samples_regions and truth_regions:
        shell('bedops -i <(sort-bed {truth_regions}) <(sort-bed {samples_regions}) > {output}')
        return output
    elif truth_regions:
        shell('sort-bed {truth_regions} > {output}')
        return output
    elif samples_regions:
        shell('sort-bed {samples_regions} > {output}')
        return output
    else:
        return None


rule cnv_to_bed_sample:
    input:
        lambda wc: config['samples'][wc.sample]
    output:
        'bed/sample_{sample}.bed'
    run:
        cnv_to_bed(input[0], output[0])

rule cnv_to_bed_truth:
    input:
        lambda wc: config['truth_variants']
    output:
        'bed/truth.bed'
    run:
        cnv_to_bed(input[0], output[0])


anno_cmd = 'annotate_bed.py {input} -o {output} --short --work-dir {params.work_dir} --short -a all' \
           ' --coding-only -g ' + config['genome']

rule annotate_bed_sample:
    input:
        rules.cnv_to_bed_sample.output[0]
    output:
        'bed_anno/sample_{sample}.bed'
    params:
        work_dir = 'work/annotate_bed/sample_{sample}'
    shell:
       anno_cmd

rule annotate_bed_truth:
    input:
        rules.cnv_to_bed_truth.output[0]
    output:
        'bed_anno/truth.bed'
    params:
        work_dir = 'work/annotate_bed/truth'
    shell:
        anno_cmd


chrom_order = get_chrom_order(config['genome'])


def _read_cn_data_by_chrom_gene(bed_fpath):
    """ Reads a BED file in form of:
        chr start     end       gene   |cn|event
        21  42818744  42820836	MX1    ||Del
        1   1         2533000	PRKCZ  |3
        Returns a dict of {"chr:gene" -> list_of_cn}
    """
    cn_data_by_chrom_gene = defaultdict(list)
    for r in BedTool(bed_fpath):
        genes = r[3].split(',')
        g_cn_event = r[4].split('|')
        cn, event = None, None
        if len(g_cn_event) == 2:
            cn = g_cn_event[1]
            cn = int(cn) if cn else None
        if len(g_cn_event) == 3:
            event = g_cn_event[2]
        if cn is not None:
            if cn == 2 and config.get('check_gt') is not True:
                continue
        for g in genes:
            if g and g != '.':
                chrom = r[0]
                cn_data_by_chrom_gene[(chrom, g)].append((cn, event))
    return cn_data_by_chrom_gene
    # sorted_gene_chroms = sorted(cns_by_gene_chrom.keys(), key=lambda k: chr_order[k[0]])
    # cns_by_gene = {c + ':' + g : cns_by_gene_chrom[(c, g)] for (c, g) in sorted_gene_chroms}
    # return cns_by_gene

# def _read_cns_by_gene(bed_fpath):
#     """ Reads a BED file in form of chr,start,end,gene,?|cn
#         Returns a dict of {"chr:gene" -> list_of_cn}
#     """
#     cns_by_gene_chrom = defaultdict(list)
#     for r in BedTool(bed_fpath):
#         genes = r[3].split(',')
#         cn = int(r[4].split('|')[1])
#         for g in genes:
#             if g and g != '.':
#                 cns_by_gene_chrom[(r[0], g)].append(cn)
#     sorted_gene_chroms = sorted(cns_by_gene_chrom.keys(), key=lambda k: chr_order[k[0]])
#     cns_by_gene = {c + ':' + g : cns_by_gene_chrom[(c, g)] for (c, g) in sorted_gene_chroms}
#     return cns_by_gene

def _data_by_gene__to__gene_cn_set(data_by_gene):
    gene_cn_set = set()
    for gene, data in data_by_gene.items():
        for cn, event in data:
            gene_cn_set.add((gene, cn))
    return gene_cn_set

def _cn_to_event(cn):
    cn = int(cn)
    if cn < 2: return 'Del'
    if cn > 2: return 'Amp'
    return 'Mix'

def _data_by_gene__to__gene_event_set(data_by_gene):
    gene_event_set = set()
    for gene, data in data_by_gene.items():
        for cn, event in data:
            if event is None and cn is not None:
                event = _cn_to_event(cn)
            gene_event_set.add((gene, event))
    return gene_event_set

# def _read_gene_set(bed_fpath):
#     return set(r[3] for r in BedTool(bed_fpath))
#
# def _read_gene_cn_set(bed_fpath):
#     return set((r[3], int(r[4].split('|')[1])) for r in BedTool(bed_fpath))

def _aggr_cns(cns):
    cns = sorted(cns, key=lambda cn: abs(cn - 2))
    return cns[-1]

def _metrics_from_sets(truth_set, sample_set):
    tp = len(sample_set & truth_set)
    fp = len(sample_set - truth_set)
    fn = len(truth_set - sample_set)
    truth = len(truth_set)
    called = len(sample_set)
    recall = tp / truth if truth else 0
    prec = tp / called if called else 0
    return truth, tp, fp, fn, recall, prec


rule table:
    input:
        sample_beds = expand('bed_anno/sample_{sample}.bed', sample=sorted(config['samples'].keys())),
        truth_bed = 'bed_anno/truth.bed'
    output:
        'eval/table.tsv'
    params:
        samples = sorted(config['samples'].keys())
    run:
        cn_by_gene_by_sname = defaultdict(dict)
        for sample_bed, sname in zip([input.truth_bed] + input.sample_beds, ['truth'] + params.samples):
            # cn_by_gene_by_sname[sname] = _read_cns_by_chrom_gene(sample_bed)
            # g, cns in _read_cns_by_chrom_gene(sample_bed).items()
            cn_by_gene_by_sname[sname] = {
                g: ', '.join(set([f'{event}:{cn}' if cn is not None else f'{event}' for cn, event in
                             [(cn, event or _cn_to_event(cn)) for (cn, event) in vals]]))
                for g, vals in _read_cn_data_by_chrom_gene(sample_bed).items()
            }
        df = pd.DataFrame(cn_by_gene_by_sname, columns=['truth'] + params.samples)
        # index = [(chr_order[c], g) for (c, g) in [c_g.split(':') for c_g in cn_by_gene_by_sname.keys()]]
        # df = df.reindex(index=[(chr_order[c], g) for (c, g) in [c_g.split(':') for c_g in df.index]])
        # print(df)
        print(df.to_string(index=True, na_rep='.'))
        with open(output[0], 'w') as out_f:
            df.to_csv(out_f, sep='\t', index=True)


# def _stats_to_df(stat_by_sname, include_cn=True):
#     idx = pd.MultiIndex.from_arrays([
#         ['Sample', 'GENE', 'GENE', 'GENE' , 'GENE'  , 'GENE', 'EVENT', 'EVENT', 'EVENT', 'EVENT' , 'EVENT'] + (['CN', 'CN', 'CN', 'CN'    , 'CN'  ] if include_cn else []),
#         [''      , 'TP'  , 'FP'  , 'FN'   , 'Recall', 'Prec', 'TP'   , 'FP'   , 'FN'   , 'Recall', 'Prec' ] + (['TP', 'FP', 'FN', 'Recall', 'Prec'] if include_cn else [])
#         ],
#         names=['1', '2'])
#
#     data = []
#     s_truth = i_truth = None
#     for sname, stats in stat_by_sname.items():
#         [g_truth, g_tp, g_fp, g_fn, g_rec, g_prec], [e_truth, e_tp, e_fp, e_fn, e_rec, e_prec] = stats[0:2]
#         if include_cn:
#             [c_truth, c_tp, c_fp, c_fn, c_rec, c_prec] = stats[2]
#
#         d = {
#             ('Sample', ''): sname,
#             ('EVENT', 'TP'): int(e_tp),
#             ('EVENT', 'FP'): int(e_fp),
#             ('EVENT', 'FN'): int(e_fn),
#             ('EVENT', 'Recall'): float(e_rec),
#             ('EVENT', 'Prec'): float(e_prec),
#             ('GENE', 'TP'): int(g_tp),
#             ('GENE', 'FP'): int(g_fp),
#             ('GENE', 'FN'): int(g_fn),
#             ('GENE', 'Recall'): float(g_rec),
#             ('GENE', 'Prec'): float(g_prec),
#         }
#         if include_cn:
#             d.extend({
#                 ('CN', 'TP'): int(c_tp),
#                 ('CN', 'FP'): int(c_fp),
#                 ('CN', 'FN'): int(c_fn),
#                 ('CN', 'Recall'): float(c_rec),
#                 ('CN', 'Prec'): float(c_prec),
#             })
#         data.append(d)
#     )
#     return pd.DataFrame(data, columns=idx)

def _stats_to_df(stat_by_sname):
    idx = ['Sample', 'TP', 'FP', 'FN', 'Recall', 'Prec']
    data = []
    s_truth = i_truth = None
    for sname, stats in stat_by_sname.items():
        truth, tp, fp, fn, rec, prec = stats
        data.append({
            ('Sample'): sname,
            ('TP'): int(tp),
            ('FP'): int(fp),
            ('FN'): int(fn),
            ('Recall'): float(rec),
            ('Prec'): float(prec)
        })

    return pd.DataFrame(data, columns=idx)

# Combine all stats to get single report:
rule report:
    input:
        # stats_files = expand('eval/{sample}_stats.tsv', sample=sorted(config['samples'].keys()))
        sample_beds = expand('bed_anno/sample_{sample}.bed', sample=sorted(config['samples'].keys())),
        truth_bed  = 'bed_anno/truth.bed',
        t = 'eval/table.tsv'   # we want to show the table before the report
    output:
        'eval/report.tsv'
    params:
        samples = sorted(config['samples'].keys())
    run:
        gene_stats_by_sname = dict()
        event_stats_by_sname = dict()
        cn_stats_by_sname = dict()
        for sample_bed, sname in zip(input.sample_beds, params.samples):
            sample_data_by_gene   = _read_cn_data_by_chrom_gene(sample_bed)
            truth_data_by_gene    = _read_cn_data_by_chrom_gene(input.truth_bed)

            sample_gene_set       = set(sample_data_by_gene.keys())
            truth_gene_set        = set(truth_data_by_gene.keys())
            gene_stats_by_sname[sname] = _metrics_from_sets(truth_gene_set, sample_gene_set)

            sample_gene_event_set = _data_by_gene__to__gene_event_set(sample_data_by_gene)
            truth_gene_event_set  = _data_by_gene__to__gene_event_set(truth_data_by_gene)
            event_stats_by_sname[sname] = _metrics_from_sets(truth_gene_event_set, sample_gene_event_set)

            sample_gene_cn_set    = _data_by_gene__to__gene_cn_set(sample_data_by_gene)
            truth_gene_cn_set     = _data_by_gene__to__gene_cn_set(truth_data_by_gene)
            include_cn = all(cn is not None for (g, cn) in truth_gene_cn_set)
            if include_cn:
                cn_stats_by_sname[sname] = _metrics_from_sets(truth_gene_cn_set, sample_gene_cn_set)

        with open(output[0], 'a') as out_fh:
            df = _stats_to_df(gene_stats_by_sname)
            print('Gene level comparison')
            dislay_stats_df(df)
            out_fh.write('Gene level comparison\n')
            df.to_csv(out_fh, sep='\t', index=False)
            out_fh.write('\n')

            df = _stats_to_df(event_stats_by_sname)
            print('\nEvent level comparison (Amp, Del)')
            dislay_stats_df(df)
            out_fh.write('\nEvent level comparison (Amp, Del)\n')
            df.to_csv(out_fh, sep='\t', index=False)
            out_fh.write('\n')

            if cn_stats_by_sname:
                df = _stats_to_df(cn_stats_by_sname)
                print('\nCN level comparison')
                dislay_stats_df(df)
                out_fh.write('\nCN level comparison\n')
                df.to_csv(out_fh, sep='\t', index=False)
                out_fh.write('\n')


# ##########################
# ######### NARROW #########
# # Extracts target regions and PASSed calls:
# rule narrow_samples_to_regions:
#     input:
#         lambda wildcards: config['samples'][wildcards.sample]
#     output:
#         'narrow/{sample}.regions.pass.vcf.gz'
#     run:
#         regions = merge_regions()
#         regions = ('-T ' + regions) if regions else ''
#         shell('bcftools view {input} {regions} -f .,PASS -Oz -o {output}')


############################
######### EVALUATE #########
# rule overlap_events:
#     input:
#         rules.annotate_bed.output[0]
#     output:
#         'eval/{sample}'
#     run:
#         pass

#
# # Count TP, FN and FN VCFs to get stats for each sample:
# rule eval:
#     input:
#         fp = rules.bcftools_isec.output.fp,
#         fn = rules.bcftools_isec.output.fn,
#         tp = rules.bcftools_isec.output.tp
#     output:
#         'eval/{sample}_stats.tsv'
#     run:
#         fp_snps, fp_inds = count_variants(input.fp)
#         fn_snps, fn_inds = count_variants(input.fn)
#         tp_snps, tp_inds = count_variants(input.tp)
#
#         with open(output[0], 'w') as f:
#             writer = csv.writer(f, delimiter='\t')
#             writer.writerow([
#                 '#SNP TP', 'SNP FP', 'SNP FN', 'SNP Recall', 'SNP Precision',
#                  'IND TP', 'IND FP', 'IND FN', 'IND Recall', 'IND Precision'
#             ])
#             # https://en.wikipedia.org/wiki/Precision_and_recall :
#             # precision                  = tp / (tp + fp)
#             # recall = sensitivity = tpr = tp / (tp + fn)
#
#             snps_truth = tp_snps + fn_snps
#             snps_recall = tp_snps / snps_truth if snps_truth else 0
#             snps_called = tp_snps + fp_snps
#             snps_prec = tp_snps / snps_called if snps_called else 0
#
#             inds_truth = tp_inds + fn_inds
#             inds_recall = tp_inds / inds_truth if inds_truth else 0
#             inds_called = tp_inds + fp_inds
#             inds_prec = tp_inds / inds_called if inds_called else 0
#
#             snps_f1 = f_measure(1, snps_prec, snps_recall)
#             snps_f2 = f_measure(2, snps_prec, snps_recall)
#             snps_f3 = f_measure(3, snps_prec, snps_recall)
#
#             inds_f1 = f_measure(1, inds_prec, inds_recall)
#             inds_f2 = f_measure(2, inds_prec, inds_recall)
#             inds_f3 = f_measure(3, inds_prec, inds_recall)
#
#             writer.writerow([
#                 snps_truth, tp_snps, fp_snps, fn_snps, snps_recall, snps_prec, snps_f1, snps_f2, snps_f3,
#                 inds_truth, tp_inds, fp_inds, fn_inds, inds_recall, inds_prec, inds_f1, inds_f2, inds_f3,
#             ])
#
# # Combine all stats to get single report:
# rule report:
#     input:
#         stats_files = expand(rules.eval.output, sample=sorted(config['samples'].keys()))
#         # sompy_files = expand(rules.sompy.output, sample=sorted(config['samples'].keys()))
#     output:
#         'eval/report.tsv'
#     params:
#         samples = sorted(config['samples'].keys())
#     run:
#         stats_by_sname = dict()
#         for stats_file, sname in zip(input.stats_files, params.samples):
#             with open(stats_file) as f:
#                 stats_by_sname[sname] = f.readlines()[1].strip().split('\t')
#         df = stats_to_df(stats_by_sname)
#
#         dislay_stats_df(df)
#
#         # Writing raw data to the TSV file
#         with open(output[0], 'w') as out_f:
#             df.to_csv(out_f, sep='\t', index=False)

