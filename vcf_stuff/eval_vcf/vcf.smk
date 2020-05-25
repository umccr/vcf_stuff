# Snakemake file for evaluation of a set of VCF files against a truth set.
# Handles multiallelic, normalization, inconsistent sample names.

# Usage:
# snakemake -p --configfile=config.yaml

import os
import gzip
import csv

from ngs_utils.file_utils import add_suffix, open_gzipsafe
from ngs_utils.vcf_utils import get_tumor_sample_name
from reference_data import api as refdata
from vcf_stuff.vcf_normalisation import make_normalise_cmd
from vcf_stuff.evaluation import dislay_stats_df, f_measure
from vcf_stuff.eval_vcf import vcf_stats_to_df


FAST = config.get('fast', False)

T_NAME = config.get('tumor_vcf_sample')
N_NAME = config.get('normal_vcf_sample')


rule all:
    input:
        'eval/report.tsv'
    output:
        'report.tsv'
    shell:
        'ln -s {input} {output}'


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

# def get_regions_input(_):
#     ret = [f for f in (truth_regions, samples_regions) if f]
#     print("ret: " + str(ret))
#     return ret

# rule prep_regions:
#     input:
#         get_regions_input
#     output:
#         'regions/regions.bed'
#     run:
#         if samples_regions and truth_regions:
#             shell('bedops -i <(sort-bed {truth_regions}) <(sort-bed {samples_regions}) > {output}')
#         elif truth_regions:
#             shell('sort-bed {truth_regions} > {output}')
#         elif samples_regions:
#             shell('sort-bed {samples_regions} > {output}')

##########################
######### NARROW #########
# Extracts target regions and PASSed calls, and removing fields from bcbio that cause errors like:
# "Error: wrong number of fields in INFO/af_exac_all at 1:26646611, expected 2, found 1"
rule narrow_samples_to_regions_and_pass:
    input:
        lambda wildcards: config['samples'][wildcards.sample]
    output:
        'narrow/{sample}.regions.pass.vcf.gz'
    run:
        regions = merge_regions()
        regions = ('-T ' + regions) if regions else ''
        rm_cmd = ''
        if not FAST:
            anno_to_remove = [
                'ac_exac_all',
                'an_exac_all',
                'ac_adj_exac_afr',
                'an_adj_exac_afr',
                'ac_adj_exac_amr',
                'an_adj_exac_amr',
                'ac_adj_exac_eas',
                'an_adj_exac_eas',
                'ac_adj_exac_fin',
                'an_adj_exac_fin',
                'ac_adj_exac_nfe',
                'an_adj_exac_nfe',
                'ac_adj_exac_oth',
                'an_adj_exac_oth',
                'ac_adj_exac_sas',
                'an_adj_exac_sas',
                'num_exac_Het',
                'num_exac_Hom',
                'rs_ids',
                'fitcons',
                'encode_consensus_gm12878',
                'encode_consensus_h1hesc',
                'encode_consensus_helas3',
                'encode_consensus_hepg2',
                'encode_consensus_huvec',
                'encode_consensus_k562',
                'rmsk',
                'hapmap1',
                'hapmap2',
                'stam_mean',
                'stam_names',
                'af_exac_all',
                'af_adj_exac_afr',
                'af_adj_exac_amr',
                'af_adj_exac_eas',
                'af_adj_exac_fin',
                'af_adj_exac_nfe',
                'af_adj_exac_oth',
                'af_adj_exac_sas',
                'max_aaf_all',
            ]
            anno_to_remove_found = []
            with open_gzipsafe(input[0]) as f:
                for l in f:
                    if l.startswith('##'):
                        for to_rm in anno_to_remove:
                            if '##INFO=<ID=' + to_rm in l:
                                anno_to_remove_found.append(to_rm)
                    if not l.startswith('#'):
                        break
            if anno_to_remove_found:
                rm_cmd = ' -Ou | bcftools annotate -x ' + ','.join(f'INFO/{to_rm}' for to_rm in anno_to_remove_found)
        shell('bcftools view {input} {regions} -f.,PASS' + rm_cmd + ' -Oz -o {output}')

prev_rule = rules.narrow_samples_to_regions_and_pass

if config.get('anno_dp_af'):
    # Propagate FORMAT fields into INFO (using NORMAL_ prefix for normal samples matches):
    rule anno_dp_af:
        input:
            rules.narrow_samples_to_regions_and_pass.output[0]
        output:
            'narrow/{sample}.regions.pass.anno.vcf.gz'
        params:
            t_name = T_NAME,
            n_name = N_NAME,
        run:
            t_name_arg = f"-tn {params.t_name}" if params.t_name else ""
            n_name_arg = f"-nn {params.n_name}" if params.n_name else ""
            shell(f'pcgr_prep {t_name_arg} {n_name_arg} {input}'
                  f' | bgzip -c > {output} && tabix -f -p vcf {output}')
    prev_rule = rules.anno_dp_af

elif config.get('remove_anno'):
    rule remove_anno:
        input:
            rules.narrow_samples_to_regions_and_pass.output[0]
        output:
            'narrow/{sample}.regions.pass.clean.vcf.gz'
        shell:
            'bcftools annotate -x INFO,FORMAT {input} | bgzip -c > {output}'
    prev_rule = rules.remove_anno

# Extract only tumor sample and tabix:
rule narrow_samples_to_tumor_sample:
    input:
        prev_rule.output[0]
    output:
        add_suffix(prev_rule.output[0], 'tumor')
    run:
        sn = get_tumor_sample_name(input[0])
        sn_opt = ''
        if sn:
            sn_opt = f'-s {sn} '
        shell('bcftools view {sn_opt} {input} -Oz -o {output} && tabix -p vcf -f {output}')

#### TRUTH: extract target, PASSed and tumor sample from the truth VCF:
rule narrow_truth_to_target:
    input:
        config['truth_variants']
    output:
        'narrow/truth_variants.vcf.gz'
    run:
        regions = merge_regions()
        regions = ('-T ' + regions) if regions else ''
        rm_cmd = ''
        if not FAST:
            rm_cmd = ' -Ou | bcftools annotate -x INFO,FORMAT'
        shell('bcftools view {input} {regions}' + rm_cmd + ' -Oz -o {output} && tabix -p vcf -f {output}')

############################
######### NORMALSE #########
# Normalise query VCFs:
prev_sample_rule = rules.narrow_samples_to_tumor_sample
prev_truth_rule = rules.narrow_truth_to_target
if not FAST:
    rule normalise_sample:
        input:
            vcf = prev_sample_rule.output[0],
            ref = config['reference_fasta']
        output:
            vcf = 'normalise/{sample}/{sample}.vcf.gz',
            tbi = 'normalise/{sample}/{sample}.vcf.gz.tbi'
        shell:
            make_normalise_cmd('{input.vcf}', '{output[0]}', '{input.ref}')
    prev_sample_rule = rules.normalise_sample

    # Normalise truth VCFs:
    rule normalise_truth:
        input:
            vcf = prev_truth_rule.output[0],
            ref = config['reference_fasta']
        output:
            vcf = 'normalise/truth_variants.vcf.gz',
            tbi = 'normalise/truth_variants.vcf.gz.tbi'
        shell:
            make_normalise_cmd('{input.vcf}', '{output[0]}', '{input.ref}')
    prev_truth_rule = rules.normalise_truth

##########################
####### Annotation #######
if config.get('anno_pon'):
    rule anno_pon_sample:
        input:
            prev_sample_rule.output[0]
        output:
            'anno_pon/{sample}/{sample}.vcf.gz'
        shell:
            'pon_anno {input} -o {output} && tabix -p vcf {output}'
    prev_sample_rule = rules.anno_pon_sample

    rule anno_pon_truth:
        input:
            prev_truth_rule.output[0]
        output:
            'anno_pon/truth.pon.vcf.gz'
        shell:
            'pon_anno {input} -o {output} && tabix -p vcf {output}'
    prev_truth_rule = rules.anno_pon_truth

if config.get('anno_tricky'):
    # Overlap normalised calls with tricky regions and annotate into INFO:
    rule prep_giab_bed:
        input:
            refdata.get_ref_file(config['genome'], ['truth_sets', 'giab', 'bed'])
        output:
            'anno_tricky_regions/giab_conf.bed.gz'
        shell:
            'cat {input} | bgzip -c > {output} && tabix -p bed {output}'

    rule prep_tricky_toml:
        input:
            tricky_bed = refdata.get_ref_file(run.genome_build, key='tricky'),
            giab_conf_bed = rules.prep_giab_bed.output[0]
        output:
            'anno_tricky_regions/tricky_vcfanno.toml'
        params:
            toml_text  = lambda wc, input, output: f'''
[[annotation]]
file="{input.tricky_bed}"
names=["TRICKY"]
columns=[4]
ops=["self"]

[[annotation]]
file="{input.giab_conf_bed}"
names=["GIAB_CONF"]
columns=[3]
ops=["flag"]
'''.replace('\n', r'\\n').replace('"', r'\"'),
        shell:
            'printf "{params.toml_text}" > {output}'

    rule anno_tricky_sample:
        input:
            vcf = prev_sample_rule.output[0],
            toml = rules.prep_tricky_toml.output
        output:
            'anno_tricky_regions/{sample}/{sample}.vcf.gz'
        shell:
            'vcfanno {input.toml} {input.vcf} | bgzip -c > {output} && tabix -p vcf -f {output}'
    prev_sample_rule = rules.anno_tricky_sample

    # Overlap normalised truth calls with tricky regions and annotate into INFO:
    rule anno_tricky_truth:
        input:
            vcf = prev_truth_rule.output[0],
            toml = rules.prep_tricky_toml.output
        output:
            'anno_tricky_regions/truth_variants.vcf.gz'
        shell:
            'vcfanno {input.toml} {input.vcf} | bgzip -c > {output} && tabix -p vcf -f {output}'
    prev_truth_rule = rules.anno_tricky_truth

if config.get('anno_gnomad'):
    rule prep_gn_toml:
        input:
            gnomad_vcf = refdata.get_ref_file(config['genome'], 'gnomad'),
            giab_conf_bed = rules.prep_giab_bed.output[0]
        output:
            'anno_gnomad/gnomad_vcfanno.toml'
        params:
            toml_text = lambda wc, input, output: f'''
[[annotation]]
file="{input.gnomad_vcf}"
fields = ["AF",]
names = ["gnomAD_AF"]
ops=["self"] 
'''.replace('\n', r'\\n').replace('"', r'\"'),
        shell:
            'printf "{params.toml_text}" > {output}'

    rule anno_gnomad_sample:
        input:
            toml = rules.prep_gn_toml.output[0],
            vcf = prev_sample_rule.output[0]
        output:
            'anno_gnomad/{sample}/{sample}.vcf.gz'
        shell:
            'vcfanno {input.toml} {input.vcf} | bgzip -c > {output} && tabix -p vcf {output}'
    prev_sample_rule = rules.anno_gnomad_sample

    rule anno_gnomad_truth:
        input:
            toml = rules.prep_gn_toml.output[0],
            vcf = prev_truth_rule.output[0]
        output:
            'anno_gnomad/truth_variants.vcf.gz'
        shell:
            'vcfanno {input.toml} {input.vcf} | bgzip -c > {output} && tabix -p vcf -f {output}'
    prev_truth_rule = rules.anno_gnomad_truth

############################
######### EVALUATE #########
# Run bcftools isec to get separate VCFs with TP, FN and FN:
rule bcftools_isec:
    input:
        sample_vcf = prev_sample_rule.output[0],
        truth_vcf = prev_truth_rule.output[0]
    params:
        output_dir = 'eval/{sample}_bcftools_isec'
    output:
        fp = 'eval/{sample}_bcftools_isec/0000.vcf',
        fn = 'eval/{sample}_bcftools_isec/0001.vcf',
        tp = 'eval/{sample}_bcftools_isec/0002.vcf'
    run:
        shell('bcftools isec {input.sample_vcf} {input.truth_vcf} -p {params.output_dir}')

def count_variants(vcf):
    snps = set()
    indels = set()
    with (gzip.open(vcf) if vcf.endswith('.gz') else open(vcf)) as f:
        for l in [l for l in f if not l.startswith('#')]:
            chrom, pos, _, ref, alt = l.split('\t')[:5]
            if len(ref) == len(alt) == 1:
                snps.add((chrom, pos, ref, alt))
            else:
                indels.add((chrom, pos, ref, alt))
    return len(snps), len(indels)

# Count TP, FN and FN VCFs to get stats for each sample:
rule eval:
    input:
        fp = rules.bcftools_isec.output.fp,
        fn = rules.bcftools_isec.output.fn,
        tp = rules.bcftools_isec.output.tp
    output:
        'eval/{sample}_stats.tsv'
    run:
        fp_snps, fp_inds = count_variants(input.fp)
        fn_snps, fn_inds = count_variants(input.fn)
        tp_snps, tp_inds = count_variants(input.tp)

        with open(output[0], 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([
                '#SNP TP', 'SNP FP', 'SNP FN', 'SNP Recall', 'SNP Precision',
                 'IND TP', 'IND FP', 'IND FN', 'IND Recall', 'IND Precision'
            ])
            # https://en.wikipedia.org/wiki/Precision_and_recall :
            # precision                  = tp / (tp + fp)
            # recall = sensitivity = tpr = tp / (tp + fn)

            snps_truth = tp_snps + fn_snps
            snps_recall = tp_snps / snps_truth if snps_truth else 0
            snps_called = tp_snps + fp_snps
            snps_prec = tp_snps / snps_called if snps_called else 0

            inds_truth = tp_inds + fn_inds
            inds_recall = tp_inds / inds_truth if inds_truth else 0
            inds_called = tp_inds + fp_inds
            inds_prec = tp_inds / inds_called if inds_called else 0

            snps_f2 = f_measure(2, snps_prec, snps_recall)
            inds_f2 = f_measure(2, inds_prec, inds_recall)

            writer.writerow([
                snps_truth, tp_snps, fp_snps, fn_snps, snps_recall, snps_prec, snps_f2,
                inds_truth, tp_inds, fp_inds, fn_inds, inds_recall, inds_prec, inds_f2,
            ])

# Combine all stats to get single report:
rule report:
    input:
        stats_files = expand(rules.eval.output, sample=sorted(config['samples'].keys()))
        # sompy_files = expand(rules.sompy.output, sample=sorted(config['samples'].keys()))
    output:
        'eval/report.tsv'
    params:
        samples = sorted(config['samples'].keys())
    run:
        stats_by_sname = dict()
        for stats_file, sname in zip(input.stats_files, params.samples):
            with open(stats_file) as f:
                stats_by_sname[sname] = f.readlines()[1].strip().split('\t')
        df = vcf_stats_to_df(stats_by_sname)

        dislay_stats_df(df)

        # Writing raw data to the TSV file
        with open(output[0], 'w') as out_f:
            df.to_csv(out_f, sep='\t', index=False)


# rule eval:
#     input:
#         rules.index_samples.output,
#         rules.index_truth.output,
#         regions = rules.prep_target.output,
#         truth_vcf = rules.narrow_truth_to_target.output,
#         sample_vcf = rules.narrow_samples_to_target.output
#     output:
#         '{sample}/{sample}.re.a/weighted_roc.tsv.gz'
#     shell:
#         '{rtgeval}/run-eval -s {sdf}'
#         ' -b {input.regions}'
#         ' {input.truth_vcf}'
#         ' {input.sample_vcf}'

# rule count_truth:
#     input:
#         truth_variants
#     output:
#         snps = 'truth.snps',
#         indels = 'truth.indels'
#     run:
#         snps, indels = count_variants(truth_variants)
#         with open(output.snps, 'w') as o:
#             o.write(snps)
#         with open(output.indels, 'w') as o:
#             o.write(indels)

# rule report:
#     input:
#         stats_files = expand(rules.eval.output, sample=config['samples'].keys()),
#         truth_snps = rules.count_truth.output.snps,
#         truth_indels = rules.count_truth.output.indels
#     output:
#         'report.tsv'
#     params:
#         samples = config['samples']
#     run:
#         truth_snps = int(open(input.truth_snps).read())
#         truth_indels = int(open(input.truth_indels).read())

#         out_lines = []
#         out_lines.append(['', 'SNP', ''  , ''  , 'INDEL', ''  , ''  ])
#         out_lines.append(['', 'TP' , 'FP', 'FN', 'TP'   , 'FP', 'FN'])

#         for stats_file, sname in zip(input.stats_files, params.samples):
#             data = defaultdict(dict)
#             with open(stats_file) as f:
#                 for l in f:
#                     if l:
#                         event_type, change_type, metric, val = l.strip().split()[:4]
#                         if event_type == 'allelic':
#                             try:
#                                 val = int(val)
#                             except ValueError:
#                                 val = float(val)
#                             data[change_type][metric] = val
#             pprint.pprint(data)
#             try:
#                 out_lines.append([sname, truth_snps   - data['SNP']['FN'],   data['SNP']['FP'],   data['SNP']['FN'],
#                                          truth_indels - data['INDEL']['FN'], data['INDEL']['FP'], data['INDEL']['FN']])
#             except KeyError:
#                 print('Some of the required data for ' + sname + ' not found in ' + fp)

#         with open(output[0], 'w') as out_f:
#             for fields in out_lines:
#                 print(fields)
#                 out_f.write('\t'.join(map(str, fields)) + '\n')



