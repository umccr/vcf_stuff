# Snakemake file for evaluation of a set of VCF files against a truth set.
# Handles multiallelic, normalization, inconsistent sample names.

# Usage:
# snakemake -p --configfile=config.yaml

import os
import gzip
import csv

from ngs_utils.file_utils import add_suffix, open_gzipsafe
from ngs_utils.vcf_utils import get_tumor_sample_name
from hpc_utils.hpc import get_loc, get_ref_file
from vcf_stuff.filtering import get_gnomad_lua
from vcf_stuff.vcf_normalisation import make_normalise_cmd
from vcf_stuff.evaluation import dislay_stats_df, f_measure
from vcf_stuff.eval_vcf import vcf_stats_to_df


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
        max_aaf_all_defined = False
        with open_gzipsafe(input[0]) as f:
            for l in f:
                if l.startswith('##'):
                    if 'max_aaf_all' in l:
                        max_aaf_all_defined = True
                        break
                if not l.startswith('#'):
                    max_aaf_all_defined = False
                    break
        rm_cmd = ''
        if max_aaf_all_defined:
            rm_cmd = ' -Ov | bcftools annotate -x INFO/ac_exac_all,INFO/an_exac_all,INFO/ac_adj_exac_afr,INFO/an_adj_exac_afr,INFO/ac_adj_exac_amr,INFO/an_adj_exac_amr,INFO/ac_adj_exac_eas,INFO/an_adj_exac_eas,INFO/ac_adj_exac_fin,INFO/an_adj_exac_fin,INFO/ac_adj_exac_nfe,INFO/an_adj_exac_nfe,INFO/ac_adj_exac_oth,INFO/an_adj_exac_oth,INFO/ac_adj_exac_sas,INFO/an_adj_exac_sas,INFO/num_exac_Het,INFO/num_exac_Hom,INFO/rs_ids,INFO/fitcons,INFO/encode_consensus_gm12878,INFO/encode_consensus_h1hesc,INFO/encode_consensus_helas3,INFO/encode_consensus_hepg2,INFO/encode_consensus_huvec,INFO/encode_consensus_k562,INFO/rmsk,INFO/hapmap1,INFO/hapmap2,INFO/stam_mean,INFO/stam_names,INFO/af_exac_all,INFO/af_adj_exac_afr,INFO/af_adj_exac_amr,INFO/af_adj_exac_eas,INFO/af_adj_exac_fin,INFO/af_adj_exac_nfe,INFO/af_adj_exac_oth,INFO/af_adj_exac_sas,INFO/max_aaf_all'
        shell('bcftools view {input} {regions} -f .,PASS ' + rm_cmd + ' -Oz -o {output}')

prev_rule = rules.narrow_samples_to_regions_and_pass

if config.get('anno_dp_af'):
    # Propagate FORMAT fields into INFO (using NORMAL_ prefix for normal samples matches):
    rule anno_dp_af:
        input:
            rules.narrow_samples_to_regions_and_pass.output[0]
        output:
            'narrow/{sample}.regions.pass.anno.vcf.gz'
        shell:
            'pcgr_prep {input} | bgzip -c > {output}'

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
        assert sn
        shell('bcftools view -s {sn} {input} -Oz -o {output} && tabix -p vcf -f {output}')

#### TRUTH: extract target, PASSed and tumor sample from the truth VCF:
rule narrow_truth_to_target:
    input:
        config['truth_variants']
    output:
        'narrow/truth_variants.vcf.gz'
    run:
        regions = merge_regions()
        regions = ('-T ' + regions) if regions else ''
        shell('bcftools view {input} {regions} -Ou | bcftools annotate -x INFO,FORMAT -Oz -o {output} && tabix -p vcf -f {output}')

############################
######### NORMALSE #########
# Normalise query VCFs:
rule normalise_sample:
    input:
        vcf = rules.narrow_samples_to_tumor_sample.output[0],
        ref = config['reference_fasta']
    output:
        vcf = 'normalise/{sample}/{sample}.vcf.gz',
        tbi = 'normalise/{sample}/{sample}.vcf.gz.tbi'
    shell:
        make_normalise_cmd('{input.vcf}', '{output[0]}', '{input.ref}')

# Normalise truth VCFs:
rule normalise_truth:
    input:
        vcf = rules.narrow_truth_to_target.output[0],
        ref = config['reference_fasta']
    output:
        vcf = 'normalise/truth_variants.vcf.gz',
        tbi = 'normalise/truth_variants.vcf.gz.tbi'
    shell:
        make_normalise_cmd('{input.vcf}', '{output[0]}', '{input.ref}')

##########################
####### Annotation #######
prev_sample_rule = rules.normalise_sample
prev_truth_rule = rules.normalise_truth
if config.get('anno_pon'):
    rule anno_pon_sample:
        input:
            rules.normalise_sample.output[0]
        output:
            'anno_pon/{sample}/{sample}.vcf.gz'
        shell:
            'pon_anno {input} -o {output} && tabix -p vcf {output}'
    prev_sample_rule = rules.anno_pon_sample

    rule anno_pon_truth:
        input:
            rules.normalise_truth.output[0]
        output:
            'anno_pon/truth.pon.vcf.gz'
        shell:
            'pon_anno {input} -o {output} && tabix -p vcf {output}'
    prev_truth_rule = rules.anno_pon_truth

if config.get('anno_tricky'):
    # Overlap normalised calls with tricky regions and annotate into INFO:
    loc = get_loc()

    rule prep_giab_bed:
        input:
            get_ref_file(config['genome'], ['truth_sets', 'giab', 'bed'])
        output:
            'anno_tricky_regions/giab_conf.bed.gz'
        shell:
            'cat {input} | bgzip -c > {output} && tabix -p bed {output}'

    rule prep_tricky_toml:
        input:
            tricky_bed = os.path.join(loc.extras, 'GRCh37_tricky.bed.gz'),
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
            gnomad_vcf = get_ref_file(config['genome'], 'gnomad'),
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
    snps = 0
    indels = 0
    with (gzip.open(vcf) if vcf.endswith('.gz') else open(vcf)) as f:
        for l in [l for l in f if not l.startswith('#')]:
            _, _, _, ref, alt = l.split('\t')[:5]
            if len(ref) == len(alt) == 1:
                snps += 1
            else:
                indels += 1
    return snps, indels

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

            snps_f1 = f_measure(1, snps_prec, snps_recall)
            snps_f2 = f_measure(2, snps_prec, snps_recall)
            snps_f3 = f_measure(3, snps_prec, snps_recall)

            inds_f1 = f_measure(1, inds_prec, inds_recall)
            inds_f2 = f_measure(2, inds_prec, inds_recall)
            inds_f3 = f_measure(3, inds_prec, inds_recall)

            writer.writerow([
                snps_truth, tp_snps, fp_snps, fn_snps, snps_recall, snps_prec, snps_f1, snps_f2, snps_f3,
                inds_truth, tp_inds, fp_inds, fn_inds, inds_recall, inds_prec, inds_f1, inds_f2, inds_f3,
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



