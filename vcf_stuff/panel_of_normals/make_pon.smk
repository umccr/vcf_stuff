import sys
import os
import glob
import re
from os.path import join
from cyvcf2 import VCF, Writer
from ngs_utils.file_utils import get_ungz_gz, splitext_plus, add_suffix, verify_file, safe_mkdir
from ngs_utils.bcbio import BcbioProject, NoConfigDirException, NoDateStampsException, MultipleDateStampsException
from ngs_utils.logger import warn
from vcf_stuff.panel_of_normals import package_path
from hpc_utils import hpc


do_merge_multiallelic = True  # True if want to compare SITE, False if ALT


if 'genomes_dir' in config:
    hpc.genomes_dir = config.get('genomes_dir')


PON_FILE = f'GRCh37/pon.vcf.gz'  # basename of the resulting file (will be generated 2: with .snps and .indels suffixes)


normals_tsv = 'normals_test.tsv'
if hpc.name == 'spartan' and config.get('test') != 'yes':
    normals_tsv = 'normals.tsv'
normals_tsv = verify_file(join(package_path(), normals_tsv), is_critical=True)


rule all:
    input:
        add_suffix(PON_FILE, 'snps'),
        add_suffix(PON_FILE, 'indels'),
        add_suffix(PON_FILE, 'snps').replace('GRCh37', 'hg38'),
        add_suffix(PON_FILE, 'indels').replace('GRCh37', 'hg38'),


checkpoint load_data:
    input:
        normals_tsv
    output:
        'work/bams.tsv'
    run:
        def _find_bcbio_run(bcbio_path):
            run = None

            def _find_handle_datestamps(bp):
                try:
                    run = BcbioProject(bp, silent=True)
                except NoDateStampsException:
                    warn(f'WARN: cannot parse bcbio run {bp} - no datestamp dir found')
                except MultipleDateStampsException:
                    warn(f'WARN: cannot parse bcbio run {bp} - multiple datestamp dirs found')
                else:
                    return run

            try:
                run = _find_handle_datestamps(bcbio_path)
            except NoConfigDirException:
                subdirs = os.listdir(bcbio_path)
                if len(subdirs) == 1:
                    bcbio_path = join(bcbio_path, subdirs[0])
                    try:
                        run = _find_handle_datestamps(bcbio_path)
                    except NoConfigDirException:
                        warn(f'WARN: cannot parse bcbio run {bcbio_path} - no config dir')
            return run

        bam_by_sample = dict()

        with open(input[0]) as f:
            # reader = csv.DictReader(f, fieldnames=f.readline().strip().split('\t'), delimiter='\t')

            total_bcbio_runs = 0
            total_samples = 0
            cannot_read_project = 0
            found_normals = 0
            not_found_bam = 0

            for line in f:
                total_bcbio_runs += 1
                bcbio_path, sample_ids = line.strip().split('\t')
                sample_ids = set(sample_ids.split(','))
                total_samples += len(sample_ids)

                bcbio = _find_bcbio_run(bcbio_path)
                if not bcbio:
                    cannot_read_project += 1
                    continue
                normals = []
                for b in bcbio.batch_by_name.values():
                    if b.normal:
                        if b.normal.name not in sample_ids:
                            print(f'WARN: {b.normal.name} not in requested normals.tsv samples for project {bcbio.final_dir}: {sample_ids}')
                        else:
                            normals.append(b.normal)
                            if b.normal.name == 'PRJ180156_1567_8073315T':
                                print('Found PRJ180156_1567_8073315T in normals in project ', bcbio.final_dir)
                if not normals:
                    warn(f'WARN: not found normals in run {bcbio.final_dir}')
                for n in normals:
                    if not n.bam or not verify_file(n.bam):
                        not_found_bam += 1
                        warn(f'WARN: not found BAM for normal {n.name}, run {bcbio.final_dir}')
                    else:
                        found_normals += 1
                        bam_by_sample[n.name] = n.bam

        print(f'Done. From {total_bcbio_runs} bcbio, found {found_normals} normal VCFs from {total_samples} samples in normals.csv. '
              f'For {cannot_read_project} bcbio run(s), could not parse folder structure. '
              f'For {not_found_bam} sample(s), not found BAM file.')

        with open(output[0], 'w') as out:
            for sn, bam in bam_by_sample.items():
                out.write(f'{sn}\t{bam}\n')


def recall_with_mutect_input(wc):
    with open(checkpoints.load_data.get().output[0]) as f_:
        bam_by_sample = dict(l.strip().split() for l in f_.readlines())
        return bam_by_sample[wc.sample]

rule recall_with_mutect:
    input:
        bam = recall_with_mutect_input,
        ref_fa = hpc.get_ref_file('GRCh37', key='fa'),
    output:
        vcf = 'work/recall/{sample}.vcf.gz'
    params:
        xms = 2000,
        xmx = 4000,
        tmp_dir = 'work/recall/tmp/{sample}'
    resources:
        mem_mb = 4000
    run:
        safe_mkdir(params.tmp_dir)
        shell(f"""
gatk --java-options "-Xms{params.xms}m -Xmx{params.xmx}m -Djava.io.tmpdir={params.tmp_dir}" \
   Mutect2 \
   -R {input.ref_fa} \
   -I {input.bam} \
   -tumor {wildcards.sample} \
   -O {output.vcf}
""")


# CreateSomaticPanelOfNormals is experimental and doesn't work:
#
#
# rule gatk_make_vcf_list:
#     input:
#         expand(rules.recall_with_mutect.output.vcf, sample=vcf_by_sample.keys()),
#     output:
#         list = 'work/gatk_vcfs.list'
#     run:
#         with open(output.list, 'w') as out:
#             for vcf in input:
#                 out.write(vcf + '\n')
#
#
# rule combine_vcfs_gatk:
#     input:
#         vcfs = expand(rules.recall_with_mutect.output.vcf, sample=vcf_by_sample.keys()),
#         list = rules.gatk_make_vcf_list.output.list
#     output:
#         vcf = 'work/gatk_merged.vcf.gz',
#         tbi = 'work/gatk_merged.vcf.gz.tbi'
#     # threads: 32
#     shell:
#         'gatk --java-options "-Xmx4g" CreateSomaticPanelOfNormals -vcfs {input.list} -O {output.vcf}'


rule clean_vcf:
    input:
        vcf = rules.recall_with_mutect.output.vcf
    output:
        vcf = 'work/clean/{sample}.clean.vcf.gz',
        tbi = 'work/clean/{sample}.clean.vcf.gz.tbi'
    shell:
        'bcftools annotate -x INFO,FORMAT {input.vcf} -Oz -o {output.vcf} '
        '&& tabix -p vcf {output.vcf}'


rule normalise_vcf:
    input:
        rules.clean_vcf.output.vcf
    output:
        vcf = 'work/clean/{sample}.clean.norm.vcf.gz',
        tbi = 'work/clean/{sample}.clean.norm.vcf.gz.tbi'
    shell:
        'norm_vcf {input} -o {output.vcf}'


def combine_vcfs_bcftools_input(wildcards):
    with open(checkpoints.load_data.get().output[0]) as f_:
        bam_by_sample = dict(l.strip().split() for l in f_.readlines())
        return expand(rules.normalise_vcf.output.vcf, sample=bam_by_sample.keys())

# Merge VCFs, make multiallelic indels as we will ignore ALT for indels
rule combine_vcfs_bcftools:
    input:
        combine_vcfs_bcftools_input
    output:
        vcf = 'work/clean_merged.vcf.gz',
        tbi = 'work/clean_merged.vcf.gz.tbi'
    threads: 32
    run:
        assert all(i.endswith('.vcf.gz') for i in input)
        shell('bcftools merge -m indels {input} --threads {threads} -Oz -o {output.vcf} '
              '&& tabix -p vcf {output.vcf}')
        # 'bcftools concat {input} -d all --allow-overlaps --threads {threads} -Oz -o {output}'


rule count_hits:
    input:
        rules.combine_vcfs_bcftools.output[0]
    output:
        get_ungz_gz(PON_FILE)[0]
    run:
        vcf = VCF(input[0])
        # vcf.add_info_to_header({'ID': 'PoN_cohorts', 'Description': 'Panel of normal hits', 'Type': 'Integer', 'Number': '1'})
        vcf.add_info_to_header({'ID': 'PoN_samples', 'Description': 'Panel of normal hits', 'Type': 'Integer', 'Number': '1'})
        w = Writer(output[0], vcf)

        written_lines = 0
        for v in vcf:
            v_samples = [s for s, gt in zip(vcf.samples, v.genotypes) if not any([g == -1 for g in gt])]
            # cohorts = set(cohort_by_sample[s] for s in v_samples)
            # v.INFO['PoN_cohorts'] = len(cohorts)
            v.INFO['PoN_samples'] = len(v_samples)
            w.write_record(v)
            written_lines += 1
        print(f'Written {written_lines} lines')
        w.close()
        # failed_w.close()
        vcf.close()


rule bgzip_and_tabix_final_vcf:
    input:
        rules.count_hits.output[0]
    output:
        vcf = PON_FILE,
        tbi = PON_FILE + '.tbi'
    shell:
        'bgzip {input} && tabix -p vcf {output.vcf}'


rule split_snps_indels:
    input:
        rules.bgzip_and_tabix_final_vcf.output[0]
    output:
        snps_vcf = add_suffix(PON_FILE, 'snps'),
        inds_vcf = add_suffix(PON_FILE, 'indels'),
        snps_tbi = add_suffix(PON_FILE, 'snps') + '.tbi',
        inds_tbi = add_suffix(PON_FILE, 'indels') + '.tbi'
    shell:
        'bcftools view -v snps {input} -Oz -o {output.snps_vcf} && tabix -p vcf {output.snps_vcf} && '
        'bcftools view -V snps {input} -Oz -o {output.inds_vcf} && tabix -p vcf {output.inds_vcf}'

rule to_hg19:  # need to do that before converting to hg38. GRCh-to-hg liftover files don't work, so we have to add `chr` prefixes manually.
    input:
        vcf = add_suffix(PON_FILE, '{type}'),
    output:
        vcf = add_suffix(PON_FILE, '{type}').replace('GRCh37', 'hg19'),
    shell: '''
gunzip -c {input.vcf} \
| py -x "x.replace('##contig=<ID=', '##contig=<ID=chr') if x.startswith('#') else 'chr' + x" \
| py -x "x.replace('chrMT', 'chrM')" \
| grep -v chrG \
| gzip -c > {output.vcf}
'''

rule to_hg38_unsorted:
    input:
        vcf = rules.to_hg19.output.vcf,
        chain = hpc.get_ref_file(key='hg19ToHg38'),
        hg38_fa = hpc.get_ref_file('hg38', key='fa')
    output:
        vcf = temp(add_suffix(PON_FILE, '{type}').replace('GRCh37', 'hg38') + '.unsorted')
    shell:
        'CrossMap.py vcf {input.chain} {input.vcf} {input.hg38_fa} {output.vcf}'

rule to_hg38:
    input:
        vcf = rules.to_hg38_unsorted.output.vcf,
        hg38_noalt_bed = join(hpc.extras, 'hg38_noalt.bed'),
    output:
        vcf = add_suffix(PON_FILE, '{type}').replace('GRCh37', 'hg38')
    shell:
        'bcftools view {input.vcf} -T {input.hg38_noalt_bed}' \
        ' | bcftools sort -Oz -o {output.vcf}' \
        ' && tabix -p vcf {output.vcf}'













