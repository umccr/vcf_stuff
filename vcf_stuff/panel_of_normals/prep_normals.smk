# Snakemake file for preparing PoN filtering

import sys
import os
import glob
import re
from os.path import join
import csv
from cyvcf2 import VCF, Writer
from ngs_utils.file_utils import get_ungz_gz, splitext_plus, add_suffix, verify_file
from ngs_utils.bcbio import BcbioProject, NoConfigDirException, NoDateStampsException, MultipleDateStampsException
from ngs_utils.logger import warn
from vcf_stuff.panel_of_normals import package_path

from hpc_utils import hpc


def str_to_lua_variable_name(name):
    name = os.path.basename(name)
    name = name.replace('.vcf.gz', '')
    name = re.sub('[^0-9a-zA-Z_]', '_', name) # Remove invalid characters
    name = re.sub('^[^a-zA-Z_]+', '_', name)  # Remove leading characters until we find a letter or underscore
    return name


do_merge_multiallelic = True  # True if want to compare SITE, False if ALT


# basename of the resulting file (will be generated 2: with .snps and .indels suffixes)
PON_FILE = f'GRCh37/panel_of_normals.vcf.gz'


def _load_data():
    vcf_by_sample = dict()

    if hpc.name == 'spartan':
        with open(verify_file(join(package_path(), 'normals.tsv'))) as f:
            # reader = csv.DictReader(f, fieldnames=f.readline().strip().split('\t'), delimiter='\t')

            total_bcbio_runs = 0
            total_samples = 0
            cannot_read_project = 0
            found_normals = 0
            not_found_vcf = 0
            found_multiple_vcfs = 0

            for line in f:
                total_bcbio_runs += 1
                bcbio_path, sample_ids = line.strip().split('\t')
                sample_ids = set(sample_ids.split(','))
                total_samples += len(sample_ids)

                try:
                    bcbio = BcbioProject(bcbio_path, silent=True)
                except NoConfigDirException:
                    warn(f'WARN: cannot parse bcbio run {bcbio_path} - no config dir')
                    cannot_read_project += 1
                except NoDateStampsException:
                    warn(f'WARN: cannot parse bcbio run {bcbio_path} - no datestamp dir found')
                    cannot_read_project += 1
                except MultipleDateStampsException:
                    warn(f'WARN: cannot parse bcbio run {bcbio_path} - multiple datestamp dirs found')
                    cannot_read_project += 1
                else:
                    normals = []
                    for b in bcbio.batch_by_name.values():
                        if b.normal:
                            if b.normal.name not in sample_ids:
                                print(f'WARN: {b.normal.name} not in requested normals.tsv samples for project {bcbio_path}: {sample_ids}')
                            else:
                                normals.append(b.normal)
                                if b.normal.name == 'PRJ180156_1567_8073315T':
                                    print('Found PRJ180156_1567_8073315T in normals in project ', bcbio_path)
                    if not normals:
                        warn(f'WARN: not found normals in run {bcbio_path}')
                    for n in normals:
                        found_vcf = glob.glob(join(bcbio.date_dir, f'{n.name}*-ensemble*.vcf.gz'))
                        if not found_vcf:
                            warn(f'WARN: not found VCF for normal {n.name}, run {bcbio_path}')
                            not_found_vcf += 1
                        elif len(found_vcf) > 1:
                            warn(f'WARN: Found multiple VCF for normal {n.name}: {found_vcf}, run {bcbio_path}')
                            found_multiple_vcfs += 1
                        else:
                            vcf_by_sample[n.name] = found_vcf[0]
                            found_normals += 1
                        # normal = b.normal if b.normal else b.tumor
                        # if normal.name == row['SampleName'] or normal.name == row['SampleID']:
                        #     found_a_sample = True
                        #     found_vcf = glob.glob(join(bcbio.date_dir, f'{normal.name}*-ensemble*.vcf.gz'))
                    # if not found_a_sample:
                    #     warn(f'WARN: Not found sample "{row["SampleName"]}"/"{row["SampleID"]}" for run {row["Results"]}')
                    #     not_found_in_run += 1

        print(f'Done. From {total_bcbio_runs} bcbio, found {found_normals} normal VCFs from {total_samples} samples in normals.csv. '
              f'For {cannot_read_project} bcbio run(s), could not parse folder structure. '
              f'For {not_found_vcf} sample(s), could not find normal VCFs; '
              f'For {found_multiple_vcfs} sample(s), found multiple normal VCFs.')

    if hpc.name == 'vlad':
        vcf_by_sample = {
            'MB_normal_50x': '/Users/vsaveliev/git/umccr/vcf_stuff/vcf_stuff/panel_of_normals/test/test-2.vcf.gz',
            'MB': '/Users/vsaveliev/git/umccr/vcf_stuff/vcf_stuff/panel_of_normals/test/test-benchmark.vcf.gz',
            'MBmore': '/Users/vsaveliev/git/umccr/vcf_stuff/vcf_stuff/panel_of_normals/test/MB-more.vcf.gz'
        }

    # cohort_by_sample = dict()
    # vcf_by_sample = dict()
    #
    # for entry in data:
    #     if isinstance(entry, str):  # VCF path
    #         sname = VCF(entry).samples[0]
    #         cohort_by_sample[sname] = sname
    #         vcf_by_sample[sname] = entry
    #     else:  # single-item dict {cohort: [VCF paths]}
    #         cohort_name, vcf_paths = list(entry.items())[0]
    #         for vcf_path in vcf_paths:
    #             sname = VCF(vcf_path).samples[0]
    #             cohort_by_sample[sname] = cohort_name
    #             vcf_by_sample[sname] = vcf_path

    return vcf_by_sample

vcf_by_sample = _load_data()


rule all:
    input:
        add_suffix(PON_FILE, 'snps'),
        add_suffix(PON_FILE, 'indels'),
        add_suffix(PON_FILE, 'snps').replace('GRCh37', 'hg38'),
        add_suffix(PON_FILE, 'indels').replace('GRCh37', 'hg38'),


rule clean_vcf:
    input:
        vcf = lambda wc: vcf_by_sample[wc.sample]
    output:
        vcf = 'work/clean/{sample}.clean.vcf.gz',
        tbi = 'work/clean/{sample}.clean.vcf.gz.tbi'
    shell:
        'if [ ! -e {input.vcf}.tbi ] ; then tabix -p vcf {input.vcf} ; fi ; ' \
        'bcftools annotate -x INFO,FORMAT {input.vcf} -Oz -o {output.vcf} && ' \
        'tabix -p vcf {output.vcf}'


rule normalise_vcf:
    input:
        rules.clean_vcf.output.vcf
    output:
        vcf = 'work/clean/{sample}.clean.norm.vcf.gz',
        tbi = 'work/clean/{sample}.clean.norm.vcf.gz.tbi'
    shell:
        'norm_vcf {input} -o {output.vcf}'


# TODO:
# Try: take PCGR `<batch>-??????-normal.pcgr.snvs_indels.tiers.tsv`, remove TIER1 and TIER2


# Merge VCFs, make multiallelic indels as we will ignore ALT for indels
rule combine_vcfs:
    input:
        expand(rules.normalise_vcf.output.vcf, sample=vcf_by_sample.keys()),
    output:
        vcf = 'work/clean_merged.vcf.gz',
        tbi = 'work/clean_merged.vcf.gz.tbi'
    threads: 32
    shell:
        'bcftools merge -m indels {input} --threads {threads} -Oz -o {output.vcf} && tabix -p vcf {output.vcf}'
        # 'bcftools concat {input} -d all --allow-overlaps --threads {threads} -Oz -o {output}'


rule count_hits:
    input:
        rules.combine_vcfs.output[0]
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

rule to_hg19:  # need to do that before converting to hg38
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
        vcf = temp(rules.to_hg38_unsorted.output.vcf),
        hg38_noalt_bed = join(hpc.extras, 'hg38_noalt.bed'),
    output:
        vcf = add_suffix(PON_FILE, '{type}').replace('GRCh37', 'hg38')
    shell:
        'bcftools view {input.vcf} -T {input.hg38_noalt_bed}' \
        ' | bcftools sort -Oz -o {output.vcf}' \
        ' && tabix -p vcf {output.vcf}'













