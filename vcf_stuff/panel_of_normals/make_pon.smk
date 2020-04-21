localrules: all, link_bams, link_bam

from ngs_utils.vcf_utils import iter_vcf
import sys
import os
import glob
import re
from os.path import join, dirname, isfile
from cyvcf2 import VCF, Writer
from ngs_utils.file_utils import get_ungz_gz, splitext_plus, add_suffix, verify_file, safe_mkdir
from ngs_utils.bcbio import BcbioProject, NoConfigDirException, NoDateStampsException, MultipleDateStampsException
from ngs_utils.logger import warn
from vcf_stuff.panel_of_normals import package_path
from hpc_utils import hpc


if 'genomes_dir' in config:
    hpc.set_genomes_dir(config.get('genomes_dir'))


bams_tsv = config.get('bams_tsv')
vcfs_tsv = config.get('vcfs_tsv')
bam_by_sample = None
vcf_by_sample = None
if bams_tsv:
    with open(bams_tsv) as f:
        bam_by_sample = dict(l.strip().split() for l in f.readlines())
elif vcfs_tsv:
    with open(vcfs_tsv) as f:
        vcf_by_sample = dict(l.strip().split() for l in f.readlines())


do_merge_multiallelic = True  # True if want to compare SITE, False if ALT
PON_FILE = f'GRCh37/pon.vcf.gz'  # basename of the resulting file (will be generated 2: with .snps and .indels suffixes)


rule all:
    input:
        add_suffix(PON_FILE, 'snps'),
        add_suffix(PON_FILE, 'indels'),
        add_suffix(PON_FILE, 'snps').replace('GRCh37', 'hg38'),
        add_suffix(PON_FILE, 'indels').replace('GRCh37', 'hg38'),


if bam_by_sample is not None:

    rule link_bams:
        input: expand('work/bam/{sample}.bam', sample=bam_by_sample.keys())


    rule link_bam:
        input:
            bam = lambda wc: bam_by_sample[wc.sample]
        output:
            bam = 'work/bam/{sample}.bam',
        shell:
            'ln -s {input.bam} {output.bam}'


    # WARNING: when running in parallel on one cluster node
    # https://gatkforums.broadinstitute.org/gatk/discussion/comment/46465/#Comment_46465
    # > It seems like this may be caused by running out of memory in the docker container. My guess is that '-Xmx32g is
    #   larger than the memory available to docker. Try reducing the -Xmx` value to less than the memory available to
    #   the virtual machine running docker. We saw a very similar problem with exit code 247 and was fixed by increasing
    #   the available memory.
    # Request more memory for the interactive session, or run with -c
    rule recall_with_mutect:
        input:
            bam = rules.link_bam.output[0],
            ref_fa = hpc.get_ref_file('GRCh37', key='fa'),
        output:
            vcf = protected('work/recall/{chrom}/{sample}.vcf.gz')
        params:
            xms = 2000,
            xmx = 4000,
            tmp_dir = 'work/recall/tmp/{chrom}/{sample}',
        resources:
            mem_mb = 4000
        run:
            safe_mkdir(params.tmp_dir)
            shell('samtools index {input.bam} && '
                  'gatk --java-options "-Xms{params.xms}m -Xmx{params.xmx}m -Djava.io.tmpdir={params.tmp_dir}" '
                  'Mutect2 '
                  '-R {input.ref_fa} '
                  '-I {input.bam} '
                  '-tumor {wildcards.sample} '
                  '-O {output.vcf} '
                  '-L {wildcards.chrom}')


rule clean_vcf:
    input:
        vcf = rules.recall_with_mutect.output.vcf if bam_by_sample else lambda wc: vcf_by_sample[wc.sample]
    output:
        vcf = 'work/clean/{chrom}/{sample}.clean.vcf.gz',
        tbi = 'work/clean/{chrom}/{sample}.clean.vcf.gz.tbi'
    shell:
        'bcftools view -f.,PASS {input.vcf} | bcftools annotate -x INFO,FORMAT -Oz -o {output.vcf} '
        '&& tabix -p vcf {output.vcf}'


rule normalise_vcf:
    input:
        rules.clean_vcf.output.vcf
    output:
        vcf = 'work/clean/{chrom}/{sample}.clean.norm.vcf.gz',
        tbi = 'work/clean/{chrom}/{sample}.clean.norm.vcf.gz.tbi'
    shell:
        'norm_vcf {input} -o {output.vcf}'


# Merge VCFs, make multiallelic indels as we will ignore ALT for indels
rule combine_vcfs_bcftools:
    input:
        expand(rules.normalise_vcf.output.vcf.replace('{chrom}', '{{chrom}}'), sample=bam_by_sample.keys()) \
        if bam_by_sample else \
        expand(rules.normalise_vcf.output.vcf, sample=vcf_by_sample.keys(), chrom=['all'])
    output:
        vcf = 'work/{chrom}/clean_merged.vcf.gz',
        tbi = 'work/{chrom}/clean_merged.vcf.gz.tbi'
    threads: 32
    resources:
        mem_mb = 10000
    benchmark:
        'benchmarks/combine_vcfs_bcftools_{chrom}.tsv'
    run:
        assert all(i.endswith('.vcf.gz') for i in input)
        shell('bcftools merge -m indels {input} --threads {threads} -Oz -o {output.vcf} '
              '&& tabix -p vcf {output.vcf}')
        # 'bcftools concat {input} -d all --allow-overlaps --threads {threads} -Oz -o {output}'


rule count_hits:
    input:
        rules.combine_vcfs_bcftools.output[0]
    output:
        'work/{chrom}/pon.vcf.gz'
        # get_ungz_gz(PON_FILE)[0]
    resources:
        mem_mb = 4000
    run:
        def hdr(vcf):
            vcf.add_info_to_header({'ID': 'PoN_samples', 'Description': 'Panel of normal hits',
                                    'Type': 'Integer', 'Number': '1'})
        def func(rec, vcf):
            v_samples = [s for s, gt in zip(vcf.samples, rec.genotypes) if not any([g == -1 for g in gt])]
            rec.INFO['PoN_samples'] = len(v_samples)
            return rec
        iter_vcf(input[0], output[0], proc_rec=func, proc_hdr=hdr)


rule split_snps_indels:
    input:
        rules.count_hits.output[0]
    output:
        # snps_vcf = add_suffix(PON_FILE, 'snps'),
        # inds_vcf = add_suffix(PON_FILE, 'indels'),
        # snps_tbi = add_suffix(PON_FILE, 'snps') + '.tbi',
        # inds_tbi = add_suffix(PON_FILE, 'indels') + '.tbi',
        snps_vcf = 'work/{chrom}/pon.snps.vcf.gz',
        inds_vcf = 'work/{chrom}/pon.indels.vcf.gz',
        snps_tbi = 'work/{chrom}/pon.snps.vcf.gz.tbi',
        inds_tbi = 'work/{chrom}/pon.indels.vcf.gz.tbi',
    resources:
        mem_mb = 4000
    shell:
        'bcftools view -v snps {input} -Oz -o {output.snps_vcf} && tabix -p vcf {output.snps_vcf} && '
        'bcftools view -V snps {input} -Oz -o {output.inds_vcf} && tabix -p vcf {output.inds_vcf}'


# split on snps and indels by chrom
if bam_by_sample:
    rule concat_chroms:
        input:
            expand('work/{chrom}/pon.{{type}}.vcf.gz', chrom = list(map(str, range(1, 22+1))) + ['X', 'Y', 'MT'])
        output:
            vcf = add_suffix(PON_FILE, '{type}')
        resources:
            mem_mb = 4000
        shell:
            'bcftools concat {input} -n -Oz -o {output.vcf} ; tabix -p vcf {output.vcf}'
else:
    rule copy_result:
        input:
            expand('work/all/pon.{{type}}.vcf.gz')
        output:
            vcf = add_suffix(PON_FILE, '{type}')
        shell:
            'cp -r {input} {output} ; cp -r {input}.tbi {output}.tbi'


#################
## Lift over concatenated GRCh37 VCF

rule to_hg19:  # need to do that before converting to hg38. GRCh-to-hg liftover files don't work,
               # so we have to add `chr` prefixes manually.
    input:
        vcf = add_suffix(PON_FILE, '{type}'),
    output:
        vcf = add_suffix(PON_FILE, '{type}').replace('GRCh37', 'hg19'),
    resources:
        mem_mb = 4000
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
        vcf = temp(add_suffix(PON_FILE, '{type}').replace('GRCh37', 'hg38') + '.unsorted'),
        # vcf = 'work/{chrom}/hg38/pon.{type}.vcf.gz.unsorted',
    resources:
        mem_mb = 4000
    shell:
        'CrossMap.py vcf {input.chain} {input.vcf} {input.hg38_fa} {output.vcf}'

rule to_hg38:
    input:
        vcf = rules.to_hg38_unsorted.output.vcf,
        hg38_noalt_bed = join(hpc.extras, 'hg38_noalt.bed'),
    output:
        vcf = add_suffix(PON_FILE, '{type}').replace('GRCh37', 'hg38'),
        # vcf = 'work/{chrom}/hg38/pon.{type}.vcf.gz',
    resources:
        mem_mb = 4000
    shell:
        'bcftools view {input.vcf} -T {input.hg38_noalt_bed}'
        ' | bcftools sort -Oz -o {output.vcf}'
        ' && tabix -p vcf {output.vcf}'












