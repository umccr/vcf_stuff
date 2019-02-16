from os.path import isfile, join, basename
import cyvcf2
import toml
import csv
from ngs_utils.file_utils import which
from ngs_utils.file_utils import get_ungz_gz
from ngs_utils.file_utils import splitext_plus
from hpc_utils.hpc import get_ref_file
from vcf_stuff import iter_vcf, count_vars, vcf_contains_field
import subprocess
from ngs_utils.reference_data import get_key_genes_set


localrules: prep_anno_toml, prep_giab_bed, annotate


#############################
#### Reading parameters #####

SAMPLE = config['sample']
GENOME = config['genome']
GENOMES_DIR = config.get('genomes_dir')
INPUT_VCF = config['input_vcf']
OUTPUT_VCF = config['output_vcf']
assert OUTPUT_VCF.endswith('.vcf.gz'), OUTPUT_VCF
assert INPUT_VCF.endswith('.vcf.gz'), INPUT_VCF

PCGR_ENV_PATH = config.get('pcgr_env_path')
conda_cmd = ''
if PCGR_ENV_PATH:
    conda_cmd = 'export PATH=' + PCGR_ENV_PATH + '/bin:$PATH; '


# REMOVE_GERMLINE = config.get('remove_germline', 'no') == 'yes'


#############################


rule all:
    input:
        vcf = OUTPUT_VCF,
        tbi = OUTPUT_VCF + '.tbi'


# Call variants with 1%
# Apply to VarDict only: (INFO/QUAL * TUMOR_AF) >= 4
# Call ensemble
# First PoN round: remove PoN_CNT>2'
# Annotate with PCGR (VEP+known cancer databases)
#   Tier 1 - variants of strong clinical significance
#   Tier 2 - variants of potential clinical significance
#   Tier 3 - variants of unknown clinical significance
#   Tier 4 - other coding variants
#   Noncoding variants
# Tier 1-3 - keep all variants
# Tier 4 and noncoding - filter with:
#   Remove gnomad_AF <=0.02
#   Remove PoN_CNT>{0 if issnp else 1}'
#   Remove indels in "bad_promoter" tricky regions
#   Remove DP<25 & AF<5% in tricky regions: gc15, gc70to75, gc75to80, gc80to85, gc85, low_complexity_51to200bp, low_complexity_gt200bp, non-GIAB confident, unless coding in cancer genes


rule prep_hmf_hotspots:
    input:
        vcf = get_ref_file(GENOME, key='hotspots', genomes_dir=GENOMES_DIR),
    output:
        vcf = f'somatic_anno/hmf_hotspot.vcf.gz',
        tbi = f'somatic_anno/hmf_hotspot.vcf.gz.tbi',
    shell:
        'bcftools filter -i "HMF=1" {input.vcf} -Oz -o {output.vcf} && tabix -p vcf {output.vcf}'

rule prep_anno_toml:
    input:
        ga4gh_dir       = join(get_ref_file(GENOME, key='problem_regions_dir', genomes_dir=GENOMES_DIR), 'GA4GH'),
        encode          = join(get_ref_file(GENOME, key='problem_regions_dir', genomes_dir=GENOMES_DIR), 'ENCODE', 'wgEncodeDacMapabilityConsensusExcludable.bed.gz'),
        lcr             = join(get_ref_file(GENOME, key='problem_regions_dir', genomes_dir=GENOMES_DIR), 'repeats', 'LCR.bed.gz'),
        segdup          = join(get_ref_file(GENOME, key='problem_regions_dir', genomes_dir=GENOMES_DIR), 'segdup.bed.gz'),
        gnomad_vcf      = get_ref_file(GENOME, key='gnomad', genomes_dir=GENOMES_DIR),
        hmf_hotspots    = rules.prep_hmf_hotspots.output.vcf,
        hmf_giab        = get_ref_file(GENOME, key='hmf_giab_conf', genomes_dir=GENOMES_DIR),
        hmf_mappability = get_ref_file(GENOME, key='hmf_mappability', genomes_dir=GENOMES_DIR),
    output:
        f'somatic_anno/tricky_vcfanno.toml'
    run:
        with open(output[0], 'w') as f:
            f.write(f"""
[[annotation]]
file="{input.gnomad_vcf}"
fields = ["AF_popmax"]
names = ["gnomAD_AF"]
ops = ["self"]

[[annotation]]
file = "{input.hmf_hotspots}"
fields = ["HMF"]
names = ["HMF_HOTSPOT"]
ops = ["flag"]

[[annotation]]
file = "{input.hmf_giab}"
names = ["HMF_GIAB_CONF"]
columns = [3]
ops = ["flag"]

[[annotation]]
file = "{input.hmf_mappability}"
names = ["HMF_MAPPABILITY_float"]
columns = [5]
ops = ["concat"]

[[annotation]]
file = "{input.lcr}"
names = ["TRICKY_LCR"]
columns = [3]
ops = ["flag"]

[[annotation]]
file = "{input.encode}"
names = ["ENCODE"]
columns = [4]
ops = ["concat"]

[[annotation]]
file = "{input.segdup}"
columns = [3]
names = ["SEGDUP"]
ops = ["flag"]

""")
        for fn in os.listdir(join(input.ga4gh_dir)):
            if fn.endswith('.bed.gz'):
                fpath = join(input.ga4gh_dir, fn)
                assert isfile(fpath), fpath
                name = splitext_plus(basename(fpath))[0]
                with open(output[0], 'a') as f:
                    f.write(f"""
[[annotation]]
file = "{fpath}"
columns = [4]
names = ["TRICKY_{name}"]
ops = ["flag"]
""")

rule somatic_vcf_regions_anno:
    input:
        vcf = INPUT_VCF,
        toml = rules.prep_anno_toml.output[0],
    output:
        vcf = f'somatic_anno/vcfanno/{SAMPLE}-somatic.vcf.gz',
        tbi = f'somatic_anno/vcfanno/{SAMPLE}-somatic.vcf.gz.tbi',
    shell:
        'vcfanno {input.toml} {input.vcf} | bgzip -c > {output.vcf} && tabix -f -p vcf {output.vcf}'

# Possibly subset VCF to avoid PCGR choking with R stuff.
# too higly mutated samples might indicate germline contamination.
rule maybe_remove_germline:
    input:
        vcf = rules.somatic_vcf_regions_anno.output.vcf,
        tbi = rules.somatic_vcf_regions_anno.output.tbi,
    output:
        vcf = f'somatic_anno/remove_gnomad/{SAMPLE}-somatic.vcf.gz',
        tbi = f'somatic_anno/remove_gnomad/{SAMPLE}-somatic.vcf.gz.tbi',
    run:
        total_vars = count_vars(input.vcf)
        if total_vars > 500_000:
            shell('bcftools filter -e "gnomAD_AF>=0.01 & HMF_HOTSPOT=0" {input.vcf} -Oz -o {output.vcf} && tabix -f -p vcf {output.vcf}')
        else:
            shell('cp {input.vcf} {output.vcf} ; cp {input.tbi} {output.tbi} ; ')

# If the noise wasn't germline, it might be artefacts/errors from FFPE or ortherwise low quality data.
# subsetting to cancer genes in this case.
rule maybe_subset_cancer_genes:
    input:
        rm_germline_vcf = rules.maybe_remove_germline.output.vcf,
        rm_germline_tbi = rules.maybe_remove_germline.output.tbi,
        full_vcf = rules.somatic_vcf_regions_anno.output.vcf,
        full_tbi = rules.somatic_vcf_regions_anno.output.tbi,
    output:
        vcf = f'somatic_anno/cancer_genes/{SAMPLE}-somatic.vcf.gz',
        tbi = f'somatic_anno/cancer_genes/{SAMPLE}-somatic.vcf.gz.tbi',
        subset_to_cancer = f'somatic_anno/subset_to_cancer_genes.flag',
    run:
        vars_left = int(subprocess.check_output(f'bcftools view -H {input.rm_germline_vcf} | wc -l', shell=True).strip())
        if vars_left > 500_000:
            genes = get_key_genes_set()
            def func(rec):
                if rec.INFO.get('ANN') is not None and rec.INFO['ANN'].split('|')[3] in genes:
                    return rec
            iter_vcf(input.vcf, output.vcf, func)
            shell('echo YES > {output.subset_to_cancer}')
        else:
            shell('cp {input.rm_germline_vcf} {output.vcf} ; cp {input.rm_germline_tbi} {output.tbi} ; ')
            shell('echo NO > {output.subset_to_cancer}')

# Removes TRICKY_ and ANN fields
rule somatic_vcf_clean_info:
    input:
        vcf = rules.maybe_subset_cancer_genes.output.vcf,
        tbi = rules.maybe_subset_cancer_genes.output.tbi,
    output:
        vcf = f'somatic_anno/clean_info/{SAMPLE}-somatic.vcf.gz',
        tbi = f'somatic_anno/clean_info/{SAMPLE}-somatic.vcf.gz.tbi',
    run:
        def proc_hdr(vcf):
            vcf.add_info_to_header({'ID': 'TRICKY', 'Description': 'Tricky regions from bcbio folders at coverage/problem_regions/GA4GH and coverage/problem_regions/LCR', 'Type': 'String', 'Number': '1'})

        def postproc_hdr(raw_hdr):
            new_hdr = []
            for l in raw_hdr.split('\n'):
                if not l.startswith('##INFO=<ID=TRICKY_') and not l.startswith('##INFO=<ID=ANN,'):
                    new_hdr.append(l)
            return '\n'.join(new_hdr)

        def func(rec):
            if rec.INFO.get('ANN') is not None:
                del rec.INFO['ANN']
            tricky_flags = [k.replace('TRICKY_', '') for k, v in rec.INFO if k.startswith('TRICKY_')]
            if tricky_flags:
                rec.INFO['TRICKY'] = '|'.join(tricky_flags)
            for f in tricky_flags:
                del rec.INFO[f'TRICKY_{f}']
            return rec
        iter_vcf(input.vcf, output.vcf, func, proc_hdr=proc_hdr, postproc_hdr=postproc_hdr)

# Preparations: annotate TUMOR_X and NORMAL_X fields for PCGR, remove non-standard chromosomes and mitochondria, remove non-PASSed calls
rule somatic_vcf_prep:
    input:
        vcf = rules.somatic_vcf_clean_info.output.vcf,
    output:
        vcf = f'somatic_anno/prep/{SAMPLE}-somatic.vcf.gz',
        tbi = f'somatic_anno/prep/{SAMPLE}-somatic.vcf.gz.tbi'
    shell:
        'pcgr_prep {input.vcf} | bgzip -c > {output.vcf} && tabix -f -p vcf {output.vcf}'

rule somatic_vcf_pon_anno:
    input:
        vcf = rules.somatic_vcf_prep.output.vcf,
        tbi = rules.somatic_vcf_prep.output.tbi,
    params:
        genome_build = GENOME,
        pon_hits = 3,
        pon_dir = get_ref_file(GENOME, 'panel_of_normals_dir', genomes_dir=GENOMES_DIR)
    output:
        vcf = f'somatic_anno/pon/{SAMPLE}-somatic.vcf.gz',
        tbi = f'somatic_anno/pon/{SAMPLE}-somatic.vcf.gz.tbi',
    shell:
        'pon_anno {input.vcf} --pon-dir {params.pon_dir} | bgzip -c > {output.vcf} && tabix -f -p vcf {output.vcf}'
        # ' | bcftools filter -e "INFO/PoN_CNT>={params.pon_hits}" --soft-filter PoN --mode + -Oz -o {output.vcf}' \
        # ' && tabix -f -p vcf {output.vcf} '

rule somatic_vcf_pcgr_round1:
    input:
        vcf = rules.somatic_vcf_pon_anno.output.vcf,
        pcgr_data = get_ref_file(key='pcgr_data', genomes_dir=GENOMES_DIR),
    output:
        tiers = f'somatic_anno/pcgr_run/{SAMPLE}-somatic.pcgr.snvs_indels.tiers.tsv',
    params:
        output_dir = f'somatic_anno/pcgr_run',
        genome = GENOME,
        sample_name = f'{SAMPLE}-somatic',
        opt='--no-docker' if not which('docker') else '',
    resources:
        mem_mb = 20000
    shell:
        conda_cmd + which('pcgr') +
        ' {input.vcf} -g {params.genome} -o {params.output_dir} -s {params.sample_name} '
        '{params.opt} --pcgr-data {input.pcgr_data}'

rule somatic_vcf_pcgr_anno:
    input:
        vcf = rules.somatic_vcf_pon_anno.output.vcf,
        tiers = rules.somatic_vcf_pcgr_round1.output.tiers,
    output:
        vcf = f'somatic_anno/pcgr_ann/{SAMPLE}-somatic.vcf.gz',
        tbi = f'somatic_anno/pcgr_ann/{SAMPLE}-somatic.vcf.gz.tbi',
    run:
        pcgr_fields_by_snp = dict()
        cosmic_by_snp = dict()
        icgc_by_snp = dict()
        print('Reading PCGR annotation from', input.tiers)
        with open(input.tiers) as f:
            reader = csv.DictReader(f, delimiter='\t', fieldnames=f.readline().strip().split('\t'))
            for row in reader:
                change = row['GENOMIC_CHANGE']
                pcgr_fields_by_snp[change] = dict()
                for k in ['SYMBOL',
                          'TIER',
                          'CONSEQUENCE',
                          'MUTATION_HOTSPOT',
                          'INTOGEN_DRIVER_MUT',
                          'TCGA_PANCANCER_COUNT',
                          'CLINVAR_CLNSIG']:
                    val = row[k]
                    val = val.replace(' ', '_')
                    val = val.replace(',', '|')
                    if val == 'TRUE':
                        val = True
                    if val == 'FALSE':
                        val = False
                    if val and val != 'NA':
                        pcgr_fields_by_snp[change][k] = val

                cosmic = row['COSMIC_MUTATION_ID']
                cosmic_by_snp[change] = len(cosmic.split('&')) if cosmic != 'NA' else 0

                icgc = row['ICGC_PCAWG_OCCURRENCE']
                icgc_by_snp[change] = sum([int(ct.split('|')[1]) for ct in icgc.split(', ')]) if icgc != 'NA' else 0

        def func_hdr(vcf):
            vcf.add_info_to_header({'ID': 'PCGR_SYMBOL', 'Description':
                'VEP gene symbol, reported by PCGR in .snvs_indels.tiers.tsv file',
                                    'Type': 'String', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_TIER', 'Description': 'TIER as reported by PCGR in .snvs_indels.tiers.tsv file. '
                                                                      '1: strong clinical significance; '
                                                                      '2: potential clinical significance; '
                                                                      '3: unknown clinical significance; '
                                                                      '4: other coding variants',                                                                                                  'Type': 'String',  'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_CONSEQUENCE',          'Description': 'VEP consequence, reported by PCGR in .snvs_indels.tiers.tsv file',                                          'Type': 'String',  'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_MUTATION_HOTSPOT',     'Description': 'Mutation hotspot, reported by PCGR in .snvs_indels.tiers.tsv file',                                         'Type': 'String',  'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_INTOGEN_DRIVER_MUT',   'Description': 'Driver mutation by Introgen, reported by PCGR in .snvs_indels.tiers.tsv file',                              'Type': 'String',  'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_TCGA_PANCANCER_COUNT', 'Description': 'Occurences in TCGR, reported by PCGR in .snvs_indels.tiers.tsv file',                                       'Type': 'Integer', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_CLINVAR_CLNSIG',       'Description': 'ClinVar clinical significance, reported by PCGR in .snvs_indels.tiers.tsv file',                            'Type': 'String',  'Number': '1'})
            vcf.add_info_to_header({'ID': 'COSMIC_CNT',                'Description': 'Hits in COSMIC, reported by PCGR in .snvs_indels.tiers.tsv file',                                           'Type': 'Integer', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'ICGC_PCAWG_HITS',           'Description': 'Hits in ICGC PanCancer Analysis of Whole Genomes (PCAWG), reported by PCGR in .snvs_indels.tiers.tsv file', 'Type': 'Integer', 'Number': '1'})
        def func(rec):
            change = f'{rec.CHROM}:g.{rec.POS}{rec.REF}>{rec.ALT[0]}'
            pcgr_d = pcgr_fields_by_snp.get(change, {})
            if pcgr_d:
                for k, v in pcgr_d.items():
                    rec.INFO[f'PCGR_{k}'] = v
            rec.INFO['COSMIC_CNT'] = cosmic_by_snp.get(change, 0)
            rec.INFO['ICGC_PCAWG_HITS'] = icgc_by_snp.get(change, 0)
            return rec
        iter_vcf(input.vcf, output.vcf, func, func_hdr)

rule annotate:
    input:
        vcf = rules.somatic_vcf_pcgr_anno.output.vcf,
        tbi = rules.somatic_vcf_pcgr_anno.output.tbi,
        subset_to_cancer = f'somatic_anno/subset_to_cancer_genes.flag',
    output:
        vcf = OUTPUT_VCF,
        tbi = OUTPUT_VCF + '.tbi',
    shell:
        'cp {input.vcf} {output.vcf} && cp {input.tbi} {output.tbi}'




