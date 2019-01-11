from os.path import isfile, join, basename
import cyvcf2
import toml
import csv
from ngs_utils.file_utils import which
from ngs_utils.file_utils import get_ungz_gz
from ngs_utils.reference_data import get_key_genes_set
from ngs_utils.file_utils import splitext_plus
from hpc_utils.hpc import get_ref_file, get_loc
from vcf_stuff import iter_vcf


localrules: prep_anno_toml, prep_giab_bed, annotate


#############################
#### Reading parameters #####

SAMPLE = config['sample']
GENOME = config['genome']
INPUT_VCF = config['input_vcf']
OUTPUT_VCF = config['output_vcf']
assert OUTPUT_VCF.endswith('.vcf.gz'), OUTPUT_VCF
assert INPUT_VCF.endswith('.vcf.gz'), INPUT_VCF

PCGR_ENV_PATH = config.get('pcgr_env_path')
conda_cmd = ''
if PCGR_ENV_PATH:
    conda_cmd = 'export PATH=' + PCGR_ENV_PATH + '/bin:$PATH; '


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


# Preparations: annotate TUMOR_X and NORMAL_X fields for PCGR, remove non-standard chromosomes and mitochondria, remove non-PASSed calls
rule somatic_vcf_prep:
    input:
        vcf = INPUT_VCF
    output:
        vcf = f'somatic_anno/{SAMPLE}-somatic-ensemble-prep.vcf.gz',
        tbi = f'somatic_anno/{SAMPLE}-somatic-ensemble-prep.vcf.gz.tbi'
    shell:
        'pcgr_prep {input.vcf} | bcftools view -f.,PASS -Oz -o {output.vcf} && tabix -f -p vcf {output.vcf}'

rule somatic_vcf_pon_anno:
    input:
        vcf = rules.somatic_vcf_prep.output.vcf,
        tbi = rules.somatic_vcf_prep.output.tbi,
    params:
        genome_build = GENOME,
        pon_hits = 3,
        pon_dir = get_ref_file(GENOME, 'panel_of_normals_dir')
    output:
        vcf = f'somatic_anno/pon/{SAMPLE}-somatic-ensemble.vcf.gz',
        tbi = f'somatic_anno/pon/{SAMPLE}-somatic-ensemble.vcf.gz.tbi',
    shell:
        'pon_anno {input.vcf} --pon-dir {params.pon_dir} | bgzip -c > {output.vcf} && tabix -f -p vcf {output.vcf}'
        # ' | bcftools filter -e "INFO/PoN_CNT>={params.pon_hits}" --soft-filter PoN --mode + -Oz -o {output.vcf}' \
        # ' && tabix -f -p vcf {output.vcf} '

rule somatic_vcf_keygenes:
    input:
        vcf = rules.somatic_vcf_pon_anno.output.vcf,
    output:
        vcf = f'somatic_anno/keygenes/{SAMPLE}-somatic-ensemble.vcf.gz',
    run:
        genes = get_key_genes_set()
        def func(rec):
            if rec.INFO.get('ANN') is not None and rec.INFO['ANN'].split('|')[3] in genes:
                return rec
        iter_vcf(input.vcf, output.vcf, func)

rule somatic_vcf_pcgr_ready:
    input:
        full_vcf = rules.somatic_vcf_pon_anno.output.vcf,
        keygenes_vcf = rules.somatic_vcf_keygenes.output.vcf,
    output:
        vcf = f'somatic_anno/pcgr_input/{SAMPLE}-somatic-ensemble.vcf.gz',
        tbi = f'somatic_anno/pcgr_input/{SAMPLE}-somatic-ensemble.vcf.gz.tbi',
    run:
        total_vars = int(subprocess.check_output(f'bcftools view -H {input.full_vcf} | wc -l', shell=True).strip())
        vcf = input.full_vcf if total_vars <= 500_000 else input.keygenes_vcf  # to avoid PCGR choking on too many variants
#        def func(rec):
#            if rec.INFO.get('ANN') is not None:
#                rec['ANN'] = None
#            return rec
#        iter_vcf(vcf, output.vcf, func)
        shell(f'bcftools annotate -x INFO/ANN {vcf} -Oz -o {output.vcf} && tabix -f -p vcf {output.vcf}')

rule somatic_vcf_pcgr_round1:
    input:
        vcf = rules.somatic_vcf_pcgr_ready.output.vcf,
        pcgr_data = get_loc().pcgr_data,
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
        vcf = rules.somatic_vcf_pcgr_ready.output.vcf,
        tiers = rules.somatic_vcf_pcgr_round1.output.tiers,
    output:
        vcf = f'somatic_anno/pcgr_ann/{SAMPLE}-somatic-ensemble.vcf.gz',
        tbi = f'somatic_anno/pcgr_ann/{SAMPLE}-somatic-ensemble.vcf.gz.tbi',
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
            vcf.add_info_to_header({'ID': 'PCGR_TIER', 'Description':
                'TIER as reported by PCGR in .snvs_indels.tiers.tsv file. '
                '1: strong clinical significance; '
                '2: potential clinical significance; '
                '3: unknown clinical significance; '
                '4: other coding variants',
                'Type': 'String', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_CONSEQUENCE', 'Description':
                'VEP consequence, reported by PCGR in .snvs_indels.tiers.tsv file',
                                    'Type': 'String', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_MUTATION_HOTSPOT', 'Description':
                'Mutation hotspot, reported by PCGR in .snvs_indels.tiers.tsv file',
                                    'Type': 'String', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_INTOGEN_DRIVER_MUT', 'Description':
                'Driver mutation by Introgen, reported by PCGR in .snvs_indels.tiers.tsv file',
                                    'Type': 'String', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_TCGA_PANCANCER_COUNT', 'Description':
                'Occurences in TCGR, reported by PCGR in .snvs_indels.tiers.tsv file',
                                    'Type': 'Integer', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_CLINVAR_CLNSIG', 'Description':
                'ClinVar clinical significance, reported by PCGR in .snvs_indels.tiers.tsv file',
                                    'Type': 'String', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'COSMIC_CNT', 'Description':
                'Hits in COSMIC, reported by PCGR in .snvs_indels.tiers.tsv file',
                                    'Type': 'Integer', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'ICGC_PCAWG_HITS', 'Description':
                'Hits in ICGC PanCancer Analysis of Whole Genomes (PCAWG), reported by PCGR in .snvs_indels.tiers.tsv file',
                                    'Type': 'Integer', 'Number': '1'})
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

rule prep_giab_bed:
    input:
        get_ref_file(GENOME, ['truth_sets', 'giab', 'bed'])
    output:
        bed = f'somatic_anno/giab_conf.bed.gz',
        tbi = f'somatic_anno/giab_conf.bed.gz.tbi'
    shell:
        'cat {input} | bgzip -c > {output.bed} && tabix -f -p bed {output.bed}'

rule prep_hmf_hotspots:
    input:
        get_ref_file(GENOME, key='hmf_hotspot'),
    output:
        vcf = f'somatic_anno/hmf_hotspot.vcf.gz',
        tbi = f'somatic_anno/hmf_hotspot.vcf.gz.tbi',
    params:
        ungz = f'somatic_anno/hmf_hotspot.vcf'
    shell: """
echo "##fileformat=VCFv4.2" >> {params.ungz} &&
echo "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO" >> {params.ungz} &&
gunzip -c {input} | py -x "print('\t'.join([x.split()[0], x.split()[1], '.', x.split()[2], x.split()[3], '.', '.', 'HS']))" >> {params.ungz} &&
bgzip {params.ungz} &&
tabix -f -p vcf {output.vcf}
"""

rule prep_anno_toml:
    input:
        ga4gh_dir       = directory(join(get_ref_file(GENOME, key='problem_regions_dir'), 'GA4GH')),
        encode          = join(get_ref_file(GENOME, key='problem_regions_dir'), 'ENCODE', 'wgEncodeDacMapabilityConsensusExcludable.bed.gz'),
        lcr             = join(get_ref_file(GENOME, key='problem_regions_dir'), 'repeats', 'LCR.bed.gz'),
        giab_conf_bed   = rules.prep_giab_bed.output.bed,
        gnomad_vcf      = get_ref_file(GENOME, key='gnomad'),
        hmf_hotspots    = rules.prep_hmf_hotspots.output.vcf,
        hmf_giab        = get_ref_file(GENOME, key='hmf_giab_conf'),
        hmf_mappability = get_ref_file(GENOME, key='hmf_mappability'),
    output:
        f'somatic_anno/tricky_vcfanno.toml'
    run:
        with open(output[0], 'w') as f:
            f.write(f"""
[[annotation]]
file = "{input.giab_conf_bed}"
names = ["GIAB_CONF"]
columns = [3]
ops = ["flag"]

[[annotation]]
file="{input.gnomad_vcf}"
fields = ["AF"]
names = ["gnomAD_AF"]
ops = ["self"]

[[annotation]]
file = "{input.hmf_hotspots}"
fields = ["HS"]
names = ["HMF_HOTSPOT"]
ops = ["flag"]

[[annotation]]
file = "{input.hmf_giab}"
names = ["HMF_GIAB_CONF"]
columns = [3]
ops = ["flag"]

[[annotation]]
file = "{input.hmf_mappability}"
names = ["HMF_MAPPABILITY"]
columns = [5]
ops = ["self"]

[[annotation]]
file = "{input.lcr}"
names = ["TRICKY_LCR"]
columns = [3]
ops = ["flag"]

[[annotation]]
file = "{input.encode}"
names = ["ENCODE"]
columns = [4]
ops = ["self"]

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
        vcf = rules.somatic_vcf_pcgr_anno.output.vcf,
        toml = rules.prep_anno_toml.output[0]
    output:
        vcf = f'somatic_anno/regions/{SAMPLE}-somatic-ensemble.vcf.gz',
        tbi = f'somatic_anno/regions/{SAMPLE}-somatic-ensemble.vcf.gz.tbi',
    shell:
        'vcfanno {input.toml} {input.vcf} | bgzip -c > {output.vcf} && tabix -f -p vcf {output.vcf}'


rule annotate:
    input:
        vcf = rules.somatic_vcf_regions_anno.output.vcf,
        tbi = rules.somatic_vcf_regions_anno.output.tbi,
    output:
        vcf = OUTPUT_VCF,
        tbi = OUTPUT_VCF + '.tbi',
    shell:
        'cp {input.vcf} {output.vcf} && cp {input.tbi} {output.tbi}'




