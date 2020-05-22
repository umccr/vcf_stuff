from cyvcf2.cyvcf2 import VCF
from os.path import isfile, join, basename, dirname
import cyvcf2
import toml
import csv
import yaml
from ngs_utils.file_utils import which, safe_mkdir
from ngs_utils.file_utils import get_ungz_gz
from ngs_utils.file_utils import splitext_plus
from ngs_utils.logger import warn
from ngs_utils.vcf_utils import iter_vcf, count_vars, vcf_contains_field
from ngs_utils.reference_data import get_key_genes_set
from reference_data import api as refdata


localrules: prep_anno_toml, annotate


#############################
#### Reading parameters #####

SAMPLE = config['sample']
GENOME = config['genome']
INPUT_VCF = config['input_vcf']
OUTPUT_VCF = config['output_vcf']
T_NAME = config.get('tumor_vcf_sample')
N_NAME = config.get('normal_vcf_sample')
assert OUTPUT_VCF.endswith('.vcf.gz'), OUTPUT_VCF
assert INPUT_VCF.endswith('.vcf.gz'), INPUT_VCF

PCGR_ENV_PATH = config.get('pcgr_env_path')
conda_cmd = ''
if PCGR_ENV_PATH:
    conda_cmd = 'export PATH=' + PCGR_ENV_PATH + '/bin:$PATH; '

if config.get('input_genomes_url'):
    refdata.find_genomes_dir(config.get('input_genomes_url'))


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
        vcf = refdata.get_ref_file(GENOME, key='hotspots'),
    output:
        vcf = f'somatic_anno/hmf_hotspot.vcf.gz',
        tbi = f'somatic_anno/hmf_hotspot.vcf.gz.tbi',
    shell:
        'bcftools filter -i "HMF=1" {input.vcf} -Oz -o {output.vcf} && tabix -p vcf {output.vcf}'

rule prep_anno_toml:
    input:
        ga4gh_dir       = join(refdata.get_ref_file(GENOME, key='problem_regions_dir'), 'GA4GH'),
        lcr             = join(refdata.get_ref_file(GENOME, key='problem_regions_dir'), 'repeats', 'LCR.bed.gz'),
        segdup          = join(refdata.get_ref_file(GENOME, key='problem_regions_dir'), 'segdup.bed.gz'),
        gnomad_vcf      = refdata.get_ref_file(GENOME, key='gnomad'),
        hmf_hotspots    = rules.prep_hmf_hotspots.output.vcf,
        hmf_giab        = refdata.get_ref_file(GENOME, key='hmf_giab_conf'),
        encode          = join(refdata.get_ref_file(GENOME, key='problem_regions_dir'), 'ENCODE',
                              {'hg38': 'encode4_unified_blacklist.bed.gz', 'GRCh37': 'blacklist.v2.bed.gz'}[GENOME])
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
file = "{input.lcr}"
names = ["TRICKY_LCR"]
columns = [3]
ops = ["flag"]

[[annotation]]
file = "{input.encode}"
names = ["ENCODE"]
columns = [3]
ops = ["flag"]

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
# If the noise wasn't germline, it might be artefacts/errors from FFPE or ortherwise low quality data.
# subsetting to cancer genes in this case.
rule maybe_subset_highly_mutated:
    input:
        vcf = rules.somatic_vcf_regions_anno.output.vcf,
        tbi = rules.somatic_vcf_regions_anno.output.tbi,
    output:
        vcf = f'somatic_anno/subset/{SAMPLE}-somatic.vcf.gz',
        tbi = f'somatic_anno/subset/{SAMPLE}-somatic.vcf.gz.tbi',
        subset_highly_mutated_stats = f'somatic_anno/subset_highly_mutated_stats.yaml',
    params:
        no_gnomad_vcf = f'somatic_anno/subset/no_gnomad/{SAMPLE}.vcf.gz',
    run:
        total_vars = count_vars(input.vcf)
        vars_no_gnomad = None
        vars_cancer_genes = None
        if total_vars > 500_000:
            warn(f'Found {total_vars}>500k somatic variants, start with removing gnomAD_AF>0.01')
            def func(rec, vcf):
                gnomad_af = rec.INFO.get('gnomAD_AF')
                if gnomad_af is not None and float(gnomad_af) >= 0.01 \
                        and not rec.INFO.get('HMF_HOTSPOT')\
                        and not rec.INFO.get('SAGE_HOTSPOT'):
                    return None
                else:
                    return rec
            safe_mkdir(dirname(params.no_gnomad_vcf))
            iter_vcf(input.vcf, params.no_gnomad_vcf, func)
            vars_no_gnomad = count_vars(params.no_gnomad_vcf)
            if vars_no_gnomad > 500_000:
                warn(f'After removing gnomAD_AF>0.01, still having {vars_no_gnomad}>500k somatic variants left. '
                     f'So _instead_ subsetting to cancer genes.')
                genes = get_key_genes_set()
                def func(rec, vcf):
                    if rec.INFO.get('ANN') is not None and rec.INFO['ANN'].split('|')[3] in genes:
                        return rec
                iter_vcf(input.vcf, output.vcf, func)
                vars_no_gnomad = None
                vars_cancer_genes = count_vars(output.vcf)
                warn(f'After subsetting to cancer genes, left with {vars_cancer_genes} variants')
            else:
                warn(f'Removing gnomAD>0.01 reduced down to {vars_no_gnomad}<=500k somatic variants, '
                     f'no need to subset to cancer genes keeping all of them for futher reporting')
                shell('cp {params.no_gnomad_vcf} {output.vcf} ; cp {params.no_gnomad_vcf}.tbi {output.tbi}')
        else:
            shell('cp {input.vcf} {output.vcf} ; cp {input.tbi} {output.tbi} ; ')

        with open(output.subset_highly_mutated_stats, 'w') as out:
            stats = dict(
                total_vars=total_vars
            )
            if vars_no_gnomad is not None:
                stats['vars_no_gnomad'] = vars_no_gnomad
            if vars_cancer_genes is not None:
                stats['vars_cancer_genes'] = vars_cancer_genes
            yaml.dump(stats, out, default_flow_style=False)

# Removes TRICKY_ and ANN fields
rule somatic_vcf_clean_info:
    input:
        vcf = rules.maybe_subset_highly_mutated.output.vcf,
        tbi = rules.maybe_subset_highly_mutated.output.tbi,
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

        def func(rec, vcf):
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
    params:
        t_name = T_NAME,
        n_name = N_NAME,
    run:
        t_name_arg = f"-tn {params.t_name}" if params.t_name else ""
        n_name_arg = f"-nn {params.n_name}" if params.n_name else ""
        shell(f'pcgr_prep {t_name_arg} {n_name_arg} {input.vcf}' 
              f' | bgzip -c > {output.vcf} && tabix -f -p vcf {output.vcf}')

rule somatic_vcf_pon_anno:
    input:
        vcf = rules.somatic_vcf_prep.output.vcf,
        tbi = rules.somatic_vcf_prep.output.tbi,
    params:
        genome_build = GENOME,
        pon_hits = 3,
        pon_dir = refdata.get_ref_file(GENOME, 'panel_of_normals_dir')
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
        pcgr_data = refdata.get_ref_file(GENOME, key='pcgr_data'),
    output:
        tiers = f'somatic_anno/pcgr_run/{SAMPLE}-somatic.pcgr.snvs_indels.tiers.tsv',
        vcf = f'somatic_anno/pcgr_run/{SAMPLE}-somatic.pcgr_ready.vep.vcf.gz',
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
        pcgr_tiers = rules.somatic_vcf_pcgr_round1.output.tiers,
        pcgr_vcf = rules.somatic_vcf_pcgr_round1.output.vcf,
    output:
        vcf = f'somatic_anno/pcgr_ann/{SAMPLE}-somatic.vcf.gz',
        tbi = f'somatic_anno/pcgr_ann/{SAMPLE}-somatic.vcf.gz.tbi',
    params:
        genome_build = GENOME
    run:
        # Reading information from PCGR tiers file
        pcgr_fields_by_snp = dict()
        cosmic_by_snp = dict()
        icgc_by_snp = dict()
        with open(input.pcgr_tiers) as f:
            reader = csv.DictReader(f, delimiter='\t', fieldnames=f.readline().strip().split('\t'))
            for row in reader:
                change = row['GENOMIC_CHANGE']
                # Fixing the chromosome name 1 -> chr1 as PCGR strips the chr prefixes
                if params.genome_build == 'hg38':
                    change = 'chr' + change.replace('MT', 'M')
                pcgr_fields_by_snp[change] = dict()
                for k in ['SYMBOL',
                          'TIER',
                          'CONSEQUENCE',
                          'MUTATION_HOTSPOT',
                          'PUTATIVE_DRIVER_MUTATION',
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

        # Reading CSQ from PCGR VCF file as the tiers file has it incomplete
        csq_by_snp = dict()
        vcf = VCF(input.pcgr_vcf)
        for rec in vcf:
            change = f'{rec.CHROM}:g.{rec.POS}{rec.REF}>{rec.ALT[0]}'
            # Fixing the chromosome name 1 -> chr1 as PCGR strips the chr prefixes
            if params.genome_build == 'hg38':
                change = 'chr' + change.replace('MT', 'M')
            try:
                csq = rec.INFO['CSQ']
            except:
                pass
            else:
                csq_by_snp[change] = csq

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
            vcf.add_info_to_header({'ID': 'PCGR_PUTATIVE_DRIVER_MUTATION',   'Description': 'Putative driver mutation, reported by PCGR in .snvs_indels.tiers.tsv file',                           'Type': 'String',  'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_TCGA_PANCANCER_COUNT', 'Description': 'Occurences in TCGR, reported by PCGR in .snvs_indels.tiers.tsv file',                                       'Type': 'Integer', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_CLINVAR_CLNSIG',       'Description': 'ClinVar clinical significance, reported by PCGR in .snvs_indels.tiers.tsv file',                            'Type': 'String',  'Number': '1'})
            vcf.add_info_to_header({'ID': 'COSMIC_CNT',                'Description': 'Hits in COSMIC, reported by PCGR in .snvs_indels.tiers.tsv file',                                           'Type': 'Integer', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'ICGC_PCAWG_HITS',           'Description': 'Hits in ICGC PanCancer Analysis of Whole Genomes (PCAWG), reported by PCGR in .snvs_indels.tiers.tsv file', 'Type': 'Integer', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'CSQ', 'Description': 'Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|RefSeq|DOMAINS|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|CHECK_REF',
                                    'Type': 'String', 'Number': '.'})
        def func(rec, vcf):
            change = f'{rec.CHROM}:g.{rec.POS}{rec.REF}>{rec.ALT[0]}'
            pcgr_d = pcgr_fields_by_snp.get(change, {})
            if pcgr_d:
                for k, v in pcgr_d.items():
                    rec.INFO[f'PCGR_{k}'] = v
            rec.INFO['COSMIC_CNT'] = cosmic_by_snp.get(change, 0)
            rec.INFO['ICGC_PCAWG_HITS'] = icgc_by_snp.get(change, 0)
            rec.INFO['CSQ'] = csq_by_snp.get(change)
            return rec
        iter_vcf(input.vcf, output.vcf, func, func_hdr)

rule annotate:
    input:
        vcf = rules.somatic_vcf_pcgr_anno.output.vcf,
        tbi = rules.somatic_vcf_pcgr_anno.output.tbi,
        subset_highly_mutated_stats = f'somatic_anno/subset_highly_mutated_stats.yaml',
    output:
        vcf = OUTPUT_VCF,
        tbi = OUTPUT_VCF + '.tbi',
    shell:
        'cp {input.vcf} {output.vcf} && cp {input.tbi} {output.tbi}'




