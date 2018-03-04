from os.path import join

samples = [
	'batch1-ensemble-annotated',
	'batch1-ensemble-annotated-bwa',
	'batch1-mutect2-annotated',
	'batch1-mutect2-annotated-bwa',
	'batch1-strelka2-annotated',
	'batch1-strelka2-annotated-bwa',
	'batch1-vardict-annotated',
	'batch1-vardict-annotated-bwa',
]

tricky_bed = '/home/vlad/bcbio/genomes/Hsapiens/GRCh37/coverage/problem_regions/GA4GH/merged.sorted.bed.gz'


rule all:
	input:
		expand('{sample}_bcftools_isec/000{k}.tricky.pcgr.vcf', sample=samples, k=[0,1,2])


TOML = f'''[[annotation]]
file="{tricky_bed}"
names=["TRICKY"]
columns=[4]
ops=["self"]'''


rule anno_tricky:
	input:
		vcf = '{sample}_bcftools_isec/000{k}.vcf'
	output:
		'{sample}_bcftools_isec/000{k}.tricky.vcf'
	params:
		toml = TOML.replace('\n', r'\n')
	shell:
		"vcfanno <(echo '{params.toml}') {input} > {output}"


rule prep_pcgr:
	input:
		vcf = '{sample}_bcftools_isec/000{k}.tricky.vcf'
	output:
		'{sample}_bcftools_isec/000{k}.tricky.pcgr.vcf'
	params:
		toml = TOML.replace('\n', r'\n')
	shell:
		"pcgr_prep {input} > {output}"
