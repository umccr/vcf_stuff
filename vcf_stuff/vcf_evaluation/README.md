Small variant validation
========================

## Workflow

`Snakefile` is a workflow file for validation of variant calling by comparing VCF files to truth set. Usage:

1. Prepare a `config.yaml` file similar to `config_eval_MB_grch37_spartan.yaml` or `config_eval_dream_grch37_local.yaml`. 
Not that `vcf_sample` should correspond to the main sample name in multi-sample VCF files.

2. Make sure all dependencies from `environment.yaml` are installed - e.g. create a new conda environment:
```
conda env create --name snakemake-validation --file environment.yaml
source activate snakemake-validation
```
3. Run snakemake in the directory where `Snakemake` file is located:
```
snakemake -p --configfile=config.yaml
``` 

## Tools

Summary of tools and datasets that can be used for variant calling comparison with truth set and visualization.

* [rtgeval](https://github.com/lh3/rtgeval) - wrapper for RTG's vcfeval.

* [Concordance between variant callers](https://github.com/vladsaveliev/concordance)

* [JS-based Venn diagrams for BED and VCF files](https://github.com/vladsaveliev/Venn)

* [bcftools isec](https://samtools.github.io/bcftools/bcftools-man.html#isec)

* Bcbio:
  * [Examples of validation](https://github.com/bcbio/bcbio_validations)
  * [Somatic variant callers validation for tumor-only samples blog post](http://bcb.io/2015/03/05/cancerval/)

## Somatic validation datasets

### GiaB NA12878/NA24385 somatic-like mixture
Source: https://github.com/hbc/projects/tree/master/giab_somatic
Paired
WGS
Depth:
Tumor purity:
Samples:
Location in the FS:

### ICGC-TCGA DREAM challange
Text: https://www.synapse.org/#!Synapse:syn312572
Datasets:
* Set 3
* Set 4
Paired
WGS
Depth:
Tumor purity:
Samples:
Location in the FS:

### COLO829
Metastatic melanoma cell line
Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4837349/, 2016
Depth: T: 99X, N: 103X

### ICGC MB
Paper: https://www.nature.com/articles/ncomms10001, 2015
Paired
WGS
Depth: T: 314x, N: 272x
Tumor purity: 95–98%
Samples:
Location in the FS:







