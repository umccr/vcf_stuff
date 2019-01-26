#!/usr/bin/env bash

set -x

VCF=MB_100vs50-vardict-annotated.vcf.gz
#VCF=test.vcf.gz

OUT=${VCF/.vcf.gz/.cyvcf2_writer.vcf}
time zsh -c "python test_vcf_parsers.py ${VCF} cyvcf2 ${OUT}"
#  130.51s user 2.27s system 99% cpu 2:13.93 total

time zsh -c "python test_vcf_parsers.py ${VCF} cyvcf2 > ${VCF/.vcf/.cyvcf2.vcf}"
#  143.69s user 2.62s system 99% cpu 2:27.67 total

time zsh -c "python test_vcf_parsers.py ${VCF} pysam > ${VCF/.vcf/.pysam.vcf}"
#  180.67s user 2.33s system 99% cpu 3:03.97 total

time zsh -c "python test_vcf_parsers.py ${VCF} python > ${VCF/.vcf/.python.vcf}"
#  110.59s user 2.03s system 99% cpu 1:53.26 total

time zsh -c "gunzip -c ${VCF} | py -x 'test_vcf_parsers.proc_line(x)' > ${VCF/.vcf/.py.vcf}"

set +x

