#!/usr/bin/env bash

set -x

VCF=MB_100vs50-vardict-annotated.vcf.gz
#VCF=test.vcf.gz

OUT=${VCF/.vcf.gz/.cyvcf2_writer.vcf}
time zsh -c "python test_vcf_parsers.py ${VCF} cyvcf2 ${OUT} && bgzip -f ${OUT}"
# 291.77s user 6.91s system 92% cpu 5:23.91 total
# 285.94s user 9.09s system 90% cpu 5:26.25 total

time zsh -c "python test_vcf_parsers.py ${VCF} cyvcf2 | bgzip -c > ${VCF/.vcf/.cyvcf2.vcf.gz}"
# 283.23s user 2.40s system 183% cpu 2:35.55 total
# 293.69s user 2.78s system 181% cpu 2:43.05 total

time zsh -c "python test_vcf_parsers.py ${VCF} pysam | bgzip -c > ${VCF/.vcf/.pysam.vcf.gz}"
# 325.09s user 2.67s system 171% cpu 3:11.47 total
# 331.54s user 2.66s system 173% cpu 3:12.14 total

time zsh -c "python test_vcf_parsers.py ${VCF} python | bgzip -c > ${VCF/.vcf/.python.vcf.gz}"
# 251.17s user 1.84s system 179% cpu 2:20.60 total
# 250.64s user 1.70s system 178% cpu 2:21.31 total

time zsh -c "gunzip -c ${VCF} | py -x 'test_vcf_parsers.proc_line(x)' | bgzip -c > ${VCF/.vcf/.py.vcf.gz}"
# 309.34s user 4.74s system 186% cpu 2:48.33 total

set +x
