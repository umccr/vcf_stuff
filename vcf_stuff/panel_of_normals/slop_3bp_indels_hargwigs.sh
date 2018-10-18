# Converts the indels panel of normals into a BED file, extending indels of at least 3 bases long by 10 base pairs. Idea from Hartwig's paper https://www.biorxiv.org/content/biorxiv/early/2018/09/20/415133.full.pdf
# > Removal of INDELS near a PON filtered INDEL - Regions of complex haplotype alterations are
#   often called as multiple long indels which can make it more difficult to construct an effective PON,
#   and sometimes we find residual artefacts at these locations. Hence we also filter inserts or
#   deletes which are 3 bases or longer where there is a PON filtered INDEL of 3 bases or longer
#   within 10 bases in the same sample.
#
# This approach will filter out any indel surrounding a 3-bp indel. Todo: any ideas how to apply it to only for 3bp long indels?

bcftools query panel_of_normals.indels.vcf.gz -f "%CHROM \t %POS \t %REF \t %ALT \t %PoN_samples \n" | bioawk -t '{ SLOP=length($3)-2 ; if(length($3)-2>=3) {SLOP=10+length($3)-2} else { if (length($4)-2>=3) {SLOP=10} } ; print $1, $2-1-SLOP, $2-1+SLOP, $5 }' > panel_of_normals.indels.indel3bp_slop10.bed
