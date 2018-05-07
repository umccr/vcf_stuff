## Allele frequency

By default, bcbio uses the minimal AF threshold of 10%. Before, bcbio applied to threshold only to VarDict, ignoring Strelka2 and Mutect2. We changed bcbio so it applies the threshold to Strelka2 and Mutect2 calls as well.
