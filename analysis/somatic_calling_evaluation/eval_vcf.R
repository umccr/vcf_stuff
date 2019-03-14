source("evaluation.R")



##############
## Parsing
dir = "/Users/vsaveliev/Analysis/snv_validation/mb/ICGC_MB/"
truth_file = "MB-benchmark.ANNO.FILT.vcf.gz"
called_file = "batch1-ensemble-annotated.ANNO.FILT.TP.SAMPLE.vcf.gz"
tumor_sample = "tumor_downsample"

dir = "/Users/vsaveliev/spa/extras/vlad/synced/umccr/vcf_stuff/vcf_stuff/panel_of_normals/tests/"
# dir = "/Users/vsaveliev/Analysis/snv_validation/mb/new_pon/"
truth_file = "MB-benchmark.ANNO.FILT.PON_NEW.vcf.gz"
called_file = "batch1-ensemble-annotated.ANNO.FILT.TP.SAMPLE.PON_NEW.vcf.gz"

dir = "/Users/vsaveliev/rjn/projects/Saveliev_SNV_Filtering/cancer-giab-na12878-na24385/final/20161212_giab-mixture/"
truth_file = "na12878-na24385-somatic-truth-NORM-ANNO-FILT.vcf.gz"
called_file = "24385-ensemble-annotated.ANNO.FILT.vcf.gz"
tumor_sample = "X24385.12878.30.200"

truth_vcf  = read.vcf(str_c(dir, truth_file ), split.info = T, split.samples = T)
called_vcf = read.vcf(str_c(dir, called_file), split.info = T)

truth_vcf$vcf <- truth_vcf$vcf %>% rename(TRUTH_DP = DP)

merged = merge_called_and_truth(truth_vcf, called_vcf, tumor_sample)



##############
## Exploring
# Disregarding known gnomAD to avoid cluttering FN:
som = merged %>% filter(is.na(FILT) | !str_detect(FILT, "gnomAD_common"))
# All called variants for comparison:
resc_all = som %>% mutate(is_passed = is_called)
# Dropping the old panel of normals:
resc_pon = som %>% rescue_if(FILT == "PoN")
# Dropping AF10 filter as well because there are a lot true calls below 10%:
resc_pon_af = resc_pon %>% rescue_if(FILT == "AF10") 
show_stats(merged, som, resc_all, resc_pon, resc_pon_af)

plot_freq(resc_pon_af, "NM")
plot_freq(resc_pon_af, "MQ")
plot_freq(resc_pon_af, "AF")
plot_freq(resc_pon_af, "VD")
# 1. False positive SNPs are predominantly on low frequencies. 
#    False negative SNPs are on every AF, though VD tend to be low.
# 2. False positive indels are at any freq but mostly only below 0.35.
#    False negative inels are on 0.4-0.5, with VD at 40x - 
#    same shape as AF, so does not depend on DP much



###################
## Exploring filters
vd8_fi = resc_pon_af %>% reject_if(VD < 8)
vd8_ki = resc_pon_af %>% reject_if(VD < 8, rescue = T)
vd8_rc = resc_pon_af %>% rescue_if(VD >= 8)
show_stats(merged, resc_pon_af, vd8_fi, vd8_ki, vd8_rc)
# Apparently just keeping all VD>=8 improves the F2 measure. However, we again have trust issues with the truth set.
# Now checking what kind of calls are rescued with rescue_if(VD >= 8). 
# Firstly, those are gnomAD and PoN calls, so we rescued them already and removed gnomAD from the truth set
# 987 tp snps if we reject VD<8, 1125 tp snps if we rescue all VD>=8. What kind of snps VD>=8 but rejected?
resc_pon_af %>% filter(is_snp & VD >= 8 & !is_passed) %>% select(-ends_with(".truth"), -ends_with(".called")) %>% View()
# 3 snps only? weird
resc_pon_af %>% filter(is_snp & VD >= 8 & !is_passed) %>% count()  # 3 snps with high VD not passed, which should be rescued
resc_pon_af %>% reject_if(VD < 8) %>% filter(is_snp & VD >= 8 & !is_passed) %>% count()
resc_pon_af %>% filter(is_passed) %>% count(VD < 8)
resc_pon_af %>% filter(is_passed == F) %>% count(VD > 8)
resc_pon_af %>% rescue_if(VD >= 8) %>% count()
resc_pon_af %>% reject_if(VD < 8) %>% count()
# 105 snps are removed after we rescue VD>8:
resc_pon_af %>% rescue_if(VD >= 8) %>% filter(is_snp) %>% filter(!is_passed) %>% select(-ends_with(".truth"), -ends_with(".called")) %>% count()  # 105
resc_pon_af %>% filter(VD < 8)     %>% filter(is_snp) %>% filter(!is_passed) %>% select(-ends_with(".truth"), -ends_with(".called")) %>% count()  # 105
# why rescuing increases TP/FP count so drastically? because previously filtered variants become passed, regardless of FILT. 
#  can we reproduce this by rescuing certain FILT fields, and then applying reject_if()? have to remove all LowVD fields too... though they are already VD<8
show_stats(passed, resc_pon_af, vd8_fi, vd8_rc, 
           resc_pon_af %>% rescue_if(FILT != "PASS"), 
           resc_pon_af %>% rescue_if(FILT != "PASS") %>% reject_if(VD < 8))
# no, not working
# what can be the difference between reject_if(VD < 8) and rescue_if(VD >= 8)? let's see. reject drops is_passed & VD<8. rescue adds !is_passed & VD>=8.
#  meaning that we can take rescue_if(VD >= 8), and filter by VD<8?
resc_pon_af %>% filter(is_snp) %>% rescue_if(VD >= 8) %>% filter(VD < 8) %>% show_stats()  # getting exactly the difference 1125-987 = 138! nice 
# what exactly are we rescuing with VD>=8?
resc_pon_af %>% filter(is_snp) %>% rescue_if(VD >= 8) %>% filter(VD < 8) %>% count_status() %>% filter(is_tp) %>% 
  select(VD, AF, DP, FILT, -ends_with(".truth"), -ends_with(".called"), everything()) %>% View()

# again, what kind of calls are rescued with rescue_if(VD >= 8)?
resc_pon_af %>% filter(!is_passed & VD >= 8) %>% 
  select(VD, AF, DP, FILT, -ends_with(".truth"), -ends_with(".called"), everything()) %>% View()
# mostly StrandBias and HP, which is fine

# bcbio filter:
bcbio_filt  = resc_pon_af %>% reject_if(VD < 6 & (MQ < 60.0 & NM > 2.0 | MQ < 55.0 & NM > 1.0))
bcbio_filt8 = resc_pon_af %>% reject_if(VD < 8 & (MQ < 60.0 & NM > 2.0 | MQ < 55.0 & NM > 1.0))
show_stats(resc_pon_af, bcbio_filt, bcbio_filt8, vd8_rc)

# 11 FN SNPs left. Exploring them:
resc_pon_af %>% filter(is_called & !is_passed & str_detect(CALLS, 'GOLD')) %>% select(VD, AF, DP, FILT, -ends_with(".truth"), -ends_with(".called"), everything()) %>% View()
# All with VD<=4 +LCR, can't do anything with that

# exploring filter DP < avg_depth + 4 sqrt(avg_depth) from  https://www.nature.com/articles/s41592-018-0054-7.epdf?author_access_token=Dx2QKivDEkvswwWUJy5h0dRgN0jAjWel9jnR3ZoTv0Piny0X_8RjtB9zl5gHqNc_0m_Th8-X_roiQOsDJ4nU6Efwua2hRwes-IteFGhqVMhIXU1t_C1CZWZG5k8XMH65XRZZOcoziGXFh0QUciDDIg%3D%3D
resc_all %>% filter(DP > 140) %>% select(VD, AF, DP, FILT, -ends_with(".truth"), -ends_with(".called"), everything()) %>% View()
# there are variants with DP ~300, but that might be just a CN gain, so no use

# compare gnomAD/PoN and low normal depth
merged %>% filter(is_called & str_detect(FILT, 'gnom')) %>% select(FILT, VD, AF, DP, MQ, NORMAL_VD, NORMAL_AF, NORMAL_DP, NORMAL_MQ, CHROM, POS, REF, ALT, CALLS)
# gnomAD and PoN appear on any normal depth, so not sure why normal didn't have these variants. 
#  TODO: report germline leakage in MultiQC or Rmd: gnomAD and germline PoN (old-style PoN) with high enough normal depth



#################
## Exploring FN
passed %>% count(HS)

passed %>% filter(HS) %>% select(CHROM, POS, REF, ALT, AF, VD, FILT, GIAB, PoN, HMF_HS, TRICKY,
                                 PCGR_TIER, CLNSIG, DRIVER, ICGC, PCGR_HS, TCGA, GENE)
passed2 %>% count(ICGC >= 3)
passed2 %>% count(TCGA >= 5)

(fn = passed %>% filter(is_fn) %>% select(HS, everything()))

# TODO: explore CALLS and TIERS.truth

# TODO: SAGE.
# - extend the set of hotspots by adding PCGR sources?
# - CACAO: compare hotspots and genes with PCGR and HMF hotspots

# TODO: new truth set? check brad's filters again?



################
## Exploring new PoN

df = resc_pon_af %>% mutate(
  PoN_NEW = nonna(PoN_NEW.called, PoN_NEW.truth),
  PoN_OLD = nonna(PoN.called, PoN.truth)
) %>% count_status()

new1 = df %>% reject_if(PoN_NEW >= 1)
new2 = df %>% reject_if(PoN_NEW >= 2)
new3 = df %>% reject_if(PoN_NEW >= 3)
new4 = df %>% reject_if(PoN_NEW >= 4)
new5 = df %>% reject_if(PoN_NEW >= 5)
new6 = df %>% reject_if(PoN_NEW >= 6)

old1 = df %>% reject_if(PoN_OLD >= 1)
old2 = df %>% reject_if(PoN_OLD >= 2)
old3 = df %>% reject_if(PoN_OLD >= 3)
old4 = df %>% reject_if(PoN_OLD >= 4)
old5 = df %>% reject_if(PoN_OLD >= 5)
old6 = df %>% reject_if(PoN_OLD >= 6)

show_stats(
  merged
  , no_gno
  , drop_pon
  , old1
  , old2
  , old3
  , old4
  , old5
  , old6
  , new1
  , new2
  , new3
  , new4
  , new5
  , new6
)
# Best result in SNPs: PoN_NEW >= 2 or 3: F2=90.3%
# Best result in indels: PoN_OLD >= 3 
# Explore on indels.
# Trying without permissive overlap:
#      Indel stats even wrose due to a higher number of FP.
#      SNP stats changed slightly to become worse as well, due to overlap with indels in the PoN - e.g. variant like 
#      indels:       1       35534250        .       AAT     A,AATAT,ATAT    0       .       PoN_samples=7
#      snps:         1       35534250        .       A       T       0       .       PoN_samples=1
#      Enabling permissive overlap for SNPs resulted +1 missed call with NEW>=1, and same result with >=2 and the OLD one.
# TODO: revisit idea of filtering surrounding indels.
# TODO: check normal DP in filtered sites
# TODO: check with somatic mixture







##################
### Brad's vardict filter http://bcb.io/2016/04/04/vardict-filtering/

# Exploring MQ and NM correlation
resc_pon_af %>% filter(!is_tn) %>% ggplot() + geom_density(aes(MQ, color = status)) + facet_wrap(~vartype, ncol = 1)
resc_pon_af %>% filter(!is_tn) %>% ggplot() + geom_density(aes(NM, color = status)) + facet_wrap(~vartype, ncol = 1)
plot_freq(resc_pon_af %>% mutate(NM_MQ = NM * MQ), NM_MQ)
resc_pon_af %>% filter(!is_tn) %>% ggplot() + geom_density(aes(NM * MQ, color = status)) + facet_wrap(~vartype, ncol = 1)
# filtering all NM*MQ < 40 is a good idea for SNPs
show_stats(resc_pon_af, 
           resc_pon_af %>% reject_if(is_snp & NM * MQ < 40), 
           resc_pon_af %>% reject_if(is_snp & NM * MQ < 40, rescue=T),
           resc_pon_af %>% rescue_if(is_snp & NM * MQ < 40))
# reject_if is the best one.
nm_mq = resc_pon_af %>% reject_if(is_snp & NM * MQ < 40)

resc_pon_af %>% filter(!is_tn) %>% ggplot() + geom_point(aes(MQ, NM)) + facet_grid(status~vartype)
# We see that indeed like in Brad's post, FPs are spread all over whereas true variants cluster in the corner

# Back to Brad's filter:
brads = resc_pon_af %>% reject_if(AF * DP < 6 &&
                                 (MQ < 55.0 && NM > 1.0 ||
                                  MQ < 60.0 && NM > 2.0 ||
                                  DP < 10 ||
                                  VD_QUAL < 45))
brads2 = resc_pon_af %>% reject_if(VD < 6)
show_stats(resc_all, resc_pon_af, nm_mq, brads, brads2)
# rescue_if and reject_if(rescue = T) are worse than just reject_if
# brads and brads_rc is better than either brads2 or brads2_rc. So we must consider NM and MQ before filtering.

brads3 = resc_pon_af %>% reject_if(AF * DP < 6 && NM * MQ < 40)
show_stats(resc_all, resc_pon_af, nm_mq, brads, brads2, brads3)
# Conclusion:
# - Brad's filter is the best. Though it's already run in bcbio.
# - NM * MQ < 40 filter is even better if applied to SNPs only.
# - TODO: explore coupled with new PoN.

show_stats(resc_pon_af, lcr_vd5, lcr_vd6)









