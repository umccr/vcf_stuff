library(stringr)
library(stringi)
library(plyranges)
library(readr)
library(purrr)
library(crayon)
library(bedr)
# devtools::install_github("r-lib/pillar")
# library(pillar)
library(patchwork)
if("dplyr" %in% (.packages())){
  detach("package:ggplot2", unload=TRUE) 
  detach("package:dplyr", unload=TRUE) 
  detach("package:plyr", unload=TRUE) 
} 
library(plyr)
library(ggplot2)
library(dplyr)
library(tidyr)

dir = "/Users/vsaveliev/Analysis/snv_validation/mb/ICGC_MB/"
truth_file = "MB-benchmark.ANNO.FILT.vcf.gz"
called_file = "batch1-ensemble-annotated.ANNO.FILT.TP.SAMPLE.vcf.gz"

dir = "/Users/vsaveliev/spa/extras/vlad/synced/umccr/vcf_stuff/vcf_stuff/panel_of_normals/test_chr1/old/compare/"
dir = "/Users/vsaveliev/tmp/"
truth_file = "MB-benchmark.ANNO.FILT.1.PON_OLD_NEW.vcf.gz"
called_file = "batch1-ensemble-annotated.ANNO.FILT.TP.SAMPLE.1.PON_OLD_NEW.vcf.gz"

truth_vcf = read.vcf(str_c(dir, truth_file), split.info = T, split.samples = T)
called_vcf = read.vcf(str_c(dir, called_file), split.info = T)
called_vcf$vcf$NM = str_split(called_vcf$vcf$tumor_downsample, ":", simplify = T)[, 10] %>% as.double()

# proc_tricky = function(vcf) {
#   for (key in names(vcf$vcf)) {
#     if (str_detect(key, "^TRICKY_|^AZ_")) {
#       print(key)
#       vcf$vcf[[key]] <- NULL
#     }
#   }
#   vcf
# }
# truth <- proc_tricky(truth)
# called <- proc_tricky(called)
# names(called$vcf)

get_type <- function(ref, alt) {
  var_types = c("SNP", "Indel", "MNP")
  vt = ifelse(str_length(ref) == str_length(alt), ifelse(str_length(ref) == 1, "SNP", "MNP"), "Indel")
  factor(vt, var_types, ordered = T)
}

count_status <- function(.data) {
  .data %>% mutate(
    is_tp = is_passed & is_true,
    is_fp = is_passed & !is_true,
    is_fn = !is_passed & is_true,
    is_tn = !is_passed & !is_true,
    status = ifelse(is_tp, 'TP', ifelse(is_fp, 'FP', ifelse(is_fn, 'FN', 'TN')))
  )
}

fix_fields <- function(data) {
  renamed <- data %>% 
    mutate(
      HMF_MAPPABILITY = map_dbl(map(str_split(HMF_MAPPABILITY, ","), as.double), min)
    ) %>% 
    select(
      # -AF,
      # -VD,
      # -DP,
      -QUAL,
      -ID
    ) %>% 
    rename(
      GENE = PCGR_SYMBOL,
      TCGA = PCGR_TCGA_PANCANCER_COUNT,
      ICGC = ICGC_PCAWG_HITS,
      DRIVER = PCGR_INTOGEN_DRIVER_MUT,
      CLNSIG = PCGR_CLINVAR_CLNSIG,
      PCGR_HS = PCGR_MUTATION_HOTSPOT,
      HMF_HS = HMF_HOTSPOT,
      GIAB = HMF_GIAB_CONF,
      PoN = PoN_CNT,
      MBL = HMF_MAPPABILITY,
      COSM = COSMIC_CNT,
      CSQ = PCGR_CONSEQUENCE
    )
}

nonna <- function(a, b) {
  ifelse(!is.na(a), a, b)
}

merged <- full_join(
    truth_vcf$vcf %>% as_tibble() %>% fix_fields(), 
    called_vcf$vcf %>% as_tibble() %>% fix_fields() %>% select(-TIERS),
    by = c('CHROM', 'POS', 'REF', 'ALT'),
    suffix = c('.truth', '.called')) %>% 
  select(
    -AF,
    -VD,
    -DP
  ) %>% 
  rename(
    AF = TUMOR_AF.called,
    VD = TUMOR_VD.called,
    DP = TUMOR_DP.called,
    FILT = FILTER.called
    # 'COSM', 'ENCODE', 'GIAB', 'MBL', 'HMF_HS', 'ICGC', 'PCGR_TIER', 'CLNSIG', 
    # 'CSQ', 'DRIVER', 'PCGR_HS', 'TCGA', 'GENE', 'TRICKY', 'PoN'
  ) %>% 
  mutate(
    vartype = get_type(REF, ALT),
    is_snp = vartype == "SNP",
    is_called = !is.na(FILT),
    umccrise_passed = is_called & FILT == 'PASS',
    is_passed = umccrise_passed,
    is_true = !is.na(TIERS),
    AF = as.double(AF),
    DP = as.integer(DP),
    VD = round(AF * DP)
  ) %>% 
  mutate(
    GENE = nonna(GENE.called, GENE.truth),
    TCGA = nonna(TCGA.called, TCGA.truth),
    ICGC = nonna(ICGC.called, ICGC.truth),
    DRIVER = nonna(DRIVER.called, DRIVER.truth),
    CLNSIG = nonna(CLNSIG.called, CLNSIG.truth),
    PCGR_HS = nonna(PCGR_HS.called, PCGR_HS.truth),
    HMF_HS = nonna(HMF_HS.called, HMF_HS.truth),
    GIAB = nonna(GIAB.called, GIAB.truth),
    PoN = nonna(PoN.called, PoN.truth),
    PoN_NEW = nonna(PoN_NEW.called, PoN_NEW.truth),
    PoN_OLD = nonna(PoN_OLD.called, PoN_OLD.truth),
    MBL = nonna(MBL.called, MBL.truth),
    COSM = nonna(COSM.called, COSM.truth),
    CSQ = nonna(CSQ.called, CSQ.truth),
    PCGR_TIER = nonna(PCGR_TIER.called, PCGR_TIER.truth),
    TRICKY = nonna(TRICKY.called, TRICKY.truth),
    HS = HMF_HS | !is.na(DRIVER) | !is.na(PCGR_HS) | 
         (!is.na(CLNSIG) & str_detect(CLNSIG, "pathogenic|uncertain")) |
         TCGA >= 5 | ICGC >= 3 | COSM >= 5 | PCGR_TIER %in% c("TIER_1", "TIER_2"),
    HS = replace_na(HS, F)
  ) %>% 
  count_status()

no_germ = merged %>% 
  filter(!str_detect(FILT, "gnomAD_common")) %>% 
  filter(!str_detect(FILT, "PoN"))

merged %>% count(HS)
merged %>% count(MQ)
merged %>% count(NM)

f_measure <- function(b, prec, recall) {
  Fmeasure <- (1 + b**2) * prec * recall / (b**2 * prec + recall)
}

build_stats <- function(.data) {
  stats <- .data %>% 
    count_status() %>% 
    group_by(vartype) %>% 
    summarise(
      passed = sum(is_passed),
      true = sum(is_true),
      TP = sum(is_tp),
      FP = sum(is_fp),
      FN = sum(is_fn),
      recall = TP / true,
      prec = TP / passed,
      F2 = f_measure(2, prec, recall)
    ) %>% 
    select(-passed, -true) %>% 
    mutate_if(is.double, function(v) {str_c(round(v * 100, v), "%")})
  
  metrics = names(stats %>% select(-vartype))
  
  # If you want samples to be columns, and metrics to be rows:
  stats_row <- stats %>%
    gather(metric, s, -vartype) %>%
    mutate(metric = factor(metric, levels = metrics)) %>%
    arrange(vartype)
    # mutate(s.called = ifelse(metric %in% c("prec", "recall", "F2"), str_c(round(s.called * 100, 2), "%"), s.called))
  
  # If you want samples to be rows, and metrics to be columns:  
  # stats_col <- stats %>%
  #   gather(metric, value, -vartype) %>%
  #   mutate(metric = factor(metric, levels = metrics)) %>%
  #   arrange(vartype) %>%
  #   # unite(type_metric, vartype, metric, sep=".") %>%
  #   mutate(type_metric = interaction(metric, vartype), vartype = NULL, metric = NULL) %>%
  #   spread(type_metric, value)
}

# For samples to be rows, and metrics to be columns:  
# print_stats <- function(df) {
#   names(df) <- map_chr(names(df), str_replace, ".SNP|.Indel", "")
#   df
# }
# called_stats %>% mutate(
#     type = "called"
#   ) %>% 
#   rbind(
#     passed_stats %>% mutate(type = "passed")
#   ) %>% 
#   select(type, everything()) %>% 
#   print_stats()

show_stats = function(...) {
  stat_dfs = map(list(...), build_stats)
  res = join_all(stat_dfs, by = c("metric", "vartype"), match = "all")  # suffix = argnames)
  argnames = sapply(substitute(list(...))[-1], deparse)  # %>% str_c(".", .)
  names(res)[c(-1, -2)] = argnames
  res
}

##############
### Exploring
called <- no_germ %>% mutate(is_passed = is_called)
passed <- no_germ %>% mutate(is_passed = umccrise_passed)
show_stats(called, passed)

plot_freq <- function(.data, metric) {
  .data %>% 
    count_status() %>% 
    filter(!is_tn) %>% 
    filter(!is.na(!!metric)) %>% 
    # mutate(AF = cut(AF, seq(0, 1, 0.05))) %>% 
    ggplot() +
    # geom_histogram(aes(get(metric), colour = status), binwidth = 0.05)
    # geom_density(aes(get(metric), colour = status))
    geom_density(aes(.data[[metric]], colour = status)) +
    # geom_freqpoly(aes(get(metric), colour = status), bins = 30) +
    facet_wrap(~vartype, nrow = 2, scales = "free_y")
}
plot_freq(passed, "NM")
plot_freq(passed, "MQ")
plot_freq(passed, "AF")
plot_freq(passed, "VD")
# 1. False positive SNPs are predominantly on low frequencies. 
#    False negative SNPs are on every AF, though VD tend to be low.
# 2. False positive indels are at any freq but mostly only below 0.35.
#    False negative inels are on 0.4-0.5, with VD at 40x - 
#    same shape as AF, so does not depend on DP much

#############
## Definding filtering functions
reject_if = function(.data, cond, rescue = F) {
  # checks passing variants and rejectes if cond is TRUE
  # if rescue = TRUE, also checks rejected variants, and rescues them if cond is FALSE
  cond_q = substitute(cond)
  # keep when is_passed & cond!=T
  mask = .data$is_passed & (eval(cond_q, .data) %in% c(F, NA))
  if (rescue) {
    # also rescue _all_ when cond=F
    mask = mask | eval(cond_q, .data) %in% c(F)
  }
  .data %>% mutate(is_passed = is_called & mask)
}

rescue_if = function(.data, cond) {
  cond_q = substitute(cond)
  mask = .data$is_passed | eval(cond_q, .data) %in% c(T)
  .data %>% mutate(is_passed = is_called & mask)
}

df = tribble(
  ~is_called, ~umccrise_passed, ~is_true, ~NM, ~vartype,
  TRUE,       TRUE,             TRUE,     1,   "snp",
  TRUE,       F,                TRUE,     1,   "snp",                
  F,          F,                TRUE,     1,   "snp",
  TRUE,       TRUE,             F,        1,   "snp",
  TRUE,       F,                F,        1,   "snp",
  TRUE,       TRUE,             TRUE,     0,   "snp",
  TRUE,       F,                TRUE,     0,   "snp",                
  F,          F,                TRUE,     0,   "snp",
  TRUE,       TRUE,             F,        0,   "snp",
  TRUE,       F,                F,        0,   "snp",
  TRUE,       TRUE,             TRUE,     NA,  "snp",
  TRUE,       F,                TRUE,     NA,  "snp",                
  F,          F,                TRUE,     NA,  "snp",
  TRUE,       TRUE,             F,        NA,  "snp",
  TRUE,       F,                F,        NA,  "snp"
) %>% 
  mutate(is_passed = is_called & umccrise_passed)

show_stats(df, df %>% reject_if(NM < 1), df %>% reject_if(NM < 1, rescue = T))
df %>% reject_if(NM < 1, rescue = F)
df %>% reject_if(NM < 1, rescue = T)
df %>% rescue_if(NM == 1)

###########
## Exploring filters
#vd8_fi = passed %>% reject_if(VD < 8)
#vd8_ki = passed %>% keep_if(VD >= 8 | is.na(VD))
allaf = passed %>% rescue_if(FILT == "AF10") 
allaf_vd8_fi = allaf %>% reject_if(VD < 8)
allaf_vd8_ki = allaf %>% reject_if(VD < 8, rescue = T)
allaf_vd8_rc = allaf %>% rescue_if(VD >= 8)
show_stats(passed, allaf, allaf_vd8_fi, allaf_vd8_ki, allaf_vd8_rc)
# Apparently just keeping all VD>=8 improves the F2 measure. However, we again have trust issues with the truth set.
# Using allaf because there are a lot true calls below 10%.

# Now checking what kind of calls are rescued with rescue_if(VD >= 8). Firstly, those are gnomAD and PoN calls. Removing them completely from the truth set as well.
show_stats(passed, allaf, allaf %>% reject_if(VD < 8), allaf %>% reject_if(VD < 8, rescue = T), allaf %>% rescue_if(VD >= 8))
# 987 tp snps if we reject VD<8, 1125 tp snps if we rescue all VD>=8. What kind of snps VD>=8 but rejected?
allaf %>% filter(is_snp & VD >= 8 & !is_passed) %>% select(-ends_with(".truth"), -ends_with(".called")) %>% View()
# 3 snps only? weird
allaf %>% filter(is_snp & VD >= 8 & !is_passed) %>% count()  # 3 snps with high VD not passed, which should be rescued
allaf %>% reject_if(VD < 8) %>% filter(is_snp & VD >= 8 & !is_passed) %>% count()
allaf %>% filter(is_passed) %>% count(VD < 8)
allaf %>% filter(is_passed == F) %>% count(VD > 8)
allaf %>% rescue_if(VD >= 8) %>% count()
allaf %>% reject_if(VD < 8) %>% count()
# 105 snps are removed after we rescue VD>8:
allaf %>% rescue_if(VD >= 8) %>% filter(is_snp) %>% filter(!is_passed) %>% select(-ends_with(".truth"), -ends_with(".called")) %>% count()  # 105
allaf %>% filter(VD < 8)     %>% filter(is_snp) %>% filter(!is_passed) %>% select(-ends_with(".truth"), -ends_with(".called")) %>% count()  # 105
# why rescuing increases TP/FP count so drastically? because previously filtered variants become passed, regardless of FILT. 
#  can we reproduce this by rescuing certain FILT fields, and then applying reject_if()? have to remove all LowVD fields too... though they are already VD<8
show_stats(passed, allaf, allaf %>% reject_if(VD < 8), allaf %>% rescue_if(VD >= 8), allaf %>% rescue_if(FILT != "PASS"), allaf %>% rescue_if(FILT != "PASS") %>% reject_if(VD < 8))
# no, not working
# what can be the difference between reject_if(VD < 8) and rescue_if(VD >= 8)? let's see. reject drops is_passed & VD<8. rescue adds !is_passed & VD>=8.
#  meaning that we can take rescue_if(VD >= 8), and filter by VD<8?
allaf %>% filter(is_snp) %>% rescue_if(VD >= 8) %>% filter(VD < 8) %>% show_stats()  # getting exactly the difference 1125-987 = 138! nice 
# what exactly are we rescuing with VD>=8?
allaf %>% filter(is_snp) %>% rescue_if(VD >= 8) %>% filter(VD < 8) %>% count_status() %>% filter(is_tp) %>% select(VD, AF, DP, FILT, -ends_with(".truth"), -ends_with(".called"), everything()) %>% View()

# again, what kind of calls are rescued with rescue_if(VD >= 8)?
allaf %>% filter(!is_passed & VD >= 8) %>% select(VD, AF, DP, FILT, -ends_with(".truth"), -ends_with(".called"), everything()) %>% View()
# mostly StrandBias and HP, which is fine

# bcbio filter:
bcbio_filt = allaf %>% reject_if(VD < 6 & (MQ < 60.0 & NM > 2.0 | MQ < 55.0 & NM > 1.0))
bcbio_filt8 = allaf %>% reject_if(VD < 8 & (MQ < 60.0 & NM > 2.0 | MQ < 55.0 & NM > 1.0))
show_stats(allaf, bcbio_filt, bcbio_filt8)

# 11 FN SNPs left. Exploring them:
allaf %>% filter(is_called & !is_passed & str_detect(CALLS, 'GOLD')) %>% select(VD, AF, DP, FILT, -ends_with(".truth"), -ends_with(".called"), everything()) %>% View()
# All with VD<=4 +LCR, can't do anything with that

# exploring filter DP < avg_depth + 4 sqrt(avg_depth) from  https://www.nature.com/articles/s41592-018-0054-7.epdf?author_access_token=Dx2QKivDEkvswwWUJy5h0dRgN0jAjWel9jnR3ZoTv0Piny0X_8RjtB9zl5gHqNc_0m_Th8-X_roiQOsDJ4nU6Efwua2hRwes-IteFGhqVMhIXU1t_C1CZWZG5k8XMH65XRZZOcoziGXFh0QUciDDIg%3D%3D
called %>% filter(DP > 140) %>% select(VD, AF, DP, FILT, -ends_with(".truth"), -ends_with(".called"), everything()) %>% View()
# there are variants with DP ~300, but that might be just a CN gain, so no use

# compare gnomAD/PoN and low normal depth
merged %>% filter(is_called & str_detect(FILT, 'gnom')) %>% select(FILT, VD, AF, DP, MQ, NORMAL_VD, NORMAL_AF, NORMAL_DP, NORMAL_MQ, CHROM, POS, REF, ALT, CALLS)
merged %>% count(FILT) %>% View()





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

drop_pon = merged %>% rescue_if(str_detect(FILT, "PoN"))
show_stats(
  merged,
  drop_pon,
  drop_pon %>% reject_if(PoN >= 1),
  drop_pon %>% reject_if(PoN_NEW >= 1),
  drop_pon %>% reject_if(PoN_OLD >= 1),
  drop_pon %>% reject_if(PoN >= 2),
  drop_pon %>% reject_if(PoN_NEW >= 2),
  drop_pon %>% reject_if(PoN_OLD >= 2)
)
# Best result in SNPs: PoN_NEW>=2: F2=84.6%, achieving max TP and min FN, with the same FP as PoN_NEW>=1
# Best result in indels: PoN_OLD>=2, however PoN_NEW>=2 is fine also. The new one removes much more FP (8 vs 20), however misses 6 calls vs 4. 
# Explore on indels: trying without permissive overlap.
po = merged %>% rescue_if(str_detect(FILT, "PoN"))
show_stats(
  drop_pon,
  drop_pon %>% reject_if(PoN >= 1),
  drop_pon %>% reject_if(PoN >= 2),
  drop_pon %>% reject_if(PoN_NEW >= 1),
  po %>% reject_if(PoN_NEW >= 1),
  drop_pon %>% reject_if(PoN_OLD >= 1),
  po %>% reject_if(PoN_OLD >= 1),
  drop_pon %>% reject_if(PoN_NEW >= 2),
  po %>% reject_if(PoN_NEW >= 2),
  drop_pon %>% reject_if(PoN_OLD >= 2),
  po %>% reject_if(PoN_OLD >= 2)
)
# Indel stats got worse (F2 67% down to 57-58%) due to a higher number of FP. TODO: revisit idea of filtering surrounding indels.
# SNP stats changed slightly to become worse as well, due to overlap with indels in the PoN (e.g. variant like 
# indels:       1       35534250        .       AAT     A,AATAT,ATAT    0       .       PoN_samples=7
# snps:         1       35534250        .       A       T       0       .       PoN_samples=1
# )
# 
# Maybe make permissive overlap for SNPs too?
po = merged %>% rescue_if(str_detect(FILT, "PoN"))
show_stats(
  drop_pon,
  drop_pon %>% reject_if(PoN >= 1),
  drop_pon %>% reject_if(PoN >= 2),
  drop_pon %>% reject_if(PoN_NEW >= 1),
  po %>% reject_if(PoN_NEW >= 1),
  drop_pon %>% reject_if(PoN_OLD >= 1),
  po %>% reject_if(PoN_OLD >= 1),
  drop_pon %>% reject_if(PoN_NEW >= 2),
  po %>% reject_if(PoN_NEW >= 2),
  drop_pon %>% reject_if(PoN_OLD >= 2),
  po %>% reject_if(PoN_OLD >= 2)
)
# +1 missed call with NEW>=1, and same result with >=2 and the OLD one.

# Back to the old approach with permissive overlap only for indels.
# Now trying >=3
show_stats(
  merged,
  drop_pon,
  drop_pon %>% reject_if(PoN_NEW >= 1),
  drop_pon %>% reject_if(PoN_OLD >= 1),
  drop_pon %>% reject_if(PoN_NEW >= 2),
  drop_pon %>% reject_if(PoN_OLD >= 2),
  drop_pon %>% reject_if(PoN_NEW >= 3),
  drop_pon %>% reject_if(PoN_OLD >= 3)
)
# The old approach getting worse (especially for indels), but the new approach stays the same.

# Now removing gnomAD first to assess how the PoN helps specifically with artefacts.
no_gno = merged %>% 
  filter(!str_detect(FILT, "gnomAD_common"))
drop_pon = no_gno %>% rescue_if(str_detect(FILT, "PoN"))

new1 = drop_pon %>% reject_if(PoN_NEW >= 1)
new2 = drop_pon %>% reject_if(PoN_NEW >= 2)
new3 = drop_pon %>% reject_if(PoN_NEW >= 3)
old1 = drop_pon %>% reject_if(PoN_OLD >= 1)
old2 = drop_pon %>% reject_if(PoN_OLD >= 2)
old3 = drop_pon %>% reject_if(PoN_OLD >= 3)
show_stats(
  merged,
  no_gno,
  drop_pon,
  new1,
  new2,
  new3,
  old1,
  old2,
  old3
)
# SNPs: the new panel whos the best result regardless of the threshold.
# Indels: NEW>=1 cuts down the max FP, however misses 1 variant compared to >=2 and >=3. In any case, much better than the old one. 
#   Also, even better than the current PoN run on all samples: removes 1 FP indel, keeping the rest stats the same.




















