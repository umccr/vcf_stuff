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

dir = "/Users/vsaveliev/Analysis/snv_validation/mb/ICGC_MB/"
truth_vcf = read.vcf(str_c(dir, "MB-benchmark.ANNO.FILT.vcf.gz"),
                     split.info = T, split.samples = T)
called_vcf = read.vcf(str_c(dir, "batch1-ensemble-annotated.ANNO.FILT.TP.SAMPLE.vcf.gz"), 
                      split.info = T)
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
called <- merged %>% mutate(is_passed = is_called)
passed <- merged %>% mutate(is_passed = umccrise_passed)
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
reject_if = function(.data, cond) {
  # Rejects additionally to existing is_passed
  cond_q <- substitute(cond)
  .data %>% mutate(is_passed = is_called & ifelse(eval(cond_q, .data) %in% c(F, NA), is_passed, F))
}

keep_if = function(.data, cond) {
  # Override existing is_passed if cond is true
  cond_q <- substitute(cond)
  .data %>% mutate(is_passed = is_called & ifelse(eval(cond_q, .data) %in% c(F, NA), is_passed, T))
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

show_stats(df, df %>% reject_if(NM < 1), df %>% keep_if(NM >= 1))
df
df %>% reject_if(NM < 1)
df %>% keep_if(NM >= 1)

###########
## Exploring filters
vd8_fi = passed %>% reject_if(VD < 8)
vd8_ki = passed %>% keep_if(VD >= 8 | is.na(VD))
allaf = passed %>% keep_if(FILT == "AF10") 
allaf_vd8_fi = allaf %>% reject_if(VD < 8)
allaf_vd8_ki = allaf %>% keep_if(VD >= 8 | is.na(VD))
show_stats(passed, vd8_fi, vd8_ki, allaf, allaf_vd8_fi, allaf_vd8_ki)
# Apparently keeping all VD>=8 improves the F2 measure. However, we again have trust issues
#   with the truth set.

# brad_filt = passed %>% filt_var(QD < 10.0 && AD[1] / (AD[1] + AD[0]) < 0.25 && ReadPosRankSum < 0.0
# bcbio filter:
bcbio_filt = allaf %>% reject_if(VD < 6 & (MQ < 60.0 & NM > 2.0 | MQ < 55.0 & NM > 1.0))
show_stats(allaf, allaf_vd8_fi, allaf_vd8_ki, bcbio_filt)

cool_filt = allaf_vd8_ki %>% reject_if(MQ < 60.0 & NM > 2.0 | MQ < 50.0 & NM > 1.0)
show_stats(allaf, allaf_vd8_fi, allaf_vd8_ki, bcbio_filt, cool_filt)







#################
## Exploring FN
passed %>% count(HS)

passed %>% filter(HS) %>% select(CHROM, POS, REF, ALT, AF, VD, FILT, GIAB, PoN, HMF_HS, TRICKY,
                                 PCGR_TIER, CLNSIG, DRIVER, ICGC, PCGR_HS, TCGA, GENE)
passed2 %>% count(ICGC >= 3)
passed2 %>% count(TCGA >= 5)

(fn = passed %>% filter(is_fn) %>% select(HS, everything()))

# TODO: explore CALLS and TIERS.truth

# SAGE. Explore how it changes in CCR180148_MH18F001P062-sage.vcf.gz (not MB unfortanately because MB doesn't have hotspots)
# - what are "inframe" hotspots?
# - add all PASS SAGE variants into the resulting VCF
# - add FILTER=SAGE_lowconf into resulting VCF if a passing variants is not confirmed by SAGE
# - extend the set of hotspots by adding PCGR sources?
# - CACAO: compare hotspots and genes with PCGR and HMF hotspots







