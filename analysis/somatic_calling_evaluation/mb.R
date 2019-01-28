library(ggplot2)
library(stringr)
library(stringi)
library(plyranges)
library(readr)
library(purrr)
library(crayon)
library(bedr)
# devtools::install_github("r-lib/pillar")
# library(pillar)
library(plyr)
library(tidyr)
library(dplyr)
library(patchwork)

dir = "/Users/vsaveliev/Analysis/snv_validation/mb/"
truth_vcf = read.vcf(str_c(dir, "MB-benchmark.ANNO.FILT.vcf.gz"), split.info = T)
called_vcf = read.vcf(str_c(dir, "batch1-ensemble-annotated.ANNO.FILT.TP.vcf.gz"), split.info = T)

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
      -GIAB_CONF,
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
    is_passed = is_called & FILT == 'PASS',
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

# TODO: why CLNSIG is missing??? Check PCGR run
merged %>% count(HS)
merged %>% count(COSM)

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
  res <- join_all(stat_dfs, by = c("metric", "vartype"), match = "all")  # suffix = argnames)
  argnames = sapply(substitute(list(...))[-1], deparse)  # %>% str_c(".", .)
  names(res)[c(-1, -2)]  <- argnames
  res
}

##############
### Exploring
called <- merged %>% mutate(is_passed = is_called)
passed <- merged
show_stats(called, passed)

# brad_filt <- merged %>% mutate(is_passed = 
# QD < 10.0 && AD[1] / (AD[1] + AD[0]) < 0.25 && ReadPosRankSum < 0.0
# merged %>% filter(!is.na(ReadPosRankSum)) %>% count(CALLERS)

plot_freq <- function(.data, metric, binwidth = NULL) {
  browser()
  .data %>% 
    count_status() %>% 
    filter(!is_tn) %>% 
    filter(!is.na(!!metric)) %>% 
    # mutate(AF = cut(AF, seq(0, 1, 0.05))) %>% 
    ggplot() +
    # geom_histogram(aes(get(metric), colour = status), binwidth = 0.05)
    # geom_density(aes(get(metric), colour = status))
    geom_freqpoly(aes(.data[[metric]], colour = status), binwidth = binwidth) +
    # geom_freqpoly(aes(get(metric), colour = status), bins = 30) +
    facet_wrap(~vartype, nrow = 2, scales = "free_y")
}

plot_freq(merged, "AF") + plot_freq(merged, "VD")
# 1. False positive SNPs are predominantly on low frequencies. 
#    False negative SNPs are on every AF, though VD tend to be low.
# 2. False positive indels are at any freq but mostly only below 0.35.
#    False negative inels are on 0.4-0.5, with VD at 40x - 
#    same shape as AF, so does not depend on DP much

plot_freq(merged, "AF", binwidth = 0.05) + plot_freq(merged %>% filter(TUMOR_VD.called < 25), "VD", binwidth = 1)

vd8 = merged %>% mutate(is_passed = (TUMOR_VD.called >= 8 | is.na(TUMOR_VD.called)))
passed_vd8 = merged %>% mutate(is_passed = (is_passed & TUMOR_VD.called >= 8 | is.na(TUMOR_VD.called)))
show_stats(called, passed, vd8, passed_vd8)
# apparently VD>=8 is the best filter. TODO: explore which filters we should remove

plot_freq(vd8, "TUMOR_AF.called", binwidth = 0.05) + 
  plot_freq(vd8 %>% filter(TUMOR_VD.called < 25), "TUMOR_VD.called", binwidth = 1)

plot_freq(passed_vd8, "TUMOR_AF.called", binwidth = 0.05) + 
  plot_freq(passed_vd8 %>% filter(TUMOR_VD.called < 25), "TUMOR_VD.called", binwidth = 1)
# VD>=8 good for indels, but not that good for SNPs



### Exploring FN
passed %>% count(HS)

passed %>% filter(HS) %>% select(CHROM, POS, REF, ALT, AF, VD, FILT, GIAB, PoN, HMF_HS, TRICKY,
                                 PCGR_TIER, CLNSIG, DRIVER, ICGC, PCGR_HS, TCGA, GENE)
passed2 %>% count(ICGC >= 3)
passed2 %>% count(TCGA >= 5)

(fn = passed %>% filter(is_fn) %>% select(HS, everything()))

# TODO: explore CALLS and TIERS.truth







