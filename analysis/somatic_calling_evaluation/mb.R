library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(stringi)
library(plyranges)
library(readr)
library(purrr)
library(crayon)
library(bedr)
# devtools::install_github("r-lib/pillar")
library(pillar)

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

merged <- full_join(truth_vcf$vcf, called_vcf$vcf,
    by = c('CHROM', 'POS', 'REF', 'ALT', 
           'COSMIC_CNT', 'ENCODE', 'GIAB_CONF', 'HMF_GIAB_CONF', 'HMF_MAPPABILITY', 'HMF_HOTSPOT', 'ICGC_PCAWG_HITS', 
           'PCGR_TIER', 'PCGR_CLINVAR_CLNSIG', 'PCGR_CONSEQUENCE', 'PCGR_INTOGEN_DRIVER_MUT', 'PCGR_MUTATION_HOTSPOT', 'PCGR_TCGA_PANCANCER_COUNT', 'PCGR_SYMBOL',
           'TRICKY', 'PoN_CNT'),
    suffix = c('.c', '.t'))

f_measure <- function(b, prec, recall) {
  Fmeasure <- (1 + b**2) * prec * recall / (b**2 * prec + recall)
}

get_type <- function(ref, alt) {
  var_types = c("SNP", "Indel", "MNP")
  vt = ifelse(str_length(ref) == str_length(alt), ifelse(str_length(ref) == 1, "SNP", "MNP"), "Indel")
  factor(vt, var_types, ordered = T)
}

build_stats <- function(data) {
  stats <- data %>% 
    mutate(
      vartype = get_type(REF, ALT),
      is_snp = vartype == "SNP",
      is_called = !is.na(FILTER.c) & is_passed,
      is_true = !is.na(TIERS.t)
    ) %>% 
    group_by(vartype) %>% 
    summarise(
      called = sum(is_called),
      true = sum(is_true),
      TP = sum(is_called & is_true),
      FP = sum(is_called & !is_true),
      FN = sum(!is_called & is_true),
      recall = TP / true,
      prec = FP / called,
      F2 = f_measure(2, prec, recall)
    ) %>% 
    select(-called, -true)
  
  metrics = names(stats %>% select(-vartype))
  
  # If you want samples to be columns, and metrics to be rows:
  stats_row <- stats %>%
    gather(metric, value, -vartype) %>%
    mutate(metric = factor(metric, levels = metrics)) %>%
    arrange(vartype)

  # If you want samples to be rows, and metrics to be columns:  
  # stats_col <- stats %>%
  #   gather(metric, value, -vartype) %>%
  #   mutate(metric = factor(metric, levels = metrics)) %>%
  #   arrange(vartype) %>%
  #   # unite(type_metric, vartype, metric, sep=".") %>%
  #   mutate(type_metric = interaction(metric, vartype), vartype = NULL, metric = NULL) %>%
  #   spread(type_metric, value)
}

called_stats <- build_stats(merged %>% mutate(is_passed = T))
passed_stats <- build_stats(merged %>% mutate(is_passed = FILTER.c == "PASS"))

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

# For samples to be columns, and metrics to be rows:  
dplyr::full_join(called_stats, passed_stats, by = c("metric", "vartype"), 
                 suffix = c(".called", ".passed"))










