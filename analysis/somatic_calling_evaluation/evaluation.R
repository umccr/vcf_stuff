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


############
## Parsing

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

fix_somatic_anno_fields <- function(data) {
  renamed <- data %>% 
    mutate(
      HMF_MAPPABILITY = map_dbl(map(str_split(HMF_MAPPABILITY, ","), as.double), min)
    ) %>% 
    select(-matches("^GENE$|^CLNSIG$")) %>% 
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

# proc_tricky = function(vcf) {
#   for (key in names(vcf$vcf)) {
#     if (str_detect(key, "^TRICKY_|^AZ_")) {
#       print(key)
#       vcf$vcf[[key]] <- NULL
#     }
#   }
#   vcf
# }

extract_fmt_field = function(format, sample_data, field) {
  ind = format %>% str_split(":") %>% map(~ which(. == field))
  val = sample_data %>% str_split(":")
  map2_dbl(val, ind, function(v, i) { possibly(parse_double, NA)(v[[i]]) })
  
  # data %>% 
  #   mutate(
  #     tmp_fld_pos = FORMAT %>% str_split(":") %>% possibly(map_int, otherwise = c(NA))(~ which(. == field))
  #   ) %>% 
  #   mutate(
  #     !!new_fld_q := str_split(eval(sn, .), ":") %>% map(~ .[tmp_fld_pos][[1]])
  #   ) %>% 
  #   select(-tmp_fld_pos)
}

# "FORMAT:TTT" %>% str_split(":") %>% possibly(map_int, otherwise = c(NA))(~ which(. == 'TT'))
#called$FORMAT %>% head(43)
# extract_fmt_field(called$FORMAT %>% head(), called[[tumor_sample]] %>% head(), "NM", new_field = "NM_VD")
# called$NM = extract_fmt_field(called_vcf$vcf$FORMAT, called_vcf$vcf[[tumor_sample]], "NM", new_field = "NM_VD")
# add_format_field(called, "NM", tumor_sample, new_field = "NM_VD")

#merge_called_and_truth(truth_vcf, called_vcf, tumor_sample)

merge_called_and_truth = function(truth_vcf, called_vcf, tumor_sample) {
#  tumor_sample = substitute(tumor_sample)
  
  truth_data = truth_vcf$vcf %>% as_tibble()
  called_data = called_vcf$vcf %>% as_tibble()
  
  # print(tumor_sample)
  # print(class(tumor_sample))
  # print(called_data[[tumor_sample]])
  # print(called_data[[eval(tumor_sample)]])
  called_data$NM      = extract_fmt_field(called_data$FORMAT, called_data[[tumor_sample]], "NM")
  called_data$VD_QUAL = extract_fmt_field(called_data$FORMAT, called_data[[tumor_sample]], "QUAL")
  called_data$SBF     = extract_fmt_field(called_data$FORMAT, called_data[[tumor_sample]], "SBF")
  
  # a = called_data %>% head()
  # a$NM = extract_fmt_field(a$FORMAT, a[[tumor_sample]], "NM")
  # called_data %>%  select(CALLERS, VD, NM, MQ) %>% mutate(NM * MQ)
  # called_data %>% count(!is.na(NM), !is.na(TUMOR_MQ))
  
  merged <- full_join(
    truth_data %>% fix_somatic_anno_fields(), 
    called_data %>% fix_somatic_anno_fields(),
    by = c('CHROM', 'POS', 'REF', 'ALT'),
    suffix = c('.truth', '.called')) %>% 
    select(-matches("^AF$|^VD$|^DP$|^MQ$")) %>% 
    rename(
      AF = TUMOR_AF.called,
      VD = TUMOR_VD.called,
      DP = TUMOR_DP.called,
      FILT = FILTER.called,
      MQ = TUMOR_MQ.called,
      # 'COSM', 'ENCODE', 'GIAB', 'MBL', 'HMF_HS', 'ICGC', 'PCGR_TIER', 'CLNSIG', 
      # 'CSQ', 'DRIVER', 'PCGR_HS', 'TCGA', 'GENE', 'TRICKY', 'PoN'
    ) %>% 
    mutate(
      vartype = get_type(REF, ALT),
      is_snp = vartype == "SNP",
      is_called = !is.na(FILT),
      umccrise_passed = is_called & FILT == 'PASS',
      is_passed = umccrise_passed,
      is_true = !is.na(FILTER.truth),
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
}


############
## Showing stats

f_measure <- function(b, prec, recall) {
  Fmeasure <- (1 + b**2) * prec * recall / (b**2 * prec + recall)
}

build_stats <- function(.data) {
  stats <- .data %>% 
    count_status() %>% 
    group_by(vartype) %>% 
    summarise(
      called = sum(is_called),
      passed = sum(is_passed),
      true = sum(is_true),
      TP = sum(is_tp),
      FP = sum(is_fp),
      FN = sum(is_fn),
      recall = TP / true,
      prec = TP / passed,
      F2 = f_measure(2, prec, recall)
    ) %>% 
    # select(-passed, -true) %>% 
    select(TP, FP, FN, F2, called, true, vartype) %>% 
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

show_stats = function(...) {
  stat_dfs = map(list(...), build_stats)
  res = join_all(stat_dfs, by = c("metric", "vartype"), match = "all")  # suffix = argnames)
  argnames = sapply(substitute(list(...))[-1], deparse)  # %>% str_c(".", .)
  names(res)[c(-1, -2)] = argnames
  res
}


#############
## Plotting
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


#############
## Filtering
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

test_rescue = function() {
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
}



