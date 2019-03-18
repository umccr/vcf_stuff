library(stringr)
library(stringi)
library(plyranges)
library(readr)
library(purrr)
library(crayon)
library(bedr)
library(tibble)
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
    # mutate(
    #   HMF_MAPPABILITY = map_dbl(map(str_split(HMF_MAPPABILITY, ","), as.double), min)
    # ) %>% 
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
      # MBL = HMF_MAPPABILITY,
      COSM = COSMIC_CNT,
      CSQ = PCGR_CONSEQUENCE
    ) %>% 
    mutate(
      HS = HMF_HS | !is.na(DRIVER) | !is.na(PCGR_HS) | 
        (!is.na(CLNSIG) & str_detect(CLNSIG, "pathogenic|uncertain")) |
        TCGA >= 5 | ICGC >= 3 | COSM >= 5 | PCGR_TIER %in% c("TIER_1", "TIER_2"),
      HS = replace_na(HS, F)
    ) %>% 
    select(-matches("^AF$|^VD$|^DP$|^MQ$")) %>% 
    rename(
      AF = TUMOR_AF,
      VD = TUMOR_VD,
      DP = TUMOR_DP,
      FILT = FILTER,
      MQ = TUMOR_MQ
    ) %>% 
    mutate(
      vartype = get_type(REF, ALT),
      is_snp = vartype == 'SNP',
      AF = as.double(AF),
      DP = as.integer(DP),
      VD = round(AF * DP)
    )
}

# nonna <- function(...) {
#   args = c(...)
#   not_na_indices = which(!is.na(args))
#   if (length(not_na_indices) == 0) {
#     return(NA)
#   } else {
#     return(args[[min(not_na_indices)]])
#   }
# }

set_tumor_field <- function(.data, field) {
  if (!field %in% names(.data)) {
    if (str_c(field, '.called') %in% names(.data)) {
      .data[[field]] = .data[[str_c(field, '.called')]]
    }
  }
  .data
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

merge_called_and_truth = function(called_vcf, truth_vcf, tumor_sample) {
#  tumor_sample = substitute(tumor_sample)

  called_data = called_vcf$vcf %>% as_tibble() %>% fix_somatic_anno_fields()
  
  if (!is.null(truth_vcf)) {
    truth_data = truth_vcf$vcf %>% as_tibble() %>% fix_somatic_anno_fields()
  } else {
    truth_data = called_data %>% head(0)
  }
  
  # called_data$NM      = extract_fmt_field(called_data$FORMAT, called_data[[tumor_sample]], "NM")
  # called_data$VD_QUAL = extract_fmt_field(called_data$FORMAT, called_data[[tumor_sample]], "QUAL")
  # called_data$SBF     = extract_fmt_field(called_data$FORMAT, called_data[[tumor_sample]], "SBF")
  
  # a = called_data %>% head()
  # a$NM = extract_fmt_field(a$FORMAT, a[[tumor_sample]], "NM")
  # called_data %>%  select(CALLERS, VD, NM, MQ) %>% mutate(NM * MQ)
  # called_data %>% count(!is.na(NM), !is.na(TUMOR_MQ))
  
  merged <- full_join(
    truth_data, 
    called_data,
    by = c('CHROM', 'POS', 'REF', 'ALT'),
    suffix = c('.truth', '.called')
  ) %>%
  set_tumor_field("vartype") %>% 
  set_tumor_field("is_snp") %>%        #     = nonna(is_snp.called    , is_snp.truth    ),
  set_tumor_field("AF") %>%        #         = nonna(AF.called        , AF.truth        ),
  set_tumor_field("VD") %>%        #         = nonna(VD.called        , VD.truth        ),
  set_tumor_field("DP") %>%        #         = nonna(DP.called        , DP.truth        ),
  set_tumor_field("MQ") %>%        #         = nonna(MQ.called        , MQ.truth        ),
  set_tumor_field("NORMAL_AF") %>%         #  = nonna(NORMAL_AF, NORMAL_AF.called , NORMAL_AF.truth ),
  set_tumor_field("NORMAL_VD") %>%         #  = nonna(NORMAL_VD, NORMAL_VD.called , NORMAL_VD.truth ),
  set_tumor_field("NORMAL_DP") %>%         #  = nonna(NORMAL_DP, NORMAL_DP.called , NORMAL_DP.truth ),
  set_tumor_field("NORMAL_MQ") %>%         #  = nonna(NORMAL_MQ, NORMAL_MQ.called , NORMAL_MQ.truth ),
  set_tumor_field("FILT") %>%        #       = nonna(FILT.called      , FILT.truth      ),
  set_tumor_field("GENE") %>%        #       = nonna(GENE.called      , GENE.truth      ),
  set_tumor_field("TCGA") %>%        #       = nonna(TCGA.called      , TCGA.truth      ),
  set_tumor_field("ICGC") %>%        #       = nonna(ICGC.called      , ICGC.truth      ),
  set_tumor_field("DRIVER") %>%        #     = nonna(DRIVER.called    , DRIVER.truth    ),
  set_tumor_field("CLNSIG") %>%        #     = nonna(CLNSIG.called    , CLNSIG.truth    ),
  set_tumor_field("PCGR_HS") %>%         #    = nonna(PCGR_HS.called   , PCGR_HS.truth   ),
  set_tumor_field("HMF_HS") %>%        #     = nonna(HMF_HS.called    , HMF_HS.truth    ),
  set_tumor_field("GIAB") %>%        #       = nonna(GIAB.called      , GIAB.truth      ),
  set_tumor_field("PoN") %>%         #        = nonna(PoN.called       , PoN.truth       ),
  # set_tumor_field("MBL") %>%         #        = nonna(MBL.called       , MBL.truth       ),
  set_tumor_field("COSM") %>%        #       = nonna(COSM.called      , COSM.truth      ),
  set_tumor_field("CSQ") %>%         #        = nonna(CSQ.called       , CSQ.truth       ),
  set_tumor_field("PCGR_TIER") %>%         #  = nonna(PCGR_TIER.called , PCGR_TIER.truth ),
  set_tumor_field("TRICKY") %>%        #     = nonna(TRICKY.called    , TRICKY.truth    ),
  set_tumor_field("HS") %>%        #         = nonna(HS.called        , HS.truth        )
  set_tumor_field("PoN") %>% 
  mutate(
    NORMAL_AF = as.double(NORMAL_AF),
    NORMAL_DP = as.integer(NORMAL_DP),
    NORMAL_VD = round(NORMAL_AF * NORMAL_DP)
  ) %>% 
  mutate(
    is_called = !is.na(FILT),
    umccrise_passed = is_called & FILT == 'PASS',
    is_passed = umccrise_passed,
    is_true = !is.na(FILT.truth)
  ) %>% 
  mutate(
    filters = str_split(FILT, ";"),
    rescued_filters = list(character())
  ) %>% 
  count_status()
}


############
## Showing stats

f_measure <- function(b, prec, recall) {
  Fmeasure <- (1 + b**2) * prec * recall / (b**2 * prec + recall)
}

summarize_stats = function(.data) {
  .data %>% 
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
    )
}

build_stats <- function(.data) {
  stats <- .data %>% 
    summarize_stats() %>% 
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

percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
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

rescue_filt = function(.data, filt) {
  new_filt_fields = c(filt)
  .data %>%
    mutate(
      rescued_filters = rescued_filters %>% map2(new_filt_fields, union) %>% map(~.[!is.na(.)])
    ) %>%
    mutate(
      filt_diff = filters %>% map(~.[!is.na(.)]) %>% map2(rescued_filters, setdiff)
      # rescued_filters = str_c(rescued_filters, sep = ';'),
      # filt_diff = str_c(filt_diff, sep = ';')
    ) %>%
    rescue_if(map_int(filt_diff, length) == 0)
}

test_resc_filt = function() {
  df = tribble(
    ~is_called, ~is_passed, ~vartype, ~FILT,
    TRUE,       TRUE,       "snp",    "", 
    TRUE,       F,          "snp",    "PoN",
    TRUE,       F,          "snp",    "gnomAD_common",
    TRUE,       F,          "snp",    "PoN;gnomAD_common",
    TRUE,       F,          "snp",    "LowVD", 
    TRUE,       F,          "snp",    "LowVD;PoN",
    TRUE,       F,          "snp",    "LowVD;gnomAD_common",
    TRUE,       F,          "snp",    "PoN;LowVD;gnomAD_common"
  ) %>% mutate(
    filters = str_split(FILT, ";"),
    rescued_filters = list(character())
  )
  
  df %>% rescue_filt("PoN") %>% rescue_filt("LowVD")
}

rescue_if = function(.data, cond) {
  cond_q = substitute(cond)
  mask = .data$is_passed | eval(cond_q, .data) %in% c(T)
  .data %>% mutate(is_passed = is_called & mask)
}

test_rescue_if = function() {
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



