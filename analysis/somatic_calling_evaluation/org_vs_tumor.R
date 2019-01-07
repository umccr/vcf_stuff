library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(vcfR)
library(plyranges)
library(readr)



####################
### Load data
load_vars <- function(fname, purity) {
  raw <- read_tsv(str_c('/Users/vsaveliev/Analysis/SNV/org_tum/', fname, '.tsv'), col_types = 'cicccllididciccccd')
  raw %>% 
    mutate(
      PoN_CNT = replace_na(PoN_CNT, 0),
      GIAB_CONF = replace_na(GIAB_CONF, F),
      HMF_GIAB_CONF = replace_na(HMF_GIAB_CONF, F),
      HMF_HOTSPOT = ifelse(HMF_HOTSPOT == '.', F, T),
      purity = mean(purity)
    )
}

# N: CCR180147_MH18B001P062
# T: CCR180148_MH18F001P062, batch: 2016_249_18_MH_P062_1, purity: 37% - 44%
# O: CCR180172_VPT-MH062,    batch: 2016_249_18_MH_P062_2, purity: 94% - 100%
t6 <- load_vars('t6', purity = c(0.37, 0.44))
o6 <- load_vars('o6', purity = c(0.94, 1.00))

merged <- full_join(o6, t6,
  by = c('CHROM', 'POS', 'REF', 'ALT', 
         'GIAB_CONF', 'HMF_GIAB_CONF', 'PoN_CNT', 'PCGR_TIER', 
         'PCGR_SYMBOL', 'HMF_HOTSPOT', 'TRICKY', 'ENCODE', 'HMF_MAPPABILITY'), 
  suffix = c('.o', '.t'))

### Load Purple
read_purple <- function(fname) {
  # #chromosome  start      end        copyNumber  bafCount  observedBAF  actualBAF 
  # 1            1          2038216    2.1435      123       0.5422       0.5327    
  # segmentStartSupport  segmentEndSupport  method              depthWindowCount  gcContent  minStart   maxStart
  # TELOMERE             DEL                BAF_WEIGHTED        579               0.5759     1          1
  raw <- read_tsv(str_c('/Users/vsaveliev/Analysis/SNV/org_tum/', fname, '.purple.cnv'), col_types = 'ciididdcccidii') %>% 
    dplyr::transmute(
      start, 
      end,
      seqnames = `#chromosome`,
      CN = round(copyNumber), 
      BAF = actualBAF,
      # seg = str_c(segmentStartSupport, '-', segmentEndSupport),
      GC = gcContent) %>%
    dplyr::select(seqnames, start, end, CN, BAF, GC) %>% 
    as_granges()
}
t6_purple <- read_purple('t6')
o6_purple <- read_purple('o6')



###################
### Processing
data <- merged %>%
  mutate(is_snp = str_length(REF) == str_length(ALT),
         CH = str_c(REF, '>', ALT),
         nAF = ifelse(is.na(NORMAL_AF.t), NORMAL_AF.o, NORMAL_AF.t),
         nDP = ifelse(is.na(NORMAL_DP.t), NORMAL_DP.o, NORMAL_DP.t)) %>% 
  dplyr::rename(
    C = CHROM,
    POSITION = POS,
    AF.t = TUMOR_AF.t,
    DP.t = TUMOR_DP.t,
    AF.o = TUMOR_AF.o,
    DP.o = TUMOR_DP.o,
    FIL.t = FILTER.t,
    FIL.o = FILTER.o) %>% 
  mutate(
    VD.t = AF.t * DP.t,
    VD.o = AF.o * DP.o,
    AF.t = AF.t * 100.0,  # convert AFs to percentages
    AF.o = AF.o * 100.0,
    nAF.n = nAF * 100.0
  ) %>% 
  select(
    C, POSITION, CH, AF.t, DP.t, VD.t, AF.o, DP.o, VD.o, nAF, nDP, FIL.t, FIL.o, 
    everything(), 
    -ALT, -REF, -NORMAL_DP.o, -NORMAL_AF.o, -NORMAL_DP.t, -NORMAL_AF.t, -nAF)

### Add purple
data_p <- data %>% 
  dplyr::mutate(seqnames = C, start = POSITION, width = 1) %>% 
  as_granges() %>% 
  join_overlap_intersect(t6_purple) %>% 
  join_overlap_intersect(o6_purple, suffix = c(".t", ".o")) %>% 
  as_tibble() %>% 
  dplyr::select(-seqnames, -start, -end, -width, -strand)



#################
### Exploring

# TP
data_p %>% 
  filter(FIL.t == "PASS" & FIL.o == "PASS")
# 5409

# FN in tumor
data_p %>% 
  filter((FIL.t != "PASS" | is.na(FIL.t)) & FIL.o == "PASS")
# 3482 (all but 118 are not called in T at all)

# FP in tumor
data_p %>% 
  filter(FIL.t == "PASS" & (FIL.o != "PASS" | is.na(FIL.o)))
# 1829 (all but 8 are not called in O at all)

# Exploring FN in tumor: what kind of FN are called in organoid with a high AF, meaning that they should show up in tumor as well?
(high_af_missed_by_tumor <- data_p %>% 
    filter(AF.o > 40) %>%                       # high AF in organoid, meaning that tumor must catch this too
    filter(is.na(FIL.t) | FIL.t != 'PASS') %>%  # however, missed by tumor
    filter(FIL.o == 'PASS') %>%                 # not just called, but passed in organoid
    filter(!is.na(PCGR_SYMBOL), is_snp) %>%     # exploring non-intergenic snps to begin with
    select(-HMF_HOTSPOT, -ENCODE)
)
# Getting 28 such variants

# Now exploring FP in tumor: low AF tumor variants that are called, but not supported by organoid at all
(low_af_in_tumor_missed_by_o <- data_p %>% 
    filter(AF.t < 10) %>%                       # low AF in tumor
    filter(is.na(FIL.o) | FIL.o != 'PASS') %>%  # must be supported by organoid, but not
    filter(FIL.t == 'PASS')                     # not just called, but passed in tumor
)
# Total 1461

# Now compare AFs of matching variants. Ratio should approximately match purity ratio
(matching <- data_p %>% 
    filter(!is.na(FIL.o), !is.na(FIL.t)) %>% 
    mutate(
      norm_AF.t = AF.t / purity.t,
      norm_AF.o = AF.o / purity.o,
      norm_AFd = norm_AF.t - norm_AF.o)
)
# The head bunch of these variants are good: all PASS, normalised AF match, CN is 2. 
# Let's rather look at those that either T or O not passing:
matching %>% 
  filter(xor(FIL.t != 'PASS', FIL.o != 'PASS')) %>% View()
# The head bunch now have only T not-passing, which is expected given lower cellularity. All rejected variants have a LowVD. Total usch 126 variants.

# Interesting how many LowVD variants are rejected correctly? Checking.
lowvd <- data_p %>% 
  filter(str_detect(FIL.t, "LowVD"))
lowvd %>% 
  filter(is.na(FIL.o) | FIL.o != 'PASS')
# 1243 variants with LowVD in filter: not called or rejected in organoid
lowvd %>% 
  filter(is.na(FIL.o))
# 1162 variants with LowVD in filter: not called at all in organoid at all. 
lowvd %>% 
  filter(FIL.o == "PASS")
# 69 good oganoid variants were rejected in tumow as LowVD.
# We can see that our LowVD filter generally does a good job.

# For comparison, the total amount of filtered tumor variants are also not called or filtered in organoind:
data_p %>% 
  filter(!is.na(FIL.t) & FIL.t != 'PASS') %>%               # called, but filtered in tumor -do we filter well? 
  filter(is.na(FIL.o) | FIL.o != 'PASS') %>%                # not called or rejected in organoid
  filter()
# Total 3637 variants with any filter.




  
  


  
  
  
  
  


  