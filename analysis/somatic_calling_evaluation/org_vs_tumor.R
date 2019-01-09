library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(stringi)
#library(vcfR)
library(plyranges)
library(readr)
library(purrr)
library(crayon)


####################
### Load data
load_vars <- function(fname, purity) {
  raw <- read_tsv(
    str_c('/Users/vsaveliev/Analysis/SNV/org_tum/', fname, '.tsv'), 
    col_types = 'cicccllididccccciiccccciicc',
    na = c("", ".", "NA"))
  raw %>% 
    mutate(
      PoN_CNT = replace_na(PoN_CNT, 0),
      GIAB_CONF = replace_na(GIAB_CONF, F),
      HMF_GIAB_CONF = replace_na(HMF_GIAB_CONF, F),
      HMF_HOTSPOT = ifelse(HMF_HOTSPOT == ".", F, T),
      purity = mean(purity),
      PCGR_INTOGEN_DRIVER_MUT = ifelse(PCGR_INTOGEN_DRIVER_MUT == "TRUE", T, F),
      PCGR_TIER = replace_na(str_extract(PCGR_TIER, '(\\d)'), 5)
    )
}


# t <- str_split(c("1,2", "3"), ",")
# m <- map(t, as.double)
# class(m)
# t6 %>% select(HMF_HOTSPOT, HMF_MAPPABILITY)


# N: CCR180147_MH18B001P062
# T: CCR180148_MH18F001P062, batch: 2016_249_18_MH_P062_1, purity: 37% - 44%
# O: CCR180172_VPT-MH062,    batch: 2016_249_18_MH_P062_2, purity: 94% - 100%
t6 <- load_vars('t6', purity = c(0.37, 0.44))
o6 <- load_vars('o6', purity = c(0.94, 1.00))

merged <- full_join(o6, t6,
  by = c('CHROM', 'POS', 'REF', 'ALT', 
         'GIAB_CONF', 'HMF_GIAB_CONF', 'PoN_CNT', 'HMF_HOTSPOT', 'TRICKY', 'ENCODE', 'HMF_MAPPABILITY',
         'PCGR_TIER', 'PCGR_CLINVAR_CLNSIG', 'PCGR_CONSEQUENCE', 'PCGR_INTOGEN_DRIVER_MUT', 'PCGR_MUTATION_HOTSPOT', 'PCGR_TCGA_PANCANCER_COUNT', 'PCGR_SYMBOL', 
         'MSI', 'MSILEN', 'LSEQ', 'RSEQ'), 
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
### Processing: fixing AF/DP/VD values
data <- merged %>%
  mutate(
    vt = ifelse(str_length(REF) == str_length(ALT), "SNP", "Indel"),
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
    FIL.o = FILTER.o,
    HP_ELEMENT = MSILEN,
    HP_REPEATED = MSI) %>% 
  mutate(
    VD.t = round(AF.t * DP.t),
    VD.o = round(AF.o * DP.o),
    AF.t = AF.t * 100.0,  # convert AFs to percentages
    AF.o = AF.o * 100.0,  # 
    nAF = nAF * 100.0,    #
    norm_AF.t = AF.t / purity.t,      # normalise AFs by sample purity
    norm_AF.o = AF.o / purity.o,      #
    norm_AFd = norm_AF.t - norm_AF.o, #
    HMF_MAPPABILITY = map(str_split(HMF_MAPPABILITY, ","), as.double),
    HS = replace_na(HMF_HOTSPOT | PCGR_INTOGEN_DRIVER_MUT | !is.na(PCGR_MUTATION_HOTSPOT) | str_detect(PCGR_CLINVAR_CLNSIG, "pathogenic") | str_detect(PCGR_CLINVAR_CLNSIG, "uncertain") | PCGR_TCGA_PANCANCER_COUNT >= 5, F),
    GIAB_CONF = HMF_GIAB_CONF
  ) %>% 
  select(
    C, POSITION, CH, AF.t, DP.t, VD.t, AF.o, DP.o, VD.o, nAF, nDP, FIL.t, FIL.o, HS, 
    everything(), 
    -NORMAL_DP.o, -NORMAL_AF.o, -NORMAL_DP.t, -NORMAL_AF.t, -nAF, 
    -HMF_HOTSPOT, -PCGR_INTOGEN_DRIVER_MUT, -PCGR_MUTATION_HOTSPOT, -PCGR_TCGA_PANCANCER_COUNT, -PCGR_CLINVAR_CLNSIG,
    -HMF_GIAB_CONF)


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
(matching <- data_p %>% 
   filter(!is.na(FIL.o), !is.na(FIL.t)))
# 6404
(matching_pass <- matching %>% 
  filter(FIL.t == "PASS" & FIL.o == "PASS"))
# 5409

# FN in tumor
(fneg <- data_p %>% 
  filter((FIL.t != "PASS" | is.na(FIL.t)) & FIL.o == "PASS"))
# 3482 (all but 118 are not called in T at all)

# FP in tumor
(fpos <- data_p %>% 
  filter(FIL.t == "PASS" & (FIL.o != "PASS" | is.na(FIL.o))))
# 1829 (all but 8 are not called in O at all)

# Exploring FN in tumor: what kind of FN are called in organoid with a high AF, meaning that they should show up in tumor as well?
(high_af_missed_by_tumor <- data_p %>% 
    filter(AF.o > 40) %>%                       # high AF in organoid, meaning that tumor must catch this too
    filter(is.na(FIL.t) | FIL.t != 'PASS') %>%  # however, missed by tumor
    filter(FIL.o == 'PASS') %>%                 # not just called, but passed in organoid
    select(-HMF_HOTSPOT, -ENCODE)
)
# Getting 140 such variants
#high_af_missed_by_tumor %>% write_tsv('high_af_missed_by_tumor.tsv')


# Now exploring FP in tumor: low AF tumor variants that are called, but not supported by organoid at all
(low_af_in_tumor_missed_by_o <- data_p %>% 
    filter(AF.t < 10) %>%                       # low AF in tumor
    filter(is.na(FIL.o) | FIL.o != 'PASS') %>%  # must be supported by organoid, but not
    filter(FIL.t == 'PASS')                     # not just called, but passed in tumor
)
# Total 1461

# Now compare AFs of matching variants. Ratio should approximately match purity ratio
(matching_with_afs <- matching_pass %>% 
    mutate(
      norm_AF.t = AF.t / purity.t,
      norm_AF.o = AF.o / purity.o,
      norm_AFd = norm_AF.t - norm_AF.o))
# The head bunch of these variants are good: all PASS, normalised AF match, CN is 2. 
# Let's rather look at those that either T or O not passing:
matching_with_afs %>% 
  filter(xor(FIL.t != 'PASS', FIL.o != 'PASS')) %>% View()
# The head bunch now have only T not-passing, which is expected given lower cellularity. All rejected variants have a LowVD. Total usch 126 variants.

# Plotting AF differences
matching_with_afs %>% 
  ggplot() +
  geom_histogram(aes(x = norm_AFd), binwidth = 1)

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
# 69 good organoid variants were rejected in tumow as LowVD.
# We can see that our LowVD filter generally does a good job.

# For comparison, the total amount of filtered tumor variants are also not called or filtered in organoind:
data_p %>% 
  filter(!is.na(FIL.t) & FIL.t != 'PASS') %>%               # called, but filtered in tumor -do we filter well? 
  filter(is.na(FIL.o) | FIL.o != 'PASS') %>%                # not called or rejected in organoid
  filter()
# Total 3637 variants with any filter.

# Count tumor with low AF vs high AF
tumor <- data_p %>% 
  filter(FIL.t == "PASS")
tumor_snps <- tumor %>% 
  filter(vt == "SNP")
tumor_snps %>% 
  count(AF.t > 10, PoN_CNT >= 1)

# If we use striter VD threshold, would it improve FP/FN rates?
tumor_snps %>% 
  count(AF.t > 10, FIL.o == "PASS")
for (t in c(2, 3, 4, 5, 6, 7)) {
  print(str_c("VD>=", t))
  print(tumor_snps %>% 
    filter(VD.t >= t) %>% 
    mutate(Organoid = ifelse(FIL.o != "PASS" | is.na(FIL.o), "no called", "PASS"),
           AF = ifelse(AF.t >= 10, ">=10%", "<10%")) %>%   
    count(AF, Organoid))
}

# AF distributions for variants missed in organoid
tumor_snps %>% 
  filter(FIL.t == "PASS") %>% 
  mutate(af_above10 = AF.t >= 10) %>% 
  mutate(Organoid = ifelse(FIL.o != "PASS" | is.na(FIL.o), "no called", "PASS")) %>% 
  # filter(!af_above10) %>% 
  ggplot() +
  geom_histogram(aes(x = AF.t), binwidth = 1) +
  facet_wrap(~Organoid, ncol = 1) +
  xlim(0, 100)


### HP
# Checking flanking sequences of MSI variants
data_p %>% 
  filter(vt == "Indel") %>% 
  select(C, POSITION, CH, ALT, REF, vt, AF.t, VD.t, AF.o, VD.o, PCGR_TIER, HS, HP_ELEMENT, HP_REPEATED, LSEQ, RSEQ) %>% 
  mutate(
    SEQ = str_c(LSEQ, "-", REF, "-", RSEQ),
    CH_len = abs(str_length(REF) - str_length(ALT)),
  ) %>% 
  select(-LSEQ, -RSEQ) %>% 
  filter(HP_REPEATED >= 5) %>%
  filter(CH_len == 1)

# Count the amount of hostpots TIER1-3 indels in homopolymers
data_p %>% 
  filter(HS, PCGR_TIER <= 3) %>% 
  filter(vt == "Indel")
# None! Yay

# We want to filter novel indels in homopolymer, which are likely artefacts  
(novel_indels <- data_p %>% 
  mutate(
    SEQ = str_c(LSEQ, "-", REF, "-", RSEQ),
    CH_len = abs(str_length(REF) - str_length(ALT)),
  ) %>%     
  filter(!is.na(FIL.t)) %>% 
  filter(!HS | PCGR_TIER >= 4) %>% 
  filter(vt == "Indel") %>% 
  select(-LSEQ, -RSEQ))

novel_indels %>% 
  mutate(
    CH_len = cut(CH_len, c(-Inf, 0, 1, 2, 3, 4, 5, Inf))
  ) %>% 
  dplyr::group_by(HP_ELEMENT, HP_REPEATED, CH_len) %>% 
  dplyr::summarise(
    n = dplyr::n(),
    VD = mean(VD.t, rm.na = T),
    AF = mean(AF.t, rm.na = T),
    AF = cut(AF, c(0, 10, 20, 30, 40, 100))
    ) %>% 
  ggplot() +
  geom_point(aes(x = HP_ELEMENT, y = HP_REPEATED, color = VD, size = n)) +
  facet_wrap(~CH_len, nrow = 1) + 
  scale_color_gradientn(colours = rainbow(5))


# Exploring CH len
novel_indels %>% 
  filter(HP_ELEMENT == 1, CH_len < 5) %>% 
  ggplot() +
  geom_bar(aes(x = cut(AF.t, c(0, 10, 20, 30, 40, 100)))) +
  facet_wrap(~CH_len)

novel_indels %>% count(HP_REPEATED)

# 3135
novel_indels %>% 
  filter(HP_REPEATED <= 5 | CH_len > 1)
# 1635
novel_indels %>%  
  filter(HP_REPEATED <= 5)
# 545

novel_indels %>% 
  filter(HP_REPEATED > 5)

# We want the change length to be a multiple of HP element length?
novel_indels %>% 
  filter(HP_REPEATED > 5) %>% 
  mutate(rem = CH_len %% HP_ELEMENT) %>% 
  select(HP_REPEATED, HP_ELEMENT, CH_len, rem, SEQ, CH, C, POSITION, AF.t, DP.t, VD.t, FIL.t, FIL.o) %>% 
  filter(rem != 0) %>% View()
  print(n = 100)
# All these changes are either a junction of 2 HP, or just a very low complexity region.

# Conclusion: filter out indels HP_REPEATED > 5
  
  
# %>% 
#   ggplot() +
#   geom_point(aes(x = HP_ELEMENT, y = HP_REPEATED, color = VD, size = n))


novel_indels %>% 
  filter(MSI >= 5) %>% 
#  filter(CH_len == 1) %>% 
  select(C, POSITION, CH, ALT, REF, vt, AF.t, VD.t, AF.o, VD.o, MSI, MSILEN, LSEQ, RSEQ, PCGR_TIER, everything()) %>% 
  mutate(SEQ = str_c(LSEQ, "-", REF, "-", RSEQ),
         CH_len = abs(str_length(REF) - str_length(ALT)),
  ) %>% 
  select(-LSEQ, -RSEQ, -ALT, -REF, -CH_len)



# Exploring different MSI lengths

data_p %>% 
group_by(MSI, MSILEN)


  
  filter(CH_len == 1) %>% 
  select(-LSEQ, -RSEQ, -ALT, -REF, -CH_len) %>% 
  filter(MSI < 5)
  

  # filter(MSI <= 4 & VD.t < 5) %>% 
  # filter(MSI <= 5 & VD.t < 6)

# MSI_FAIL if:
# msi <=  2 and af < 0.005,
# msi <=  4 and af < 0.01,
# msi <=  7 and af < 0.03,
# msi ==  8 and af < 0.06,
# msi ==  9 and af < 0.125,
# msi == 10 and af < 0.175,
# msi == 11 and af < 0.25,
# msi == 12 and af < 0.3,
# msi >  12 and af < 0.35])  
# change_len == 3 and msi >= 5 and af < 0.1:  # ignore low AF in 3nt MSI region




data_p %>% filter(POSITION > 32698781 & POSITION < 32698821)
  

?View



















  