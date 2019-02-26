library(tidyverse)

df <- read_csv('/Users/vsaveliev/Downloads/UMCCR Samples - Sheet1.csv')
(normals <- df %>%
  select(sname = 'SampleName', type = 'Type', phenotype = 'Phenotype', bcbio_path = 'Results', notes = 'Notes') %>% 
  filter(!is.na(bcbio_path) & bcbio_path != 'n/a') %>% 
  filter(phenotype == "normal", type == "WGS") %>% 
  filter(is.na(notes)) %>% 
  mutate(sname = str_replace(sname, "17MHP031Bld-CCR170089_S1", "17MHP031Bld")) %>% 
  mutate(sname = str_replace(sname, "CCR170011_MH17B001P006_S6", "p006_normal")) %>% 
  mutate(sname = str_replace(sname, "CCR170094_MH17B001P032_S3", "p032_normal")) %>% 
  mutate(sname = str_replace(sname, "PRJ170195_GBB10_B", "GBB10_B")) %>% 
  mutate(sname = str_replace(sname, "kConFab-Mother-Blood", "kconfab_blood")) %>% 
  mutate(sname = str_replace(sname, "PRJ170198_SFRC01059_B", "diploid_blood")) %>% 
  mutate(sname = str_replace(sname, "PRJ170197_CUP_SC932_blood_S2", "cup_normal")) %>% 
  mutate(sname = str_replace(sname, "PRJ180514_CMN_N", "CMN_Normal")) %>% 
  mutate(sname = str_replace(sname, "2984Nor-PRJ170090", "IPMN2984_normal")) %>% 
  mutate(sname = str_replace(sname, "PRJ170062-2219-8073022ddnm_S2", "IPMN2219_normal")) %>% 
  mutate(sname = str_replace(sname, "3541Nor-PRJ170104", "IPMN3541_normal")) %>% 
  mutate(sname = str_replace(sname, "PRJ170052_IPMN1957_N_S4", "IPMN1957")) %>% 
  select(bcbio_path, sname) %>% 
  group_by(bcbio_path) %>% 
  summarise(samples = str_c(sname, collapse = ","))
)
  
normals %>% write_tsv('/Users/vsaveliev/git/umccr/vcf_stuff/vcf_stuff/panel_of_normals/normals.tsv', col_names = F)

