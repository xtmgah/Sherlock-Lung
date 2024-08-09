
# Table S1: Basic information ---------------------------------------------
tmp1 <- readxl::read_xlsx('../../SmokingVar/Final_0_Public_WGS.xlsx',col_names = T) %>% mutate(WGS_Data_Source=paste0(Study_Temp_Name,' (', Deposition_ID,')')) %>% select(Tumor_Barcode,WGS_Data_Source) %>% unique()

tmp1 <- tibble(Tumor_Barcode=c('LU-FF76_tumor','D02326109','17MH0013_P1'),WGS_Data_Source=c('Lee-2019 (EGA: EGAS00001002801)','Imielinski (phs000488.v1.p1)','Imielinski (phs000488.v1.p1)')) %>% 
  bind_rows(tmp1)

load('../sampleset.RData')
load('../../BBsolution_final3_short.RData')
load('../../covdata0.RData')


saminfo <- covdata0 %>% filter(Histology == 'Adenocarcinoma', Tumor_Barcode %in% hq_samples) %>% select(Tumor_Barcode) %>% left_join(sp_group_data2) %>% select(Tumor_Barcode,Group=SP_Group_New) %>% 
  left_join(
    wgs_groups_info %>% select(Study,Tumor_Barcode)
  ) %>% 
  mutate(Source=if_else(str_detect(Study,'EAGLE'),'EAGLE-WGS',if_else(Study=='Public-WGS','Public WGS','Sherlock-Lung WGS'))) %>% 
  select(Source,Group,Tumor_Barcode) %>% 
  arrange(Group,desc(Tumor_Barcode))


tibble(Tumor_Barcode=sherlock_nonsmoker)%>% 
  left_join(tmp1) %>% 
  left_join(wgs_groups_info %>% select(Tumor_Barcode,Study)) %>% 
  mutate(Source=if_else(str_detect(Study,'EAGLE'),'EAGLE-WGS',if_else(Study=='Public-WGS','Public WGS','Sherlock-Lung WGS'))) %>% 
  mutate(WGS_Data_Source=if_else(Source=='Sherlock-Lung WGS','This study (phs001697)',WGS_Data_Source)) %>% 
  mutate(WGS_Data_Source=if_else(Source=='EAGLE-WGS','This study (phs002992)',WGS_Data_Source)) %>%
  mutate(WGS_Data_Source=if_else(is.na(WGS_Data_Source) & str_starts(Tumor_Barcode,'^TCGA'),'TCGA-new (phs000178)',WGS_Data_Source)) %>% 
  mutate(WGS_Data_Source=if_else(is.na(WGS_Data_Source) & str_starts(WGS_Data_Source,'LU-'),'Lee-2019 (EGA: EGAS00001002801)',WGS_Data_Source)) %>%
  mutate(WGS_Data_Source=if_else(is.na(WGS_Data_Source) & str_starts(WGS_Data_Source,'EBUS'),'Leong-2019 (EGA:EGAD00001005287)',WGS_Data_Source)) %>%
  left_join(covdata0) %>%
  select(-Smoking) %>% 
  mutate(Age=sprintf('%.2f',Age)) %>% 
  mutate(Tumor_Purity=sprintf('%.2f',Tumor_Purity)) %>% 
  #left_join(BBsolution4 %>% select(Tumor_Barcode,Tumor_Ploid=BB_Ploidy,WGD_Status,NRPCC)) %>% 
  #left_join(sherlock_PGA %>% select(-PGA) %>% rename(PGA=PGA_WGD)) %>% 
  #left_join(rnaseq_sample) %>% 
  #left_join(methylation_sample) %>% 
  write_csv('Supplementary_Tables_samplesource.csv',col_names = T,na = '')
