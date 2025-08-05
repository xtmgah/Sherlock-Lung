# ------------------------------------------------------------------------------
# Script: Figure 2 - Tumor Evolution Analysis (LUAD cohort)
# Description: This script reproduces all main analyses and figures for tumor 
# evolution in high-quality LUAD samples.
# ------------------------------------------------------------------------------

# --- Load Required Libraries and Set Plotting Styles ---------------------------
# (You may need to install some packages if missing)
library(tidyverse)
library(ggplot2)
library(ggsankey)
library(scales)
library(hrbrthemes)   # For theme_ipsum_rc
library(cowplot)
library(forcats)
library(ggrepel)
library(ggnewscale)
library(data.table)
library(ggasym)
library(ggpmisc)
library(broom)
library(ggsci)

# --- Define Helper Functions (if any custom) -----------------------------------
# Please source your custom functions here if needed, or place them in ./functions/
# source('./functions/your_custom_functions.R')

# --- Load Input Datasets ------------------------------------------------------
# All data files should be placed in the appropriate subfolders.
# Example folder structure:
#   ./data/           - for main input files
#   ./functions/      - for custom R functions
#   ./output/         - for saving plots and results

# Replace the filenames below with your actual dataset filenames.
load('./data/BBsolution_final3_short.RData')
load('./data/sp_group_data.RData')
load('./data/covdata0.RData')
load('./data/clinical_data.RData')
load('./data/sherlock_data_all.RData')
load('./data/sherlock_variable.RData')
load('./data/sp_group_data.RData')
load('./data/id2data.RData')
load('./data/tedata.RData')
load('./data/RNASeq_Exp.RData',verbose = T)
load('./data/suvdata.RData',verbose = T)


# load function -----------------------------------------------------------
load('./data/ZTW_functions.RData')

# load analysis related data set
load('./data/Chronological_timing.RData',verbose = T)
tmp <- colnames(MRCAdata)
MRCAdata <- MRCAdata %>% select(-SP_Group) %>% left_join(sp_group_data2) %>% select(-SP_Group_New) %>% select(one_of(tmp))

## analysis limited to luad only
hq_samples2 <- covdata0 %>% filter(Histology == 'Adenocarcinoma', Tumor_Barcode %in% hq_samples) %>% pull(Tumor_Barcode)
rm(list=c('hq_samples'))

# =======================Original codes for performing association analysis for all genomic variables with latency variables ========================

#  Latency difference among all genome alteration or features
tdata <- MRCAdata %>% filter(!is.na(Latency)) %>% select(Tumor_Barcode,Acc=acceleration,Latency) #Histology == "Adenocarcinoma"; acceleration == '1x'
#tdata <- MRCAdata %>% filter(!is.na(Latency),acceleration == '1x',SP_Group %in% c('N_A','N_U')) %>% select(Tumor_Barcode,Latency) #Histology == "Adenocarcinoma"

## wicox test ## 
tdata <- sherlock_data_full %>% 
  filter(!(Type %in% c('Signature_Denovo'))) %>% 
  filter(Tumor_Barcode %in% tdata$Tumor_Barcode) %>% 
  mutate(Gene=paste0(Gene,"|",Type)) %>% 
  select(Tumor_Barcode,Gene,Alteration) %>% 
  left_join(tdata) %>% 
  #left_join(freq_mrca %>% dplyr::rename(Gene=name))  %>% 
  left_join(sp_group_data2)
#filter(Freq>0.03)

tdata_all <- bind_rows(
  tdata,
  tdata %>% mutate(SP_Group='ALL',SP_Group_New='ALL'),
  tdata %>% filter(SP_Group_New %in% c('EU_N','AS_N')) %>% mutate(SP_Group='Non-smokers',SP_Group_New='Non-smokers')
) 

#tdata <- tdata %>% filter(Freq > 0.03)

tresult_all <- tdata_all %>%
  group_by(SP_Group_New,Acc,Gene) %>%
  do(tresult = safely(wilcox.test)(Latency ~ Alteration, data=. )) %>%
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']]))) %>%
  select(SP_Group_New,Acc,Gene,fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup() %>%
  arrange(p.value) %>%
  mutate(TMP=Gene) %>%
  separate(col = 'TMP',into = c('Gene_short','Type'),sep = '\\|') %>%
  ungroup()

# add diff
tdata_all <- tdata_all %>% mutate(Alteration = case_when(
  Alteration %in% c("Female","Y","Yes","Non-Smoker","WGD") ~ "Yes",
  Alteration %in% c("Male","N","No","Smoker","nWGD") ~ "No",
  TRUE ~ NA_character_
)) 
tmp <- tdata_all %>% filter(!is.na(Alteration)) %>% group_by(SP_Group_New,Acc,Gene,Alteration) %>% summarise(value=median(Latency,na.rm = T)) %>% ungroup() %>% pivot_wider(names_from = 'Alteration',values_from = value) %>% mutate(diff=Yes-No)

# add freq
tmpsize <- tdata_all %>% drop_na() %>% select(SP_Group_New,Acc,Tumor_Barcode,Gene) %>% unique() %>% group_by(SP_Group_New,Acc) %>% count(Gene) %>% dplyr::rename(size=n)
freq_mrca <- tdata_all %>% drop_na()  %>% count(SP_Group_New,Acc,Gene,Alteration) %>% filter(!is.na(Alteration)) %>% left_join(tmpsize)%>% mutate(Freq=n/size) %>% group_by(SP_Group_New,Acc,Gene) %>% arrange(Freq) %>% dplyr::slice(1) %>% ungroup() %>% select(SP_Group_New,Acc,Gene,Freq) %>% mutate(Freq=if_else(Freq==1,0,Freq))

# tmpsize <- mdata %>% select(Tumor_Barcode,name) %>% unique()  %>% count(name) 
# freq_mrca <- mdata %>% count(name,value) %>% filter(!is.na(value)) %>% left_join(tmpsize)%>% mutate(Freq=n/size) %>% group_by(name) %>% arrange(Freq) %>% dplyr::slice(1) %>% ungroup() %>% select(name,Freq) %>% mutate(Freq=if_else(Freq==1,0,Freq))

tresult_all <- tresult_all %>% left_join(tmpsize) %>% left_join(tmp) %>% left_join(freq_mrca)

#saveRDS(tresult_all,file='latency_tresult.RDS')
tdata_all <- tdata_all %>% 
  mutate(TMP=Gene) %>%
  separate(col = 'TMP',into = c('Gene_short','Type'),sep = '\\|')

save(tresult_all,tdata_all,file='latency_tresult.RData')


# ======================== Fig. 2d - Latency association with mutational signatures ========================
# Signatures Only
load('latency_tresult.RData')
resulttmp <- tresult_all %>% 
  filter(Type=='Signature_Cosmic_final',!str_detect(Gene,'APOBEC')) %>% 
  filter(Acc=='1x') %>% 
  filter(SP_Group_New=='ALL') %>% 
  filter(Freq>0.01,Gene_short!='SBS5') %>% 
  mutate(Sig=if_else(str_detect(Gene,'^SBS'),'SBS',if_else(str_detect(Gene,'^DBS'),'DBS',if_else(str_detect(Gene,'^ID'),'ID',NA_character_)))) %>% 
  mutate(Sig=factor(Sig,levels=c('SBS','ID','DBS'))) %>% 
  group_by(Sig) %>% 
  mutate(FDR=p.adjust(p.value)) %>% 
  ungroup()

#filter(p.value<0.1)
#filter(str_detect(Gene_short,'^SBS'))

datatmp <- tdata_all %>% 
  filter(Type=='Signature_Cosmic_final') %>% 
  filter(Acc=='1x') %>% 
  filter(SP_Group_New=='ALL') %>% 
  filter(Gene_short %in% resulttmp$Gene_short) %>% 
  mutate(Sig=if_else(str_detect(Gene,'^SBS'),'SBS',if_else(str_detect(Gene,'^DBS'),'DBS',if_else(str_detect(Gene,'^ID'),'ID',NA_character_)))) %>% 
  mutate(Sig=factor(Sig,levels=c('SBS','ID','DBS')))

#tmp <- datatmp %>% count(Gene,Gene_short,Type,Alteration) %>% filter(n>5) %>% count(Gene,Gene_short,Type) %>% filter(n==2) %>% pull(Gene)
#datatmp <- datatmp %>% filter(Gene %in% tmp)
#resulttmp <- resulttmp %>% filter(Gene %in% tmp)

tmplevel <- datatmp %>% 
  filter(Alteration == "Yes") %>% 
  group_by(Gene_short) %>% 
  summarise(value=median(Latency)) %>% 
  arrange(value) %>% 
  pull(Gene_short) 

resulttmp <- resulttmp %>% mutate(Gene_short = factor(Gene_short,levels = tmplevel)) 
datatmp <- datatmp %>% mutate(Gene_short = factor(Gene_short,levels = tmplevel)) 

# barplot 
p1 <- datatmp %>% 
  #filter(Gene_short=='ID-Novel-C') %>% 
  ggplot(aes(Gene_short,Latency,fill=Alteration))+
  geom_point(pch=21,size=1,position=position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),color="black")+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25)+
  scale_fill_manual(values = c("Yes" = "#4daf4a","No" = "#cccccc"))+
  #scale_fill_manual(values = sigcol)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = "Y",grid_col = "gray90",ticks = T)+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),legend.position = 'top',panel.spacing.x = unit(0.1,'cm'))+
  labs(x = "", y = 'Latency, years before diagnosis')+
  scale_y_continuous(breaks = pretty_breaks())+ #,limits = c(-1.5,32)
  #guides(fill="none")+
  facet_grid(.~Sig,scales = 'free',space='free')+
  panel_border(color = 'black',linetype = 1)

p2 <- resulttmp %>%
  ggplot(aes(Gene_short,"1",fill=-log10(p.value)))+
  geom_tile()+
  facet_grid(.~Sig,scales = 'free',space='free')+
  scale_fill_viridis_c(option = "D")+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = FALSE,grid_col = "gray90",ticks = F)+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),legend.position = 'top',legend.key.width = unit(1,"cm"),legend.key.height = unit(0.3,"cm"),panel.spacing.x = unit(0.1,'cm'))+
  labs(x = "", y = '',fill="Wilcox test, -log10(P)\n")
#guides(fill="none")+
# panel_border(color = 'black',linetype = 1)

plot_grid(p2,p1,axis = 'v',align = 'lr',ncol = 1,rel_heights = c(1.1,4))

ggsave(filename = './output/latency_signatures.pdf',width = 16,height = 10,device = cairo_pdf)

## vacanol plot

resulttmp %>% 
  mutate(Sig=as.factor(Sig)) %>% 
  ggplot(aes((diff),-log10(FDR),fill=Sig))+
  geom_hline(yintercept = -log10(0.05),linetype=2,col='red',size=0.5)+
  geom_vline(xintercept = c(-5,5),linetype=2,col='blue',size=0.5)+
  geom_point(aes(size=Freq),pch=21,stroke=0.2)+
  scale_size_binned(labels = percent_format())+
  scale_fill_manual(values = ncicolpal[c(1,3,4)])+
  scale_x_continuous(breaks = pretty_breaks())+
  ggrepel::geom_text_repel(data=resulttmp,aes(label=Gene_short),max.overlaps = 30,size=3)+
  labs(x='Latency Change (Years)',y='-log10(FDR)',size='Frequency',fill='Signature')+
  guides(fill = guide_legend(override.aes = list(size=3.5)))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid='XY',ticks = T)+
  panel_border(color = 'black',size = 0.5)+
  coord_cartesian(clip = 'off')


ggsave(filename = './output/latency_fdr_signatures.pdf',width = 7,height = 5.5,device = cairo_pdf)


# ======================== Fig.  2a - Latency association with driver genes ========================
# Driver Mutations
#load('latency_tresult.RData')
drglist <- readRDS('./data/drivers_intogene.RDS') %>% pull(symbol)

resulttmp <- tresult_all %>% 
  filter(Type=='Mutation_Driver',Gene_short %in% drglist) %>% 
  filter(Acc=='1x') %>% 
  filter(SP_Group_New=='ALL') %>% 
  filter(Freq>0.01) %>% 
  mutate(FDR=p.adjust(p.value)) 

#filter(p.value<0.1)
#filter(str_detect(Gene_short,'^SBS'))

datatmp <- tdata_all %>% 
  filter(Type=='Mutation_Driver',Gene_short %in% drglist) %>% 
  filter(Acc=='1x') %>% 
  filter(SP_Group_New=='ALL') %>% 
  filter(Gene_short %in% resulttmp$Gene_short)

#tmp <- datatmp %>% count(Gene,Gene_short,Type,Alteration) %>% filter(n>5) %>% count(Gene,Gene_short,Type) %>% filter(n==2) %>% pull(Gene)
#datatmp <- datatmp %>% filter(Gene %in% tmp)
#resulttmp <- resulttmp %>% filter(Gene %in% tmp)

tmplevel <- datatmp %>% 
  filter(Alteration == "Yes") %>% 
  group_by(Gene_short) %>% 
  summarise(value=median(Latency)) %>% 
  arrange(value) %>% 
  pull(Gene_short) 

resulttmp <- resulttmp %>% mutate(Gene_short = factor(Gene_short,levels = tmplevel)) 
datatmp <- datatmp %>% mutate(Gene_short = factor(Gene_short,levels = tmplevel)) 
## vacanol plot

resulttmp %>% 
  ggplot(aes((diff),-log10(FDR)))+
  geom_hline(yintercept = -log10(0.05),linetype=2,col='red',size=0.5)+
  geom_vline(xintercept = c(-5,5),linetype=2,col='blue',size=0.5)+
  geom_point(aes(size=Freq),pch=21,stroke=0.1,fill='#3D4551',col='white')+
  scale_size_binned(labels = percent_format())+
  scale_x_continuous(breaks = pretty_breaks())+
  ggrepel::geom_text_repel(data=resulttmp %>% filter(FDR<0.2),aes(label=Gene_short),max.overlaps = 30,size=4)+
  labs(x='Latency Change (Years)',y='-log10(FDR)',size='Mutation frequency',fill='Signature')+
  guides(fill = guide_legend(override.aes = list(size=3.5)))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid='XY',ticks = T)+
  theme(legend.position = 'top',legend.key.width = unit(1,'cm'))+
  panel_border(color = 'black',size = 0.5)+
  coord_cartesian(clip = 'off')


resulttmp %>% 
  ggplot(aes((diff),-log10(p.value)))+
  geom_hline(yintercept = -log10(0.001),linetype=2,col='red',size=0.5)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col='orange',size=0.5)+
  geom_vline(xintercept = c(-5,5),linetype=2,col='blue',size=0.5)+
  geom_point(aes(size=Freq),pch=21,stroke=0.1,fill='#3D4551',col='white')+
  scale_size_binned(labels = percent_format())+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_y_continuous(breaks = pretty_breaks())+
  ggrepel::geom_text_repel(data=resulttmp %>% filter(p.value<0.05),aes(label=Gene_short),max.overlaps = 30,size=3.5)+
  labs(x='Latency Change (Years)',y='-log10(FDR)',size='Mutation frequency',fill='Signature')+
  guides(fill = guide_legend(override.aes = list(size=3.5)))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid='XY',ticks = T)+
  panel_border(color = 'black',size = 0.5)+
  coord_cartesian(clip = 'off')

ggsave(filename = './output/latency_driverGene_signatures.pdf',width = 5,height = 5.5,device = cairo_pdf)

## for EGFR
tmp <- resulttmp <- tresult_all %>% 
  filter(Type=='Mutation_Driver',Gene_short %in% 'EGFR') %>% 
  filter(Acc=='1x') %>% 
  filter(SP_Group_New %in% c('AS_N','EU_N','EU_S',"Others")) %>% 
  mutate(label=paste0(SP_Group_New,'\nP=',scientific_format()(p.value))) %>%  select(SP_Group_New,label)

tdata_all %>% 
  filter(Type=='Mutation_Driver',Gene_short %in% 'EGFR') %>% 
  filter(Acc=='1x') %>% 
  filter(SP_Group_New %in% c('AS_N','EU_N','EU_S',"Others")) %>% 
  mutate(Aternation=factor(Alteration,levels=c('Yes','No'))) %>% 
  left_join(tmp) %>% 
  ggplot(aes(Alteration,Latency,fill=Alteration))+
  geom_point(pch=21,size=2,position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.5),color="black")+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25)+
  scale_fill_manual(values = c("Yes" = "#4daf4a","No" = "#cccccc"))+
  #scale_fill_manual(values = sigcol)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = "Y",grid_col = "gray90",ticks = T)+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),legend.position = 'bottom',panel.spacing.x = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 2))+
  labs(x = "", y = 'Latency, years before diagnosis',fill='EGFR Mutation')+
  scale_y_continuous(breaks = pretty_breaks())+ #,limits = c(-1.5,32)
  #guides(fill="none")+
  facet_grid(.~label,scales = 'free',space='free')+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = './output/latency_driverGene_signatures2.pdf',width = 6,height = 6,device = cairo_pdf)


## for KRAS
tmp <- resulttmp <- tresult_all %>% 
  filter(Type=='Mutation_Driver',Gene_short %in% 'KRAS') %>% 
  filter(Acc=='1x') %>% 
  filter(SP_Group_New %in% c('AS_N','EU_N','EU_S',"Others")) %>% 
  mutate(label=paste0(SP_Group_New,'\nP=',scientific_format()(p.value))) %>%  select(SP_Group_New,label)


tdata_all %>% 
  filter(Type=='Mutation_Driver',Gene_short %in% 'KRAS') %>% 
  filter(Acc=='1x') %>% 
  filter(SP_Group_New %in% c('AS_N','EU_N','EU_S','Others')) %>% 
  mutate(Aternation=factor(Alteration,levels=c('Yes','No'))) %>% 
  left_join(tmp) %>% 
  ggplot(aes(Alteration,Latency,fill=Alteration))+
  geom_point(pch=21,size=2,position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.5),color="black")+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25)+
  scale_fill_manual(values = c("Yes" = "#4daf4a","No" = "#cccccc"))+
  #scale_fill_manual(values = sigcol)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = "Y",grid_col = "gray90",ticks = T)+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),legend.position = 'bottom',panel.spacing.x = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 2))+
  labs(x = "", y = 'Latency, years before diagnosis',fill='KRAS Mutation')+
  scale_y_continuous(breaks = pretty_breaks())+ #,limits = c(-1.5,32)
  #guides(fill="none")+
  facet_grid(.~label,scales = 'free',space='free')+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = './output/latency_driverGene_signatures_kras.pdf',width = 6,height = 6,device = cairo_pdf)

# # example of regression
# tdata_all %>% 
#   filter(Type=='Mutation_Driver',Gene_short %in% 'EGFR') %>% 
#   filter(Acc=='1x') %>% 
#   left_join(covdata0) %>% 
#   filter(SP_Group_New %in% c('AS_N','EU_N','EU_S')) %>%
#   group_by(SP_Group_New) %>% 
#   do(tidy(lm(Latency~Alteration+Gender+Histology+Tumor_Purity,data=.))) %>% 
#   filter(term=='AlterationYes') %>% 
#   arrange(p.value)
# #do(tidy(t.test(Latency~Alteration,data=.)))


# ======================== Fig.  2b-c - Latency association with EGFR, Gender, Sex ========================
tmp <- sherlock_data_full %>% 
  filter(Type=='Mutation_Driver',Gene=='EGFR') %>% 
  select(Tumor_Barcode,EGFR=Alteration)

tdata <- MRCAdata %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  filter(acceleration=='1x', Tumor_Barcode %in% hq_samples2) %>%
  left_join(sp_group_data2) %>% 
  left_join(covdata0 %>% select(Tumor_Barcode,Gender,Tumor_Purity)) %>% 
  filter(SP_Group_New %in% c('AS_N','EU_N','EU_S','Others')) %>% 
  left_join(tmp)

tdata %>% do(tidy(wilcox.test(Latency~Gender,data=.)))
tmp <- tdata %>% 
  group_by(SP_Group_New) %>% do(tidy(wilcox.test(Latency~Gender,data=.))) %>% 
  mutate(label=paste0(SP_Group_New,'\nP=',scientific_format()(p.value))) %>%  select(SP_Group_New,label)

tdata %>% 
  left_join(tmp) %>% 
  ggplot(aes(Gender,Latency,fill=Gender))+
  geom_point(pch=21,size=2,position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.5),color="black")+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25)+
  scale_fill_manual(values = c("Female" = "#4daf4a","Male" = "#cccccc"))+
  #scale_fill_manual(values = sigcol)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = "Y",grid_col = "gray90",ticks = T)+
  #theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),legend.position = 'top',panel.spacing.x = unit(0.1,'cm'))+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),legend.position = 'bottom',panel.spacing.x = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 2))+
  labs(x = "", y = 'Latency, years before diagnosis',fill='Sex')+
  scale_y_continuous(breaks = pretty_breaks())+ #,limits = c(-1.5,32)
  facet_grid(.~label,scales = 'free',space='free')+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = './output/latency_gender_signatures.pdf',width = 6,height = 5.8,device = cairo_pdf)


p1 <- tdata %>% 
  filter(SP_Group_New == 'EU_N') %>% 
  ggplot(aes(Gender,Latency,fill=Gender))+
  geom_point(pch=21,size=2,position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.5),color="black")+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25)+
  scale_fill_manual(values = c("Female" = "#984ea3","Male" = "#cccccc"))+
  #scale_fill_manual(values = sigcol)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = "Y",grid_col = "gray90",ticks = T)+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),legend.position = 'top',panel.spacing.x = unit(0.1,'cm'))+
  labs(x = "", y = 'Latency, years before diagnosis',fill='Sex')+
  scale_y_continuous(breaks = pretty_breaks())+ #,limits = c(-1.5,32)
  facet_grid(.~EGFR,scales = 'free',space='free')+
  panel_border(color = 'black',linetype = 1)


p2 <- tdata %>% 
  filter(SP_Group_New == 'EU_N') %>% 
  ggplot(aes(EGFR,Latency,fill=EGFR))+
  geom_point(pch=21,size=2,position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.5),color="black")+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25)+
  scale_fill_manual(values = c("Yes" = "#4daf4a","No" = "#cccccc"))+
  #scale_fill_manual(values = sigcol)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = "Y",grid_col = "gray90",ticks = T)+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),legend.position = 'top',panel.spacing.x = unit(0.1,'cm'))+
  labs(x = "", y = 'Latency, years before diagnosis',fill='EGFR Mutation')+
  scale_y_continuous(breaks = pretty_breaks())+ #,limits = c(-1.5,32)
  facet_grid(.~Gender,scales = 'free',space='free')+
  panel_border(color = 'black',linetype = 1)

plot_grid(p1,p2,align = 'v',axis = 'lr')


ggsave(filename = './output/latency_gender_EGFR_signatures2.pdf',width = 7.5,height = 5.5,device = cairo_pdf)

tdata %>% 
  filter(SP_Group_New=='EU_N') %>%
  group_by(Gender) %>% 
  do(tidy(wilcox.test(Latency~EGFR,data=.)))

tdata %>% 
  filter(SP_Group_New=='EU_N') %>%
  group_by(EGFR) %>% 
  do(tidy(wilcox.test(Latency~Gender,data=.)))


tdata %>% 
  group_by(SP_Group_New) %>% 
  do(tidy(wilcox.test(Latency~Gender,data=.)))

tmpresult <- tdata %>% 
  group_by(SP_Group_New) %>% 
  do(tidy(lm(Latency~Gender+EGFR+Tumor_Purity,data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term!='(Intercept)') %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  mutate(term=str_remove(term,"EGFR")) %>% 
  mutate(term=str_remove(term,"Histology")) %>% 
  mutate(term=str_remove(term,"Gender"))


tmpresult <- tmpresult %>% mutate(label=paste0('β = ',round(estimate,2), ', p = ',scientific_format(digits = 3)(p.value))) 
# %>%
#   bind_rows(
#     tibble(term=c("No","Male"),label = "",conf.low = 0, conf.high=0,estimate=0,p.value = 1) %>% mutate(SP_Group_New='AS_N'),
#     tibble(term=c("No","Male"),label = "",conf.low = 0, conf.high=0,estimate=0,p.value = 1) %>% mutate(SP_Group_New='EU_N'),
#     tibble(term=c("No","Male"),label = "",conf.low = 0, conf.high=0,estimate=0,p.value = 1) %>% mutate(SP_Group_New='EU_S')
#   )


# tmplevels <- c('Tumor_Purity','Male','Female','No','Yes')
# tmplables <- c('Tumor purity', 'Male','Female','EGFR wildtype','EGFR mutant')

tmplevels <- c('Tumor_Purity','Female','Yes')
tmplables <- c('Tumor purity', 'Female\nRef=Male','EGFR Mutant\nRef=Wildtype')

tmpresult %>% 
  mutate(term=factor(term,levels=tmplevels,labels=tmplables)) %>% 
  ggplot(aes(estimate, term, xmin = conf.low, xmax = conf.high, height = 0,color=p.value<0.05)) +
  geom_point(pch=19,size=3) +
  geom_vline(xintercept = 0, lty = 4) +
  geom_errorbarh(height = .3)+
  scale_color_manual(values = c("FALSE"='black',"TRUE"='red'))+
  ggrepel::geom_text_repel(aes(label=label),size=3.5,nudge_y = 0.2)+
  facet_grid(~SP_Group_New)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(panel.spacing = unit(0.2,"cm"))+
  labs(x = "Regression coefficient for tumor latency", y = NULL)+
  guides(color="none")+
  scale_x_continuous(breaks = pretty_breaks())+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = './output/latency_gender_EGFR_signatures.pdf',width = 14,height = 3,device = cairo_pdf)


# ======================== Fig.  2e - Latency association (Multivarible regression) ========================

tmp1 <- sherlock_data_full %>% 
  filter(Type=='Mutation_Driver',Gene=='EGFR') %>% 
  select(Tumor_Barcode,EGFR=Alteration)

tmp2 <- sherlock_data_full %>% 
  filter(Type=='Mutation_Driver',Gene=='KRAS') %>% 
  select(Tumor_Barcode,KRAS=Alteration)

tmp3 <- sherlock_data_full %>% 
  filter(Type=='Mutation_Driver',Gene=='TP53') %>% 
  select(Tumor_Barcode,TP53=Alteration)

tdata <- MRCAdata %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  filter(acceleration=='1x', Tumor_Barcode %in% hq_samples2) %>%
  left_join(sp_group_data2) %>% 
  left_join(covdata0 %>% select(Tumor_Barcode,Gender,Smoking,Assigned_Population,Tumor_Purity)) %>% 
  left_join(id2data %>% select(Tumor_Barcode, ID2_Present)) %>% 
  #filter(SP_Group_New %in% c('AS_N','EU_N','EU_S')) %>% 
  left_join(tmp1) %>% 
  left_join(tmp2) %>% 
  left_join(tmp3) %>% 
  mutate(Smoking = if_else(Smoking == 'Unknown',NA_character_,Smoking)) %>% 
  mutate(Smoking = factor(Smoking, levels=c('Smoker','Non-Smoker')))


tmpresult <- tdata %>% 
  do(tidy(lm(Latency~Gender+Assigned_Population+Smoking+ID2_Present+EGFR+KRAS+Tumor_Purity,data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term!='(Intercept)') %>% 
  ungroup() %>% 
  arrange(desc(estimate)) %>% 
  #mutate(term=str_remove(term,"EGFR")) %>% 
  mutate(term=str_remove(term,"Histology")) %>% 
  mutate(term=str_remove(term,"Gender")) %>% 
  mutate(label=paste0('β = ',round(estimate,2), ', p = ',scientific_format(digits = 3)(p.value)))


# tmpresult <- tmpresult %>% mutate(label=paste0('β = ',round(estimate,2), ', p = ',scientific_format(digits = 3)(p.value))) %>%
#   bind_rows(
#     tibble(term=c("No","Male"),label = "",conf.low = 0, conf.high=0,estimate=0,p.value = 1) %>% mutate(SP_Group_New='AS_N'),
#     tibble(term=c("No","Male"),label = "",conf.low = 0, conf.high=0,estimate=0,p.value = 1) %>% mutate(SP_Group_New='EU_N'),
#     tibble(term=c("No","Male"),label = "",conf.low = 0, conf.high=0,estimate=0,p.value = 1) %>% mutate(SP_Group_New='EU_S')
#   )


# tmplevels <- c('Tumor_Purity','Male','Female','No','Yes')
# 
# tmplables <- c('Tumor purity', 'Male','Female','EGFR wildtype','EGFR mutant')

tmplevels <-  tmpresult$term
tmplables <- tmpresult$term

#tmplables <- c('EGFR Mutant\nRef=Wildtype','Female\nRef=Male','Smokers\nRef=Nonsmokers','Tumor Purity','Ancestry EAS\nRef=EUR','Ancestry Others\nRef=EUR','KRAS Mutant\nRef=Wildtype','ID2 Signature Present\nRef=Absent')

tmplables <- c('EGFR Mutant\nRef=Wildtype','Female\nRef=Male','Tumor Purity','Ancestry EAS\nRef=EUR','Nonsmokers\nRef=Smokers','Ancestry Others\nRef=EUR','KRAS Mutant\nRef=Wildtype','ID2 Signature Present\nRef=Absent')


tmpresult %>% 
  mutate(p.value = if_else(p.value>0.05,0.5,p.value)) %>% 
  mutate(term=factor(term,levels=tmplevels,labels=tmplables)) %>% 
  ggplot(aes(estimate, term, xmin = conf.low, xmax = conf.high, height = 0,color=-log10(p.value))) +
  geom_vline(xintercept = 0, lty = 5) +
  geom_errorbarh(height = .3)+
  geom_point(pch=19,size=4) +
  #scale_color_manual(values = c("FALSE"='black',"TRUE"='red'))+
  scale_color_material(palette = 'red')+
  ggrepel::geom_text_repel(aes(label=label),size=3.5,nudge_y = 0.4)+
  #facet_grid(~SP_Group_New)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid_col = 'gray85')+
  theme(panel.spacing = unit(0.2,"cm"))+
  labs(x = "Regression coefficient", y = NULL)+
  guides(color="none")+
  scale_x_continuous(breaks = pretty_breaks())+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = './output/latency_multivariables.pdf',width = 8,height = 5,device = cairo_pdf)



summary(lm(Latency~Gender+Assigned_Population+Smoking+Tumor_Purity+Age,data=tdata))$r.squared

summary(lm(Latency~Gender+Assigned_Population+Smoking+ID2_Present+Tumor_Purity+Age,data=tdata))$r.squared

summary(lm(Latency~Gender+Assigned_Population+Smoking+EGFR+KRAS+Tumor_Purity+Age,data=tdata))$r.squared


summary(lm(Latency~Gender+Assigned_Population+Smoking+ID2_Present+EGFR+KRAS+Tumor_Purity+Age,data=tdata))$r.squared


tdata <- tdata %>% filter(Smoking == 'Non-Smoker',Assigned_Population %in% c('EAS','EUR'))
#summary(lm(Latency~SP_Group_New+ID2_Present + KRAS +EGFR+SP_Group_New:EGFR + Tumor_Purity + Gender,data=tdata)) %>% tidy() %>% arrange(p.value)
summary(lm(Latency~Assigned_Population+EGFR+EGFR:Assigned_Population + Tumor_Purity + Gender,data=tdata)) %>% tidy() %>% arrange(p.value)




# Additional Supplementary Figures ----------------------------------------

# Supplementary Fig. 12 ---------------------------------------------------

#Supplementary Fig. 12a, overall 
mrate %>% 
  left_join(sp_group_data) %>%
  #filter(Tumor_Barcode %in% hq_samples2) %>% 
  #filter(SP_Group!="Others") %>% 
  drop_na() %>% 
  ggplot(aes(age,log2(CpG.TpG.Gb)))+
  geom_point(pch=21,size=2,fill=ncicolpal[14]) +
  geom_smooth(method = 'lm')+
  #facet_wrap(~SP_Group_New,nrow = 1,scales = 'free_x')+
  scale_fill_manual(values = sp_group_color_new)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  labs(x = "Age at diagnosis (years)", y = "SNVs/Gb (log2)")+
  guides(fill="none")+
  scale_x_continuous(breaks = pretty_breaks())+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = './output/mutation_rate_age_all.pdf',width = 4,height = 3.6,device = cairo_pdf)

# Supplementary Fig. 12b, seperated by the group
mrate %>% 
  left_join(sp_group_data) %>% 
  #filter(SP_Group!="Others") %>% 
  #filter(Tumor_Barcode %in% hq_samples2) %>% 
  drop_na() %>% 
  ggplot(aes(age,log2(CpG.TpG.Gb),fill=SP_Group_New))+
  geom_point(pch=21,size=2) +
  geom_smooth(method = 'lm')+
  facet_wrap(~SP_Group_New,nrow = 1,scales = 'free_x')+
  scale_fill_manual(values = sp_group_color_new)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  labs(x = "Age at diagnosis (years)", y = "SNVs/Gb (log2)")+
  guides(fill="none")+
  scale_x_continuous(breaks = pretty_breaks())+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = './output/mutation_rate_age.pdf',width = 12,height = 3.5,device = cairo_pdf)

mrate %>% do(tidy(lm(log2(CpG.TpG.Gb)~age,data=.))) %>% filter(term!="(Intercept)") 

mrate %>% group_by(SP_Group) %>% do(tidy(lm(log2(CpG.TpG.Gb)~age,data=.))) %>% filter(term!="(Intercept)") %>% ungroup() %>% arrange(p.value)
#mrate %>% left_join(clinical_data) %>% filter(Histology == 'Adenocarcinoma') %>% group_by(SP_Group) %>% do(tidy(lm(log2(CpG.TpG.Gb)~age,data=.))) %>% filter(term!="(Intercept)") %>% ungroup() %>% arrange(p.value)
#sorder <- read_rds('../SmokingVar/sorder.RDS')
#mrate %>% left_join(sorder) %>% filter(!is.na(APOBEC_Signature)) %>% group_by(APOBEC_Signature) %>% do(tidy(lm(CpG.TpG.Gb~age,data=.))) %>% filter(term!="(Intercept)") %>% ungroup() %>% arrange(p.value)


# Supplementary Fi. 13a ---------------------------------------------------
# average rate as barplot
as.data.frame(qRateDeam) %>% 
  rownames_to_column() %>% 
  pivot_longer(cols = -rowname) %>% 
  mutate(rowname = paste0("d",str_remove(rowname,"%"))) %>% 
  pivot_wider(names_from = rowname,values_from = value) %>% 
  filter(name!="Others") %>% 
  dplyr::rename(SP_Group=name) %>% 
  left_join(sp_group_data) %>% 
  ggplot(aes(SP_Group_New,d50,fill=SP_Group_New))+
  geom_col(width = 0.7)+
  scale_fill_manual(values = sp_group_color_new)+
  guides(fill="none")+
  geom_errorbar(aes(ymin=d0,ymax=d100),width=0.1,size=0.7)+
  labs(y = "CpG>TpG rate [SNVs/Gb/yr]", x= "")+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = './output/mutation_rate_age_average.pdf',width = 3,height = 4,device = cairo_pdf)



# Latency Plot by differnt accelration rate ------------------------------------------------------------
u <- base::setdiff(names(finalSnv), remove)
guessAccel <- sapply(subclonesTimeAbs, function(x) "5x")
qSubclone <- sapply(subclonesTimeAbs, function(x) apply(x[,"hat",][rownames(x)%in%u,,drop=FALSE], 2, quantile, c(0.05,0.25,0.5,0.75,0.95), na.rm=TRUE), simplify='array')
subclonesTimeAbsType <- sapply(names(subclonesTimeAbs), function(n) {x <- subclonesTimeAbs[[n]]; x[,,guessAccel[n]][setdiff(rownames(x),remove), 1:3, drop=FALSE]})
nSubclones <- sapply(subclonesTimeAbsType, function(x) sum(!is.na(x[,1])))

qSubclone <- qSubclone[,,-2]

m <- qSubclone["50%","5x",]#t[1,3,]
names(m) <- dimnames(qSubclone)[[3]]
m[nSubclones < 5] <- NA
o <- rev(order(m, na.last=NA))
n <- dimnames(qSubclone)[[3]]

cairo_pdf(filename = 'Acceleration_alterative.pdf',width = 3,height = 5)
par(mar=c(3,3,1,1), mgp=c(2,.5,0), tcl=0.25,cex=1, bty="L", xpd=FALSE, las=1)
plot(NA,NA, xlim=c(0.5,length(m[o])+0.5), ylab="Latency [yr]", xlab="", xaxt="n", yaxs="i", ylim=c(0,18))
abline(h=seq(10,20,10), col="#DDDDDD", lty=3)
x <- seq_along(m[o])
mg14::rotatedLabel(x, labels=as.character(sp_group_lift[names(rev(sort(m)))]))
b <- .3
rect(seq_along(o)-b,qSubclone["50%","1x",o],seq_along(o)+b,qSubclone["50%","10x",o], col=alpha(as.character(sp_group_color[n[o]]),0.4), border=1)
rect(seq_along(o)-b,qSubclone["50%","2.5x",o],seq_along(o)+b,qSubclone["50%","7.5x",o], col=alpha(as.character(sp_group_color[n[o]]),0.9), border=1)
rect(seq_along(o)-b,qSubclone["50%","20x",o],seq_along(o)+b,qSubclone["50%","10x",o], col=alpha(as.character(sp_group_color[n[o]]),0.2), border=1)
segments(seq_along(o)-b,qSubclone["50%","5x",o],seq_along(o)+b,qSubclone["50%","5x",o])
dev.off()



# Tumor latency vs Grade --------------------------------------------------
tmp <- readxl::read_xlsx('../Clinical/Data_harmonization/Delivery 2024-09-09/sherlock_20240909.xlsx') %>% clean_names()

grade_info <- tmp %>% select(sherlock_pid,subject_id,study_site,grade_differentiation,grade_diff_later_cancer) %>% filter(!is.na(grade_differentiation),!is.na(sherlock_pid)) %>% 
  mutate(sherlock_pid = str_remove(sherlock_pid,'^NSLC-')) %>% 
  left_join(
    wgs_groups_info %>% select(sherlock_pid=Sherlock_PID,Tumor_Barcode)
  ) %>% 
  filter(!is.na(Tumor_Barcode)) %>% 
  select(Tumor_Barcode,grade_differentiation) %>% 
  filter(grade_differentiation %in% c(1,2,3,4)) %>%
  mutate(Grade=if_else(grade_differentiation %in% c(1,2),'Low grade (1&2)','High grade (3&4)')) 


tdata <- MRCAdata %>% 
  filter(acceleration == "1x") %>% 
  mutate(Group=if_else(Latency>8,'High','Low')) %>%
  left_join(grade_info) %>% 
  filter(grade_differentiation %in% c(1,2,3,4)) %>%
  left_join(sp_group_data2)  %>% 
  left_join(wgs_groups_info %>% select(Study,Tumor_Barcode)) 
  #filter(Study %in% c('Hong_Kong','Toronto'))

tdata %>% select(Grade,Group) %>% table() %>% fisher.test()

tdata %>% 
  count(Group,Grade) %>% 
  mutate(Grade=fct_rev(Grade),Group=fct_rev(Group)) %>% 
  ggplot(aes(Group,n,fill=Grade,Latency))+
  geom_col(position = "fill",width = 0.85)+
  scale_fill_jama()+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16)+
  theme(plot.margin = margin(4,4,4,4),plot.title = element_text(hjust = 0.5,face = 'plain'))+
  labs(x = "Tumor latency", y = 'Proportion (%)',fill='Tumor grade',title = 'P=0.0088')+
  #guides(fill="none")+
  panel_border(color = 'black',linetype = 1)

#ggsave(filename = './output/latency_grade_enrichment_hongkong_toronto.pdf',width = 3.6,height = 6.5,device = cairo_pdf)


tdata <- tdata %>% 
  filter(Study %in% c('Hong_Kong','Toronto')) %>% 
  mutate(Study='Study site (Hong Kong + Toronto)') %>% 
  bind_rows(
    tdata %>% mutate(Study="ALL")
  )
  
my_comparisons <- list(c('Low grade (1&2)','High grade (3&4)'))
tdata%>% 
  mutate(Grade=fct_rev(Grade)) %>% 
  ggplot(aes(Grade,Latency,fill=Grade))+
  geom_boxplot(width=0.6,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1.5,position = position_jitter(width = 0.15,height = 0.05),color="white",stroke=0.02)+
  scale_fill_jama()+
  facet_wrap(~Study,nrow = 1)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 4),axis.text.x = element_text(angle = 30,hjust = 1))+
  labs(x = "", y = 'Tumor latency years')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = './output/latency_grade_enrichment_boxplot.pdf',width = 5,height = 7,device = cairo_pdf)






# ID2/L1 vs tumor grade ---------------------------------------------------

tdata <- id2data %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  left_join(grade_info) %>% 
  filter(!is.na(Grade)) %>% 
  left_join(wgs_groups_info %>% select(Study,Tumor_Barcode))
  #filter(Study %in% c('Hong_Kong','Toronto'))

tdata %>% select(ID2_Present,Grade) %>% table() %>% fisher.test()

tdata %>% 
  count(ID2_Present,Grade) %>% 
  mutate(Grade=fct_rev(Grade),ID2_Present=fct_rev(ID2_Present)) %>% 
  ggplot(aes(ID2_Present,n,fill=Grade))+
  geom_col(position = "fill",width = 0.75)+
  scale_fill_jama()+
  #facet_wrap(~Study)+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16)+
  theme(plot.margin = margin(4,4,4,4),plot.title = element_text(hjust = 0.5,face = 'plain'))+
  labs(x = "Mutational signature ID2", y = 'Proportion (%)',fill='Tumor grade',title = 'P=0.00013')+
  #guides(fill="none")+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = './output/ID2_grade_enrichment.pdf',width = 4,height = 6.5,device = cairo_pdf)



tdata <- tedata %>% 
  mutate(L1_Present=if_else(Total_L1>0,'Present','Absent')) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  left_join(grade_info) %>% 
  filter(!is.na(Grade)) %>% 
  left_join(wgs_groups_info %>% select(Study,Tumor_Barcode))
#filter(Study %in% c('Hong_Kong','Toronto'))

tdata %>% select(L1_Present,Grade) %>% table() %>% fisher.test()

tdata %>% 
  count(L1_Present,Grade) %>% 
  mutate(Grade=fct_rev(Grade),L1_Present=fct_rev(L1_Present)) %>% 
  ggplot(aes(L1_Present,n,fill=Grade))+
  geom_col(position = "fill",width = 0.75)+
  scale_fill_jama()+
  #facet_wrap(~Study)+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16)+
  theme(plot.margin = margin(4,4,4,4),plot.title = element_text(hjust = 0.5,face = 'plain'))+
  labs(x = "Somatic L1 insertion", y = 'Proportion (%)',fill='Tumor grade',title = 'P=0.026')+
  #guides(fill="none")+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = './output/L1_grade_enrichment.pdf',width = 4,height = 6.5,device = cairo_pdf)



# Tumor latency vs proliferation maker ---------------------
tmp <- rdata1 %>% 
  filter(RNAseq_Type == 'Tumor', Gene %in% c('MYBL2','BUB1','PLK1','CCNE1','CCNB1','BUB1','FOXM1','TOP2A','MKI67'))

tdata <- MRCAdata %>% 
  filter(acceleration == "1x") %>% 
  mutate(Group=if_else(Latency>8,'High','Low')) %>% 
  left_join(tmp) %>% 
  left_join(sp_group_data2) %>% 
  filter(!is.na(Exp),!is.na(Latency))


tdata %>% group_by(Gene) %>% do(tidy(wilcox.test(Exp~Group,data=.))) %>% ungroup() %>%  mutate(FDR=p.adjust(p.value))

my_comparisons <- list(c("High",'Low'))
tdata%>% 
  mutate(Group=fct_rev(Group)) %>% 
  ggplot(aes(Group,Exp,fill=Group))+
  geom_boxplot(width=0.7,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1.5,position = position_jitter(width = 0.15,height = 0.05),color="white",stroke=0.02)+
  scale_fill_jama()+
  facet_wrap(~Gene,nrow = 1)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 4))+
  labs(x = "Tumor latency", y = 'RNA-Seq expression log2(CPM)')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = './output/latency_proliferation_expression.pdf',width = 8,height = 6,device = cairo_pdf)



# TP53 vs Tumor latenc/MRCA age

tdata <- MRCAdata %>% 
  filter(acceleration == "1x") %>% 
  mutate(Group=if_else(Latency>8,'High','Low')) %>% 
  left_join(
    sherlock_data_full %>% filter(Gene=='TP53',Type=='Mutation_Driver') %>% select(Tumor_Barcode,TP53=Alteration)
  ) %>% 
  left_join(sp_group_data2)

tdata <- tdata %>% mutate(SP_Group_New = 'ALL') %>% bind_rows(tdata)

my_comparisons <- list(c("Yes",'No'))

tdata %>% 
  ggplot(aes(TP53,MRCA_age,fill=TP53))+
  geom_boxplot(width=0.7,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1.5,position = position_jitter(width = 0.15,height = 0.05),color="white",stroke=0.02)+
  scale_fill_jama()+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  facet_wrap(~SP_Group_New,nrow = 1)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 4))+
  labs(x = "TP53 mutation", y = 'Tumor MRCA years')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = './output/TP53_MRCA_years.pdf',width = 5,height = 5,device = cairo_pdf)


tmp <- rdata1 %>% 
  filter(RNAseq_Type == 'Tumor', Gene %in% c('MYBL2','BUB1','PLK1','CCNE1','CCNB1','BUB1','FOXM1','TOP2A','MKI67'))

sherlock_data_full %>%
  filter(Gene=='TP53',Type=='Mutation_Driver') %>% 
  select(Tumor_Barcode,TP53=Alteration) %>%
  left_join(tmp) %>% 
  left_join(sp_group_data2) %>% 
  filter(Tumor_Barcode %in% hq_samples2,!is.na(Exp)) %>% 
  ggplot(aes(TP53,Exp,fill=TP53))+
  geom_boxplot(width=0.7,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1.5,position = position_jitter(width = 0.15,height = 0.05),color="white",stroke=0.02)+
  scale_fill_jama()+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  facet_wrap(~Gene,nrow = 1)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 4))+
  labs(x = "TP53 mutation", y = 'RNA-Seq expression log2(CPM)')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)
  #stat_compare_means(comparisons = my_comparisons)

ggsave(filename = './output/TP53_proliferation_expression.pdf',width = 9,height = 6,device = cairo_pdf)



# Latency difference between ID2 present vs ID2 absent in KRAS mutant group --------------------------
my_comparisons <- list(c("Absent",'Present'))
MRCAdata %>% 
  filter(acceleration == '1x') %>% 
  left_join(id2data) %>% 
  left_join(
    sherlock_data_full %>%
      filter(Gene=='KRAS',Type=='Mutation_Driver') %>% 
      select(Tumor_Barcode,KRAS=Alteration)
  ) %>% 
  left_join(sp_group_data2) %>% 
  filter(KRAS=='Yes',Tumor_Barcode %in% hq_samples2) %>% 
  filter(SP_Group_New %in% c('EU_N','EU_S')) %>% 
  ggplot(aes(ID2_Present,Latency,fill=ID2_Present))+
  geom_boxplot(width=0.6,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=3,position = position_jitter(width = 0.15,height = 0),color="black",stroke=0.02)+
  scale_fill_manual(values = id2color)+
  facet_wrap(~SP_Group_New)+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 4))+
  labs(x = "Mutational signature ID2", y = 'Tumor latency years')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons,method = 't.test',method.args = list(alternative='greater'))


ggsave(filename = './output/KRAS_ID2_present_vs_absent_latency.pdf',width = 5,height = 6,device = cairo_pdf)



# Latency difference between AS_N and EU_N in the EGFR mutant group --------------------------
my_comparisons <- list(c("AS_N",'EU_N'))
MRCAdata %>% 
  filter(acceleration == '1x') %>% 
  left_join(id2data) %>% 
  left_join(
    sherlock_data_full %>%
      filter(Gene=='EGFR',Type=='Mutation_Driver') %>% 
      select(Tumor_Barcode,EGFR=Alteration) %>% 
      mutate(EGFR=if_else(EGFR=='Yes','EGFR mutant','EGFR wildtype'))
  ) %>% 
  left_join(sp_group_data2) %>% 
  filter(Tumor_Barcode %in% hq_samples2, SP_Group_New %in% c('AS_N','EU_N'),!is.na(Latency)) %>% 
  ggplot(aes(SP_Group_New,Latency,fill=SP_Group_New))+
  geom_boxplot(width=0.6,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=3,position = position_jitter(width = 0.15,height = 0),color="black",stroke=0.02)+
  facet_wrap(~EGFR)+
  scale_fill_manual(values = sp_group_color_new)+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 4),axis.title.x = element_blank())+
  labs(y = 'Tumor latency years')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = './output/EGFR_SP_Group_latency.pdf',width = 4,height = 6,device = cairo_pdf)


# SP_group Latency difference----------------------------------------------------------------
library(ggpubr)
my_comparisons <- list( c("AS_N", "EU_N"), c("AS_N", "EU_S"), c("EU_N", "EU_S") )
my_comparisons <- list( c("AS_N", "EU_N"), c("AS_N", "EU_S"), c("EU_N", "EU_S"),c("AS_N", "Others"), c("Others", "EU_S"), c("EU_N", "Others") )


MRCAdata %>% 
  filter(acceleration == "1x") %>% 
  group_by(SP_Group) %>% 
  summarise(Latency =median(Latency,na.rm = T),Age =median(Age,na.rm = T),MRCA_age =median(MRCA_age,na.rm = T))

MRCAdata %>% 
  left_join(sp_group_data) %>% 
  filter(acceleration == "1x") %>% 
  filter(!is.na(Latency)) %>% 
  #filter(SP_Group!="Others") %>% 
  ggplot(aes(SP_Group_New,Latency,fill=SP_Group_New))+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1.5,position = position_jitter(width = 0.15,height = 0.05),color="black")+
  scale_fill_manual(values = sp_group_color_new)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  labs(x = "", y = 'Latency, years before diagnosis')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = './output/latency_spgroup_1x.pdf',width = 4,height = 6,device = cairo_pdf)



MRCAdata %>% 
  left_join(sp_group_data) %>% 
  filter(acceleration == "1x") %>% 
  #filter((SP_Group %in% c('AS_N','EU_N') & acceleration == "1x")|(SP_Group %in% c('EU_S') & acceleration == "7.5x")) %>% 
  filter(!is.na(Latency)) %>% 
  #filter(SP_Group!="Others") %>% 
  ggplot(aes(SP_Group_New,MRCA_age,fill=SP_Group_New))+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1.5,position = position_jitter(width = 0.15,height = 0),color="black")+
  scale_fill_manual(values = sp_group_color_new)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  labs(x = "", y = 'MRCA Age (years)')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)#+
#stat_compare_means(label.y = 100) 

ggsave(filename = './output/MRCA_spgroup.pdf',width = 4,height = 6,device = cairo_pdf)



MRCAdata %>% 
  left_join(sp_group_data) %>% 
  filter(acceleration == "1x") %>% 
  #filter((SP_Group %in% c('AS_N','EU_N') & acceleration == "1x")|(SP_Group %in% c('EU_S') & acceleration == "7.5x")) %>% 
  filter(!is.na(Latency)) %>% 
  #filter(SP_Group!="Others") %>% 
  ggplot(aes(SP_Group_New,Age,fill=SP_Group_New))+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1.5,position = position_jitter(width = 0.15,height = 0),color="black")+
  scale_fill_manual(values = sp_group_color_new)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  labs(x = "", y = 'Age at diagnosis (years)')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = './output/MRCA_spgroup_age.pdf',width = 4,height = 6,device = cairo_pdf)

MRCAdata %>% 
  left_join(sp_group_data) %>% 
  filter(acceleration == "1x") %>% 
  #filter((SP_Group %in% c('AS_N','EU_N') & acceleration == "1x")|(SP_Group %in% c('EU_S') & acceleration == "5x")) %>% 
  filter(!is.na(Latency)) %>% 
  filter(SP_Group!="Others") %>% 
  filter(SP_Group!='N_U') %>% 
  do(tidy(wilcox.test(MRCA_age~SP_Group,data=.)))




# WGD latency -------------------------------------------------------------
library(ggpubr)
my_comparisons <- list( c("WGD", "nWGD"))

MRCAdata %>% 
  left_join(sp_group_data) %>% 
  filter(acceleration == "1x") %>% 
  filter(!is.na(Latency)) %>% 
  filter(SP_Group!="Others") %>% 
  left_join(BBsolution4 %>% select(Tumor_Barcode,WGD_Status)) %>% 
  ggplot(aes(WGD_Status,Latency,fill=WGD_Status))+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1.5,position = position_jitter(width = 0.15,height = 0.05),color="black")+
  facet_wrap(~SP_Group_New)+
  scale_fill_jama()+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(panel.spacing = unit(0.1,"cm"))+
  labs(x = "", y = 'Latency, years before diagnosis')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = './output/latency_WGD_spgroup.pdf',width = 5,height = 6,device = cairo_pdf)



# De-novo group Latency difference----------------------------------------------------------------
load('../RDS/sherlock_data_all.RData')

tmp <- MRCAdata %>% 
  filter(acceleration == "1x") %>% 
  left_join(sherlock_data %>% filter(Type == 'Overall_Feature',Gene=='DN_Group') %>% select(Tumor_Barcode,DN_Group=Alteration)) 

tmp %>% 
  group_by(DN_Group) %>% 
  summarise(Latency =median(Latency,na.rm = T),Age =median(Age,na.rm = T),MRCA_age =median(MRCA_age,na.rm = T))

tmp %>% 
  filter(acceleration == "1x") %>% 
  filter(!is.na(Latency),!is.na(DN_Group)) %>% 
  ggplot(aes(DN_Group,Latency,fill=DN_Group))+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1.5,position = position_jitter(width = 0.15,height = 0.05),color="black")+
  scale_fill_manual(values = ncicolpal)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  labs(x = "", y = 'Latency, years before diagnosis')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = './output/latency_dngroup.pdf',width = 4,height = 5,device = cairo_pdf)

MRCAdata %>% 
  filter(acceleration == "1x") %>% 
  filter(!is.na(Latency)) %>% 
  filter(SP_Group!="Others") %>% 
  ggplot(aes(SP_Group,MRCA_age,fill=SP_Group))+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1.5,position = position_jitter(width = 0.15,height = 0.05),color="black")+
  scale_fill_manual(values = sp_group_color)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  labs(x = "", y = 'MRCA Age')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)








 

# OLD visualization for all groups and target genes -----------------------


# tresult <- tresult %>% filter(Freq>0.03)
## for trended on only major alterations
tmp <- tresult %>% filter(Yes>5|Yes<2|(p.value<0.05 & Freq>0.2))
tresult %>% 
  ggplot(aes((Yes),-log10(p.value),fill=Type))+
  geom_hline(yintercept = -log10(0.05),linetype=2,col='red',size=0.5)+
  geom_point(aes(size=Freq),pch=21,stroke=0.2)+
  scale_fill_manual(values = sherlock_type_colors[sort(unique(tresult$Type))])+
  scale_size_binned()+
  #ggrepel::geom_text_repel(data=tmp,aes(label=Gene_short),max.overlaps = 30)+
  labs(x='Latency, years before diagnosis',y='Wilcox test, -log10(pvalue)')+
  guides(fill = guide_legend(override.aes = list(size=3.5)))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid="XY",ticks = T)+
  panel_border(color = 'black',size = 0.5)+
  coord_cartesian(clip = 'off')

tmp <- tresult %>% filter(abs(diff)>5|(p.value<0.05 & Freq>0.2)|p.value<0.01)
tresult %>% 
  ggplot(aes((diff),-log10(p.value),fill=Type))+
  geom_hline(yintercept = -log10(0.05),linetype=2,col='red',size=0.5)+
  geom_point(aes(size=Freq),pch=21,stroke=0.2)+
  scale_fill_manual(values = sherlock_type_colors[sort(unique(tresult$Type))])+
  scale_size_binned()+
  ggrepel::geom_text_repel(data=tmp,aes(label=Gene_short),max.overlaps = 30,size=3)+
  labs(x='Latency Change (Year)',y='Wilcox test, -log10(pvalue)')+
  guides(fill = guide_legend(override.aes = list(size=3.5)))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid="XY",ticks = T)+
  panel_border(color = 'black',size = 0.5)+
  coord_cartesian(clip = 'off')

ggsave(filename = './output/latency_sherlock_all.pdf',width = 12,height = 8,device = cairo_pdf)


# barplot # 
tmp <- tresult %>% 
  filter(qvalue<0.12) %>%
  select(Gene,p.value,Gene_short,Type) %>% 
  left_join(tdata) %>% 
  mutate(Alteration = if_else(Alteration == "Female", "Yes", if_else(Alteration == "Male","No",Alteration))) %>%
  mutate(Gene_short = if_else(Gene_short == "Gender","Gender_Female",Gene_short)) %>% 
  mutate(Gene_short = paste0(Gene_short,"\n",Type))

tmplevel <- tmp %>% 
  filter(Alteration == "Yes") %>% 
  group_by(Gene_short) %>% 
  summarise(value=median(Latency)) %>% 
  arrange(value) %>% 
  pull(Gene_short)

tmp <- tmp %>% mutate(Gene_short = factor(Gene_short,levels = tmplevel)) 
p1 <- tmp %>% 
  left_join(wgs_groups_info %>% select(Tumor_Barcode,SP_Group)) %>% 
  left_join(sp_group_data) %>% 
  filter(Alteration == 'No') %>% 
  ggplot(aes(Gene_short,Latency,fill=SP_Group_New))+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1,position=position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),color="black")+
  #scale_fill_manual(values = c("Yes" = "#008B45FF","No" = "#cccccc"))+
  scale_fill_manual(values =sp_group_color_new)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = "Y",grid_col = "gray90",ticks = T)+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),legend.position = 'top')+
  labs(x = "", y = 'Latency, years before diagnosis')+
  scale_y_continuous(breaks = pretty_breaks(),limits = c(-1.5,32))+
  #guides(fill="none")+
  panel_border(color = 'black',linetype = 1)

p2 <- tmp %>% 
  filter(Alteration == "Yes") %>% 
  ggplot(aes(Gene_short,"1",fill=-log10(p.value)))+
  geom_tile()+
  scale_fill_viridis_c(option = "D")+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = FALSE,grid_col = "gray90",ticks = F)+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),legend.position = 'top',legend.key.width = unit(1,"cm"))+
  labs(x = "", y = '',fill="Wilcox test, -log10(pvalue)\n")
#guides(fill="none")+
# panel_border(color = 'black',linetype = 1)

plot_grid(p2,p1,axis = 'v',align = 'lr',ncol = 1,rel_heights = c(1.1,4))

ggsave(filename = './output/latency_fdr01.pdf',width = 16,height = 9,device = cairo_pdf)


# targeable alterations
genelist <- c('TP53','EGFR','KRAS','HER2',"BRAF","ALK","RET","ROS1","NTRK","MET","ERBB2","NF1","MAPK2K1","KRAS","FRFRG2","3p21.31",'6q22.31','9p21.3','12p11.22','13q21.33','16q24.1','17p12','18q22.3','22q12.2','5p15.33','8q24.21','10p11.1','12q15','14q21.1','19q13.2','NRAS','RIT1','FGFR1')

tmp <- tresult %>% 
  filter(Type %in% c('Mutation_Driver','SCNA_Focal_Cytoband','Fusion')) %>% 
  filter(Gene_short %in% genelist) %>% 
  select(Gene,p.value,Gene_short,Type) %>% 
  left_join(tdata) %>% 
  mutate(Alteration = if_else(Alteration == "Female", "Yes", if_else(Alteration == "Male","No",Alteration))) %>%
  mutate(Gene_short = if_else(Gene_short == "Gender","Gender_Female",Gene_short)) %>% 
  mutate(Gene_short = paste0(Gene_short,"\n",Type))

exclude_tmp <- tmp %>% filter(Alteration == 'Yes')%>% count(Gene) %>% filter(n<5) %>% pull(Gene)
tmp <- tmp %>% filter(!(Gene %in% exclude_tmp))

tmplevel <- tmp %>% 
  filter(Alteration == "Yes") %>% 
  group_by(Gene_short) %>% 
  summarise(value=median(Latency,na.rm = T)) %>% 
  arrange(value) %>% 
  pull(Gene_short)

tmp <- tmp %>% mutate(Gene_short = factor(Gene_short,levels = tmplevel)) 
p1 <- tmp %>% 
  ggplot(aes(Gene_short,Latency,fill=Alteration))+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1,position=position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),color="black")+
  scale_fill_manual(values = c("Yes" = "#FF7F0EFF","No" = "#cccccc"))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = "Y",grid_col = "gray90",ticks = T)+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),legend.position = 'top')+
  labs(x = "", y = 'Latency, years before diagnosis')+
  scale_y_continuous(breaks = pretty_breaks())+
  #guides(fill="none")+
  panel_border(color = 'black',linetype = 1)

p2 <- tmp %>% 
  filter(Alteration == "Yes") %>% 
  ggplot(aes(Gene_short,"1",fill=-log10(p.value)))+
  geom_tile()+
  scale_fill_viridis_c(option = "D")+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = FALSE,grid_col = "gray90",ticks = F)+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),legend.position = 'top',legend.key.width = unit(1,"cm"))+
  labs(x = "", y = '',fill="Wilcox test, -log10(pvalue)\n")
#guides(fill="none")+
# panel_border(color = 'black',linetype = 1)

plot_grid(p2,p1,axis = 'v',align = 'lr',ncol = 1,rel_heights = c(1.1,4))

ggsave(filename = './output/latency_targetable.pdf',width = 16,height = 9,device = cairo_pdf)



# SBS17b ------------------------------------------------------------------
sherlock_data_full %>% 
  filter(Gene == 'SBS17b') %>% 
  left_join(clinical_data) %>% 
  ggplot(aes(Alteration,Age))+geom_boxplot()

sherlock_data_full %>% 
  filter(Gene == 'SBS17b') %>% 
  left_join(clinical_data) %>% 
  do(tidy(wilcox.test(Age~Alteration,data=.)))

load('../RDS/sherlock_variable.RData')

sherlock_variable %>% 
  filter(name == 'SBS17b') %>%
  left_join(clinical_data) %>% 
  filter(value>0) %>% 
  ggplot(aes(log2(value),Age))+geom_point()+geom_smooth(method = 'lm')

sherlock_variable %>% 
  filter(name == 'SBS17b') %>%
  left_join(clinical_data) %>% 
  filter(value>0) %>% 
  do(tidy(cor.test(.$Age,.$value)))


# SBS17b ------------------------------------------------------------------



# Germline L1 -------------------------------------------------------------
load('Chronological_timing_short.RData')
load('../RDS/sherlock_data_all.RData')
load('../BBsolution_final3_short.RData')

sherlock_data_full

tdata <- MRCAdata %>% 
  filter(acceleration=='1x') %>% 
  left_join(
    sherlock_data_full %>% filter(str_detect(Gene,'All_L1')) %>% select(Tumor_Barcode,Gene,Alteration)
  ) %>% 
  filter(Tumor_Barcode %in% hq_samples2)

tdata %>% group_by(Gene) %>% do(tidy(wilcox.test(Latency~Alteration,data=.)))

tdata %>% 
  ggplot(aes(Gene,Latency,fill=Alteration))+
  geom_point(pch=21,size=2,position=position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),stroke=0.05,col='white')+
  geom_boxplot(width=0.5,pch=21,outlier.shape = NA,alpha=0.3)+
  scale_fill_manual(values = id2color)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  labs(x = "", y = 'Latency, years before diagnosis')+
  panel_border(color = 'black',linetype = 1)

ggsave(file='./output/L1_latency_diff.pdf',width = 6,height = 5,device = cairo_pdf())





