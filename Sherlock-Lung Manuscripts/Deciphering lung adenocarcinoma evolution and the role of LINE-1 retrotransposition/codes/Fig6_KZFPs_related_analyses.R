# ------------------------------------------------------------------------------
# Script: Figure 6 - Tumor Evolution Analysis (LUAD cohort)
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
library(ggpubr)

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
# load Sherlock-lung data -------------------------------------------------
load('./data/BBsolution_final3_short.RData',verbose = T)
load('./data/covdata0.RData',verbose = T)
load('./data/sherlock_data_all.RData',verbose = T)
load('./data/sherlock_variable.RData',verbose = T)
load('./data/sp_group_data.RData',verbose = T)
load('./data/id2data.RData',verbose = T)
load('./data/tedata.RData',verbose = T)
load('./data/annotable_ens90.RData',verbose = T)
protein_coding_genes <- grch38 %>% filter(biotype=='protein_coding') %>% pull(symbol) %>% unique()

# load function -----------------------------------------------------------
load('./data/ZTW_functions.RData')

# load analysis related data set
load('./data/kzfp_genes.RData')
load('./data/RNASeq_Exp.RData',verbose = T)

# --- Subset to LUAD Only ------------------------------------------------------
# Filter to high-quality LUAD samples only
hq_samples2 <- covdata0 %>%
  filter(Histology == 'Adenocarcinoma', Tumor_Barcode %in% hq_samples) %>%
  pull(Tumor_Barcode)
rm(list = c('hq_samples'))


id2data <- id2data %>% filter(Tumor_Barcode %in% hq_samples2)
tedata <- tedata %>% filter(Tumor_Barcode %in% hq_samples2)
rdata1 <- rdata1 %>% filter(Tumor_Barcode %in% hq_samples2, Gene %in% kzfp_genes,Gene %in% protein_coding_genes)


# Fig. 6a: Expression difference of KZFP protein coding genes between IDs Present and Absent in tumors  --------
#Analysis of differentially expressed KZFP protein coding genes between tumors with ID2 signatures and those without. Horizontal dashed lines represent significance thresholds (FDR < 0.05 in orange and FDR < 0.01 in red). The top 20 significant genes are annotated with gene names. 

rdata1 <- rdata1 %>% left_join(id2data %>% select(Tumor_Barcode,ID2,ID2_Present)) %>% filter(!is.na(ID2_Present))
rdata1 <- rdata1 %>% 
  left_join(
    sherlock_data_full %>% filter(Gene == 'TP53', Type=='Mutation_Driver') %>% select(Tumor_Barcode,TP53=Alteration)
  ) %>% 
  left_join(
    sherlock_data_full %>% filter(Gene == 'KRAS', Type=='Mutation_Driver') %>% select(Tumor_Barcode,KRAS=Alteration)
  ) %>% 
  left_join(
    sherlock_data_full %>% filter(Gene == 'EGFR', Type=='Mutation_Driver') %>% select(Tumor_Barcode,EGFR=Alteration)
  )

rdata1 <- rdata1 %>% left_join(covdata0 %>% select(Tumor_Barcode,Smoking,Gender,Age))

fcdata <- rdata1 %>% 
  filter(RNAseq_Type == 'Tumor') %>% 
  group_by(Gene,ID2_Present) %>% 
  summarise(Exp=median(Exp)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = ID2_Present,values_from = Exp) %>% 
  mutate(FC=log2(2^(Present)/(2^`Absent`)))

## Wilcox test ##
# tresult <- rdata1 %>%
#   filter(RNAseq_Type == 'Tumor') %>%
#   group_by(Gene) %>%
#   do(tidy(wilcox.test(Exp~ID2_Present, data=.))) %>%
#   ungroup() %>%
#   arrange(p.value) %>%
#   select(-method,-alternative) %>%
#   mutate(FDR=p.adjust(p.value,method="BH")) %>%
#   left_join(fcdata) %>%
#   left_join(kzfp %>% select(Gene=Label,chrom))

## logistic regression using glm
tresult <- rdata1 %>% 
  filter(RNAseq_Type == 'Tumor') %>% 
  mutate(ID2_Present=as.factor(ID2_Present)) %>%
  group_by(Gene) %>% 
  do(tidy(glm(ID2_Present ~ Exp + Smoking + Age + Gender + TP53 + KRAS + EGFR, family='binomial',data=.,),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term=='Exp') %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  mutate(FDR=p.adjust(p.value,method="BH")) %>% 
  left_join(fcdata) %>% 
  left_join(kzfp %>% select(Gene=Label,chrom))


tresult %>% 
  ggplot(aes(FC,-log10(FDR),fill=chrom))+
  geom_vline(xintercept = c(-0.5,0,0.5),linetype=2,col='#cccccc')+
  geom_hline(yintercept = -log10(0.01),linetype=2,col='red',size=0.5)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col='orange',size=0.5)+
  geom_point(pch=21,stroke=0.5,col='black',size=3)+
  scale_x_continuous(breaks = pretty_breaks(),limits = c(-2.2,2.2))+
  scale_y_continuous(breaks = pretty_breaks(n=7))+
  scale_fill_manual(values = rev(ncicolpal))+
  ggrepel::geom_text_repel(data=tresult %>% slice(1:20), aes(label=Gene),force = 10,size=3.5)+
  #ggrepel::geom_text_repel(data=tresult %>%filter(Gene %in% targets), aes(label=Gene),force = 10,size=4,col='red')+
  labs(x='Fold Change (log2)',y='-log10(FDR)')+
  guides(fill = guide_legend(override.aes = list(size=4)))+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid=FALSE,ticks = T)+
  panel_border(color = 'black',size = 0.5)

#ggsave(filename = 'KZFPs_DEGs_tumor_ID2_present_absent.pdf',width = 7,height = 6,device = cairo_pdf)
ggsave(filename = './output/KZFPs_DEGs_tumor_ID2_present_absent_multivariable.pdf',width = 7,height = 6,device = cairo_pdf)




# Fig. 6b: ZNF695 expression difference among ID2 present and absent group-----------------------------------------------------------------
#Box plots illustrate the differential expression of ZNF695 among normal tissue or blood, tumors without ID2 signatures, and tumors with ID2 signatures. 

my_comparisons <- list(c("Normal", "Tumor-ID2 absent"),c("Normal", "Tumor-ID2 present"),c("Tumor-ID2 absent", "Tumor-ID2 present"))

rdata1 %>% 
  mutate(Group=if_else(RNAseq_Type == 'Normal','Normal',if_else(ID2_Present == 'Present','Tumor-ID2 present','Tumor-ID2 absent'))) %>% 
  filter(Gene=='ZNF695') %>% 
  ggplot(aes(Group,Exp))+
  geom_point(aes(fill=Group),pch=21,size=2,position = position_jitter(width = 0.15,height = 0),color="black")+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25,fill=NA)+
  scale_fill_manual(values = pal_jama()(3)[3:1])+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(axis.text.x = element_text(angle = 30,hjust = 0.5,vjust = 0.5))+
  labs(x = "", y = 'ZNF695 RNA-Seq expression log2(CPM)')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = './output/KZFPs_DEGs_tumor_ID2_present_absent_example.pdf',width = 4,height = 7,device = cairo_pdf)

# adjusted for the cofounders e.g, smoking status, TP53, KRAS mutant status
tmpresult <- rdata1 %>% 
  mutate(Smoking=if_else(Smoking == 'Unknown',NA_character_,Smoking)) %>% 
  filter(Gene=='ZNF695',RNAseq_Type == 'Tumor') %>% 
  mutate(ID2_Present=as.factor(ID2_Present)) %>%
  do(tidy(glm(ID2_Present ~ Exp + Smoking + Age + Gender + TP53 + KRAS + EGFR, family='binomial',data=.,),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term!='(Intercept)') %>% 
  ungroup() %>% 
  arrange(desc(estimate)) %>% 
  #mutate(term=str_remove(term,"EGFR")) %>% 
  mutate(term=str_remove(term,"Smoker")) %>% 
  mutate(term=str_remove(term,"Yes")) %>% 
  mutate(label=paste0('P = ',scientific_format(digits = 3)(p.value)))


tmpresult %>% 
  mutate(p.value2 = if_else(p.value>0.05,NA,-log10(p.value))) %>% 
  #mutate(conf.high = if_else(conf.high>10,NA,conf.high)) %>% 
  #mutate(term=factor(term,levels=tmplevels,labels=tmplables)) %>% 
  mutate(term=fct_inorder(term)) %>% 
  ggplot(aes(estimate, term, xmin = conf.low, xmax = conf.high, height = 0,color=p.value2)) +
  geom_vline(xintercept = 0, lty = 5) +
  geom_errorbarh(height = .3)+
  geom_point(pch=19,size=4) +
  #scale_color_manual(values = c("FALSE"='black',"TRUE"='red'))+
  scale_color_gradient(low = "#EE5250FF",high = "#B71B1BFF",na.value = 'gray20')+
  ggrepel::geom_text_repel(aes(label=label),size=3.5,nudge_y = 0.4)+
  #facet_grid(~SP_Group_New)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid_col = 'gray85')+
  theme(panel.spacing = unit(0.2,"cm"))+
  labs(x = "Logistic regression coefficient", y = NULL)+
  guides(color="none")+
  scale_x_continuous(breaks = pretty_breaks())+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = './output/KZFPs_DEGs_tumor_ID2_present_absent_example2.pdf',width = 6,height = 4,device = cairo_pdf)


# Fig. 6c: Correlation with ID2 ----------------------------------------------
#Pearson correlations between KZFP protein-coding gene expression and indels attributed to mutational signature ID2. Horizontal dashed lines indicate significance thresholds (FDR < 0.05 in orange and FDR < 0.01 in red). The top 20 significant genes are annotated with gene names.

# pearson correlation
# tresult <- rdata1 %>%
#   filter(RNAseq_Type == 'Tumor') %>%
#   left_join(id2data %>% select(Tumor_Barcode,ID2)) %>%
#   filter(ID2>0) %>%
#   group_by(Gene) %>%
#   do(tidy(cor.test(.$Exp,log2(.$ID2+1)))) %>%
#   ungroup() %>%
#   arrange(p.value) %>%
#   select(-method,-alternative) %>%
#   mutate(FDR=p.adjust(p.value,method="BH")) %>%
#   left_join(kzfp %>% select(Gene=Label,chrom))

# adjust for the covariables
tresult <- rdata1 %>% 
  filter(RNAseq_Type == 'Tumor') %>% 
  filter(ID2>0) %>% 
  group_by(Gene) %>% 
  do(tidy(lm(log2(ID2) ~ Exp + Smoking + Age + Gender + TP53 + KRAS + EGFR, data=.,),conf.int = TRUE, conf.level = 0.95)) %>% 
  #do(tidy(cor.test(.$Exp,log2(.$ID2+1)))) %>% 
  ungroup() %>% 
  filter(term=='Exp') %>% 
  arrange(p.value) %>% 
  mutate(FDR=p.adjust(p.value,method="BH")) %>% 
  left_join(kzfp %>% select(Gene=Label,chrom))


tresult %>% 
  ggplot(aes(estimate,-log10(FDR),fill=chrom))+
  geom_vline(xintercept = 0,linetype=2)+
  geom_hline(yintercept = -log10(0.01),linetype=2,col='red',size=0.5)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col='orange',size=0.5)+
  geom_point(pch=21,stroke=0.5,col='black',size=3)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  scale_fill_manual(values = rev(ncicolpal))+
  ggrepel::geom_text_repel(data=tresult %>% slice(1:20), aes(label=Gene),force = 10,size=3.5)+
  #ggrepel::geom_text_repel(data=tresult %>%filter(Gene %in% targets), aes(label=Gene),force = 10,size=4,col='red')+
  labs(x="Linear regression coefficient",y='-log10(FDR)')+
  guides(fill = guide_legend(override.aes = list(size=4)))+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid=FALSE,ticks = T)+
  panel_border(color = 'black',size = 0.5)


#ggsave(filename = 'KZFPs_tumor_ID2_correlation.pdf',width = 7,height = 6,device = cairo_pdf)
ggsave(filename = './output/KZFPs_tumor_ID2_correlation_mutivariable.pdf',width = 7,height = 6,device = cairo_pdf)

# Fig. 6d: individual correction ------------------------------------------
#Correlation between ZNF695 expression and deletions attributed to mutational signature ID2. Pearson correlation coefficients and corresponding p-values are displayed above the plot.

rdata1 %>% 
  filter(RNAseq_Type == 'Tumor') %>% 
  left_join(id2data %>% select(Tumor_Barcode,ID2)) %>% 
  filter(Gene=='ZNF695') %>% 
  #filter(Gene=='ZNF658B') %>% 
  filter(ID2>0) %>% 
  ggplot(aes(Exp,log2(ID2)))+
  geom_point(pch=21,stroke=0.5,fill=ncicolpal[2],size=2.5)+
  geom_smooth(method = 'lm')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'Y',strip_text_face = 'bold')+
  labs(y='Number of ID2 deletions (log2)',x='ZNF695 RNA-Seq expression log2(CPM)')+
  theme(panel.spacing = unit(0.2,'lines'),strip.text.x = element_text(face = 'plain',hjust = 0.5))+
  coord_cartesian(clip = 'off')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  guides(fill="none")

ggsave(filename = './output/KZFPs_tumor_ID2_correlation_example.pdf',width = 6,height = 5,device = cairo_pdf)

# adjust for the covariables
tmpresult <- rdata1 %>% 
  filter(Gene=='ZNF695',RNAseq_Type == 'Tumor') %>% 
  left_join(id2data %>% select(Tumor_Barcode,ID2)) %>% 
  filter(ID2>0) %>% 
  do(tidy(lm(log2(ID2) ~ Exp + Smoking + Age + Gender + TP53 + KRAS + EGFR, data=.,),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term!='(Intercept)') %>% 
  ungroup() %>% 
  arrange(desc(estimate)) %>% 
  #mutate(term=str_remove(term,"EGFR")) %>% 
  mutate(term=str_remove(term,"Smoker")) %>% 
  mutate(term=str_remove(term,"Yes")) %>% 
  mutate(label=paste0('P = ',scientific_format(digits = 3)(p.value)))


tmpresult %>% 
  mutate(p.value2 = if_else(p.value>0.05,NA,-log10(p.value))) %>% 
  #mutate(term=factor(term,levels=tmplevels,labels=tmplables)) %>% 
  mutate(term=fct_inorder(term)) %>% 
  ggplot(aes(estimate, term, xmin = conf.low, xmax = conf.high, height = 0,color=p.value2)) +
  geom_vline(xintercept = 0, lty = 5) +
  geom_errorbarh(height = .3)+
  geom_point(pch=19,size=4) +
  #scale_color_manual(values = c("FALSE"='black',"TRUE"='red'))+
  scale_color_gradient(low = "#EE5250FF",high = "#B71B1BFF",na.value = 'gray20')+
  ggrepel::geom_text_repel(aes(label=label),size=3.5,nudge_y = 0.4)+
  #facet_grid(~SP_Group_New)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid_col = 'gray85')+
  theme(panel.spacing = unit(0.2,"cm"))+
  labs(x = "Linear regression coefficient", y = NULL)+
  guides(color="none")+
  scale_x_continuous(breaks = pretty_breaks())+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = './output/KZFPs_tumor_ID2_correlation_example2.pdf',width = 6,height = 4,device = cairo_pdf)



# Fig. 6e ----------------------------------------------------------------
#Correlation between ZNF695 RNA-Seq expression and the median of DNA methylation levels across the genome locations of the CpG island on chr22q12.1. Pearson correlation coefficients and corresponding p-values are displayed. 
# codes generated Fig. 6e can be found in the Fig5_LINE1_Methylation.R 



# Fig. 6f -----------------------------------------------------------------
#Differentially expressed ZNF695 target genes identified between tumors with and without mutational signature ID2.

ZNF695_targets <- unique(c("UBR5", "PPP2R1A", "PPP2R1B", "PSMB3", "PPP4C", "HECTD3", "USP9X", "BRE", "HECTD1", "TIMM50", "PPP2CB", "PSMC3", "PSMD3", "PSMD6", "PSMA6", "PSMA3", "PDCD5", "PSMD11", "DDI2", "PSMC1", "PSMC2", "PSMD7", "PSMC5", "DNAJA1", "PSMA2", "PSMC6","MEOX2", "TRIM28", "HECTD1", "H2BC21","MER4D", "MER4D1", "HAL1-2a_MD"))
length(ZNF695_targets)

load('./data/tresult_wilcox_ID2.RData',verbose = T)

tresult <- tresult %>% 
  filter(Gene %in% ZNF695_targets) %>% 
  mutate(FDR=p.adjust(p.value,method='BH'))


tresult %>% 
  ggplot(aes(FC,-log10(FDR)))+
  geom_vline(xintercept = c(-1,0,1),linetype=2,col='#cccccc')+
  geom_hline(yintercept = -log10(0.01),linetype=2,col='red',size=0.5)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col='orange',size=0.5)+
  geom_point(pch=21,col='black',size=3,stroke = 0.1,fill=ncicolpal[1])+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_y_continuous(breaks = pretty_breaks(n=7))+
  #scale_fill_manual(values = c(ncicolpal[c(2,3,14)],'#cccccc'),na.value = '#cccccc')+
  ggrepel::geom_text_repel(data=tresult %>% arrange(FDR) %>% slice(1:50), aes(label=Gene),force = 10,size=3)+
  ggrepel::geom_text_repel(data=tresult %>% filter(Gene=='ZNF695'), aes(label=Gene),force = 10,size=4,col='red')+
  labs(x='Fold Change (log2)',y='-log10(FDR)')+
  guides(fill = guide_legend(override.aes = list(size=5)))+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid=FALSE,ticks = T)+
  theme(legend.position = 'top')+
  panel_border(color = 'black',size = 0.5)

ggsave(filename = './output/tresult_wilcox_ID2_ZNF695_targets.pdf',width = 6,height = 6,device = cairo_pdf)

geneset1 <- tresult$Gene

# general target genes for kzfps
kzfps_targets <- readxl::read_xlsx('./data/ID2 Paper - KRAB-ZFP Target Genes Master List.xlsx',skip = 1,col_names = T) %>% pull(`Target Genes`)

length(kzfps_targets)

load('./data/tresult_wilcox_ID2.RData',verbose = T)

tresult <- tresult %>% 
  filter(Gene %in% kzfps_targets) %>% 
  mutate(FDR=p.adjust(p.value,method='BH'))


tresult %>% 
  ggplot(aes(FC,-log10(FDR)))+
  geom_vline(xintercept = c(-1,0,1),linetype=2,col='#cccccc')+
  geom_hline(yintercept = -log10(0.01),linetype=2,col='red',size=0.5)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col='orange',size=0.5)+
  geom_point(pch=21,col='black',size=3,stroke = 0.1,fill=ncicolpal[1])+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_y_continuous(breaks = pretty_breaks(n=7))+
  #scale_fill_manual(values = c(ncicolpal[c(2,3,14)],'#cccccc'),na.value = '#cccccc')+
  ggrepel::geom_text_repel(data=tresult %>% arrange(FDR) %>% slice(1:50), aes(label=Gene),force = 10,size=3)+
  ggrepel::geom_text_repel(data=tresult %>% filter(Gene=='ZNF695'), aes(label=Gene),force = 10,size=4,col='red')+
  labs(x='Fold Change (log2)',y='-log10(FDR)')+
  guides(fill = guide_legend(override.aes = list(size=5)))+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid=FALSE,ticks = T)+
  theme(legend.position = 'top')+
  panel_border(color = 'black',size = 0.5)

ggsave(filename = './output/tresult_wilcox_ID2_KZFPs_targets.pdf',width = 6,height = 6,device = cairo_pdf)


