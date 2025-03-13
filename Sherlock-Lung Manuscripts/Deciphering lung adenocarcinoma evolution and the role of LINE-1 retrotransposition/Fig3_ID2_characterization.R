set_wd()
libztw()
pdfhr2()
myggstyle()

# load Sherlock-lung data -------------------------------------------------
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/BBsolution_final3_short.RData')
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/covdata0.RData')
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/Clinical/clinical_data.RData')
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/RDS/sherlock_data_all.RData')
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/RDS/sherlock_variable.RData')

load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/RDS/sherlock_variable.RData')

load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/sp_group_data.RData')

load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/ID2_TE/id2data.RData')
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/ID2_TE/tedata.RData')

load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/RNASeq/RNASeq_Exp.RData',verbose = T)
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/Signature_ludmil3/Signature_Lumidl_CN.RData',verbose = T)
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/Signature_ludmil3/Signature_Lumidl_SV.RData',verbose = T)
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/Signature_ludmil3/sherlock_profiles.RData',verbose = T)


# load function -----------------------------------------------------------
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/ZTW_functions.RData')
#source('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/ZTW_functions.R')
source('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/Sherlock_functions.R')


# load analysis related data set
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/MutationTimeR_HQ_LUAD/Chronological_timing_short.RData',verbose = T)
tmp <- colnames(MRCAdata)
MRCAdata <- MRCAdata %>% select(-SP_Group) %>% left_join(sp_group_data2) %>% select(-SP_Group_New) %>% select(one_of(tmp))

load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/Survival/suvdata.RData',verbose = T)
## analysis limited to luad only
hq_samples2 <- covdata0 %>% filter(Histology == 'Adenocarcinoma', Tumor_Barcode %in% hq_samples) %>% pull(Tumor_Barcode)
rm(list=c('hq_samples'))





# Fig. 3a-b -----------------------------------------------------------------


# ID2 with all proliferation markers
## new limitation to hqsamples for luad only
load('../covdata0.RData')
hq_samples2 <- covdata0 %>% filter(Histology == 'Adenocarcinoma', Tumor_Barcode %in% hq_samples) %>% pull(Tumor_Barcode)
lq_samples2 <- covdata0 %>% filter(Histology == 'Adenocarcinoma', !(Tumor_Barcode %in% hq_samples)) %>% pull(Tumor_Barcode)

rm(list=c('hq_samples'))


genelist <- c('MKI67','TOP2A','MYBL2','BUB1','PLK1','CCNE1','CCNB1','BUB1','FOXM1')
load('../RDS/sherlock_data_all.RData')
tdata <- rdata1 %>% filter(Gene %in% genelist)
load('../ID2_TE/id2data.RData')

tdata <- id2data %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  left_join(tdata) %>% 
  filter(!is.na(Exp))


fcdata <- tdata %>% 
  group_by(RNAseq_Type,Gene,ID2_Present) %>% 
  summarise(mvalue=median(Exp)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = ID2_Present,values_from = mvalue) %>% 
  mutate(FC=Present/Absent)


tdata %>% 
  group_by(RNAseq_Type,Gene) %>% 
  do(tidy(wilcox.test(Exp~ID2_Present, data=.))) %>% 
  left_join(fcdata) %>% 
  group_by(RNAseq_Type) %>% 
  mutate(FDR=p.adjust(p.value)) %>% 
  ggplot(aes(log2(FC),-log10(FDR),fill=RNAseq_Type))+
  geom_hline(yintercept = -log10(0.05),linetype=2,col='red',size=0.5)+
  # geom_vline(xintercept = c(-5,5),linetype=2,col='blue',size=0.5)+
  geom_point(pch=21,stroke=0.4,size=3)+
  scale_fill_nejm()+
  scale_x_continuous(breaks = pretty_breaks(n=6))+
  scale_y_continuous(breaks = pretty_breaks(n = 6))+
  ggrepel::geom_text_repel(aes(label=Gene),size=3,max.overlaps = 20,force = 20)+
  labs(x='log2(Fold change)',y='-log10(FDR)',fill='RNA-Seq data')+
  guides(fill = guide_legend(override.aes = list(size=3.5)))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid='XY',ticks = T)+
  theme(legend.position = 'top')+
  panel_border(color = 'black',size = 0.5)+
  coord_cartesian(clip = 'off')

ggsave(filename = 'ID2_Proliferation_Markers.pdf',width = 5,height = 5,device = cairo_pdf())


tdata %>% 
  left_join(covdata0) %>% 
  group_by(Gene,RNAseq_Type) %>% 
  do(tidy(lm(Exp~ID2_Present + Assigned_Population + Gender + Smoking + Tumor_Purity + Age, data=.))) %>% 
  filter(term=='ID2_PresentPresent') %>% 
  arrange((p.value))

# association

tdata2 <- tdata %>% 
  filter(ID2>0, RNAseq_Type=='Tumor') %>% 
  select(Tumor_Barcode,ID2,Gene,Exp) %>% 
  pivot_wider(names_from = Gene,values_from = Exp) %>% 
  mutate(ID2=log2(ID2))

library(ggcorrplot)

corr <- round(cor(tdata2[,-1]), 1)
p.mat <- cor_pmat(tdata2[,-1])

ggcorrplot(corr)
ggcorrplot(corr,method = "circle",insig = "blank",
           hc.order = TRUE,
           type = "lower",
           p.mat = p.mat,
           lab = TRUE)

ggsave(filename = 'ID2_Proliferation_Markers2.pdf',width = 5,height = 5,device = cairo_pdf())



# box plot
tdata %>% 
  ggplot(aes(ID2_Present,Exp))+
  geom_point(aes(fill=ID2_Present),pch=21,size=2.5,position = position_jitter(width = 0.2,height = 0))+
  geom_boxplot(width=0.7,outlier.shape = NA,fill=NA,size=1)+
  facet_grid(RNAseq_Type~Gene,scales ='free') +
  scale_fill_manual(values = id2color)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'Y')+
  guides(fill='none')+#fill = guide_legend(override.aes=list(shape=21)),)+
  labs(x='ID2 Mutational Signature',y='RNA-Seq gene expression log2(CPM)')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),panel.spacing = unit(0.1,'cm'))+
  panel_border(color = 'black',linetype = 1)

ggsave('ID2_Proliferation_RNASeq.pdf',width = 10,height = 9,device = cairo_pdf)

tdata %>% 
  group_by(RNAseq_Type,Gene) %>% 
  do(tidy(wilcox.test(Exp~ID2_Present,data=.))) 

# scatter plot
tmp <- tdata %>% 
  filter(ID2>0,RNAseq_Type=='Tumor') %>% 
  group_by(Gene) %>% 
  do(tidy(cor.test(log2(.$ID2),.$Exp))) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  mutate(label=paste0(Gene,'\nR=',round(estimate,2),'; P=',scientific(p.value))) %>% 
  select(Gene,label) %>% 
  mutate(label=fct_inorder(label))

tdata %>% 
  filter(ID2>0,RNAseq_Type=='Tumor') %>% 
  left_join(tmp) %>% 
  ggplot(aes(log2(ID2),Exp))+
  geom_point(pch=21,stroke=0.5,fill=id2color[2],size=2)+
  facet_wrap(~label,nrow = 2)+
  geom_smooth(method = 'lm')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'Y',strip_text_face = 'bold')+
  labs(x='Number of ID2 deletions (log2)',y='RNA-Seq expression log2(CPM)')+
  theme(panel.spacing = unit(0.2,'lines'),strip.text.x = element_text(face = 'plain',hjust = 0.5))+
  coord_cartesian(clip = 'off')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  guides(fill="none")

ggsave('ID2_Proliferation_RNASeq_tumor.pdf',width = 8,height = 5.5,device = cairo_pdf)


tmp <- tdata %>% 
  filter(ID2>0,RNAseq_Type=='Normal') %>% 
  group_by(Gene) %>% 
  do(tidy(cor.test(log2(.$ID2),.$Exp))) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  mutate(label=paste0(Gene,'\nR=',round(estimate,2),'; P=',round(p.value,2))) %>% 
  select(Gene,label) %>% 
  mutate(label=fct_inorder(label))

tdata %>% 
  filter(ID2>0,RNAseq_Type=='Normal') %>% 
  left_join(tmp) %>% 
  ggplot(aes(log2(ID2),Exp))+
  geom_point(pch=21,stroke=0.5,fill=id2color[2],size=2)+
  facet_wrap(~label,nrow = 2)+
  geom_smooth(method = 'lm')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'Y',strip_text_face = 'bold')+
  labs(x='Number of ID2 deletions (log2)',y='RNA-Seq expression log2(CPM)')+
  theme(panel.spacing = unit(0.2,'lines'),strip.text.x = element_text(face = 'plain',hjust = 0.5))+
  coord_cartesian(clip = 'off')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  guides(fill="none")

ggsave('ID2_Proliferation_RNASeq_normal.pdf',width = 8,height = 5.5,device = cairo_pdf)


# Fig. 3c-------------------------------------------------------------------
# ID2 survival 
load('../covdata0.RData')
source('../Survival/Survival_function.R')

suvdata <- suvdata %>% left_join(wgs_groups_info %>% select(Tumor_Barcode,SP_Group))
suvdata <- suvdata %>% 
  left_join(
    sherlock_data_full %>% filter(Gene=='TP53',Type=='Mutation_Driver') %>% select(Tumor_Barcode,TP53_Status=Alteration) %>% mutate(TP53_Status = factor(TP53_Status, levels=c('No','Yes')))
  ) 

suvdata_tmp <- id2data %>% 
  rename(Key = ID2_Present) %>% 
  left_join(suvdata) %>% 
  mutate(Key = factor(Key, levels = c('Absent','Present'))) %>% 
  filter(Tumor_Barcode %in% hq_samples2) 
#mutate(Survival_Month = if_else(Survival_Month > 60, NA, Survival_Month))

SurvZTWms(suvdata_input = suvdata_tmp,plot = TRUE,keyname='Mutational Signature ID2: ',gcolors = as.character(id2color),pvalsize = 4,width = 5,height = 5,filename = 'ID2-survival.pdf')
SurvZTWms_tp53(suvdata_input = suvdata_tmp,plot = TRUE,keyname='Mutational Signature ID2: ',gcolors = as.character(id2color),pvalsize = 4,width = 5,height = 4,filename = 'ID2-survival-tp53.pdf')



# Fig. 3d ------------------------------------------------------------------
#load('hq_samples2.RData')
clinical_data %>% select(Tumor_Barcode,Metastasis_old)

tdata <- id2data %>% 
  left_join(clinical_data %>% select(Tumor_Barcode,Metastasis=Metastasis_after_diagnosis) %>% mutate(Metastasis = factor(Metastasis, levels = c('No','Yes')))
  ) %>% 
  filter(!is.na(Metastasis)) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  left_join(sp_group_data2) %>% 
  left_join(
    sherlock_data_full %>% filter(Gene=='TP53',Type=='Mutation_Driver') %>% select(Tumor_Barcode,TP53_Status=Alteration) %>% mutate(TP53_Status = factor(TP53_Status, levels=c('No','Yes'),labels=c('TP53 wildtype','TP53 mutant')))
  )

tdata %>% do(tidy(fisher.test(.$ID2_Present,.$Metastasis)))
tdata %>% group_by(TP53_Status) %>%  do(tidy(fisher.test(.$ID2_Present,.$Metastasis)))
tdata %>% group_by(SP_Group_New) %>% do(tidy(fisher.test(.$ID2_Present,.$Metastasis)))
tdata %>% filter(ID2>0) %>%  do(tidy(wilcox.test(ID2~Metastasis,data=.)))

my_comparisons = list(c('Yes','No'))
myggstyle()
tdata %>% 
  filter(ID2>0) %>% 
  ggplot(aes(Metastasis,log2(ID2)))+
  geom_point(aes(fill=Metastasis),pch=21,size=3,position = position_jitter(width = 0.2,height = 0),stroke=0.2)+
  geom_boxplot(width=0.6,outlier.shape = NA,fill=NA,col=ncicolpal[1],size=0.6)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),panel.spacing = unit(0.2, "lines"),strip.background = element_rect(fill = 'gray95'),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,4,4))+
  labs(x = "Tumor Metastasis", y = 'ID2 deletions (log2)')+
  scale_fill_manual(values = id2color)+
  guides(fill="none")+
  #facet_wrap(~TP53_Status,nrow = 1)+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(file='ID2_metastasis2.pdf',width = 3,height = 6,device = cairo_pdf)
ggsave(file='ID2_metastasis.pdf',width = 2,height = 6,device = cairo_pdf)


barplot_fisher(tdata,'ID2_Present','Metastasis',var1lab = 'Mutational signature ID2',var2lab = 'Tumor metastasis',filename = 'ID2_metastasis_enrchment.pdf')

barplot_fisher(tdata %>% filter(TP53_Status=='TP53 wildtype'),'ID2_Present','Metastasis',var1lab = 'Mutational signature ID2',var2lab = 'Tumor metastasis',filename = 'ID2_metastasis_enrchment_tp53_wildtype.pdf')

barplot_fisher(tdata %>% filter(TP53_Status=='TP53 mutant'),'ID2_Present','Metastasis',var1lab = 'Mutational signature ID2',var2lab = 'Tumor metastasis',filename = 'ID2_metastasis_enrchment_tp53_mutant.pdf')

tdata %>% 
  left_join(covdata0) %>%
  filter(Smoking=='Non-Smoker') %>% 
  do(tidy(glm(Metastasis ~ ID2_Present + TP53_Status + Tumor_Purity + Age + Gender, family='binomial',data=.)))



# Fig. 3e -----------------------------------------------------------------
# Association with ID2 tumors ---------------------------------------------
load('id2data.RData')
load('../RDS/sherlock_data_all.RData')
load('../BBsolution_final3_short.RData')
load('../covdata0.RData')
drglist <- readRDS('../../../Collaborators/Nuria/Update2/drivers_intogene.RDS') %>% pull(symbol)
drglist <- unique(c(drglist,c('RET','ALK','AXL','NRG1','MET','FGFR2','ROS1','RB1','ERBB4','EGFR')))

drglist <- c(drglist,paste0(drglist,'-Amp'),paste0(drglist,'-Del'))

features=c('WGD_Status','Distal_Proximal','A3A_or_A3B','Kataegis','HRDetect_Status','Forte','All_L1')

tdata <- sherlock_data_full %>% 
  mutate(Type=if_else(Gene=='All_L1','Overall_Feature',Type)) %>% 
  filter(Tumor_Barcode %in% hq_samples2,Gene %in% c(features,drglist),!(Type %in% c('Gene_Alterations','Multiple_Mutations','Pathogenic_Germline'))) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  mutate(Alteration = if_else(Alteration == 'WGD','Yes',if_else(Alteration == 'nWGD','No',Alteration))) %>% 
  mutate(Alteration = if_else(Alteration == 'A3A','Yes',if_else(Alteration == 'A3B','No',Alteration))) %>% 
  left_join(id2data) %>% 
  mutate(Type=factor(Type,levels=c('Overall_Feature','Mutation_Driver','SCNA_Focal_Gene','Fusion')))

tresult <- tdata %>% 
  left_join(covdata0) %>% 
  group_by(Type,Gene) %>% 
  mutate(ID2_Present= as.integer(as.factor(ID2_Present))-1) %>% 
  do(tresult = safely(glm)(ID2_Present ~ Alteration + Assigned_Population + Gender + Age + Smoking + Tumor_Purity, family = "binomial", data=. )) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']]))) %>% 
  select(Type,Gene,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term!='(Intercept)', term=='AlterationYes') %>% 
  filter(abs(estimate)<10) %>% 
  mutate(FDR=p.adjust(p.value))


tresult %>% 
  filter(!is.na(Type)) %>% 
  ggplot(aes((estimate),-log10(FDR)))+
  geom_hline(yintercept = -log10(0.01),linetype=2,col='red',size=0.5)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col='orange',size=0.5)+
  geom_point(aes(fill=Type),pch=21,stroke=0.1,col='black',size=4)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_y_continuous(breaks = pretty_breaks())+
  ggrepel::geom_text_repel(data=tresult %>%filter(FDR<0.05), aes(label=Gene),force = 10,size=4)+
  scale_fill_manual(values = ncicolpal)+
  labs(x='Regression coefficient for tumor latency',y='-log10(FDR)',fill='Alterations or features')+
  guides(fill = guide_legend(override.aes = list(size=5)))+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid='XY',ticks = T)+
  panel_border(color = 'black',size = 0.5)

ggsave(filename = 'ID2_hq_tumor_genomic_association.pdf',width = 7.5,height = 5,device = cairo_pdf)





# Fig. 3f -----------------------------------------------------------------

# ID2 and hypoxia 
library(ggdist)
tresult <- sherlock_variable %>% 
  left_join(id2data) %>% 
  group_by(name) %>% 
  do(tidy(wilcox.test(value~ID2_Present,data=.)))


tresult %>% 
  filter(str_detect(name,'hypoxia_score'))


tdata <- sherlock_variable %>% 
  filter(str_detect(name,'hypoxia_score')) %>% 
  filter(name=='hypoxia_score_Combined_GeneSet') %>% 
  left_join(id2data) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  left_join(sp_group_data2) %>% 
  ungroup()

tdata <- tdata %>% mutate(SP_Group_New='ALL') %>% bind_rows(tdata) #%>% filter(SP_Group_New!='Others')


tdata %>% 
  ggplot(aes(SP_Group_New,value,fill=ID2_Present))+
  stat_eye(side = "left",position = position_dodge(width = 0.2),col='black')+
  scale_fill_manual(values = id2color)+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 18,axis_title_just = 'm',axis_title_size = 22)+
  labs(x = "", y = 'Hypoxia scorer',fill='ID2 mutational signature')+
  #guides(fill="none")+
  theme(legend.position = 'top')+
  panel_border(color = 'black',linetype = 1)

ggsave(file='ID2_hypoxia_score.pdf',width = 9,height = 9,device = cairo_pdf)



tdata %>% group_by(SP_Group_New) %>% do(tidy(wilcox.test(value~ID2_Present,data=.)))





# Others ------------------------------------------------------------------


# ATR -  CHECK1  vs ID2 -------------------------------------------------------------------------
load('../ID2_TE/id2data.RData')
tdata <- rdata1 %>% filter(Gene %in% c('CHEK1','CHEK2'))

tdata <- id2data %>% 
  left_join(tdata) %>% 
  filter(!is.na(Exp),!is.na(ID2)) %>% 
  left_join(wgs_groups_info)

tdata %>% 
  filter(RNAseq_Type=='Tumor') %>% 
  ggplot(aes(ID2_Present,Exp))+
  geom_point(aes(fill=ID2_Present),pch=21,size=2.5,position = position_jitter(width = 0.15,height = 0))+
  geom_boxplot(width=0.7,outlier.shape = NA,fill=NA)+
  facet_grid(RNAseq_Type~Gene,scales = 'free') +
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'Y')+
  guides(fill = guide_legend(override.aes=list(shape=21)))+
  labs(x='ID2 mutational signature',y='RNA-Seq Expression log2(CPM)')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),panel.spacing = unit(0.25, "cm"))+
  panel_border(color = 'black',linetype = 1)+
  scale_fill_manual(values = c("#01665e","#ff7f00"))+
  guides(fill="none")

ggsave('ID2_Checkpoint_RNASeq.pdf',width = 5,height = 8,device = cairo_pdf)

tdata %>% 
  group_by(RNAseq_Type,Gene) %>% 
  do(tidy(wilcox.test(Exp~ID2_Present,data=.))) %>% 
  arrange(p.value)


tdata %>% 
  filter(ID2>0) %>% 
  filter(RNAseq_Type=='Tumor') %>% 
  ggplot(aes(log2(ID2),Exp))+
  geom_point(aes(fill=RNAseq_Type),pch=21,stroke=0.2,size=2.5)+
  geom_smooth(method = 'lm')+
  facet_wrap(~Gene,scales = 'free',nrow = 2) +
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'Y',strip_text_face = 'bold')+
  labs(y='RNA-Seq Expression log2(CPM)',x='ID2 mutations (log2)')+
  theme(panel.spacing = unit(0.5,'lines'))+
  scale_fill_manual(values = c("#ff7f00"))+
  panel_border(color = 'black',linetype = 1)+
  guides(fill="none")

ggsave('ID2_Checkpoint_RNASeq2.pdf',width = 4.5,height = 8,device = cairo_pdf)


tdata %>% 
  group_by(RNAseq_Type,Gene) %>% 
  filter(ID2>0) %>% 
  do(tidy(cor.test(.$Exp,log2(.$ID2),data=.))) %>% 
  arrange(p.value)

# tdata %>% 
#   mutate(Exp=2^Exp) %>% 
#   group_by(RNAseq_Type,ID2) %>% 
#   summarise(mvalue=median(Exp))
# 







# Figure xxx: ID2 association with SV count and signature -----------------------------
my_comparisons <- list(c("Absent",'Present'))

tdata <- sherlock_variable %>% 
  filter(str_detect(name,'SV')) %>% 
  left_join(id2data) %>% 
  filter(Tumor_Barcode %in% hq_samples2)

tdata %>% group_by(name) %>% do(tidy(wilcox.test(log2(value+1) ~ ID2_Present, data=.)))

tdata %>% 
  ggplot(aes(ID2_Present,log2(value+1),fill=ID2_Present))+
  geom_boxplot(width=0.6,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=2.5,position = position_jitter(width = 0.15,height = 0),color="white",stroke=0.05)+
  scale_fill_manual(values = id2color)+
  facet_wrap(~name)+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 4))+
  labs(x = "Mutational signature ID2", y = 'log2(SV count + 1)')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = 'SV_ID2_present_vs_absent.pdf',width = 5,height = 6,device = cairo_pdf)

# SV profiles

tdata <- sherlock_sv32_profile %>% 
  separate(col = MutationType,into = c('a','b','c'),sep = '_') %>% group_by(Tumor_Barcode,a,b) %>% summarise(Contribution = sum(Contribution)) %>% ungroup() %>% 
  left_join(id2data) %>% 
  filter(Tumor_Barcode %in% hq_samples2)

tdata %>% group_by(a,b) %>% do(tidy(wilcox.test(log2(Contribution+1) ~ ID2_Present, data=.))) %>% ungroup()%>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value)) %>% View()

tdata %>% 
  ggplot(aes(ID2_Present,log2(Contribution+1),fill=ID2_Present))+
  geom_boxplot(width=0.6,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=2,position = position_jitter(width = 0.15,height = 0),color="white",stroke=0.05)+
  scale_fill_manual(values = id2color)+
  facet_grid(a~b)+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 4))+
  labs(x = "Mutational signature ID2", y = 'log2(SV count + 1)')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = 'SV_subtypes_ID2_present_vs_absent.pdf',width = 6,height = 8,device = cairo_pdf)

my_comparisons <- list(c("Absent",'Present'))

tdata <- sherlock_variable %>% 
  filter(str_detect(name,'SV')) %>% 
  left_join(id2data) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  left_join(covdata0 %>% select(Tumor_Barcode, Smoking)) %>% 
  filter(Smoking!='Unknown') %>% 
  filter(name=='SV')

tdata %>% group_by(Smoking,name) %>% do(tidy(wilcox.test(log2(value+1) ~ ID2_Present, data=.)))

tdata %>% 
  ggplot(aes(ID2_Present,log2(value+1),fill=ID2_Present))+
  geom_boxplot(width=0.6,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=2.5,position = position_jitter(width = 0.15,height = 0),color="white",stroke=0.05)+
  scale_fill_manual(values = id2color)+
  facet_wrap(~Smoking)+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 4))+
  labs(x = "Mutational signature ID2", y = 'log2(SV count + 1)')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = 'SV_ID2_present_vs_absent_Smoking.pdf',width = 4,height = 6,device = cairo_pdf)

# # Figure xxx: ID2 association with PGA and CNV signature ----------------
my_comparisons <- list(c("Absent",'Present'))
load('../sherlock_PGA.RData',verbose = T)

tdata <- sherlock_PGA %>% 
  left_join(id2data) %>% 
  filter(Tumor_Barcode %in% hq_samples2)

tdata %>% do(tidy(wilcox.test(PGA_WGD ~ ID2_Present, data=.)))

tdata %>% 
  ggplot(aes(ID2_Present,PGA_WGD,fill=ID2_Present))+
  geom_boxplot(width=0.6,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=2.5,position = position_jitter(width = 0.15,height = 0),color="white",stroke=0.05)+
  scale_fill_manual(values = id2color)+
  scale_y_continuous(breaks = pretty_breaks(n = 7),labels = percent_format())+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 4))+
  labs(x = "Mutational signature ID2", y = 'Percent genome altered by copy number')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = 'PGA_ID2_present_vs_absent.pdf',width = 2.7,height = 6,device = cairo_pdf)


my_comparisons <- list(c("Absent",'Present'))
load('../sherlock_PGA.RData',verbose = T)

tdata <- sherlock_PGA %>% 
  left_join(id2data) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  left_join(covdata0 %>% select(Tumor_Barcode, Smoking)) %>% 
  filter(Smoking!='Unknown')

tdata %>% do(tidy(wilcox.test(PGA_WGD ~ ID2_Present, data=.)))

tdata %>% 
  ggplot(aes(ID2_Present,PGA_WGD,fill=ID2_Present))+
  geom_boxplot(width=0.6,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=2.5,position = position_jitter(width = 0.15,height = 0),color="white",stroke=0.05)+
  facet_wrap(~Smoking)+
  scale_fill_manual(values = id2color)+
  scale_y_continuous(breaks = pretty_breaks(n = 7),labels = percent_format())+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 4))+
  labs(x = "Mutational signature ID2", y = 'Percent genome altered by copy number')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = 'PGA_ID2_present_vs_absent_smoking.pdf',width = 4,height = 6,device = cairo_pdf)


# Signature
# Figure 2c-d ---------------------------------------------------------------
nmutation_input <-  0
conflicts_prefer(dplyr::lag)
source('/Users/zhangt8/NIH-Work/MutationSignature/mSigPortal/CBIIT/mSigPortal/Codes/Sigvisualfunc.R')
load('/Users/zhangt8/NIH-Work/MutationSignature/mSigPortal/CBIIT/mSigPortal/Codes/Sigvisualfunc.RData')

sigcol <- ncicolpal[1:length(colnames(cn68_activity)[-1])]
names(sigcol) <- colnames(cn68_activity)[-1]



tmp <- id2data %>% filter(Tumor_Barcode %in% hq_samples2,ID2_Present == "Present")%>% pull(Tumor_Barcode)
sigdata <- cn68_activity %>% rename(Samples=Tumor_Barcode) %>% filter(Samples %in% tmp)
data1 <- prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,output_plot = paste0('ID2_present_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))

tmp <- id2data %>% filter(Tumor_Barcode %in% hq_samples2,ID2_Present == "Absent")%>% pull(Tumor_Barcode)
sigdata <- cn68_activity %>% rename(Samples=Tumor_Barcode) %>% filter(Samples %in% tmp)
data2 <- prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,output_plot = paste0('ID2_absent_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))

tdata <- bind_rows(
  data1$freq_data %>% mutate(Group='Present'),
  data2$freq_data %>% mutate(Group='Absent')
)

myggstyle()

plotdata <- tdata %>% 
  filter(Type == 'Prevalence by samples') %>% 
  arrange(desc(Value)) %>% 
  mutate(Signature = fct_inorder(Signature)) %>% 
  mutate(Value2=if_else(Value<0.1,0.1,Value)) %>% 
  mutate(Value = if_else(Group=='Absent',-1*Value,Value)) %>% 
  mutate(Value2 = if_else(Group=='Absent',-1*Value2,Value2)) 

p1 <- plotdata %>% 
  ggplot(aes(Value, y=Signature,fill=Signature))+
  geom_col(width = 0.75,col='gray10',linewidth=0.3)+
  geom_vline(xintercept = 0,col='gray10')+
  geom_text(data = plotdata %>% filter(Group=='Present'),aes(label=Percentage),vjust=0.5,hjust=2)+
  geom_text(data = plotdata %>% filter(Group=='Absent'),aes(label=Percentage),vjust=0.5,hjust=-1)+
  scale_fill_manual(values = sigcol)+
  scale_x_continuous(limits = c(-1,1),breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),labels = c("100%","75%","50%","25%","0%","25%","50%","75%","100%"))+
  labs(x='Prevalence by samples (%)',y=NULL,fill='Group')+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 12,grid = 'Y',ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  panel_border(color = 'gray10',size = 1)+
  theme(legend.position = 'none',panel.grid.major.y = element_line(linetype = 2,colour = 'gray90'))


testdata <- cn68_activity %>%
  pivot_longer(-Tumor_Barcode) %>% 
  mutate(value=if_else(value>0,'Present','Absent')) %>% 
  mutate(value=factor(value,levels=c('Absent','Present'))) %>% 
  left_join(id2data) %>% 
  filter(Tumor_Barcode %in% hq_samples2,ID2_Present %in% c('Absent','Present')) %>% 
  mutate(ID2_Present = factor(ID2_Present, levels=c('Absent','Present'))) %>% 
  group_by(name) %>% 
  mutate(value=as.factor(value)) %>% 
  do(tresult = safely(stats::fisher.test)(.$value,.$ID2_Present)) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']]))) %>% 
  select(name,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  mutate(FDR=p.adjust(p.value,method='BH'))

p2 <- testdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR),fill=name))+
  geom_hline(yintercept = -log10(0.05),col=ncicolpal[2],linetype=2)+
  geom_hline(yintercept = -log10(0.01),col=ncicolpal[1],linetype=2)+
  geom_vline(xintercept = 0,col="gray10",linetype=1,linewidth=0.5)+
  geom_point(pch=21,size=3.5,col='black',stroke=0.2)+
  ggrepel::geom_text_repel(data=testdata %>% filter(FDR<0.05),aes(label=name))+
  scale_fill_manual(values = sigcol)+
  scale_x_continuous(breaks = pretty_breaks(n = 7),limits = c(-4.5,4.5),position = "top")+
  scale_y_continuous(breaks = pretty_breaks(n=7))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 12,grid = F,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  labs(x = 'Odd ratio (log2)', y = '-log10(FDR)\n')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  coord_cartesian(clip = 'off')

plot_grid(p2,p1,align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2,2))

ggsave(filename = 'CNV_signature_ID2_present_vs_absent.pdf',width = 6,height = 5.5,device = cairo_pdf)



