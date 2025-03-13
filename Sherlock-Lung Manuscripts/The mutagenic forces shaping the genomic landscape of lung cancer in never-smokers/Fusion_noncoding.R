set_wd()
libztw()
pdfhr2()
pdfhr()
myggstyle()

library(gtsummary)
sci_notation <- function(x) ifelse(is.na(x), NA, format(x, digits = 3, scientific = TRUE))




# Fusion ------------------------------------------------------------------
load('../../BBsolution_final3_short.RData')
load('../sampleset.RData',verbose = T)
load('../../RDS/sherlock_data_all.RData')
load('../../RDS/sherlock_variable.RData')
load('../../sherlock_PGA.RData')
load('../../covdata0.RData')

conflicted::conflicts_prefer(ggpubr::get_legend)

# fusion landscape

# 9c ----------------------------------------------------------------------
library(data.table)
source('../../RDS/Oncoplots_functions.R')
tmp <- read_csv('../../RDS/oncoplot_colors.csv')
landscape_colors <- tmp$Color
names(landscape_colors) <- tmp$Name

drglist <- readRDS('../../../../Collaborators/Nuria/Update2/drivers_intogene.RDS') %>% pull(symbol)

# drive gene
data_top <- sherlock_data %>% 
  filter(Type=='Mutation_Driver',Gene %in% drglist) %>% 
  mutate(Alteration = str_remove(Alteration,":.*")) %>% 
  unique()

# fusion
genelist <- c('RET','ALK','AXL','NRG1','MET','FGFR2','ROS1','RB1','ERBB4','EGFR','PTPRT','PTPRB','PTPRD')
data_fusion <- sherlock_data %>% filter(Type=='Fusion',Gene %in% genelist)

# focal CNV 
genelist <- c('BAP1','FAT1','ROS1','CDKN2A','KRAS','BRCA2','RB1','B2M','FANCA','TP53','SMAD4','CHEK2','TERT','IRF4','GRHL2','MYC','RET','KDM5A','MDM2','NKX2-1','GRIN2A','GRB2','KEAP1')
data_focal <- sherlock_data %>% filter(Type=='SCNA_Focal_Gene',Gene %in% genelist)

# CNV ARM
data_arm <- sherlock_data %>% filter(Type=='SCNA_Arm')

# catagoly value
data_feature <- sherlock_overall %>% 
  select(Subject,Tumor_Barcode,Assigned_Population,Gender,Smoking,Histology,WGD_Status,SCNA_Group,MMR_Status,HR_Status,HRD_Type,HLA_LOH,HRDetect_Status,Kataegis,EBV) %>% 
  pivot_longer(cols = -c(Subject,Tumor_Barcode)) %>% 
  filter(!value %in% c('No','None','N')) %>% 
  mutate(Type='Genomic_Feature') %>% 
  select(Subject,Tumor_Barcode,Gene=name,Alteration = value,Type) 

data_feature <- data_feature %>% filter(Gene %in% c('Histology','Gender','Assigned_Population')) %>% mutate(Gene = if_else(Gene=='Gender','Sex',Gene)) %>% mutate(Gene = if_else(Gene=='Assigned_Population','Ancestry',Gene))


data_feature %>% filter(Gene=='Histology') %>% pull(Alteration) %>% unique()->tmp 
landscape_colors[tmp] <- c("#F28833","#AD2741",'#357ABA',"#69AF58","#2C655F")
# continues value 
data_feature2 <- sherlock_overall %>% 
  select(Subject,Tumor_Barcode,Tumor_Purity,Subclonal_Mutation_Ratio,NRPCC) %>% 
  left_join(sherlock_PGA %>% select(-PGA) %>% rename(PGA=PGA_WGD)) %>% 
  pivot_longer(cols = -c(Subject,Tumor_Barcode)) %>% 
  mutate(Type=name) %>% 
  select(Subject,Tumor_Barcode,Gene=name,Alteration= value,Type) 

# TMB
data_tmb <- sherlock_overall %>% 
  select(Subject,Tumor_Barcode,SNV,INDEL,SV) %>% 
  pivot_longer(cols = -c(Subject,Tumor_Barcode)) %>% 
  mutate(Type='Mutations') %>% 
  select(Subject,Tumor_Barcode,Gene=name,Alteration= value,Type) 

data_feature2 <- data_feature2 %>% filter(Gene=='Tumor_Purity')

# Define the sample set
sample_level0 <-  data_fusion %>% filter(Tumor_Barcode %in% sherlock_nonsmoker) %>% pull(Tumor_Barcode) %>% unique()

result_top <- oncoplot(data = data_top,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_scna <- oncoplot(data = data_focal,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_arm <- oncoplot(data = data_arm,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_feature <- oncoplot(data = data_feature,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)

fgene_level <- c(result_fusion$gene_level[!result_fusion$gene_level %in% c('MET','FGFR2','EGFR', 'ERBB4', 'ROS1', 'RET', 'ALK')],result_fusion$gene_level[result_fusion$gene_level %in% c('MET','FGFR2','EGFR', 'ERBB4', 'ROS1', 'RET', 'ALK')]
)

result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,gene_level = fgene_level,GeneSortOnly = TRUE)


sample_new_level <- result_top$sample_level %>% 
  left_join(
    result_fusion$sample_level
  ) %>% left_join(
    result_scna$sample_level
  ) %>% left_join(
    result_arm$sample_level
  ) %>% left_join(
    result_feature$sample_level
  ) %>% arrange(Fusion,Mutation_Driver,SCNA_Focal_Gene,Genomic_Feature) %>% 
  pull(Tumor_Barcode)

result_top <- oncoplot(data = data_top,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar = -0.2,bmar = -0.05)
result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,gene_level = fgene_level)
result_scna <- oncoplot(data = data_focal,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,)
#result_arm <- oncoplot(data = data_arm,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,)
result_feature <- oncoplot(data = data_feature,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE,removeAlterbyColor=TRUE)

result_purity <- oncoplot2(data = data_feature2 %>% filter(Gene == 'Tumor_Purity'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(),tmar=0,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
# result_ratio <- oncoplot2(data = data_feature2 %>% filter(Gene=='Subclonal_Mutation_Ratio'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(option = "C"),tmar=-0.2,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
# result_nrpcc <- oncoplot2(data = data_feature2 %>% filter(Gene =='NRPCC'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material(palette = 'green'),tmar=-0.2,bmar =0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
# result_cnvratio <- oncoplot2(data = data_feature2 %>% filter(Gene =='PGA'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material(palette = 'blue'),tmar=-0.2,bmar =0,p2_axis_hidden = TRUE,p2_hidden = TRUE)

result_tmb <- oncoplot3(data = data_tmb,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar=0.2,bmar=-0.2,height = 6)

oncoplot_final <- oncoplot_combined(result_tmb,result_fusion,result_top,result_feature,result_purity) #result_scna,result_ratio,result_nrpcc,result_cnvratio

save_plot(filename = 'Genome_landscape_fusion_tmp.pdf',plot = oncoplot_final,base_height = 10,base_width = 10,device=cairo_pdf)





## comparison with or without fusion
genelist <- c('RET','ALK','AXL','NRG1','MET','FGFR2','ROS1','RB1','ERBB4','EGFR','PTPRT','PTPRB','PTPRD')
rtk_genes <- c('ALK','RET','ROS1','ERBB4','EGFR','FGFR2','MET')
data_fusion <- sherlock_data_full %>% filter(Type=='Fusion',Gene %in% genelist, Tumor_Barcode %in% sherlock_nonsmoker)

# data_fusion <- data_fusion %>%
#   filter(Gene %in% rtk_genes) %>%
#   mutate(Gene='RTK-RAS Fusions') %>% 
#   group_by(Tumor_Barcode) %>%
#   arrange(Tumor_Barcode,desc(Alteration)) %>%
#   slice(1) %>%
#   ungroup() %>%
#   bind_rows(data_fusion)

tmp1 <- sherlock_data_full %>% filter(Gene=='EGFR',Type=='Mutation_Driver',Tumor_Barcode %in% sherlock_nonsmoker) %>% select(Tumor_Barcode,Var=Alteration) %>% mutate(Group='EGFR Mutation (Wildtype <-> Mutant)')
tmp2 <- sherlock_data_full %>% filter(Gene=='TP53',Type=='Mutation_Driver',Tumor_Barcode %in% sherlock_nonsmoker) %>% select(Tumor_Barcode,Var=Alteration) %>% mutate(Group='TP53 Mutation (Wildtype <-> Mutant)')
tmp3 <- sherlock_data_full %>% filter(Gene=='SETD2',Type=='Mutation_Driver',Tumor_Barcode %in% sherlock_nonsmoker) %>% select(Tumor_Barcode,Var=Alteration) %>% mutate(Group='SETD2 Mutation (Wildtype <-> Mutant)')
tmp4 <- covdata0 %>% filter(Tumor_Barcode %in% sherlock_nonsmoker) %>% mutate(Gender = factor(Gender,levels=c('Male','Female'),labels=c('No','Yes'))) %>% select(Tumor_Barcode,Var=Gender) %>% mutate(Group='Sex (Male <-> Female)')
tmp5 <- covdata0 %>% filter(Tumor_Barcode %in% sherlock_nonsmoker) %>% mutate(Assigned_Population = factor(Assigned_Population,levels=c('EUR','EAS'),labels=c('No','Yes'))) %>% select(Tumor_Barcode,Var=Assigned_Population) %>% mutate(Group='Ancestry (EUR <-> EAS)')

tmp <- bind_rows(tmp1,tmp2,tmp3,tmp4,tmp5) 

tdata <- data_fusion %>% left_join(tmp) %>% left_join(covdata0)

pdata <- tdata %>% 
  group_by(Gene,Group) %>% 
  do(tresult = safely(fisher.test)(.$Alteration,.$Var )) %>% 
  #do(tresult = safely(glm)(as.factor(Alteration) ~ Var + Gender + Age + Tumor_Purity, family='binomial',data=. )) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']]))) %>% 
  select(Gene,Group,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  #filter(term=='VarYes') %>% 
  arrange(p.value) %>% 
  group_by(Group) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  ungroup() %>% 
  mutate(Group = fct_inorder(Group))

pdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR),fill=Group))+
  geom_hline(yintercept = -log10(0.05),col=ncicolpal[2],linetype=2)+
  geom_hline(yintercept = -log10(0.01),col=ncicolpal[1],linetype=2)+
  geom_vline(xintercept = 0,col="gray10",linetype=1,linewidth=0.5)+
  geom_point(pch=21,size=3.5,col='white',stroke=0.2)+
  ggrepel::geom_text_repel(data=pdata %>% filter(FDR<0.05),aes(label=Gene))+
  facet_wrap(~Group,nrow = 1)+
  scale_fill_nejm()+
  # scale_x_continuous(breaks = pretty_breaks(n = 7),limits = c(-4.05,4.05))+
  scale_y_continuous(breaks = pretty_breaks(n=5))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 12,grid = F,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing.x = unit('0.25','cm'),strip.text.x = element_text(face = 'bold',hjust = 0.5),legend.position = 'bottom')+
  labs(x = 'Odds ratio (log2)', y = '-log10(FDR)',fill='Enrichment')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  coord_cartesian(clip = 'off')+
  guides(fill='none')

ggsave(filename = 'WGS_fusion_association_tmp.pdf',width = 14,height = 3,device = cairo_pdf)



## signature ##
#load('../Signature_Lumidl.RData')
tmp <- ludmil_activity_all_obs %>% 
  pivot_longer(-Tumor_Barcode,names_to = 'Group',values_to = 'Var') %>% 
  mutate(Var=if_else(Var,'Yes','No')) 
  

tdata <- data_fusion %>% left_join(tmp) %>% left_join(covdata0) %>% 
  mutate(Alteration = factor(Alteration, levels = c('No','Yes'))) %>% 
  mutate(Var = factor(Var, levels = c('No','Yes'))) 

pdata <- tdata %>% 
  group_by(Gene,Group) %>% 
  #do(tresult = safely(fisher.test)(.$Alteration,.$Var )) %>% 
  do(tresult = safely(glm)(as.factor(Alteration) ~ Var + Gender + Age + Tumor_Purity + Histology, family='binomial',data=. )) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']]))) %>% 
  select(Gene,Group,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  filter(term=='VarYes') %>% 
  arrange(p.value) %>% 
  group_by(Group) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  ungroup() %>% 
  mutate(Group = fct_inorder(Group)) %>% 
  filter(abs(estimate)<10)

pdata %>% 
  ggplot(aes((estimate),-log10(FDR),fill=Gene))+
  geom_hline(yintercept = -log10(0.05),col=ncicolpal[2],linetype=2)+
  geom_hline(yintercept = -log10(0.01),col=ncicolpal[1],linetype=2)+
  geom_vline(xintercept = 0,col="gray10",linetype=1,linewidth=0.5)+
  geom_point(pch=21,size=3.5,col='white',stroke=0.2)+
  ggrepel::geom_text_repel(data=pdata %>% filter(FDR<0.05),aes(label=Group))+
  facet_wrap(~Gene,nrow = 2)+
  scale_fill_manual(values = ncicolpal)+
  # scale_x_continuous(breaks = pretty_breaks(n = 7),limits = c(-4.05,4.05))+
  scale_y_continuous(breaks = pretty_breaks(n=5))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 12,grid = F,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing.x = unit('0.25','cm'),strip.text.x = element_text(face = 'bold',hjust = 0.5),legend.position = 'bottom')+
  labs(x = "Logistic regression coefficient", y = '-log10(FDR)',fill='Enrichment')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  coord_cartesian(clip = 'off')+
  guides(fill='none')

ggsave(filename = 'WGS_fusion_association2_tmp.pdf',width = 14,height = 5,device = cairo_pdf)



## Age ## 
genelist <- c('RET','ALK','AXL','NRG1','MET','FGFR2','ROS1','RB1','ERBB4','EGFR','PTPRT','PTPRB','PTPRD')
rtk_genes <- c('ALK','RET','ROS1','ERBB4','EGFR','FGFR2')
data_fusion <- sherlock_data %>% filter(Type=='Fusion',Gene %in% genelist, Tumor_Barcode %in% sherlock_nonsmoker)

genecol <- c('#cccccc',pal_jama()(6))
names(genecol) <- c('Non-Fusion',rtk_genes)
# data_fusion <- data_fusion %>%
#   filter(Gene %in% rtk_genes) %>%
#   mutate(Gene='RTK-RAS Fusions') %>% 
#   group_by(Tumor_Barcode) %>%
#   arrange(Tumor_Barcode,desc(Alteration)) %>%
#   slice(1) %>%
#   ungroup() %>%
#   bind_rows(data_fusion)

pdata <- tibble(Tumor_Barcode=sherlock_nonsmoker) %>% filter(!(Tumor_Barcode %in% data_fusion$Tumor_Barcode)) %>% mutate(Gene='Non-Fusion') %>% 
  bind_rows(data_fusion %>% select(Tumor_Barcode,Gene)) %>% 
  left_join(
    covdata0 %>% select(Tumor_Barcode,Age)
  )

pdata <- pdata %>% 
  filter(Gene %in% c(rtk_genes,'Non-Fusion','RTK-RAS Fusions'))

pdata <- pdata %>% mutate(Gene=fct_reorder(Gene,Age))



#my_comparisons <-  list(c('Non-Fusion','ALK'),c('Non-Fusion','EGFR'),c('Non-Fusion','ERBB4'),c('Non-Fusion','FGFR2'),c('Non-Fusion','NRG1'),c('Non-Fusion','PTPRB'),c('Non-Fusion','PTPRD'),c('Non-Fusion','PTPRT'),c('Non-Fusion','RB1'),c('Non-Fusion','RET'),c('Non-Fusion','ROS1'),c('Non-Fusion','RTK-RAS Fusions'))

my_comparisons <-  list(c('Non-Fusion','ALK'),c('Non-Fusion','EGFR'),c('Non-Fusion','ERBB4'),c('Non-Fusion','FGFR2'),c('Non-Fusion','RET'),c('Non-Fusion','ROS1')) #,c('Non-Fusion','RTK-RAS Fusions')

stat.test <- pdata %>% 
  rstatix::wilcox_test(Age~ Gene, comparisons = my_comparisons) %>%
  rstatix::add_xy_position(x = "Gene") %>%
  mutate(myformatted.p = if_else(p<0.001,sprintf("P = %.2e",p),sprintf("P = %.2f",p))) %>% 
  mutate(Gene = my_comparisons[[1]][1]) %>% 
  mutate(col=if_else(p<0.06,ncicolpal[1],'black')) 
  #filter(p<0.06) 

stat.test

pdata %>% 
  ggplot(aes(Gene,Age,fill=Gene))+
  #geom_boxplot()+
  #facet_wrap(~name,scales = 'free_y',nrow = 1)+
  #geom_violin(trim=T,size=0.2)+
  #geom_boxplot(width=0.15, fill="gray95",color='black',outlier.shape = 21,outlier.size = 0.8)+
  ggbeeswarm::geom_quasirandom(pch=21,size=2,width = 0.2,color="black",stroke=0.2)+
  geom_boxplot(width=0.7,outlier.shape = NA,alpha =0.6,size=0.4)+
  scale_fill_manual(values = genecol)+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  labs(x=NULL,y='Age at diagnosis',color='Group')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = 'Y')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,4,4))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+  panel_border(color = 'black',size = 0.5)+
  theme(legend.position = 'none')+
  stat_pvalue_manual(data = stat.test, col='col',label = "myformatted.p")+
  scale_color_manual(values = c(ncicolpal[1],'black'))

ggsave(filename = 'WGS_fusion_association_age_tmp.pdf',width = 4,height =7,device = cairo_pdf)  

### RTK-Pathway fusion P=0.030



## TMB ## 
genelist <- c('RET','ALK','AXL','NRG1','MET','FGFR2','ROS1','RB1','ERBB4','EGFR','PTPRT','PTPRB','PTPRD')
rtk_genes <- c('ALK','RET','ROS1','ERBB4','EGFR','FGFR2')
data_fusion <- sherlock_data %>% filter(Type=='Fusion',Gene %in% genelist, Tumor_Barcode %in% sherlock_nonsmoker)

# data_fusion <- data_fusion %>%
#   filter(Gene %in% rtk_genes) %>%
#   mutate(Gene='RTK-RAS Fusions') %>%
#   group_by(Tumor_Barcode) %>%
#   arrange(Tumor_Barcode,desc(Alteration)) %>%
#   slice(1) %>%
#   ungroup() %>%
#   bind_rows(data_fusion)

pdata <- tibble(Tumor_Barcode=sherlock_nonsmoker) %>% filter(!(Tumor_Barcode %in% data_fusion$Tumor_Barcode)) %>% mutate(Gene='Non-Fusion') %>% 
  bind_rows(data_fusion %>% select(Tumor_Barcode,Gene)) %>% 
  left_join(
    sherlock_variable %>% filter(name=='TMB') %>% mutate(TMB=log2(value)) %>% select(Tumor_Barcode,TMB)
    )

pdata <- pdata %>% 
  filter(Gene %in% c(rtk_genes,'Non-Fusion','RTK-RAS Fusions'))

pdata <- pdata %>% mutate(Gene=fct_reorder(Gene,TMB))



#my_comparisons <-  list(c('Non-Fusion','ALK'),c('Non-Fusion','EGFR'),c('Non-Fusion','ERBB4'),c('Non-Fusion','FGFR2'),c('Non-Fusion','NRG1'),c('Non-Fusion','PTPRB'),c('Non-Fusion','PTPRD'),c('Non-Fusion','PTPRT'),c('Non-Fusion','RB1'),c('Non-Fusion','RET'),c('Non-Fusion','ROS1'),c('Non-Fusion','RTK-RAS Fusions'))

my_comparisons <-  list(c('Non-Fusion','ALK'),c('Non-Fusion','EGFR'),c('Non-Fusion','ERBB4'),c('Non-Fusion','FGFR2'),c('Non-Fusion','RET'),c('Non-Fusion','ROS1')) #,c('Non-Fusion','RTK-RAS Fusions')

stat.test <- pdata %>% 
  rstatix::wilcox_test(TMB~ Gene, comparisons = my_comparisons) %>%
  rstatix::add_xy_position(x = "Gene") %>%
  mutate(myformatted.p = if_else(p<0.001,sprintf("P = %.2e",p),sprintf("P = %.2f",p))) %>% 
  mutate(Gene = my_comparisons[[1]][1]) %>% 
  mutate(col=if_else(p<0.06,ncicolpal[1],'black')) 
#filter(p<0.06) 

stat.test

pdata %>% 
  ggplot(aes(Gene,TMB,fill=Gene))+
  #geom_boxplot()+
  #facet_wrap(~name,scales = 'free_y',nrow = 1)+
  #geom_violin(trim=T,size=0.2)+
  #geom_boxplot(width=0.15, fill="gray95",color='black',outlier.shape = 21,outlier.size = 0.8)+
  ggbeeswarm::geom_quasirandom(pch=21,size=2,width = 0.2,color="black",stroke=0.2)+
  geom_boxplot(width=0.7,outlier.shape = NA,alpha =0.6,size=0.4)+
  scale_fill_manual(values = genecol)+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  labs(x=NULL,y='log2(Number of mutations per megabase)',color='Group')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = 'Y')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),panel.spacing = unit(0.4, "lines"),strip.background = element_blank(),strip.text = element_text(hjust = 0.5),plot.margin = margin(4,4,4,4))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+  panel_border(color = 'black',size = 0.5)+
  theme(legend.position = 'none')+
  stat_pvalue_manual(data = stat.test, col='col',label = "myformatted.p")+
  scale_color_manual(values = c(ncicolpal[1],'black'))

ggsave(filename = 'WGS_fusion_association_tmb_tmp.pdf',width = 4,height =7,device = cairo_pdf)  

### RTK-Pathway fusion P=0.030

### parterner
data_fusion <-  read_delim('/Volumes/data/NSLC2/SV/Meerkat_Analysis/all_fusion_sample.txt',delim = '\t',col_names = T)
data_fusion <- data_fusion %>% filter(type1 == "gene-gene",str_detect(type3,'frame'), geneA %in% genelist | geneB %in% genelist,barcode %in% sherlock_nonsmoker) %>% 
  select(Tumor_Barcode=barcode,geneA,geneB)

known_fusion <- readxl::read_xlsx('~/NIH-Work/EAGLE_NSLC/Collaborators/Lixing/Fusion/table.S4.driver.fusions.xlsx')
known_fusion <- known_fusion %>% select(Tumor_Barcode=Sample_ID,geneA=Partner_5prime,geneB=Partner_3prime)


data_fusion <- bind_rows(data_fusion,known_fusion) %>% unique()

tmp <- data_fusion %>% 
  filter(geneA == 'ALK' | geneB =='ALK') %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  unique() %>% 
  count(value,sort=T) %>% 
  filter(value!='ALK')

tmp <- data_fusion %>% 
  filter(geneA == 'ALK' | geneB =='ALK') %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  filter(value!='ALK') %>% 
  mutate(value=factor(value,levels=tmp$value)) %>% 
  group_by(Tumor_Barcode) %>% 
  arrange(value) %>% 
  slice(1) %>% 
  ungroup() %>% 
  count(value,sort=T)


source('~/NIH-Work/R/ZTW_function/ztw.R')
PieDonut_ztw(tmp,aes(pies=value,count=n),mainCol = pal_bmj()(6),start=3*pi/2,showNum = T,showRatioPie = F,showRatioThreshold = 0.02,showDonutName = T,pieLabelSize=3.5,titlesize = 5,r0=0.3,title='AS_N', showPieName = FALSE,family = 'Roboto Condensed', labelpositionThreshold =0.02)
ggsave(filename = 'ALK_fusion.pdf',width = 2.5,height = 2.5,device = cairo_pdf)


tmp <- data_fusion %>% 
  filter(geneA == 'RET' | geneB =='RET') %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  unique() %>% 
  count(value,sort=T) %>% 
  filter(value!='RET')

tmp <- data_fusion %>% 
  filter(geneA == 'RET' | geneB =='RET') %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  filter(value!='RET') %>% 
  mutate(value=factor(value,levels=tmp$value)) %>% 
  group_by(Tumor_Barcode) %>% 
  arrange(value) %>% 
  slice(1) %>% 
  ungroup() %>% 
  count(value,sort=T)

source('~/NIH-Work/R/ZTW_function/ztw.R')
PieDonut_ztw(tmp,aes(pies=value,count=n),mainCol = ncicolpal,start=3*pi/2,showNum = T,showRatioPie = F,showRatioThreshold = 0.02,showDonutName = T,pieLabelSize=3.5,titlesize = 5,r0=0.3,title='AS_N', showPieName = FALSE,family = 'Roboto Condensed', labelpositionThreshold =0.02,title)
ggsave(filename = 'RET_fusion.pdf',width = 2.5,height = 2.5,device = cairo_pdf)


tmp <- data_fusion %>% 
  filter(geneA == 'ROS1' | geneB =='ROS1') %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  unique() %>% 
  count(value,sort=T) %>% 
  filter(value!='ROS1')

tmp <- data_fusion %>% 
  filter(geneA == 'ROS1' | geneB =='ROS1') %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  filter(value!='ROS1') %>% 
  mutate(value=factor(value,levels=tmp$value)) %>% 
  group_by(Tumor_Barcode) %>% 
  arrange(value) %>% 
  slice(1) %>% 
  ungroup() %>% 
  count(value,sort=T)

source('~/NIH-Work/R/ZTW_function/ztw.R')
PieDonut_ztw(tmp,aes(pies=value,count=n),mainCol = ncicolpal,start=3*pi/2,showNum = T,showRatioPie = F,showRatioThreshold = 0.02,showDonutName = T,pieLabelSize=3.5,titlesize = 5,r0=0.3,title='AS_N', showPieName = FALSE,family = 'Roboto Condensed', labelpositionThreshold =0.02,title)
ggsave(filename = 'ROS1_fusion.pdf',width = 2.5,height = 2.5,device = cairo_pdf)

tmp <- data_fusion %>% 
  filter(geneA == 'NRG1' | geneB =='NRG1') %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  unique() %>% 
  count(value,sort=T) %>% 
  filter(value!='NRG1')

tmp <- data_fusion %>% 
  filter(geneA == 'NRG1' | geneB =='NRG1') %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  filter(value!='NRG1') %>% 
  mutate(value=factor(value,levels=tmp$value)) %>% 
  group_by(Tumor_Barcode) %>% 
  arrange(value) %>% 
  slice(1) %>% 
  ungroup() %>% 
  count(value,sort=T)

source('~/NIH-Work/R/ZTW_function/ztw.R')
PieDonut_ztw(tmp,aes(pies=value,count=n),mainCol = ncicolpal,start=3*pi/2,showNum = T,showRatioPie = F,showRatioThreshold = 0.02,showDonutName = T,pieLabelSize=3.5,titlesize = 5,r0=0.3,title='AS_N', showPieName = FALSE,family = 'Roboto Condensed', labelpositionThreshold =0.02,title)
ggsave(filename = 'NRG1_fusion.pdf',width = 4.5,height = 4.5,device = cairo_pdf)


tmp <- data_fusion %>% 
  filter(geneA == 'ERBB4' | geneB =='ERBB4') %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  unique() %>% 
  count(value,sort=T) %>% 
  filter(value!='ERBB4')

tmp <- data_fusion %>% 
  filter(geneA == 'ERBB4' | geneB =='ERBB4') %>% 
  pivot_longer(cols = -Tumor_Barcode) %>% 
  filter(value!='ERBB4') %>% 
  mutate(value=factor(value,levels=tmp$value)) %>% 
  group_by(Tumor_Barcode) %>% 
  arrange(value) %>% 
  slice(1) %>% 
  ungroup() %>% 
  count(value,sort=T)

source('~/NIH-Work/R/ZTW_function/ztw.R')
PieDonut_ztw(tmp,aes(pies=value,count=n),mainCol = ncicolpal,start=3*pi/2,showNum = T,showRatioPie = F,showRatioThreshold = 0.02,showDonutName = T,pieLabelSize=3.5,titlesize = 5,r0=0.3,title='AS_N', showPieName = FALSE,family = 'Roboto Condensed', labelpositionThreshold =0.02,title)
ggsave(filename = 'ERBB4_fusion.pdf',width = 4.5,height = 4.5,device = cairo_pdf)







# Noncoding ---------------------------------------------------------------
load('/Volumes/data/NSLC2/ActiveDriverWGS/sherlock_1217/sherlock_nocoding.RData',verbose = T)
load('/Volumes/data/NSLC2/ActiveDriverWGS/sherlock_1217/varinfo.RData',verbose = T)

tmp <- sherlock_nocoding_detail %>%
  filter(Tumor_Barcode %in% sherlock_nonsmoker) %>% 
  select(Tumor_Barcode,encode_id) %>%
  unique() %>%
  count(encode_id,sort=T) %>%
  mutate(freq_nonsmoker=n/871) %>% 
  rename(id=encode_id)



data <- sherlock_nocoding %>% left_join(tmp)
data %>%
  mutate(gene=if_else(gene=='ADGRG6','GPR126',gene)) %>% 
  filter(fdr_element<0.01) %>%
  mutate(gene = if_else(freq_nonsmoker<0.01,"",gene)) %>% 
  ggplot(aes(freq_nonsmoker ,(-log10(pp_element)),fill=element_type))+
  geom_point(pch=21,size=3)+
  geom_text_repel(aes(label=gene,col=element_type),size=3,force = 5,max.overlaps = 50)+
  labs(x="Mutation Frequency", y="Significance, -log10(pvalue)",fill="",col="")+
  scale_fill_aaas()+
  scale_color_aaas()+
  scale_x_continuous(breaks = pretty_breaks(),expand = c(0.01,0),,labels=percent_format())+
  scale_y_continuous(breaks = pretty_breaks(n=5))+
  theme_ipsum_rc(base_size = 16,axis_title_just = "m",axis_title_size = 16,axis = TRUE)+
  theme(legend.position = "right")+
  panel_border(color = "black",size = 0.4)+
  coord_cartesian(clip = 'off')

ggsave(filename = 'nonsmokers_nonCoding_drivers3.pdf',width = 9.5,height = 5,device = cairo_pdf )





query_gene='EGFR'

idtmp <- sherlock_nocoding %>% filter(gene==query_gene) %>% dplyr::slice(1) %>% pull(id)
etmp <- sherlock_nocoding_detail %>% filter(id==idtmp) %>% dplyr::slice(1)

dtmp <- sherlock_nocoding_detail %>%
  filter(id==idtmp,Tumor_Barcode %in% sherlock_nonsmoker) %>%
  #mutate(mutation_type=paste0(ref.mut,"->", alt.mut)) %>%
  count(chrom,start.mut ,ref.mut,alt.mut) %>%
  arrange(desc(n)) %>%
  mutate(chrom=str_remove(chrom,'chr')) %>%
  select(CHROM=chrom,POS=start.mut,REF=ref.mut,ALT=alt.mut,n) %>%
  left_join(varinfo)

pmax <- dtmp$n[1]
pmax <- if_else(pmax<2,2L,pmax)
# attachContext(mutData = dtmp %>% filter(str_length(ref.mut)==1,str_length(alt.mut)==1) %>% mutate(chrom=str_remove(chrom,'chr')), chr_colName = 'chrom',start_colName = 'start.mut',end_colName = 'start.mut',nucl_contextN = 3,BSGenomeDb = hg38)

COLORS6 = c(  "#2EBAED", "#000000", "#DE1C14",  "#D4D2D2", "#ADCC54", "#F0D0CE",'#9467BDFF','#FF7F0EFF')
names(COLORS6) <- c('C>A','C>G','C>T','T>A','T>C','T>G','DEL','INS')

dtmp %>% ggplot(aes(POS,n))+geom_linerange(aes(x=POS,ymin=0, ymax=n),color="#aaaaaa",size=0.3)+geom_jitter(aes(fill=mutType2),pch=21,width = 0,height = 0.05,size=2)+theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 16,grid = "yY",axis = "y",ticks = T)+scale_fill_manual(values = COLORS6)+guides(fill = guide_legend(override.aes=list(shape=21)))+labs(x=paste0("\nChromosome ",str_remove(etmp$chrom,'chr'), "(bp)"),y="Number of Mutations\n",fill='Mutation Type')+scale_x_comma(limits =c(etmp$start,etmp$end),expand = c(0.01,0.01) )+scale_y_continuous(expand = c(0,0),limits = c(0,pmax*1.15))+
  #geom_rect(xmin=etmp$start, xmax=etmp$end, ymin=pmax*1.05, ymax=pmax*1.1, fill='#cccccc', color="black", alpha=0.5)+
  geom_text(x=(etmp$start+etmp$end)/2,y=pmax*1.075,label=paste0(idtmp,'/GPR126'),size=5)

ggsave(file = "mutation_plot_ADGRG6.pdf",heigh=5,width=8,device=cairo_pdf)

ggsave(file = "mutation_plot_EGFR.pdf",heigh=5,width=8,device=cairo_pdf)



