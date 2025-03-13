
# Fig. 5a -----------------------------------------------------------------
# figure modified based a publication by Levin and Mora (ref 58), the mechanism depicts non-LTR retrotransposons mobilizing through target-site-primed reverse transcription (TPRT)


libztw()
set_wd()

# load the data -----------------------------------------------------------
load('../../BBsolution_final3_short.RData')

tmp <- readxl::read_xlsx('email/Sherlock data MC (Methylation levels){BE}.xlsx',sheet = 1,col_names = T,n_max = 1)
tmp <- tmp %>% 
  pivot_longer(cols = -c(`Region Chr22:`, `...2` )) %>%
  select(Pos=name,CpG_ID=value)

tmp2 <- readxl::read_xlsx('email/Sherlock data MC (Methylation levels){BE}.xlsx',sheet = 1,col_names = T,skip = 1)
mdata <- tmp2 %>% 
  pivot_longer(cols = -c(`CpG n.`,`Tumor/Normal/Control` )) %>% 
  rename(Barcode=`CpG n.`, Type=`Tumor/Normal/Control`,CpG_ID=name,Beta=value) %>% 
  filter(str_detect(Barcode,'_')) %>% 
  left_join(tmp)

mdata <- mdata %>% 
  mutate(Sherlock_PID=str_remove(Barcode,'^[0-9]*_')) %>% 
  mutate(Sherlock_PID=str_remove(Sherlock_PID,'_[TN]$')) %>% 
  filter(!CpG_ID %in% c('Tum met cluster','Normal Tissue Met cluster')) %>% 
  mutate(Sherlock_PID = if_else(Sherlock_PID =='HOKO-153','HOKO-0153',Sherlock_PID)) %>% 
  mutate(Sherlock_PID = if_else(Sherlock_PID =='ITLU_0462','ITLU-0462',Sherlock_PID)) %>% 
  mutate(Sherlock_PID = if_else(Sherlock_PID =='TWAIN-0016','TWAN-0016',Sherlock_PID)) %>% 
  mutate(Sherlock_PID = if_else(Sherlock_PID =='TWAIN-0042','TWAN-0042',Sherlock_PID)) %>% 
  left_join(wgs_groups_info %>% select(Sherlock_PID,Subject)) %>% 
  mutate(Subject = if_else(Barcode=='64_100_PERCENT','100% methylated',Subject)) %>% 
  mutate(Subject = if_else(Barcode=='63_66.6_PERCENT','66.6% methylated_3',Subject)) %>% 
  mutate(Subject = if_else(Barcode=='62_66.6_PERCENT','66.6% methylated_2',Subject)) %>% 
  mutate(Subject = if_else(Barcode=='61_66.6_PERCENT','66.6% methylated_1',Subject)) %>% 
  mutate(Subject = if_else(Barcode=='60_33.3_PERCENT','33.3% methylated',Subject)) %>% 
  mutate(Subject = if_else(Barcode=='59_UNMETHYLATED','0% methylated',Subject)) %>% 
  mutate(Beta=as.numeric(Beta)) %>% 
  mutate(Pos=as.integer(Pos)) %>% 
  mutate(Type=if_else(Type=='T','Tumor',Type)) %>% 
  mutate(Type=if_else(Type=='N','Normal',Type)) %>% 
  mutate(Type=if_else(Type=='C','Control',Type))

mdata_control <- mdata %>% filter(Type=='Control')

tmp3 <- readxl::read_xlsx('email/Sherlock data MC 2024-2023 (Methylation levels){BE}.xlsx',sheet = 2,col_names = T,skip = 2)
mdata <- tmp3 %>% 
  pivot_longer(cols = -c(`Year of analysis`,`Sample name`,`Tumor/Normal/Control` )) %>% 
  select(-`Year of analysis`) %>% 
  rename(Barcode=`Sample name`, Type=`Tumor/Normal/Control`,Pos=name,Beta=value)

mdata <- mdata %>% 
  mutate(Sherlock_PID=str_remove(Barcode,'^[0-9]*_')) %>% 
  mutate(Sherlock_PID=str_remove(Sherlock_PID,'_[TN]$')) %>% 
  mutate(Sherlock_PID=str_remove(Sherlock_PID,'[TN]$')) %>% 
  mutate(Sherlock_PID=str_replace(Sherlock_PID,'_','-')) %>% 
  #filter(!CpG_ID %in% c('Tum met cluster','Normal Tissue Met cluster')) %>% 
  mutate(Sherlock_PID = if_else(Sherlock_PID =='HOKO-153','HOKO-0153',Sherlock_PID)) %>% 
  mutate(Sherlock_PID = if_else(Sherlock_PID =='HOKO-004','HOKO-0004',Sherlock_PID)) %>%
  mutate(Sherlock_PID = if_else(Sherlock_PID =='HOKO-009','HOKO-0009',Sherlock_PID)) %>%
  mutate(Sherlock_PID = if_else(Sherlock_PID =='ITLU_0462','ITLU-0462',Sherlock_PID)) %>% 
  mutate(Sherlock_PID = if_else(Sherlock_PID =='TWAIN-0016','TWAN-0016',Sherlock_PID)) %>% 
  mutate(Sherlock_PID = if_else(Sherlock_PID =='TWAIN-0042','TWAN-0042',Sherlock_PID)) %>% 
  left_join(wgs_groups_info %>% select(Sherlock_PID,Subject)) %>% 
  # mutate(Subject = if_else(Barcode=='64_100_PERCENT','100% methylated',Subject)) %>% 
  # mutate(Subject = if_else(Barcode=='63_66.6_PERCENT','66.6% methylated_3',Subject)) %>% 
  # mutate(Subject = if_else(Barcode=='62_66.6_PERCENT','66.6% methylated_2',Subject)) %>% 
  # mutate(Subject = if_else(Barcode=='61_66.6_PERCENT','66.6% methylated_1',Subject)) %>% 
  # mutate(Subject = if_else(Barcode=='60_33.3_PERCENT','33.3% methylated',Subject)) %>% 
  # mutate(Subject = if_else(Barcode=='59_UNMETHYLATED','0% methylated',Subject)) %>% 
  mutate(Beta=as.numeric(Beta)) %>% 
  mutate(Pos=as.integer(Pos)) %>% 
  mutate(Type=if_else(Type=='T','Tumor',Type)) %>% 
  mutate(Type=if_else(Type=='N','Normal',Type)) %>% 
  mutate(Type=if_else(Type=='C','Control',Type))

mdata <- mdata %>% left_join(
  mdata_control %>% select(CpG_ID,Pos) %>% unique()
)

mdata <- bind_rows(mdata_control,mdata)

save(mdata,file='CGW_data.RData')



## new limitation to hqsamples for luad only
load('../../BBsolution_final3_short.RData')
load('../../covdata0.RData')
hq_samples2 <- covdata0 %>% filter(Histology == 'Adenocarcinoma', Tumor_Barcode %in% hq_samples) %>% select(Tumor_Barcode) %>% left_join(wgs_groups_info) %>% pull(Subject)
rm(list=c('hq_samples'))

load('../../ID2_TE/id2data.RData')
load('../../sp_group_data.RData')

tmpdata <- wgs_groups_info %>% 
  left_join(sp_group_data2) %>% 
  left_join(id2data) %>% 
  select(Tumor_Barcode,Normal_Barcode,Subject,SP_Group_New,ID2,ID2_ratio,ID2_Present)

mdata <- mdata %>% left_join(tmpdata) %>% mutate(Sample=if_else(Type=='Tumor',Tumor_Barcode,if_else(Type=='Normal',Normal_Barcode,Subject)))


# 
# ## Absent 6 and Present 23 have been selected before, no we need to select addtional samples for the sequecing
# tmp <- tmpdata %>% filter(Subject %in% hq_samples2, !(Subject %in% mdata$Subject))
# 
# tmp %>% filter(ID2_Present=="Present") %>% arrange(desc(ID2)) %>% select(Tumor_Barcode,ID2,ID2_ratio,ID2_Present)->select1
# select1 <- select1 %>% left_join(wgs_groups_info) %>%
#   filter(Study!='Public-WGS') %>% 
#   select(Subject,Subject_ID,Sherlock_PID,Tumor_Barcode,Tumor_Sample_ID,Normal_Barcode,Normal_Sample_ID,contains("ID2"))
# 
# tmp %>% filter(ID2_Present=="Absent") %>% arrange((ID2)) %>% select(Tumor_Barcode,ID2,ID2_ratio,ID2_Present)->select2
# select2 <- select2 %>% left_join(wgs_groups_info) %>%
#   filter(Study!='Public-WGS') %>% 
#   select(Subject,Subject_ID,Sherlock_PID,Tumor_Barcode,Tumor_Sample_ID,Normal_Barcode,Normal_Sample_ID,contains("ID2"))
# 
# bind_rows(select1,select2) %>%
#   write_csv('additional_selelection.csv',col_names = T)


# Analysis ----------------------------------------------------------------

# phase plot

# tmp <- mdata %>% 
#   filter(Type=='Normal') %>%
#   select(Sample,Pos,Beta) %>% 
#   pivot_wider(names_from = Pos,values_from = Beta)
# 
# tmpm <- as.matrix(tmp[,-1])
# rownames(tmpm) <- tmp$Sample
# tmpclust <- hclust(dist(tmpm))
# order1 <- rev(rownames(tmpm)[tmpclust$order])


tmp <- mdata %>% 
  filter(Type=='Tumor',ID2_Present=='Absent') %>%
  select(Sample,Pos,Beta) %>% 
  pivot_wider(names_from = Pos,values_from = Beta)

tmpm <- as.matrix(tmp[,-1])
rownames(tmpm) <- tmp$Sample
tmpclust <- hclust(dist(tmpm))
order2 <- (rownames(tmpm)[tmpclust$order])


tmp <- mdata %>% 
  filter(Type=='Tumor',ID2_Present=='Present') %>%
  select(Sample,Pos,Beta) %>% 
  pivot_wider(names_from = Pos,values_from = Beta)

tmpm <- as.matrix(tmp[,-1])
rownames(tmpm) <- tmp$Sample
tmpclust <- hclust(dist(tmpm))
order3 <- (rownames(tmpm)[tmpclust$order])

order1 <- tibble(Tumor_Barcode=c(order3,order2)) %>% left_join(wgs_groups_info) %>% pull(Normal_Barcode)

samleves <- c(rev(c('100% methylated','66.6% methylated_1','66.6% methylated_2','66.6% methylated_3','33.3% methylated','0% methylated')),order1,order3,order2)

mdata <- mdata %>% mutate(Sample=factor(Sample,levels=samleves)) %>% arrange(Sample)

mdata <- mdata %>% 
  mutate(ID2=if_else(Type=='Normal',NA_integer_,ID2)) %>% 
  mutate(ID2_Present=if_else(Type=='Normal',NA_character_,ID2_Present))

myggstyle()

mdata %>% 
  #filter(Type=='Tumor') %>% 
  ggplot(aes(Pos,Sample,fill=Beta))+
  geom_point(pch=21,size=2)+
  scale_fill_viridis_c(limits = c(0, 1))+
  scale_x_continuous(breaks = sort(unique(mdata$Pos)),labels = sort(unique(mdata$Pos)),expand = c(0,0))+
  facet_grid(Type~.,scales = 'free_y',space='free')+
  theme_ipsum_rc(axis_title_just = 'm',ticks = T,axis = 'Y')+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1),axis.text.y = element_text(size=7),panel.spacing = unit(0.5,'cm'),legend.position = 'top',legend.key.width = unit(2.5,'cm'), panel.grid.minor.x = element_blank())+
  labs(x='',y='',fill='Methylation beta value\n')+
  ggnewscale::new_scale_fill()+
  geom_tile(aes(x=29059215,Sample,fill=log2(ID2+1),width=5))+
  scale_fill_viridis_c(option = 'A')+
  ggnewscale::new_scale_fill()+
  geom_tile(aes(x=29059210,Sample,fill=SP_Group_New,width=5))+
  scale_fill_manual(values = sp_group_color_new)+
  coord_cartesian(clip = 'off')

ggsave(filename = 'methylation_haplotype_tmp.pdf',width = 22,height = 14,device = cairo_pdf)




# Fig. 5c -----------------------------------------------------------------


# Difference among difference groups --------------------------------------

# tumor vs normal

tmpx <- mdata %>% 
  filter(Type!='Control') %>% 
  select(Subject,Type,Pos,Beta) %>% 
  group_by(Pos) %>% 
  do(tidy(wilcox.test(Beta~Type,data=.))) %>% 
  ungroup() %>% 
  select(Pos,p.value) %>% 
  mutate(Type="Tumor",Beta=0)


p1 <- mdata %>% 
  filter(Type!='Control') %>% 
  select(Subject,Type,Pos,Beta) %>% 
  ggplot(aes(as.factor(Pos),Beta,fill=Type))+
  geom_point(pch=21,size=2,position=position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),stroke=0.5,col='black')+
  geom_boxplot(width=0.8,pch=21,outlier.shape = NA,alpha=0.3)+
  scale_fill_d3()+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'Xx')+
  labs(x = "Chromsome 22 coordinate", y = 'Methylation beta value')+
  panel_border(color = 'black',linetype = 1)+
  theme(axis.text.x = element_text(angle = 90,hjust = 0,vjust = 0.5))

p2 <- tmpx %>% 
  ggplot(aes(as.factor(Pos),0,fill=-log10(p.value)))+
  geom_tile(height=0.5)+
  scale_fill_viridis_c()+
  theme_minimal(base_family = 'Roboto Condensed')+
  theme(axis.title = element_blank(),axis.text = element_blank(),legend.position = 'top',panel.grid = element_blank(),legend.key.width = unit(2.5,'cm'))+
  labs(fill='Tumor vs. Normal, -log10(P)\n')

plot_grid(p1,p2,align = 'v',axis = 'lr',rel_heights = c(6,1),ncol = 1)

ggsave(filename = 'methylation_haplotype_tumor_normal.pdf',width = 16,height = 7,device = cairo_pdf)


tmpx <- mdata %>% 
  filter(Type=='Tumor') %>% 
  select(Subject,ID2_Present,Pos,Beta) %>% 
  group_by(Pos) %>% 
  do(tidy(wilcox.test(Beta~ID2_Present,data=.))) %>% 
  ungroup() %>% 
  select(Pos,p.value) %>% 
  arrange(p.value)


p1 <- mdata %>% 
  filter(Type=='Tumor') %>% 
  select(Subject,ID2_Present,Pos,Beta) %>% 
  ggplot(aes(as.factor(Pos),Beta,fill=ID2_Present))+
  geom_point(pch=21,size=2,position=position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),stroke=0.5,col='black')+
  geom_boxplot(width=0.8,pch=21,outlier.shape = NA,alpha=0.3)+
  scale_fill_nejm()+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'Xx')+
  labs(x = "Chromsome 22 coordinate", y = 'Methylation beta value')+
  panel_border(color = 'black',linetype = 1)+
  theme(axis.text.x = element_text(angle = 90,hjust = 0,vjust = 0.5))

p2 <- tmpx %>% 
  ggplot(aes(as.factor(Pos),0,fill=-log10(p.value)))+
  geom_tile(height=0.5)+
  scale_fill_viridis_c()+
  theme_minimal(base_family = 'Roboto Condensed')+
  theme(axis.title = element_blank(),axis.text = element_blank(),legend.position = 'top',panel.grid = element_blank(),legend.key.width = unit(2.5,'cm'))+
  labs(fill='Tumor vs. Normal, -log10(P)\n')

plot_grid(p1,p2,align = 'v',axis = 'lr',rel_heights = c(6,1),ncol = 1)

ggsave(filename = 'methylation_haplotype_tumor_ID2.pdf',width = 16,height = 7,device = cairo_pdf)




## boxplots for the group
my_comparisons <-  list(c('Normal','Tumor (ID2 Present)'),c('Normal','Tumor (ID2 Absent)'),c('Tumor (ID2 Absent)','Tumor (ID2 Present)'))
mdata %>% 
  mutate(Group=case_when(
    Type=='Normal' ~ 'Normal',
    Type=='Tumor' & ID2_Present == 'Present' ~ 'Tumor (ID2 Present)',
    Type=='Tumor' & ID2_Present == 'Absent' ~ 'Tumor (ID2 Absent)',
    TRUE ~ NA
  )) %>% 
  filter(!is.na(Group)) %>% 
  select(Sample,Group,Beta,CpG_ID,Pos) %>% 
  group_by(Group,Sample,) %>% 
  summarise(Beta=median(Beta,na.rm=T)) %>% 
  ungroup() %>% 
  ggplot(aes(Group,Beta,fill=Group))+
  geom_point(aes(fill=Group),pch=21,size=2,position = position_jitter(width = 0.15,height = 0),color="black",stroke = 0.1)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.25,fill=NA)+
  scale_fill_manual(values =pal_jama()(3)[c(3,1,2)])+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+
  labs(x = "", y = 'Methylation beta value')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = 'methylation_group_comparaison.pdf',width = 2.6,height = 7,device = cairo_pdf)

# linear relationship
mdata %>% 
  filter(Type=='Tumor', ID2>0) %>% 
  select(Sample,ID2,Beta,CpG_ID) %>% 
  group_by(Sample,ID2) %>% 
  summarise(Beta=median(Beta,na.rm=T)) %>% 
  ungroup() %>% 
  ggplot(aes(Beta,log2(ID2)))+
  geom_point(pch=21,stroke=0.5,fill=ncicolpal[2],size=2.5)+
  geom_smooth(method = 'lm')+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = 'Y',strip_text_face = 'bold',ticks = T)+
  labs(x='Methylation beta value',y='Number of ID2 deletions (log2)')+
  theme(panel.spacing = unit(0.2,'lines'),strip.text.x = element_text(face = 'plain',hjust = 0.5))+
  coord_cartesian(clip = 'off')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  guides(fill="none")

ggsave(filename = 'methylation_ID2_correlation.pdf',width = 4.5,height = 4,device = cairo_pdf)

mdata %>% 
  filter(Type=='Tumor', ID2>0) %>% 
  select(Sample,ID2,Beta,CpG_ID) %>% 
  group_by(Sample,ID2) %>% 
  summarise(Beta=median(Beta,na.rm=T)) %>% 
  ungroup() %>% 
  do(tidy(cor.test(.$Beta,log2(.$ID2),data=.))) %>% 
  ungroup() %>% 
  arrange(p.value)







# Fig. 5b -----------------------------------------------------------------

# Polar plots 
pdata <- mdata %>% 
  mutate(Group=case_when(
    Type=='Normal' ~ 'Normal',
    Type=='Tumor' & ID2_Present == 'Present' ~ 'Tumor (ID2 Present)',
    Type=='Tumor' & ID2_Present == 'Absent' ~ 'Tumor (ID2 Absent)',
    TRUE ~ NA
  )) %>% 
  filter(!is.na(Group)) %>% 
  select(Sample,Group,Beta,CpG_ID,Pos) %>% 
  group_by(Group,Pos) %>% 
  summarise(Beta=median(Beta,na.rm=T)) %>% 
  ungroup() %>% 
  arrange(Pos) %>% 
  mutate(Pos=fct_inorder(as.factor(comma_format()(Pos))))


pdata %>% 
  ggplot(aes(x = factor(Pos), y = Beta, group = Group,color=Group)) + 
  ylim(0, NA)+
  geom_point(stat = 'identity',size=2) + 
  geom_line()+
  scale_y_continuous(limits = c(0.016,1),trans='log2')+
  coord_polar(start = - pi * 1/24)+
  theme_ipsum_rc(base_size = 14,ticks = FALSE,grid = "Xx")+
  scale_color_manual(values = pal_jama()(3)[c(3,1,2)])+
  theme(axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid = element_line(linetype = 2),axis.text.x = element_text(size = 10),plot.margin = margin(4,4,4,4),legend.position = 'top')+
  geom_hline(yintercept = seq(0,1,by=0.2), color = "gray90",size=0.15) +
  annotate('text', x = 1, y = c(0.0243,0.217,0.475,0.971),label=c('0%','33%','66%','100%'),size=4,col='#542788')
#annotate('text', x = 1, y =seq(0,1,by=0.2),label=seq(0,1,by=0.2),size=4,col='#542788')+
#guides(color="none")

ggsave(filename = 'methylation_haplotype_polar.pdf',width = 5.5,height = 5.5,device = cairo_pdf)






# Fig. 6e -----------------------------------------------------------------

# Methylation data and ZNF695 expression ----------------------------------
load('../../RNASeq/RNASeq_Exp.RData',verbose = T)
load('../../ID2_TE/hq_samples2.RData')

pdata <- mdata %>% 
  filter(Type!='Control') %>% 
  select(Subject,Type,Beta,CpG_ID,Pos) %>% 
  group_by(Subject,Type) %>% 
  summarise(Beta=median(Beta,na.rm=T)) %>% 
  ungroup() %>% 
  left_join(wgs_groups_info %>% select(Subject,Tumor_Barcode)) %>% 
  left_join(
    rdata1 %>% filter(Gene=='ZNF695') %>% select(Tumor_Barcode,Exp,Type=RNAseq_Type)
  )


pdata %>%
  #filter(Tumor_Barcode %in% hq_samples2) %>% 
  group_by(Type) %>% 
  do(tidy(cor.test(.$Exp, .$Beta))) %>% 
  arrange(p.value)

pdata %>% 
  ggplot(aes(Exp,Beta,fill=Type))+
  geom_point(pch=21,size=3)+
  stat_smooth(aes(fill = Type, color = Type), method = "lm") +
  stat_cor(size=5,col=ncicolpal[1]) +
  facet_wrap(~Type,scales = 'free')+
  #stat_regline_equation(label.y = c(60,64),label.x = c(4,4)) +
  scale_fill_manual(values = pal_jama()(3)[c(3,2)]) +
  scale_color_manual(values = pal_jama()(3)[c(3,2)]) +
  scale_y_continuous(breaks = pretty_breaks(n = 5))+
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  #ggh4x::facet_grid2(TL~Subtype,scale='free',independent = 'y')+
  labs(x='ZNF695 RNA-Seq expression log2(CPM)',y='22q12.1 CpG probes\nmedian methylation beta value')+
  theme_ipsum_rc(axis_text_size = 14,axis_title_just = 'm',axis_title_size = 16,grid = 'XY',ticks = T)+
  panel_border(color = 'black')+
  theme(#axis.text.x = element_blank(),
    panel.spacing = unit(0.4,'cm'),
    #axis.ticks.x = element_blank(),
    plot.margin = margin(4,4,4,4),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text = element_text(hjust = 0.5,face = 'bold',size = 13)
  )+
  guides(fill='none',col='none')

ggsave(filename = 'chr22_cpg_methyaltion_znf695_exp.pdf',width = 10,height =4,device = cairo_pdf)  





# Fig. 5d ------------------------------------------------------------------



# Read L1EM result --------------------------------------------------------
l1em <- read_delim('filter_L1HS_FPM.txt',delim = '\t',col_names = T)
l1em <- l1em %>% 
  mutate(tmp=family.category.locus.strand) %>% 
  separate(col = tmp,into = c('family','category','locus','strand'),sep = '\\.') %>% 
  mutate(strand=str_sub(strand,1,1))

# Cytoband 

cytoband <- read_delim('cytoBand.txt.gz',delim = '\t',col_names = F)
colnames(cytoband) <- c('chrom','start','end','cytoband','region')
cytoband <- cytoband %>%mutate(cytoband=paste0(str_remove(chrom,'chr'),cytoband))

tmp <- l1em %>% 
  select(locus) %>% 
  unique() %>% 
  mutate(tmp=locus) %>% 
  separate(col = tmp,into = c('chrom','start','end'),convert = T) %>% 
  select(chrom:end,locus)

tmp <- bed_intersect(tmp,cytoband) %>% 
  select(locus=locus.x,cytoband=cytoband.y) %>% 
  group_by(locus) %>% 
  summarise(cytoband=paste0(cytoband,collapse = '-'))

l1em <- l1em %>% left_join(tmp)

l1em_all <- l1em %>% mutate(L1exp=only + `3prunon`, only_ratio=only/L1exp) %>% left_join(cdata)
l1em <- l1em_all %>% filter(L1exp>2,only_ratio>0.9) 

save(l1em,l1em_all,file='L1EM.RData')


# Landscape of L1EM in tumor and normal in all available dataset-----------------------------------
library(tidyHeatmap)
library(circlize)

tdata <- l1em %>% 
  group_by(Sample, cytoband,dataset,RNAseq_Type) %>% 
  summarise(L1exp=sum(L1exp)) %>% 
  mutate(L1exp=log2(L1exp)) %>% 
  ungroup()

cytoband_list <- tdata %>% count(cytoband,sort=T) %>% filter(n>20) %>% pull(cytoband)

tdata <- tdata %>% filter(cytoband %in% cytoband_list) %>% complete(nesting(Sample,dataset,RNAseq_Type),cytoband,fill = list(L1exp=0))


cairo_pdf(file = 'L1-Exp-RNAseq-all-tumor.pdf',width = 9,height = 10,family = 'Roboto Condensed')
cairo_pdf(file = 'L1-Exp-RNAseq-all-normal.pdf',width = 9,height = 4,family = 'Roboto Condensed')

tdata %>% 
  #filter(RNAseq_Type =='Tumor') %>% 
  filter(RNAseq_Type =='Normal') %>% 
  group_by(dataset) %>% 
  heatmap(Sample,cytoband,L1exp,
          #palette_value = circlize::colorRamp2(c(seq(0,20,4),40,60),c('white',viridis::viridis(7))),
          palette_value = circlize::colorRamp2(seq(0,6),c('white',viridis::viridis(6))),
          rect_gp = grid::gpar(col = "#161616", lwd = 0.1),
          column_names_gp = gpar(fontsize = 10),
          show_row_names=FALSE,column_title='',row_title='',
          heatmap_legend_param=list(at=seq(1,6),labels=2^(seq(1,6)),title= 'LINE-1 RNA (FPM)', title_position = "leftcenter-rot",legend_height=unit(6,"cm")),
          palette_grouping = list(ncicolpal[c(1:2,4)]) 
  )  #annotation_tile(dataset)

dev.off()


# Landscape of 1217 sample only ---------------------------------------------
tdata <- l1em %>% 
  select(-dataset,-RNAseq_Type) %>% 
  filter(Sample %in% cdata3$Sample) %>% 
  left_join(cdata3) %>% 
  left_join(wgs_groups_info %>% select(Tumor_Barcode,SP_Group) %>% left_join(sp_group_data2) %>% select(-SP_Group)) %>% 
  group_by(Sample, cytoband,SP_Group_New,RNAseq_Type) %>% 
  summarise(L1exp=sum(L1exp)) %>% 
  mutate(L1exp=log2(L1exp)) %>% 
  ungroup()

cytoband_list <- tdata %>% count(cytoband,sort=T) %>% filter(n>15) %>% pull(cytoband)

tdata <- tdata %>% filter(cytoband %in% cytoband_list) %>% complete(nesting(Sample,SP_Group_New,RNAseq_Type),cytoband,fill = list(L1exp=0))

tmpsel <- cdata3 %>% filter(Tumor_Barcode %in% hq_samples2) %>% pull(Sample)
library(tidyHeatmap)
library(grid)

#cairo_pdf(file = 'L1-Exp-RNAseq-hq-luad-tumor.pdf',width = 5,height = 8,family = 'Roboto Condensed')
cairo_pdf(file = 'L1-Exp-RNAseq-hq-luad-normal.pdf',width = 5,height = 6,family = 'Roboto Condensed')

tdata %>% 
  #filter(RNAseq_Type =='Tumor') %>% 
  filter(RNAseq_Type =='Normal') %>%
  filter(Sample %in% tmpsel) %>% 
  group_by(SP_Group_New) %>% 
  heatmap(Sample,cytoband,L1exp,
          #palette_value = circlize::colorRamp2(c(seq(0,20,4),40,60),c('white',viridis::viridis(7))),
          palette_value = circlize::colorRamp2(seq(0,6),c('white',viridis::viridis(6))),
          rect_gp = grid::gpar(col = "#161616", lwd = 0.1),
          column_names_gp = gpar(fontsize = 10),
          show_row_names=FALSE,column_title='',row_title='',
          heatmap_legend_param=list(at=seq(1,6),labels=2^(seq(1,6)),title= 'LINE-1 RNA (FPM)', title_position = "leftcenter-rot",legend_height=unit(6,"cm")),
          palette_grouping = list(sp_group_color) 
  )  #annotation_tile(dataset)

dev.off()

## overall 
tdata <- l1em %>% 
  select(-dataset,-RNAseq_Type) %>% 
  filter(Sample %in% cdata3$Sample) %>% 
  left_join(cdata3) %>% 
  left_join(wgs_groups_info %>% select(Tumor_Barcode,SP_Group) %>% left_join(sp_group_data2) %>% select(-SP_Group)) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  group_by(Sample,Tumor_Barcode,SP_Group_New,RNAseq_Type) %>% 
  summarise(L1exp=sum(L1exp)) %>% 
  #mutate(L1exp=log2(L1exp)) %>% 
  ungroup()

tdata %>%
  filter(SP_Group_New != 'Others') %>% 
  ggplot(aes(SP_Group_New,log2(L1exp),fill=RNAseq_Type))+
  geom_boxplot(width=0.6,pch=21,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1.5,position=position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),color="black")+
  scale_fill_manual(values = ncicolpal[c(20,1)])+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  labs(x = "", y = 'LINE-1 RNA (FPM)',fill='Sample Type')+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = 'L1-Exp-RNAseq-hq-luad-overall.pdf',width = 5,height = 5,device = cairo_pdf)

tdata %>% group_by(SP_Group_New) %>% do(tidy(wilcox.test(L1exp~RNAseq_Type,data=.))) %>% ungroup() %>% filter(SP_Group_New != 'Others') %>% select(SP_Group_New,p.value) %>% arrange(p.value) %>% write_clip()

tdata %>% filter(SP_Group_New %in% c('EU_S','AS_N')) %>% group_by(RNAseq_Type) %>% do(tidy(wilcox.test(L1exp~SP_Group_New,data=.))) %>% ungroup()  %>% select(RNAseq_Type,p.value) %>% arrange(p.value) %>% write_clip()
tdata %>% filter(SP_Group_New %in% c('EU_S','EU_N')) %>% group_by(RNAseq_Type) %>% do(tidy(wilcox.test(L1exp~SP_Group_New,data=.))) %>% ungroup()  %>% select(RNAseq_Type,p.value) %>% arrange(p.value) %>% write_clip()


tdata %>% mutate(SP_Group_New=if_else(SP_Group_New == 'AS_N','EU_N',SP_Group_New)) %>% filter(SP_Group_New %in% c('EU_S','EU_N')) %>% group_by(RNAseq_Type) %>% do(tidy(wilcox.test(L1exp~SP_Group_New,data=.))) %>% ungroup()  %>% select(RNAseq_Type,p.value) %>% arrange(p.value) %>% write_clip()


load('../SmokingVar/sdata.RData')

my_comparisons <- list(c("Smoker", "Non-Smoker"))
tdata %>%
  left_join(covdata0) %>% 
  left_join(
    sdata %>% select(Tumor_Barcode,CIGT_CURRENT_STATUS) %>% mutate(CIGT_CURRENT_STATUS = if_else(CIGT_CURRENT_STATUS ==1,'Current Smoker','Former Smoker'))
  ) %>% 
  filter(Smoking != 'Unknown') %>% 
  ggplot(aes(Smoking,log2(L1exp),fill=Smoking))+
  geom_boxplot(width=0.6,pch=21,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1.5,position=position_jitterdodge(jitter.width = 0.5,dodge.width = 0.5),color="black")+
  #scale_fill_manual(values = ncicolpal[c(20,1)])+
  facet_wrap(~RNAseq_Type)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  labs(x = "", y = 'LINE-1 RNA (FPM)',fill='Sample Type')+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)


ggsave(filename = 'L1-Exp-RNAseq-hq-luad-overall_smoking.pdf',width = 6.5,height = 5,device = cairo_pdf)


my_comparisons <- list(c("Current Smoker", "Non-Smoker"),c("Former Smoker", "Non-Smoker"),c("Current Smoker", "Former Smoker"))

tdata %>%
  left_join(covdata0) %>% 
  left_join(
    sdata %>% select(Tumor_Barcode,CIGT_CURRENT_STATUS)
  ) %>% 
  mutate(CIGT_CURRENT_STATUS = if_else(is.na(CIGT_CURRENT_STATUS),Smoking,if_else(CIGT_CURRENT_STATUS ==1,'Current Smoker','Former Smoker'))) %>% 
  filter(Smoking != 'Unknown',CIGT_CURRENT_STATUS != 'Smoker') %>% 
  ggplot(aes(CIGT_CURRENT_STATUS,log2(L1exp),fill=CIGT_CURRENT_STATUS))+
  geom_boxplot(width=0.6,pch=21,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1.5,position=position_jitterdodge(jitter.width = 0.3,dodge.width = 0.5),color="black")+
  #scale_fill_manual(values = ncicolpal[c(20,1)])+
  facet_wrap(~RNAseq_Type)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  labs(x = "", y = 'LINE-1 RNA (FPM)',fill='Sample type')+
  panel_border(color = 'black',linetype = 1)+
  theme(axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1))+
  stat_compare_means(comparisons = my_comparisons)


ggsave(filename = 'L1-Exp-RNAseq-hq-luad-overall_smoking2.pdf',width = 8,height = 6,device = cairo_pdf)


# Fig. 5e -----------------------------------------------------------------


# ID2 
load('id2data.RData')

tdata <- l1em %>% 
  select(-dataset,-RNAseq_Type) %>% 
  filter(Sample %in% cdata3$Sample) %>% 
  left_join(cdata3) %>% 
  left_join(wgs_groups_info %>% select(Tumor_Barcode,SP_Group) %>% left_join(sp_group_data) %>% select(-SP_Group)) %>% 
  group_by(Tumor_Barcode,SP_Group_New,RNAseq_Type) %>% 
  summarise(L1exp=sum(L1exp)) %>% 
  #mutate(L1exp=log2(L1exp)) %>% 
  ungroup() %>% 
  filter(Tumor_Barcode %in% hq_samples2)

tdata <- tdata %>% left_join(id2data) 

tdata %>% filter(SP_Group_New != "Others") %>% group_by(SP_Group_New) %>% do(tidy(wilcox.test(L1exp~ID2_Present,data=.))) %>% ungroup() %>% select(SP_Group_New,p.value) %>% write_clip()
#%>% select(RNAseq_Type,p.value) %>% arrange(p.value) %>% write_clip()

tdata %>% 
  filter(SP_Group_New!='Others') %>% 
  ggplot(aes(SP_Group_New,log2(L1exp),fill=ID2_Present))+
  geom_boxplot(width=0.6,pch=21,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=1.5,position=position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),color="black")+
  scale_fill_manual(values = ncicolpal[c(2:1)])+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  labs(x = "", y = 'LINE-1 RNA (FPM)',fill='ID2 Signature')+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = 'ID2-L1-Exp-RNAseq.pdf',width = 5,height = 6,device = cairo_pdf) 


