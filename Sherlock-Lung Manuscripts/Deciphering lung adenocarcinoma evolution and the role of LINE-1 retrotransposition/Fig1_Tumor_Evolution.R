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
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/sp_group_data.RData')

load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/ID2_TE/id2data.RData')

# load analysis related data set
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/MutationTimeR_HQ_LUAD/Chronological_timing_short.RData',verbose = T)

load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/MutationTimeR_HQ_LUAD/clsdata.RData',verbose = T)

# load function -----------------------------------------------------------
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/ZTW_functions.RData')

## analysis limited to luad only
hq_samples2 <- covdata0 %>% filter(Histology == 'Adenocarcinoma', Tumor_Barcode %in% hq_samples) %>% pull(Tumor_Barcode)
rm(list=c('hq_samples'))




# Fig. 1a ---------------------------------------------------------------

library(ggsankey)
tdata0 <- covdata0 %>% 
  left_join(sp_group_data2) %>% 
  mutate(Data=if_else(Tumor_Barcode %in% hq_samples,'High-quality tumors','Others tumors')) %>% 
  select(Tumor_Barcode,Data,Histology,SP_Group_New) %>% 
  filter(Histology == "Adenocarcinoma", Tumor_Barcode %in% hq_samples)

tdata0 %>% count(Data,Histology,SP_Group_New)

tdata <- tdata0 %>% make_long(Data,SP_Group_New)
tdata_nr <- 
  tdata %>% 
  filter(!is.na(node)) %>% 
  group_by(x, node)%>% 
  summarise(count = n())

tdata <- tdata %>% left_join(tdata_nr)

tdata$node <- factor(tdata$node,levels = c('High-quality tumors','AS_N','EU_N','EU_S','Others'))

ggplot(tdata, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6, node.color = "gray30",) +
  geom_sankey_text(size = 4, color = "black") +
  geom_sankey_text(aes(label = count), size = 3.5, check_overlap = TRUE) +
  scale_fill_manual(values = c( "#ff7f00","#E64B35FF","#f781bf","#14315C","#4DBBD5FF","#00A087FF", "#947100", "#4BBFC6", "#FACE00", "#3C5488FF")) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme_ipsum_rc(axis = FALSE,axis_text_size = 16)+
  theme(axis.text.y = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

ggsave(filename = 'Sherlock_hq_LUAD_WGS_dataset.pdf',width = 7,height = 5,device = cairo_pdf)




# Fig. 1b -----------------------------------------------------------------

# WGD Frequency 
tmp <- BBsolution4 %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  left_join(wgs_groups_info) %>% 
  count(SP_Group,WGD_Status) %>% 
  group_by(SP_Group) %>% 
  mutate(Freq = n/sum(n)) %>% 
  ungroup()

tmp <- BBsolution4 %>%
  left_join(wgs_groups_info) %>%
  left_join(covdata0 %>% select(Tumor_Barcode,Histology)) %>%
  filter(Tumor_Barcode %in% hq_samples, Histology == 'Adenocarcinoma') %>%
  count(SP_Group,WGD_Status) %>%
  group_by(SP_Group) %>%
  mutate(Freq = n/sum(n)) %>%
  ungroup()
# 

tmp %>%
  left_join(sp_group_data) %>% 
  #filter(SP_Group != "Others") %>% 
  ggplot(aes(x=SP_Group_New, y=n, fill = SP_Group_New,alpha=WGD_Status))+
  geom_bar(stat="identity", position ="fill",width = 0.8)+
  geom_text(aes(label=paste0(n,'\n',sprintf("%1.1f", Freq*100),"%")), position=position_fill(vjust=0.5), colour="black",size=3.5)+
  scale_y_continuous(breaks = pretty_breaks(),labels = percent_format(),expand = c(0,0))+
  scale_fill_manual(values = sp_group_color_new)+
  scale_alpha_manual(values = c(0.5,1))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 16,axis_text_size = 12,grid = 'Yy',ticks = FALSE)+
  #panel_border(color = 'black')+
  labs(fill="Group",x="", y = 'Percentage')+
  guides(fill="none")

ggsave(filename = 'wgd_frequency_hq_luad.pdf',width = 5,height = 6,device = cairo_pdf)



# Fig. 1c -----------------------------------------------------------------
load('../DP_info_data.RData')
load('../ZTW_functions.RData')
load('../covdata0.RData')

# all mutations
# tmp <- DP_info_data %>% filter(no.chrs.bearing.mut == 2) %>% count(Tumor_Barcode,name = "CNV2") 
# 
# tmp <- DP_info_data %>% 
#   count(Tumor_Barcode,name='Total') %>% 
#   left_join(tmp)
# tmp <-  tmp %>% replace_na(list(CNV2=0))

# # CpG only
# tmp <- DP_info_data %>% filter(str_detect(mutType,'\\[C>T\\]G'),no.chrs.bearing.mut == 2) %>% count(Tumor_Barcode,name = "CNV2") 
# 
# tmp <- DP_info_data %>% 
#   filter(str_detect(mutType,'\\[C>T\\]G')) %>% 
#   count(Tumor_Barcode,name='Total') %>% 
#   left_join(tmp)

## SBS1 + SBS5
load('../Signature_ludmil3/Mutation_Signature_Probability_SBS.RData',verbose = T)
tmp <- Mutation_Signature_Probability_SBS %>% left_join(DP_info_data %>% select(ID,no.chrs.bearing.mut)) %>% filter(no.chrs.bearing.mut == 2) %>% group_by(Tumor_Barcode) %>% summarise(CNV2=sum(SBS1+SBS5,na.rm=T)) 

tmp <- Mutation_Signature_Probability_SBS %>%
  left_join(DP_info_data %>% select(ID,no.chrs.bearing.mut)) %>% 
  group_by(Tumor_Barcode) %>% 
  summarise(Total=sum(SBS1+SBS5,na.rm=T)) %>% 
  left_join(tmp)



tmp <-  tmp %>% replace_na(list(CNV2=0))

tmp <- BBsolution4 %>%
  select(Tumor_Barcode,WGD_Status) %>% 
  left_join(wgs_groups_info) %>% 
  left_join(tmp) %>% 
  mutate(Ratio=CNV2/Total)

tmp %>% filter(WGD_Status=='WGD') %>% 
  left_join(sp_group_data) %>% 
  filter(SP_Group!="Others") %>% 
  ggplot(aes(Ratio,fill=SP_Group_New,group=SP_Group_New))+
  geom_density(alpha=0.8)+
  scale_x_continuous(breaks = pretty_breaks(),labels = percent_format(),limits = c(0,1))+
  scale_y_continuous(breaks = pretty_breaks())+
  scale_fill_manual(values = sp_group_color_new,breaks = sp_group_major_new)+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid = FALSE,ticks = TRUE)+
  panel_border(color = 'black')+
  labs(fill='Group',x='% mutations with copy number = 2', y = 'Density')

ggsave(filename = 'WGD_CNV2_all.pdf',width = 6,height = 4,device = cairo_pdf )


tmp %>% left_join(covdata0 %>% select(Tumor_Barcode,Histology)) %>% 
  filter(WGD_Status=='WGD',Tumor_Barcode %in% hq_samples,Histology=='Adenocarcinoma') %>% 
  # tmp %>% filter(WGD_Status=='WGD',Tumor_Barcode %in% hq_samples) %>% 
  left_join(sp_group_data) %>% 
  #filter(SP_Group!="Others") %>% 
  ggplot(aes(Ratio,fill=SP_Group_New,group=SP_Group_New))+
  geom_density(alpha=0.8)+
  scale_x_continuous(breaks = pretty_breaks(),labels = percent_format(),limits = c(0,1))+
  scale_y_continuous(breaks = pretty_breaks())+
  scale_fill_manual(values = sp_group_color_new,breaks = sp_group_major_new)+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid = FALSE,ticks = TRUE)+
  panel_border(color = 'black')+
  labs(fill='Group',x='% mutations with copy number = 2', y = 'Density')

tmp %>% left_join(covdata0 %>% select(Tumor_Barcode,Histology)) %>% 
  filter(WGD_Status=='WGD',Tumor_Barcode %in% hq_samples,Histology=='Adenocarcinoma') %>% 
  # tmp %>% filter(WGD_Status=='WGD',Tumor_Barcode %in% hq_samples) %>% 
  left_join(sp_group_data) %>% 
  #filter(SP_Group!="Others") %>% 
  ggplot(aes(Ratio,col=SP_Group_New))+
  geom_density(alpha=1,size=1.6)+
  scale_x_continuous(breaks = pretty_breaks(),labels = percent_format(),limits = c(0,1))+
  scale_y_continuous(breaks = pretty_breaks())+
  #scale_fill_manual(values = sp_group_color_new,breaks = sp_group_major_new)+
  scale_color_manual(values = sp_group_color_new,breaks = names(sp_group_color_new))+
  labs(col='Group')+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid = FALSE,ticks = TRUE)+
  panel_border(color = 'black')+
  labs(fill='Group',x='% mutations with copy number = 2', y = 'Density')

ggsave(filename = 'WGD_CNV2_hq_luad.pdf',width = 6,height = 4,device = cairo_pdf )


# Fig.1d -------------------------------------------------------------------
load('../MutationTimeR_HQ_LUAD/Chronological_timing_short.RData')
clsdata %>%
  filter(!is.na(CLS)) %>%
  select(-Freq) %>%
  pivot_wider(names_from = CLS,values_from = n) %>% 
  mutate(Ratio = (`clonal [early]` + `clonal [late]`)/(`clonal [early]` + `clonal [late]` + `clonal [NA]`)) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  # left_join(sp_group_data2) %>% 
  # group_by(SP_Group_New) %>% 
  summarise(Ratio = median(Ratio, na.rm=T))
  

clsdata %>%
  left_join(wgs_groups_info) %>% 
  left_join(sp_group_data) %>% 
  #filter(SP_Group!="Others") %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  filter(!is.na(CLS)) %>% 
  ggplot(aes(CLS,Freq,fill=SP_Group_New))+
  geom_boxplot()+
  labs(x="", y = "Proportion of total mutations",fill='Group')+
  scale_fill_manual(values = sp_group_color_new,breaks = names(sp_group_color_new))+
  theme_ipsum_rc(axis_text_size = 12,axis_title_just = 'm',axis_title_size = 14)+
  theme(axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1))+
  panel_border(color = 'black')

ggsave(filename = 'CLS_SP_Group_hq_luad.pdf',width = 6,height = 6,device = cairo_pdf)



# for SBS1 only
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/MutationTimeR_HQ_LUAD/CLS2Signature.RData',verbose = T)

tdata <- tdata %>% 
  select(Tumor_Barcode:SBS1) %>% 
  filter(SBS1>0,Tumor_Barcode %in% hq_samples2) %>% 
  group_by(Tumor_Barcode,CLS) %>% 
  summarise(n=sum(SBS1)) %>% 
  ungroup()

tdata <- tdata %>% 
  group_by(Tumor_Barcode) %>% 
  mutate(Freq=n/sum(n))

tdata %>% 
  left_join(sp_group_data2) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  filter(!is.na(CLS)) %>% 
  ggplot(aes(CLS,Freq,fill=SP_Group_New))+
  geom_boxplot()+
  labs(x="", y = "Proportion of total mutations",fill='Group')+
  scale_fill_manual(values = sp_group_color_new,breaks = names(sp_group_color_new))+
  theme_ipsum_rc(axis_text_size = 12,axis_title_just = 'm',axis_title_size = 14)+
  theme(axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1))+
  panel_border(color = 'black')

ggsave(filename = 'CLS_SP_Group_hq_luad_SBS1.pdf',width = 6,height = 6,device = cairo_pdf)



# Fig. 1e -----------------------------------------------------------------

# Driver Gene Frequency by Clonality
load('../MAFtools/sherlock_maf.RData')
load('../MAFtools/sherlock_driver_mutations.RData')
genelist <- read_rds('../../../Collaborators/Nuria/Update2/drivers_intogene.RDS') %>% pull(symbol)

tmp <- sherlock_mtvcf %>% select(Tumor_Barcode, MutationID,CLS) 

#tmp <- sherlock_maf %>% 
tmp <- sherlock_driver_mutations %>% 
  mutate(Chromosome=str_remove(Chromosome,"^chr")) %>% 
  mutate(MutationID=paste(Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2,sep = ':')) %>% 
  select(Tumor_Barcode,MutationID,Hugo_Symbol,Variant_Classification) %>% 
  left_join(tmp)

tmpdata <- tmp %>% 
  left_join(wgs_groups_info %>% select(Tumor_Barcode,SP_Group)) %>% 
  filter(Hugo_Symbol %in% genelist, Tumor_Barcode %in% hq_samples2) 

tmp <- tmpdata %>% 
  count(SP_Group,Hugo_Symbol,CLS)
#filter(SP_Group != 'Others')


tmp <- tmp %>% 
  left_join(
    wgs_groups_info %>% filter(Tumor_Barcode %in% hq_samples2) %>% count(SP_Group,name='size')
  ) %>% 
  mutate(Freq=n/size) %>% 
  arrange(desc(Freq))

genelist <- tmp %>% group_by(Hugo_Symbol) %>% summarise(n=sum(n)) %>% arrange(desc(n)) %>% filter(n>8) %>% pull(Hugo_Symbol)

tmp %>% 
  left_join(sp_group_data) %>% 
  filter(Hugo_Symbol %in% genelist) %>% 
  mutate(Hugo_Symbol=factor(Hugo_Symbol,levels = genelist)) %>% 
  ggplot(aes(Hugo_Symbol,Freq,fill=CLS))+
  geom_bar(position = 'stack',stat = 'identity',size=0.1)+
  facet_wrap(~SP_Group_New,ncol = 1,scales = 'free')+
  scale_fill_manual(values = clscolors)+
  scale_y_continuous(breaks = pretty_breaks(),labels = percent_format())+
  labs(x="",y="Driver mutation freqeuncy\n",fill="")+
  theme_ipsum_rc(axis_text_size = 12,axis_title_just = 'm',axis_title_size = 16)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),panel.spacing = unit(0.1, "lines"))

ggsave(filename = 'Driver_Gene_Clonality.pdf',width = 9,height = 12,device = cairo_pdf)

tmp %>% 
  left_join(sp_group_data) %>% 
  filter(Hugo_Symbol %in% genelist) %>% 
  mutate(Hugo_Symbol=factor(Hugo_Symbol,levels = genelist)) %>% 
  ggplot(aes(Hugo_Symbol,Freq,fill=CLS))+
  geom_bar(position = 'fill',stat = 'identity',size=0.1)+
  facet_wrap(~SP_Group_New,ncol = 1,scales = 'free')+
  scale_fill_manual(values = clscolors)+
  scale_y_continuous(breaks = pretty_breaks(),labels = percent_format())+
  labs(x="",y="Proportion\n",fill="")+
  theme_ipsum_rc(axis_text_size = 12,axis_title_just = 'm',axis_title_size = 16)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),panel.spacing = unit(0.1, "lines"))

ggsave(filename = 'Driver_Gene_Clonality2.pdf',width = 9,height = 12,device = cairo_pdf)


# proportion of clonal mutations vs other mutation by fisher exact test
tmpdata %>% 
  filter(SP_Group %in% c('N_A','N_U')) %>% 
  mutate(CLS=if_else(CLS=='clonal [early]','Yes','No')) %>% 
  select(Tumor_Barcode,Hugo_Symbol,CLS,SP_Group) %>% 
  # unique() %>% 
  # count(Tumor_Barcode,Hugo_Symbol,sort=T) %>% 
  mutate(CLS=as.factor(CLS),SP_Group=as.factor(SP_Group)) %>% 
  group_by(Hugo_Symbol) %>% 
  do(tidy(fisher.test(.$CLS,.$SP_Group))) %>% 
  arrange(p.value)

tmpdata %>% 
  filter(SP_Group %in% c('N_U','S_U')) %>% 
  mutate(CLS=if_else(CLS=='clonal [early]','Yes','No')) %>% 
  select(Tumor_Barcode,Hugo_Symbol,CLS,SP_Group) %>% 
  # unique() %>% 
  # count(Tumor_Barcode,Hugo_Symbol,sort=T) %>% 
  mutate(CLS=as.factor(CLS),SP_Group=as.factor(SP_Group)) %>% 
  group_by(Hugo_Symbol) %>% 
  do(tidy(fisher.test(.$CLS,.$SP_Group))) %>% 
  arrange(p.value)

tresult <- tmpdata %>% 
  filter(SP_Group != 'Others') %>% 
  mutate(SP_Group=if_else(SP_Group %in% c('N_U','N_A'),'Yes','No')) %>% 
  mutate(CLS=if_else(CLS=='clonal [early]','Yes','No')) %>% 
  select(Tumor_Barcode,Hugo_Symbol,CLS,SP_Group) %>% 
  # unique() %>% 
  # count(Tumor_Barcode,Hugo_Symbol,sort=T) %>% 
  mutate(CLS=as.factor(CLS),SP_Group=as.factor(SP_Group)) %>% 
  group_by(Hugo_Symbol) %>% 
  do(tidy(fisher.test(.$CLS,.$SP_Group))) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(Hugo_Symbol %in% genelist) %>% 
  mutate(FDR=p.adjust(p.value))

p.adjust(c(0.001725,tresult$p.value))
myggstyle()
tresult %>% 
  left_join(
    tmp %>% group_by(Hugo_Symbol,CLS) %>% summarise(n=sum(n)) %>% mutate(Ratio=n/sum(n)) %>% filter(CLS=='clonal [early]') %>% arrange(Ratio) %>% ungroup() %>% select(Hugo_Symbol,Ratio)
  ) %>% 
  ggplot(aes(log2(estimate),-log10(p.value)))+
  geom_hline(yintercept = -log10(0.001725),linetype=2,col='#BB0E3D')+
  geom_hline(yintercept = -log10(0.05),linetype=2,col="#4daf4a")+
  geom_vline(xintercept = 0,col='#cccccc',size=0.25)+
  geom_point(aes(size=Ratio,fill=estimate>1),pch=21,stroke=0.5)+#result %>% arrange(FDR) %>% filter(FDR<0.05) %>% tail(n=1) %>% pull(p.value)aes(size=Freq),
  #scale_shape_manual(values = c(21,22,24)) +
  ggrepel::geom_text_repel(data = tresult %>% filter(p.value<0.05),aes(label=Hugo_Symbol))+
  theme_ipsum_rc(base_size = 15,axis_title_just = 'm',axis_title_size = 16,ticks = T)+
  theme(legend.position = 'top',panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),legend.key.width = unit(0.7,'cm'))+
  coord_cartesian(clip = 'off')+
  #scale_fill_manual(values= ncicolpal)+
  scale_fill_aaas()+
  scale_x_continuous(breaks = pretty_breaks(n = 7),limits = c(-5.1,5.1))+
  scale_y_continuous(breaks = pretty_breaks())+
  scale_size_binned(nice.breaks = T,n.breaks = 6)+
  panel_border(color = 'black',size = 0.3)+
  labs(x='log2(Odds ratio)',y='-log10(P-value)',fill='Enrichment',size='Proportion')+ #fill='Early clonal mutation enriched in',size='Proportion of early clonal mutation'
  guides(fill=guide_legend(override.aes = list(size=3.5)))

ggsave(filename = 'Driver_Gene_earlyclone_enrichment.pdf',width = 6,height = 6,device = cairo_pdf)



# Fig. 1f -----------------------------------------------------------------
load('ASCETIC_ID2_group_pattern.RData',verbose = T)

## 
gene_freq1 <- tdata0 %>% select(Tumor_Barcode,GENE_ID,SP_Group_New) %>% unique() %>% count(SP_Group_New,GENE_ID) %>% left_join(
  sp_group_data2 %>% filter(Tumor_Barcode %in% hq_samples2) %>% count(SP_Group_New,name = 'Total')
) %>% 
  mutate(freq=n/Total) %>% 
  #mutate(freq=n) %>% 
  select(gene=GENE_ID,name=SP_Group_New,value=n) %>% 
  pivot_wider(values_fill = 0)


gene_freq2 <- tdata0 %>% select(Tumor_Barcode,GENE_ID,Smoking) %>% unique() %>% count(Smoking,GENE_ID) %>% left_join(
  wgs_groups_info %>% filter(Tumor_Barcode %in% hq_samples2) %>% count(Smoking,name = 'Total')
) %>% 
  mutate(freq=n/Total) %>% 
  #mutate(freq=n) %>% 
  select(gene=GENE_ID,name=Smoking,value=n) %>% 
  pivot_wider(values_fill = 0)

gene_freq3 <- tdata0 %>% select(Tumor_Barcode,GENE_ID) %>% unique() %>% count(GENE_ID) %>% mutate(Total=length(hq_samples2)) %>% 
  mutate(freq=n/Total) %>% 
  #mutate(freq=n) %>% 
  select(gene=GENE_ID,ALL=n) 

gene_freq <- left_join(gene_freq1,gene_freq2) %>% left_join(gene_freq3)


# curve plot for AS_N 
# define the gene coordinate

datatmp <- df_asn %>% slice(1:5) %>% filter(from!='ELF3')

datatmp$from_rank="1"
datatmp$to_rank="2"

# add wildtype
datatmp <- datatmp %>% 
  select(gene=from,rank=from_rank) %>% 
  left_join(gene_freq %>% select(gene,AS_N)) %>% 
  mutate(from='WT',from_rank="0") %>% 
  select(from,to=gene,value=AS_N,from_rank,to_rank=rank) %>% 
  bind_rows(datatmp)



rankdata <- bind_rows(
  datatmp %>% select(gene=from,rank=from_rank),
  datatmp %>% select(gene=to,rank=to_rank)
) %>% 
  unique() %>% 
  mutate(rank=as.integer(rank)) %>% 
  arrange((rank),desc(gene)) %>% 
  mutate(rank=as.integer(as.factor(rank))) %>% 
  arrange(desc(rank),desc(gene))


# define x
datax <- rankdata %>% select(gene,x=rank) %>% mutate(x=as.integer(x))

# define y 
max_gene_rank <- rankdata %>% count(rank,sort=T) %>% slice(1) %>% pull(n)

datay <- rankdata %>% 
  left_join(
    rankdata %>% count(rank,sort=T)  
  ) %>% 
  left_join(gene_freq %>% select(gene,AS_N)) %>% 
  group_by(rank) %>%
  arrange((AS_N)) %>% 
  mutate(seq=seq_along(gene)) %>% 
  ungroup() %>% 
  mutate(y=seq*max_gene_rank/(n+1)) %>% 
  select(gene,y)

datav <- left_join(datax,datay) %>% left_join(gene_freq %>% select(gene,freq='AS_N'))

dataz <- datatmp %>% 
  arrange(value) %>% 
  left_join(datav %>% select(from=gene,x,y)) %>% 
  left_join(datav %>% select(to=gene,x2=x,y2=y)) %>% mutate(x2=0.98*x2,y2=y2) 

myggstyle()

datav %>% 
  ggplot(aes(x,y))+
  geom_segment(data = dataz,aes(x = x,y=y,xend=x2,yend=y2,size=value,col=value),arrow = arrow(length = unit(0.04,'npc')))+
  geom_point(fill='gray',col='black',pch=21,size=6)+
  geom_text(aes(label=gene),nudge_y=-0.08)+
  scale_color_viridis_c(alpha = 0.9,direction = -1)+
  theme_ipsum_rc(grid = FALSE,axis = FALSE)+
  labs(col='Number of Sample\n') +
  guides(size='none')+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),legend.position = 'top',legend.key.width = unit(1.2,'cm'),legend.key.height = unit(0.4,'cm'))

ggsave(filename = 'AS_S_ascetic_infer_new.pdf',width = 4.5,height = 4,device = cairo_pdf)


# curve plot for EU_N 
# define the gene coordinate
datatmp <- df_eun %>% slice(1:5)
datatmp$to_rank[datatmp$to_rank==4]=3

# add wildtype
datatmp <- datatmp %>% 
  select(gene=from,rank=from_rank) %>% 
  left_join(gene_freq %>% select(gene,EU_N)) %>% 
  mutate(from='WT',from_rank="0") %>% 
  select(from,to=gene,value=EU_N,from_rank,to_rank=rank) %>% 
  bind_rows(datatmp)


rankdata <- bind_rows(
  datatmp %>% select(gene=from,rank=from_rank),
  datatmp %>% select(gene=to,rank=to_rank)
) %>% 
  unique() %>% 
  mutate(rank=as.integer(rank)) %>% 
  arrange((rank),desc(gene)) %>% 
  mutate(rank=as.integer(as.factor(rank))) %>% 
  arrange(desc(rank),desc(gene))


# define x
datax <- rankdata %>% select(gene,x=rank) %>% mutate(x=as.integer(x))

# define y 
max_gene_rank <- rankdata %>% count(rank,sort=T) %>% slice(1) %>% pull(n)

datay <- rankdata %>% 
  left_join(
    rankdata %>% count(rank,sort=T)  
  ) %>% 
  left_join(gene_freq %>% select(gene,EU_N)) %>% 
  group_by(rank) %>%
  arrange((EU_N)) %>% 
  mutate(seq=seq_along(gene)) %>% 
  ungroup() %>% 
  mutate(y=seq*max_gene_rank/(n+1)) %>% 
  select(gene,y)

datav <- left_join(datax,datay) %>% left_join(gene_freq %>% select(gene,freq='EU_N'))

dataz <- datatmp %>% 
  arrange(value) %>% 
  left_join(datav %>% select(from=gene,x,y)) %>% 
  left_join(datav %>% select(to=gene,x2=x,y2=y)) %>% mutate(x2=0.975*x2,y2=y2) 

myggstyle()

datav %>% 
  ggplot(aes(x,y))+
  geom_segment(data = dataz,aes(x = x,y=y,xend=0.98*x2,yend=0.98*y2,size=value,col=value),arrow = arrow(length = unit(0.04,'npc')))+
  geom_point(fill='gray',col='black',pch=21,size=6)+
  geom_text(aes(label=gene),nudge_y=-0.2)+
  scale_color_viridis_c(alpha = 0.9,direction = -1)+
  theme_ipsum_rc(grid = FALSE,axis = FALSE)+
  labs(col='Number of Sample\n') +
  guides(size='none')+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),legend.position = 'top',legend.key.width = unit(1.2,'cm'),legend.key.height = unit(0.4,'cm'))

ggsave(filename = 'EU_N_ascetic_infer_new.pdf',width = 4.5,height = 4,device = cairo_pdf)


# curve plot for EUS 
# 
# define the gene coordinate

datatmp <- df_eus %>% filter(value>4)

datatmp$from_rank="1"
datatmp$to_rank="2"

datatmp$from_rank[5]="2"
datatmp$to_rank[5]="3"

# add wildtype
datatmp <- datatmp %>% 
  select(gene=from,rank=from_rank) %>% 
  left_join(gene_freq %>% select(gene,EU_S)) %>% 
  mutate(from='WT',from_rank="0") %>% 
  select(from,to=gene,value=EU_S,from_rank,to_rank=rank) %>% 
  bind_rows(datatmp)


rankdata <- bind_rows(
  datatmp %>% select(gene=from,rank=from_rank),
  datatmp %>% select(gene=to,rank=to_rank)
) %>% 
  unique() %>% 
  arrange(desc(rank),desc(gene))

# define x
datax <- rankdata %>% select(gene,x=rank) %>% mutate(x=as.integer(x))

# define y 
max_gene_rank <- rankdata %>% count(rank,sort=T) %>% slice(1) %>% pull(n)

datay <- rankdata %>% 
  left_join(
    rankdata %>% count(rank,sort=T)  
  ) %>% 
  left_join(gene_freq %>% select(gene,EU_S)) %>% 
  group_by(rank) %>%
  arrange((EU_S)) %>% 
  mutate(seq=seq_along(gene)) %>% 
  ungroup() %>% 
  mutate(y=seq*max_gene_rank/(n+1)) %>% 
  select(gene,y)

datav <- left_join(datax,datay) %>% left_join(gene_freq %>% select(gene,freq='EU_S'))

dataz <- datatmp %>% 
  arrange(value) %>% 
  left_join(datav %>% select(from=gene,x,y)) %>% 
  left_join(datav %>% select(to=gene,x2=x,y2=y)) %>% mutate(x2=0.98*x2,y2=y2) 

myggstyle()

datav %>% 
  ggplot(aes(x,y))+
  geom_segment(data = dataz,aes(x = x,y=y,xend=x2,yend=y2,size=value,col=value),arrow = arrow(length = unit(0.04,'npc')))+
  geom_point(fill='gray',col='black',pch=21,size=6)+
  geom_text(aes(label=gene),nudge_y=-0.15)+
  scale_color_viridis_c(alpha = 0.9,direction = -1)+
  theme_ipsum_rc(grid = FALSE,axis = FALSE)+
  labs(col='Number of Sample\n') +
  guides(size='none')+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),legend.position = 'top',legend.key.width = unit(1.2,'cm'),legend.key.height = unit(0.4,'cm'))

ggsave(filename = 'EU_S_ascetic_infer_new.pdf',width = 4.5,height = 5,device = cairo_pdf)



# curve plot for Others 
# define the gene coordinate
datatmp <- df_others %>% slice(1:5)

datatmp$from_rank="1"
datatmp$to_rank="2"
datatmp$from_rank[3]="2"
datatmp$to_rank[3]="3"

# add wildtype
datatmp <- datatmp %>% 
  select(gene=from,rank=from_rank) %>% 
  left_join(gene_freq %>% select(gene,Others)) %>% 
  mutate(from='WT',from_rank="0") %>% 
  select(from,to=gene,value=Others,from_rank,to_rank=rank) %>% 
  bind_rows(datatmp)

rankdata <- bind_rows(
  datatmp %>% select(gene=from,rank=from_rank),
  datatmp %>% select(gene=to,rank=to_rank)
) %>% 
  unique() %>% 
  mutate(rank=as.integer(rank)) %>% 
  arrange((rank),desc(gene)) %>% 
  mutate(rank=as.integer(as.factor(rank))) %>% 
  arrange(desc(rank),desc(gene))


# define x
datax <- rankdata %>% select(gene,x=rank) %>% mutate(x=as.integer(x))

# define y 
max_gene_rank <- rankdata %>% count(rank,sort=T) %>% slice(1) %>% pull(n)

datay <- rankdata %>% 
  left_join(
    rankdata %>% count(rank,sort=T)  
  ) %>% 
  left_join(gene_freq %>% select(gene,Others)) %>% 
  group_by(rank) %>%
  arrange((Others)) %>% 
  mutate(seq=seq_along(gene)) %>% 
  ungroup() %>% 
  mutate(y=seq*max_gene_rank/(n+1)) %>% 
  select(gene,y)

datav <- left_join(datax,datay) %>% left_join(gene_freq %>% select(gene,freq='Others'))

dataz <- datatmp %>% 
  arrange(value) %>% 
  left_join(datav %>% select(from=gene,x,y)) %>% 
  left_join(datav %>% select(to=gene,x2=x,y2=y)) %>% mutate(x2=0.97*x2,y2=y2) 

myggstyle()

datav %>% 
  ggplot(aes(x,y))+
  geom_segment(data = dataz,aes(x = x,y=y,xend=x2,yend=y2,size=value,col=value),arrow = arrow(length = unit(0.04,'npc')))+
  geom_point(fill='gray',col='black',pch=21,size=6)+
  geom_text(aes(label=gene),nudge_y=-0.2)+
  scale_color_viridis_c(alpha = 0.9,direction = -1)+
  theme_ipsum_rc(grid = FALSE,axis = FALSE)+
  labs(col='Number of Sample\n') +
  guides(size='none')+
  #scale_y_continuous(limits = c(0,1.1))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),legend.position = 'top',legend.key.width = unit(1.2,'cm'),legend.key.height = unit(0.4,'cm'))

ggsave(filename = 'Others_ascetic_infer_new.pdf',width = 4.5,height = 4,device = cairo_pdf)




# Fig. 1g -----------------------------------------------------------------
load('../Signature_ludmil3/Signature_Lumidl.RData')
load('../Signature_ludmil3/Mutation_Signature_Probability_SBS.RData')
load('../Signature_ludmil3/sigcol.RData')

## new limitation to hqsamples for luad only
load('../covdata0.RData')
hq_samples2 <- covdata0 %>% filter(Histology == 'Adenocarcinoma', Tumor_Barcode %in% hq_samples) %>% pull(Tumor_Barcode)
rm(list=c('hq_samples'))


## using ludmil signature analysis

apobec_ratio <- sherlock_variable %>%
  #filter(str_starts(name,'SBS_')) %>%
  #mutate(name=str_remove(name,'^ludmil_')) %>%
  filter(name %in% c('SBS2_ratio','SBS13_ratio')) %>% 
  pivot_wider(names_from = 'name',values_from = value) %>% 
  mutate(APOBEC=SBS2_ratio + SBS13_ratio) %>% 
  select(Tumor_Barcode,APOBEC_Ratio=APOBEC)

# APOBEC + Clonality
# load('../DP_info_data.RData')
# sigprobality <- read_delim('/Volumes/data/NSLC2/Mutations/Signature_1217/result_spextractor_cosmic3.2/NSLC.SBS288.all/SBS288/Suggested_Solution/COSMIC_SBS288_Decomposed_Solution/Activities/Decomposed_Mutation_Probabilities.txt',delim = '\t',col_names = T)
# colnames(sigprobality)[1] <- 'Tumor_Barcode'
# sigprobality <- sigprobality %>% mutate(Tumor_Barcode=str_replace(Tumor_Barcode,"_TPD_01$","_TPD_01.01")) #%>% rename(Tumor_Barcode=`Sample Names`)
# 
# apobec_clone_data <- DP_info_data %>% 
#   select(Tumor_Barcode,ID,MutationTypes=mutType,Clone) %>% 
#   left_join(
#     sigprobality %>% mutate(APOBEC=SBS2+SBS12) %>% select(Tumor_Barcode,MutationTypes,APOBEC)
#   )

apobec_clone_data <- Mutation_Signature_Probability_SBS %>%  mutate(APOBEC=SBS2+SBS12) %>% select(Tumor_Barcode,ID,MutationType,Clone,APOBEC)


apobec_clone_data <- apobec_clone_data %>% 
  group_by(Tumor_Barcode,Clone) %>% 
  summarise(APOBEC=sum(APOBEC,na.rm = T),Total=n_distinct(ID)) %>% 
  ungroup()


tmp <- apobec_clone_data %>% 
  select(-Total) %>% 
  pivot_wider(names_from = 'Clone',values_from = 'APOBEC') %>%
  replace_na(replace = list(`N`=0,Y=0,`NA`=0)) %>% 
  mutate(Subclone_Ratio=N/(N+Y),Enrichment=N/Y,total=N+Y+`NA`) %>% 
  filter(total>100) %>% 
  left_join(wgs_groups_info %>% select(Tumor_Barcode,SP_Group)) %>% 
  #filter(SP_Group!="Others") %>% 
  left_join(apobec_ratio)

# hq sample set
tmp <- tmp %>% filter(Tumor_Barcode %in% hq_samples2) %>% left_join(sp_group_data)


# tmp %>% 
#   select(SP_Group,Tumor_Barcode,`APOBEC Enrichment (Subclonal/Clonal)`=Enrichment,`APOBEC subclonal_Ratio`=Subclone_Ratio) %>% 
#   pivot_longer(cols = -c(SP_Group,Tumor_Barcode)) %>% 
#   ggplot(aes(value,col=SP_Group))+
#   geom_density(size=1)+
#   facet_wrap(~name,scales = "free")+
#   scale_x_continuous(breaks = pretty_breaks(),labels = label_comma())+
#   scale_y_continuous(breaks = pretty_breaks(),labels = label_comma())+
#   theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
#   labs(x="\nAPOBEC subclonal mutation ratio", y = "Density") +
#   panel_border(color = 'black') 
# 
# ggsave(filename = 'APOBEC_Subclonality.pdf',width = 9,height = 5.5,device = cairo_pdf)


tmp %>% 
  #filter(Enrichment>0) %>% 
  #mutate(Enrichment=if_else(Enrichment>5,5,Enrichment)) %>% 
  ggplot(aes(SP_Group_New,log2(Enrichment+1),fill=(APOBEC_Ratio)))+
  #annotation_logticks(sides="l")+
  geom_jitter(position=position_jitter(w=0.2,h=0),size=2.5,shape=21,col="white",stroke=0.2)+ 
  geom_boxplot(fill = NA,outlier.shape=NA,width=0.6)+
  scale_y_continuous(breaks = pretty_breaks(),expand = expansion(mult = c(0.02,0.02)))+
  labs(x="",y="log2(FC +1), FC= APOBEC Subclonal/APOBEC Clonal",fill="APOBEC Ratio")+
  theme_ipsum_rc(base_size = 14,axis_title_size = 14,axis_title_just = 'm',grid = "Y",ticks = TRUE,axis = "")+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  panel_border(size = 0.5,color = "black")+
  scale_fill_viridis_c()+
  geom_hline(yintercept = 1,linetype=2)
#ggsave(filename = 'APOBEC_Subclonality1.pdf',width = 6,height = 9,device = cairo_pdf)
ggsave(filename = 'APOBEC_Subclonality1_hq_ID2paper_tmp.pdf',width = 6,height = 9,device = cairo_pdf)


tmp %>% 
  filter(N>0) %>% 
  #filter(Enrichment>0) %>% 
  #mutate(Enrichment=if_else(Enrichment>5,5,Enrichment)) %>% 
  ggplot(aes(SP_Group_New,log2(Enrichment),fill=(APOBEC_Ratio)))+
  #annotation_logticks(sides="l")+
  geom_jitter(position=position_jitter(w=0.15,h=0.05),size=2.5,shape=21,col="white",stroke=0.2)+ 
  geom_boxplot(fill = NA,outlier.shape=NA,width=0.6)+
  scale_y_continuous(breaks = pretty_breaks(),expand = expansion(mult = c(0.02,0.02)))+
  labs(x="",y="log2(FC), FC=APOBEC Subclonal/APOBEC Clonal",fill="APOBEC Ratio")+
  theme_ipsum_rc(base_size = 14,axis_title_size = 14,axis_title_just = 'm',grid = "Y",ticks = TRUE,axis = "")+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  panel_border(size = 0.5,color = "black")+
  scale_fill_viridis_c()+
  geom_hline(yintercept = 0,linetype=2)

#ggsave(filename = 'APOBEC_Subclonality2.pdf',width = 6,height = 9,device = cairo_pdf)
ggsave(filename = 'APOBEC_Subclonality2_hq_ID2paper_tmp.pdf',width = 6,height = 9,device = cairo_pdf)

# Subclonality by the Ratio
## add signature ratio
mutation_ratio <- sherlock_variable %>%
  # filter(str_starts(name,'ludmil_')) %>%
  # mutate(name=str_remove(name,'^ludmil_')) %>%
  filter(str_starts(name,'SBS'),str_ends(name,'ratio')) %>% 
  pivot_wider(names_from = 'name',values_from = value) %>% 
  rename(APOBEC_ratio = SBS_APOBEC_ratio) %>% 
  pivot_longer(-Tumor_Barcode) %>% 
  mutate(name = str_remove(name,'_ratio')) %>% 
  rename(Mutation_Ratio = value)

tmp <- Mutation_Signature_Probability_SBS %>%  mutate(SBS_APOBEC=SBS2+SBS12) %>% select(Tumor_Barcode,ID,MutationType,Clone,starts_with('SBS'))

# tmp <- DP_info_data %>% 
#   select(Tumor_Barcode,ID,MutationTypes=mutType,Clone) %>% 
#   left_join(
#     sigprobality %>% mutate(SBS_APOBEC=SBS2+SBS12) %>% select(Tumor_Barcode,MutationTypes,contains('SBS'))
#   )
tmp <- tmp %>% filter(!is.na(Clone)) %>% group_by(Tumor_Barcode,Clone) %>% summarise_at(vars(contains('SBS')),sum,na.rm=TRUE) %>% ungroup()

tmp <- tmp %>% pivot_longer(cols = -c(Tumor_Barcode,Clone)) 

## clone vs suclone mutational signature
#excludeSBS <- c('SBS7a','SBS7b','SBS7c','SBS7d','SBS17a','SBS17b','SBS28','SBS44','SBS39','SBS21')
excludeSBS <- c('SBS21','SBS44')
tmp %>% group_by(name) %>% do(tidy(wilcox.test(value~Clone,data=.))) %>% ungroup() %>% arrange(p.value) %>% mutate(FDR=p.adjust(p.value,method = 'BH'))

tmp <- tmp %>% mutate(name=str_remove(name,'SBS_'))
tmp0 <- tmp

tmp <- tmp %>% filter(!(name %in% excludeSBS))

## for the high quality data 
tmp <- tmp %>% filter(Tumor_Barcode %in% hq_samples2)
# tmp2 <- tmp %>% group_by(Tumor_Barcode,Clone) %>% summarise(total=sum(value))

tmp2 <- tmp %>% 
  select(Tumor_Barcode,Clone,name,value) %>%
  group_by(Tumor_Barcode,Clone) %>% 
  summarise(total=sum(value)) %>% 
  ungroup() %>% 
  pivot_wider(names_from="Clone",values_from="total",values_fill = list(value = 0)) %>% 
  dplyr::rename(Ytotal=Y,Ntotal=N)

tmp3 <- tmp %>% select(Tumor_Barcode,Clone,name,value) %>% 
  pivot_wider(names_from="Clone",values_from="value",values_fill = list(value = 0)) %>% 
  left_join(tmp2) %>% 
  left_join(mutation_ratio) %>% 
  # left_join(sherlock_sis_numbers %>% select(Tumor_Barcode,SNV_Num)) %>% 
  filter(Ytotal>100,Ntotal>100,Y>20,N>20,Mutation_Ratio>0.02) %>% 
  mutate(Nratio=(N/Ntotal),Yratio=(Y/Ytotal),Fold=Nratio/Yratio) %>% 
  mutate(name=fct_reorder(name,Fold)) %>% 
  filter(name!='APOBEC')

tmpvalue <- tmp3 %>% group_by(name) %>% do(tidy(wilcox.test(.$Nratio,.$Yratio))) %>% arrange(p.value) %>% ungroup()

source("/Users/zhangt8/NIH-Work/MutationSignature/mSigPortal/CBIIT/mSigPortal/Codes/Sigvisualfunc.R")

SBScolor <- sigcol
# SBScolor['SBS36'] <- 'yellow4'
# SBScolor['APOBEC'] <- "#f781bf"
# SBScolor['SBS-Novel-A'] <- "#984ea3"
# SBScolor['SBS-Novel-C'] <- "#3D4551"

myggstyle()
library(ggnewscale)
tmp3 %>%
  #mutate(Fold=if_else(Fold>4,4,Fold)) %>% 
  ggplot(aes(name,Fold,fill=name))+
  geom_hline(yintercept = 1,linetype=2)+
  #scale_y_continuous(trans = 'log10',limits = c(0.3,5),expand = c(0,0))+
  #annotation_logticks(sides="l")+
  geom_jitter(aes(size=Mutation_Ratio),position=position_jitter(w=0.15,h=0.01),shape=21,col="gray50",stroke=0.4)+ 
  geom_boxplot(fill = NA,outlier.shape=NA,width=0.6)+
  #scale_y_continuous(breaks = pretty_breaks(),expand = expansion(mult = c(0.02,0.02)))+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels=label_log(),
    limits=c(0.1,12)
  ) +
  annotation_logticks(sides = "l",outside = T)+
  labs(x="",y="Fold Change (Subclonal/Clonal)")+
  theme_ipsum_rc(base_size = 14,axis_title_size = 14,axis_title_just = 'm',grid = "Y",ticks = TRUE,axis = "")+
  theme(axis.text.x = element_text(angle = 60,vjust = 1,hjust = 1))+
  panel_border(size = 0.5,color = "black")+
  scale_fill_manual("Signatures",values = SBScolor,breaks = levels(tmp3$name))+
  scale_size_binned(n.breaks = 6,range = c(1,5))+
  guides(fill ="none",size = guide_legend(override.aes = list(fill = 'black')) )+ #guide_legend(override.aes = list(size = 3))
  ggnewscale::new_scale_fill()+
  geom_tile(data=tmpvalue %>% mutate(Fold=0.12,p.value=if_else(p.value>0.05,NA_real_,p.value)),aes(fill=-log10(p.value)),height=.05)+
  scale_fill_viridis_c()+
  coord_cartesian(clip = "off")


#ggsave(filename = 'SBS_enrich_subclone.pdf',width = 10,height = 7,device = cairo_pdf)
ggsave(filename = 'SBS_enrich_subclone_hq_ID2paper_tmp.pdf',width = 9,height = 8,device = cairo_pdf)

load('../sp_group_data.RData')
tmp3 %>%
  left_join(wgs_groups_info %>% select(Tumor_Barcode,SP_Group)) %>% 
  #filter(SP_Group != "Others") %>% 
  left_join(sp_group_data2) %>% 
  #mutate(Fold=if_else(Fold>4,4,Fold)) %>% 
  ggplot(aes(name,log2(Fold+1),fill=name))+
  geom_hline(yintercept = 1,linetype=2)+
  #scale_y_continuous(trans = 'log10',limits = c(0.3,5),expand = c(0,0))+
  #annotation_logticks(sides="l")+
  geom_jitter(position=position_jitter(w=0.15,h=0.01),size=2,shape=21,col="#bbbbbb",stroke=0.2)+ 
  geom_boxplot(fill = NA,outlier.shape=NA,width=0.6)+
  facet_wrap(~SP_Group_New,scales = 'free_y',ncol = 1)+
  scale_y_continuous(breaks = pretty_breaks(),expand = expansion(mult = c(0.02,0.02)))+
  labs(x="",y="Fold Change (Subclonal/Clonal)")+
  theme_ipsum_rc(base_size = 14,axis_title_size = 14,axis_title_just = 'm',grid = "Y",ticks = TRUE,axis = "")+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  panel_border(size = 0.5,color = "black")+
  scale_fill_manual("Signatures",values = SBScolor,breaks = levels(tmp3$name))

#ggsave(filename = 'SBS_enrich_subclone_group.pdf',width = 18,height = 7,device = cairo_pdf)
ggsave(filename = 'SBS_enrich_subclone_group_hq_ID2paper.pdf',width = 9,height = 12,device = cairo_pdf)





# Others figures related to response ----------------------------------------------



# Fig xxx, Mutual exclusivity patterns -------------------------------------
library(data.table)
conflicts_prefer(cowplot::get_legend)
source('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/RDS/Oncoplots_functions.R')
tmp <- read_csv('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/RDS/oncoplot_colors.csv')
landscape_colors <- tmp$Color
names(landscape_colors) <- tmp$Name
landscape_colors <- c(landscape_colors,sp_group_color_new)


genelist <- c('EGFR','RBM10','KRAS','SETD1','TP53','NKX2-1','STK11','CSMD3','CDH10','LRP1B','ERBB4','MED12','RB1','CTNNB1','GNAS','KEAP1')

data_top <- sherlock_data %>% 
  filter(Gene %in% genelist, Type=='Mutation_Driver') %>% 
  mutate(Alteration = str_remove(Alteration,":.*")) %>% 
  filter(Tumor_Barcode %in% hq_samples2)

data_feature <- sp_group_data2 %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  mutate(Type='Genomic_Feature',Gene='Group') %>% 
  left_join(wgs_groups_info) %>% 
  select(Subject,Tumor_Barcode,Gene,Alteration = SP_Group_New,Type) 

result_top <- oncoplot(data = data_top,landscape_colors = landscape_colors,sample_level0 = hq_samples2,GeneSortOnly = TRUE)
result_feature <- oncoplot(data = data_feature,landscape_colors = landscape_colors,sample_level0 = hq_samples2,GeneSortOnly = TRUE)

sample_new_level <- result_top$sample_level %>% 
  left_join(
    result_feature$sample_level
  ) %>% 
  left_join(sp_group_data2) %>% 
  arrange(SP_Group_New,Genomic_Feature,Mutation_Driver) %>% 
  pull(Tumor_Barcode)

genelevel0 <- result_top$gene_level

result_top <- oncoplot(data = data_top,landscape_colors = landscape_colors,sample_level0 = hq_samples2,sample_level = sample_new_level,tmar = -0.2,bmar = -0.05)
result_feature <- oncoplot(data = data_feature,landscape_colors = landscape_colors,sample_level0 = hq_samples2,sample_level = sample_new_level,GeneSortOnly = TRUE,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE,removeAlterbyColor=TRUE)

oncoplot_final <- oncoplot_combined(result_top,result_feature)

save_plot(filename = 'Genome_landscape_group_drivers_tmp.pdf',plot = oncoplot_final,base_height = 3,base_width = 10,device=cairo_pdf)

tmpbarcode <- sp_group_data2 %>% filter(Tumor_Barcode %in% hq_samples2) %>% filter(SP_Group_New == 'Others') %>% pull(Tumor_Barcode)
n=length(tmpbarcode)
result_top <- oncoplot(data = data_top %>% filter(Tumor_Barcode %in% tmpbarcode),landscape_colors = landscape_colors,sample_level0 = tmpbarcode,sample_level = sample_new_level[sample_new_level %in% tmpbarcode],gene_level =genelevel0, tmar = -0.2,bmar = -0.05)
result_feature <- oncoplot(data = data_feature,landscape_colors = landscape_colors,sample_level0 = tmpbarcode,sample_level = sample_new_level[sample_new_level %in% tmpbarcode],gene_level = genelevel0,GeneSortOnly = TRUE,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE,removeAlterbyColor=TRUE)
oncoplot_final_tmp <- oncoplot_combined(result_top,result_feature)

save_plot(filename = 'Genome_landscape_group_drivers_tmp1.pdf',plot = oncoplot_final_tmp,base_height = 3,base_width = 0.025*n,device=cairo_pdf)
save_plot(filename = 'Genome_landscape_group_drivers_tmp2.pdf',plot = oncoplot_final_tmp,base_height = 3,base_width = 0.025*n,device=cairo_pdf)
save_plot(filename = 'Genome_landscape_group_drivers_tmp3.pdf',plot = oncoplot_final_tmp,base_height = 3,base_width = 0.027*n,device=cairo_pdf)
save_plot(filename = 'Genome_landscape_group_drivers_tmp4.pdf',plot = oncoplot_final_tmp,base_height = 3,base_width = 2,device=cairo_pdf)


# mutational exclusivity test
library(ggasym)

tdata <- sherlock_data_full %>% 
  filter(Gene %in% genelist, Type=='Mutation_Driver') %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  left_join(sp_group_data2) %>% 
  #filter(SP_Group_New=='Others') %>% 
  select(Tumor_Barcode, Gene,Alteration) %>% 
  pivot_wider(Tumor_Barcode,names_from = Gene,values_from = Alteration)

# Convert your mutation data to binary format (Yes -> 1, No -> 0)
tdata_binary <- tdata %>%
  mutate(across(-Tumor_Barcode, ~ ifelse(. == "Yes", 1, 0)))

# Get the list of genes (excluding Tumor_Barcode column)
genes <- colnames(tdata_binary)[-1]

# Create a function to perform Fisher's Exact Test and extract p-value and OR
fisher_test <- function(gene1, gene2, data) {
  tryCatch({
    # Create a contingency table for the two genes
    contingency_table <- table(data[[gene1]], data[[gene2]])
    
    # Perform Fisher's Exact Test
    fisher_result <- fisher.test(contingency_table)
    
    # Return the p-value and odds ratio as a tibble
    tibble(
      #gene1 = gene1,
      #gene2 = gene2,
      p_value = fisher_result$p.value,
      odds_ratio = fisher_result$estimate
    )
  }, error = function(e) {
    # In case of error, return NA for both p-value and odds ratio
    tibble(
      #gene1 = gene1,
      #gene2 = gene2,
      p_value = NA,
      odds_ratio = NA
    )
  })
}
fisher_results <- NULL
# Use purrr to apply the fisher_test function to all unique pairs of genes
fisher_results <- expand.grid(gene1 = genes, gene2 = genes) %>%
  filter(gene1 != gene2) %>%             # Exclude comparisons of the same gene
  distinct() %>%
  as_tibble() %>% 
  mutate(result = map2(gene1, gene2, ~ fisher_test(.x, .y, tdata_binary))) %>%
  unnest_wider(result)                   # Use unnest_wider to expand the 'result' column

# View the results with p-values and odds ratios
fisher_results <- fisher_results %>% asymmetrise(gene1,gene2)

fisher_results <- fisher_results %>% mutate(odds_ratio = if_else(p_value>0.05,NA_real_,odds_ratio)) %>% mutate(p_value = if_else(p_value < 0.05,p_value,NA_real_))%>% mutate(odds_ratio = if_else(odds_ratio==0,1/64,odds_ratio)) %>%  mutate(odds_ratio = if_else(is.infinite(odds_ratio),64,odds_ratio))

ggplot(fisher_results, aes(x = gene1, y = gene2)) +
  geom_asymmat(aes(fill_tl = -log10(p_value), fill_br = log2(odds_ratio)),col='white') +
  scale_fill_tl_gradient2(low = "white", high = "tomato",na.value = 'gray90',breaks=pretty_breaks(n=5),limits=c(0,16)) +
  scale_fill_br_gradient2(low = 'purple3',high='orange3',midpoint = 0,mid = 'white',na.value = 'gray90',breaks=pretty_breaks(n=5),limits=c(-6,6)) +
  labs(
    title = "Others",
    fill_tl = "Fisher's exact test\n-log10(p-value)",
    fill_br = "Fisher's exact test\nlog2(OR)"
  ) +
  theme_ipsum_rc()+
  theme(
    #panel.background = element_rect(fill = "grey75"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))+
  panel_border(color = 'black',size = 0.5)

ggsave(filename = 'tmp4.pdf',width = 6.2,height = 5,device = cairo_pdf)







# Response: NRPCC and tumor ploidy, ID2 ,and latency --------------------------------------------------

BBsolution

BBsolution4 %>% 
  filter(NRPCC>10) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  ggplot(aes(NRPCC, BB_Ploidy))+
  geom_point()+
  geom_smooth(method = 'lm')


BBsolution4 %>% 
  filter(NRPCC>10) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  ggplot(aes(WGD_Status, NRPCC))+
  geom_boxplot()

my_comparisons <- list(c("AS_N",'EU_N'),c("AS_N",'EU_S'),c("AS_N",'Others'),c("EU_N",'EU_S'),c("EU_N",'Others'),c("EU_S",'Others'))

BBsolution4 %>% 
  left_join(sp_group_data2) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  ggplot(aes(SP_Group_New,NRPCC,fill=SP_Group_New))+
  geom_boxplot(width=0.6,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=2.5,position = position_jitter(width = 0.15,height = 0),color="white",stroke=0.05)+
  scale_fill_manual(values = sp_group_color_new)+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  facet_wrap(~WGD_Status)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 4),axis.title.x = element_blank())+
  labs(y = 'NRPCC')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = 'NRPCC_groups2.pdf',width = 6,height = 8,device = cairo_pdf)

BBsolution4 %>% 
  left_join(sp_group_data2) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  ggplot(aes(SP_Group_New,BB_Ploidy,fill=SP_Group_New))+
  geom_boxplot(width=0.6,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=2.5,position = position_jitter(width = 0.15,height = 0),color="white",stroke=0.05)+
  scale_fill_manual(values = sp_group_color_new)+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  facet_wrap(~WGD_Status)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 4),axis.title.x = element_blank())+
  labs(y = 'Tumor ploidy')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = 'tumor_ploidy_groups.pdf',width = 4,height = 8,device = cairo_pdf)


# ID2
my_comparisons <- list(c("Present",'Absent'))

BBsolution4 %>% 
  left_join(id2data) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  ggplot(aes(ID2_Present, NRPCC, fill=ID2_Present))+
  geom_boxplot(width=0.6,outlier.shape = NA,alpha =0.25)+
  geom_point(pch=21,size=2.5,position = position_jitter(width = 0.15,height = 0),color="white",stroke=0.05)+
  scale_fill_manual(values = id2color)+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)+
  theme(plot.margin = margin(4,4,4,4),panel.spacing = unit(0.1,'cm'),strip.text.x = element_text(hjust = 0.5,face = 4))+
  labs(x="Mutational signature ID2",y = 'NRPCC')+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_compare_means(comparisons = my_comparisons)

ggsave(filename = 'NRPCC_ID2group.pdf',width = 2.5,height = 6,device = cairo_pdf)

# latency
library(ggpmisc)

BBsolution4 %>% 
  left_join(MRCAdata %>% filter(acceleration == '1x') %>% select(Tumor_Barcode, Latency)) %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  ggplot(aes(NRPCC, Latency))+
  geom_point(pch=21,stroke=0.2,col='white',size=3,fill='gray40')+
  stat_poly_line() +
  stat_poly_eq(use_label(c("R2","p.value")),label.x=0.95,label.y = 0.95)+
  scale_x_continuous(breaks = pretty_breaks(n = 5))+
  scale_y_continuous(breaks = pretty_breaks(n = 5))+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 15,grid="XY",ticks = T)+
  labs(x='NRPCC',y='Tumor latency years')+
  theme(legend.position = 'right', panel.spacing = unit(0.1,"cm"),plot.margin = margin(9,9,4,4),strip.text.x = element_text(hjust = 0.5),strip.text.y = element_text(vjust = 0.5))+
  panel_border(color = 'black',size = 0.5)+coord_cartesian(clip = 'off')

ggsave(filename = 'NRPCC_Latency_correlation.pdf',width = 5,height = 4,device = cairo_pdf)


