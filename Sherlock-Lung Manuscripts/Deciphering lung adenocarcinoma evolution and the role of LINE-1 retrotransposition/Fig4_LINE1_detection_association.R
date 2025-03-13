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
source('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/Sherlock_functions.R')


# load analysis related data set
load('/Users/zhangt8/NIH-Work/EAGLE_NSLC/SecondPaper/Biowulf/trafic.RData')

## analysis limited to luad only
hq_samples2 <- covdata0 %>% filter(Histology == 'Adenocarcinoma', Tumor_Barcode %in% hq_samples) %>% pull(Tumor_Barcode)
rm(list=c('hq_samples'))




# Fig. 4a-b ---------------------------------------------------------------


# example of germline dominat sample
trafic %>% filter(Tumor_Barcode == 'NSLC-0622-T01') %>% count(SRCTYPE,SRCID)
trafic_data <- trafic %>% filter(Tumor_Barcode == 'NSLC-0622-T01') 

#trdata <- trafic_data %>% filter(SRC!=".") %>% select(Tumor_Barcode,SRC,SRCTYPE,SRCID) %>% separate(SRC,c("chrom","start","end")) %>% mutate(start=as.integer(start),end=as.integer(end)) %>% select(chrom,start,end,Tumor_Barcode,SRCTYPE,SRCID)
bedall <- trafic_data  %>% select(CHROM,POS,SRC,SRCTYPE,Tumor_Barcode,SRCID) %>% separate(SRC,c("chrom","start","end"),convert = T) #%>% filter(SRC!='.')
bedall <- bedall %>% mutate(chrom = if_else(chrom=="",CHROM,chrom),start=if_na(start,POS),end=if_na(end,POS))
bed1 <- bedall %>% select(chr=CHROM,start=POS,value1=SRCTYPE) %>% mutate(end=start,start=start-1000000,end=end+1000000) %>% select(chr,start,end,value1) %>% as.data.frame()
bed2 <- bedall %>% select(chr=chrom,start,end,value2=SRCTYPE) %>% mutate(end=end,start=start-1000000,end=end+1000000) %>% select(chr,start,end,value2)%>% as.data.frame()

#scols <- rand_color(2, transparency = 0)
scols <- c(pal_d3(palette = 'category20')(20)[c(1,19,6)],"gray40")
names(scols) <- c('22q12.1','11q21','5q14.1','.')
scols <- scols[bedall$SRCID]

pdf('NSLC-0622-T01_Somatic_transductions.pdf',height = 5,width = 5,family = 'Roboto Condensed')
circos.initializeWithIdeogram()
circos.genomicLink(bed1, bed2, col = scols , border = NA,directional = -1,arr.length=0.3)
dev.off()

trafic_data %>%
  filter(SRC=='chr22_28663284_28669316') %>% 
  mutate(ID=as.character(seq_along(TDC))) %>%
  separate(TDC,into = c('chr','pos1','pos2'),convert = T) %>%
  mutate(pos=(pos1+pos2)/2) %>% 
  filter(pos1<28672000) %>% 
  ggplot(aes(x=pos,y=0.08,fill=ID))+
  geom_point(pch=25,size=5)+
  geom_hline(yintercept = 0,col='#777777')+
  scale_fill_manual(values = c(ncicolpal,pal_nejm()(7),pal_jama()(7)))+
  geom_segment(x=28663284,xend=28669316,y=0,yend=0,linewidth=6)+
  xlim(c(28669279,28670100))+
  ylim(c(0,2))+
  theme_ipsum_rc(axis = F,grid = F, ticks = T)+
  theme(axis.title.x = element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),legend.position = 'none')

ggsave(filename = 'tmp.pdf',width = 8,height = 5)

# example of somatic L1 driven-tumor 
trafic %>% filter(Tumor_Barcode == 'NSLC-0832-T01') %>% count(SRCTYPE,SRCID)
trafic_data <- trafic %>% filter(Tumor_Barcode == 'NSLC-0832-T01') 

#trdata <- trafic_data %>% filter(SRC!=".") %>% select(Tumor_Barcode,SRC,SRCTYPE,SRCID) %>% separate(SRC,c("chrom","start","end")) %>% mutate(start=as.integer(start),end=as.integer(end)) %>% select(chrom,start,end,Tumor_Barcode,SRCTYPE,SRCID)
bedall <- trafic_data  %>% select(CHROM,POS,SRC,SRCTYPE,Tumor_Barcode,SRCID) %>% separate(SRC,c("chrom","start","end"),convert = T) #%>% filter(SRC!='.')
bedall <- bedall %>% mutate(chrom = if_else(chrom=="",CHROM,chrom),start=if_na(start,POS),end=if_na(end,POS))
bedall <- bedall %>% mutate(SRCID = if_else(SRCTYPE == 'SOMATIC','SOMATIC',SRCID))
bed1 <- bedall %>% select(chr=CHROM,start=POS,value1=SRCTYPE) %>% mutate(end=start,start=start-1000000,end=end+1000000) %>% select(chr,start,end,value1) %>% as.data.frame()
bed2 <- bedall %>% select(chr=chrom,start,end,value2=SRCTYPE) %>% mutate(end=end,start=start-1000000,end=end+1000000) %>% select(chr,start,end,value2)%>% as.data.frame()

#scols <- rand_color(2, transparency = 0)
scols <- c(pal_d3(palette = 'category20')(20)[c(4,17)],"#01665e","gray40")
names(scols) <- c('14q23.1','7p14.3','SOMATIC','.')
scols <- scols[bedall$SRCID]

pdf('NSLC-0832-T01_Somatic_transductions.pdf',height = 5,width = 5,family = 'Roboto Condensed')
circos.initializeWithIdeogram()
#circos.initializeWithIdeogram(chromosome.index = c("chr11", "chr3"))
circos.genomicLink(bed1, bed2, col = scols , border = NA,directional = -1,arr.length=0.3)
dev.off()

# for somatic only 
bedall <- bedall %>% filter(SRCTYPE == 'SOMATIC')
bed1 <- bedall %>% select(chr=CHROM,start=POS,value1=SRCTYPE) %>% mutate(end=start,start=start-100000,end=end+100000) %>% select(chr,start,end,value1) %>% as.data.frame()
bed2 <- bedall %>% select(chr=chrom,start,end,value2=SRCTYPE) %>% mutate(end=end,start=start-100000,end=end+100000) %>% select(chr,start,end,value2)%>% as.data.frame()
pdf('NSLC-0832-T01_Somatic_transductions-somatic.pdf',height = 5,width = 5,family = 'Roboto Condensed')
#circos.initializeWithIdeogram()
circos.clear()
circos.initializeWithIdeogram(chromosome.index = c("chr3", "chr11"))
circos.genomicLink(bed1, bed2, col = '#01665e' , border = NA,directional = -1,arr.length=0.2,)
dev.off()




# Fig. 4c -----------------------------------------------------------------

# L1 insertion hotspot

trafic <- trafic %>% mutate(ID=paste(CHROM,POS,REF,sep=':'))
trdata <- trafic %>% select(ID,Tumor_Barcode,SRC,SRCTYPE,SRCID,SRCGENE,STRAND) %>% separate(SRC,c("chrom","start","end"),sep = '_') %>% mutate(start=as.integer(start),end=as.integer(end)) %>% select(ID,Chromosome=chrom,Position1=start,Position2=end,Strand=STRAND,Associated_gene=SRCGENE,Status=SRCTYPE,SRCID,Tumor_Barcode) %>% unique() %>% 
  left_join(sp_group_data2) %>% 
  mutate(Status=if_else(Status=='.','Unknow',Status))


src_teinsertion <- trdata %>% count(Chromosome,Position1,Position2,SRCID,Status,sort=T) %>% filter(Chromosome!="")
#save(src_teinsertion, file='../src_teinsertion.RData')

l1data <- trdata %>% count(Tumor_Barcode,Status) %>% pivot_wider(names_from = Status,values_from = n,values_fill = 0)

# tedata <- BBsolution4 %>% 
#   select(Tumor_Barcode,`Total L1`,`Total Alu`,`Total Others`,`Total Insertion`) %>% 
#   left_join(l1data) %>% 
#   replace_na(list(SOMATIC=0, Unknow=0, GERMLINE=0)) %>% 
#   rename(Total_L1=`Total L1`,Total_Alu=`Total Alu`,Total_Others=`Total Others`,Total_Insertion=`Total Insertion`)
# 
# save(tedata,file='tedata.RData')

trdata <- trdata %>% filter(Tumor_Barcode %in% hq_samples2)

# barplot 

trdatatmp <- bind_rows(
  trdata %>% mutate(SP_Group_New='ALL'),
  trdata #%>% filter(SP_Group_New!='Others')
)

trdatatmp <- trdatatmp %>% count(SP_Group_New,Status)

source('~/NIH-Work/R/ZTW_function/ztw.R')
pdfhr2()
myggstyle()
PieDonut_ztw(trdatatmp %>% filter(SP_Group_New=='AS_N'),aes(pies=Status,count=n),mainCol = ncicolpal[c(3,1,8)],showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.05,pieLabelSize=5,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = 'AS_N_TE_insertions_hq_luad.pdf',width = 4,height = 4)

PieDonut_ztw(trdatatmp %>% filter(SP_Group_New=='EU_N'),aes(pies=Status,count=n),mainCol = ncicolpal[c(3,1,8)],showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.01,pieLabelSize=5,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = 'EU_N_TE_insertions_hq_luad.pdf',width = 4,height = 4)

PieDonut_ztw(trdatatmp %>% filter(SP_Group_New=='EU_S'),aes(pies=Status,count=n),mainCol = ncicolpal[c(3,1,8)],showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.01,pieLabelSize=5,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = 'EU_S_TE_insertions_hq_luad.pdf',width = 4,height = 4)

PieDonut_ztw(trdatatmp %>% filter(SP_Group_New=='Others'),aes(pies=Status,count=n),mainCol = ncicolpal[c(3,1,8)],showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.01,pieLabelSize=5,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = 'Others_TE_insertions_hq_luad.pdf',width = 4,height = 4)

PieDonut_ztw(trdatatmp %>% filter(SP_Group_New=='ALL'),aes(pies=Status,count=n),mainCol = ncicolpal[c(3,1,8)],showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.01,pieLabelSize=5,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = 'ALL_TE_insertions_hq_luad.pdf',width = 4,height = 4)


tmplevs <- trdata %>% 
  filter(Status == 'GERMLINE') %>% 
  count(SRCID,sort=T) %>% 
  mutate(Seq=seq_along(SRCID)) %>% 
  mutate(ID=if_else(Seq>10,"Others",SRCID)) %>% 
  mutate(ID=fct_rev(fct_inorder(ID))) %>% 
  group_by(ID) %>% 
  summarise(n=sum(n)) %>% 
  pull(ID) %>% 
  levels()

trdata <- bind_rows(
  trdata %>% filter(Status == 'GERMLINE') %>% mutate(SP_Group_New='ALL'),
  trdata %>% filter(Status == 'GERMLINE') #%>% filter(SP_Group_New!='Others')
)


trdata %>% 
  mutate(ID=if_else(SRCID %in% tmplevs,SRCID,"Others")) %>% 
  mutate(ID=factor(ID,levels=tmplevs)) %>% 
  count(SP_Group_New,ID,sort=T) %>% 
  mutate(SP_Group_New=factor(SP_Group_New,levels=rev(c('AS_N','EU_N','Others','EU_S','ALL')))) %>% 
  ggplot(aes(SP_Group_New,n,fill=ID))+
  geom_col(width = 0.8)+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  labs(x="",y="Number of TE insertions",fill='Germline L1 Master')+
  scale_fill_manual(values = (c(pal_d3()(10),'#cccccc')),breaks = rev(tmplevs))+
  theme_ipsum_rc(grid = FALSE,ticks = TRUE,axis_title_just = 'm',axis_title_size = 14)+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(nrow = 2))+
  coord_flip()

ggsave('Germline_L1_master_hq_luad_barplots.pdf',width = 7,height = 4,device = cairo_pdf)


trdatatmp <- trdata %>% 
  mutate(ID=if_else(SRCID %in% tmplevs,SRCID,"Others")) %>% 
  mutate(ID=factor(ID,levels=tmplevs)) %>% 
  count(SP_Group_New,ID,sort=T) %>% 
  mutate(SP_Group_New=factor(SP_Group_New,levels=rev(c('AS_N','EU_N','Others','EU_S','ALL'))))

tmpcolor <- c(pal_d3()(10),'#cccccc')
names(tmpcolor) <- rev(tmplevs)


trdatatmp0 <- trdatatmp %>% filter(SP_Group_New=='AS_N')
PieDonut_ztw(trdatatmp0,aes(pies=ID,count=n),mainCol = tmpcolor,color = 'black',showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.1,pieLabelSize=6,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = 'tmp1.pdf',width = 4,height = 4)

trdatatmp0 <- trdatatmp %>% filter(SP_Group_New=='EU_N')
PieDonut_ztw(trdatatmp0,aes(pies=ID,count=n),mainCol = tmpcolor,color = 'black',showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.1,pieLabelSize=6,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = 'tmp2.pdf',width = 4,height = 4)

trdatatmp0 <- trdatatmp %>% filter(SP_Group_New=='EU_S')
PieDonut_ztw(trdatatmp0,aes(pies=ID,count=n),mainCol = tmpcolor,color = 'black',showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.1,pieLabelSize=6,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = 'tmp3.pdf',width = 4,height = 4)

trdatatmp0 <- trdatatmp %>% filter(SP_Group_New=='Others')
PieDonut_ztw(trdatatmp0,aes(pies=ID,count=n),mainCol = tmpcolor,color = 'black',showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.1,pieLabelSize=6,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = 'tmp4.pdf',width = 4,height = 4)



# Fig. 4d -----------------------------------------------------------------
# L1 insertion vs all signautres
tdata <- sherlock_data_full %>% 
  filter(Type=='Signature_Cosmic_final') %>% 
  select(Tumor_Barcode,Gene,Signature=Alteration) %>% 
  left_join(
    sherlock_data_full %>% 
      filter(Gene=='All_L1_Germline') %>% 
      select(Tumor_Barcode,L1=Alteration) 
  ) %>% 
  filter(Tumor_Barcode %in% hq_samples2,!str_detect(Gene,'APOBEC'))

# 
# tresult <- tdata %>% 
#   group_by(Gene) %>% 
#   do(tresult = safely(fisher.test)(.$L1,.$Signature)) %>% 
#   mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
#   filter(!tresult_null) %>% 
#   mutate(fit = list(tidy(tresult[['result']]))) %>% 
#   select(Gene,fit) %>% 
#   unnest(cols = c(fit)) %>% 
#   ungroup() %>% 
#   arrange(p.value) %>% 
#   mutate(FDR=p.adjust(p.value))


tresult <- tdata %>%
  left_join(covdata0) %>%
  mutate(Signature=as.factor(Signature),L1=as.factor(L1)) %>%
  group_by(Gene) %>%
  do(tresult = safely(glm)(Signature ~ L1+Smoking+Age+Gender+Assigned_Population+Tumor_Purity,family='binomial',data=.)) %>% #Histology
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>%
  filter(!tresult_null) %>%
  mutate(fit = list(tidy(tresult[['result']], exponentiate =F))) %>%
  select(Gene,fit) %>%
  unnest(cols = c(fit)) %>%
  ungroup() %>%
  arrange(p.value) %>%
  filter(term=='L1Yes') %>%
  filter(abs(estimate)<10) %>% 
  mutate(FDR=p.adjust(p.value))


tresult %>% 
  #filter(abs(estimate)<1) %>% 
  mutate(Profile=if_else(str_detect(Gene,'SBS'),'SBS',if_else(str_detect(Gene,'ID'),'ID',if_else(str_detect(Gene,'DBS'),'DBS',NA_character_)))) %>% 
  mutate(Profile=factor(Profile,levels=c('SBS','ID','DBS'))) %>% 
  ggplot(aes((estimate),-log10(FDR)))+
  geom_hline(yintercept = -log10(0.01),linetype=2,col='red',size=0.5)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col='orange',size=0.5)+
  # geom_vline(xintercept = c(-5,5),linetype=2,col='blue',size=0.5)+
  geom_point(aes(fill=Profile),pch=21,stroke=0.1,col='black',size=4)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_y_continuous(breaks = pretty_breaks())+
  scale_fill_manual(values = ncicolpal[c(1,3,4)])+
  ggrepel::geom_text_repel(data=tresult %>% filter(FDR<0.2),aes(label=Gene),max.overlaps = 30,size=4)+
  labs(x='Regression coefficient',y='-log10(FDR)',fill='Signature')+
  guides(fill = guide_legend(override.aes = list(size=5)))+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid='XY',ticks = T)+
  panel_border(color = 'black',size = 0.5)
#coord_cartesian(clip = 'off')

ggsave(filename = 'L1_regression_signatures_germline.pdf',width = 6.5,height = 5,device = cairo_pdf())
ggsave(filename = 'L1_regression_signatures_somatic.pdf',width = 6.5,height = 5,device = cairo_pdf())


# Fig. 4e -----------------------------------------------------------------

# ID2_ID1 correlation
load('../Signature_ludmil3/Signature_Lumidl.RData')
load('hq_samples2.RData')

#myggstyle()
ludmil_activity_all %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  filter(ID2>0,ID1>0) %>% 
  ggplot(aes(log2(ID2),log2(ID1)))+
  geom_point(pch=21,fill=ncicolpal[8],size=3,col='white',stroke=0.1)+
  geom_smooth(method = 'lm')+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_y_continuous(breaks = pretty_breaks())+
  labs(x='ID2 mutations (log2)',y='ID1 mutations (log2)')+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,ticks = T,grid = FALSE,base_family = 'Helvetica') +
  panel_border(color = 'black')

ggsave(file='ID2_ID1_correlation.pdf',width = 5,height = 4)

ludmil_activity_all %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  filter(ID2>0,ID1>0) %>% 
  do(tidy(cor.test(log2(.$ID2),log2(.$ID1))))


# Fig. 4f-g -----------------------------------------------------------------
# ID2 sequences motif 
load('../Signature_ludmil/Mutation_Signature_Probability_ID.RData')
load('../Signature_ludmil/Mutation_Signature_Probability.RData')
load('../Signature_ludmil/Mutation_Signature_Probability_DBS.RData')

tdata <- Mutation_Signature_Probability_ID %>%
  filter(ID2==1,MutationsType=='1:Del:T:5') %>% #,str_detect(ID,"T:")
  select(Tumor_Barcode,ID,Gene_Name,Variant_Classification,AAChange,VAF,CCF,CLASS,MutationsType,ID2) 


tdata %>% 
  filter(Tumor_Barcode %in% hq_samples) %>% 
  select(Tumor_Barcode,ID) %>% 
  separate(ID,into = c('chr','pos','ref','alt'),sep = ':') %>% 
  write_delim('Motif/ID2_indels.txt',delim = '\t',col_names = F)


tdata <- Mutation_Signature_Probability_ID %>%
  filter(ID1==1,MutationsType=='1:Ins:T:5') %>%
  select(Tumor_Barcode,ID,Gene_Name,Variant_Classification,AAChange,VAF,CCF,CLASS,MutationsType,ID1) 


tdata %>% 
  filter(Tumor_Barcode %in% hq_samples) %>% 
  select(Tumor_Barcode,ID) %>% 
  separate(ID,into = c('chr','pos','ref','alt'),sep = ':') %>% 
  write_delim('Motif/ID1_indels.txt',delim = '\t',col_names = F)



tdata <- Mutation_Signature_Probability_DBS %>%
  filter(DBS2==1,MutationsType=='CC>AA') %>%
  select(Tumor_Barcode,ID,MutationsType,DBS2) 

tdata %>% 
  filter(Tumor_Barcode %in% hq_samples) %>% 
  select(Tumor_Barcode,ID) %>% 
  separate(ID,into = c('chr','pos'),sep = ':') %>% 
  write_delim('Motif/DBS2_indels.txt',delim = '\t',col_names = F)


# seqlog by ggseqlogo

id1seq <- read_delim(file = 'Motif/ID1_seq15.txt',delim = '\t',col_names = F) %>% pull(X1)

ggplot()+
  geom_logo(id1seq,font = 'roboto_bold',rev_stack_order = T)+
  theme_logo()+
  theme(axis.text.x = element_text(size = 8))

ggsave(file='id1_motif.pdf',width = 6,height = 2.5,device = cairo_pdf())

id2seq <- read_delim(file = 'Motif/ID2_seq15.txt',delim = '\t',col_names = F) %>% pull(X1)

ggplot()+
  geom_logo(id2seq,font = 'roboto_bold')+
  theme_logo()

ggsave(file='id2_motif.pdf',width = 6,height = 2.5,device = cairo_pdf())


dbs2seq <- read_delim(file = 'Motif/DBS2_seq15.txt',delim = '\t',col_names = F) %>% pull(X1)

ggplot()+
  geom_logo(dbs2seq,font = 'roboto_bold')+
  theme_logo()



# Fig. 4h -----------------------------------------------------------------
# L1_ID2 
load('id2data.RData')
load('../RDS/sherlock_variable.RData')
load('../MutationTimeR/Chronological_timing_short.RData')
load('../sp_group_data.RData')

load('../trafic.RData')
load('../sp_group_data.RData')
trafic <- trafic %>% filter(CLASS=='L1')%>% mutate(ID=paste(CHROM,POS,REF,sep=':')) %>% mutate(SRCTYPE=if_else(SRCTYPE =='.','Unknow',SRCTYPE))
trdata <- trafic %>% select(ID,Tumor_Barcode,SRC,SRCTYPE,SRCID,SRCGENE,STRAND) %>% separate(SRC,c("chrom","start","end")) %>% mutate(start=as.integer(start),end=as.integer(end)) %>% select(ID,Chromosome=chrom,Position1=start,Position2=end,Strand=STRAND,Associated_gene=SRCGENE,Status=SRCTYPE,SRCID,Tumor_Barcode) %>% unique() %>% 
  left_join(sp_group_data2)

trdata <- trdata %>% filter(Tumor_Barcode %in% hq_samples2)

tmp <- trdata %>% count(Tumor_Barcode,Status) %>% pivot_wider(names_from = Status,values_from = n,values_fill = 0)



# tdata <- id2data %>% 
#   left_join(
#     sherlock_variable %>% filter(name == 'Total L1') %>% pivot_wider()
#   ) %>% 
#   left_join(sp_group_data2) %>% 
#   left_join(tmp) %>% 
#   filter(Tumor_Barcode %in% hq_samples)
#filter(ID2>0,`Total L1`>0)

tdata <- trdata %>% mutate(Status='ALL') %>% 
  bind_rows(trdata) %>% 
  count(Tumor_Barcode,Status) %>% 
  pivot_wider(names_from = Status,values_from = n,values_fill = 0) %>% 
  left_join(id2data) %>% 
  left_join(sp_group_data2)


tdata %>% 
  filter(ID2>0,`ALL`>0) %>% 
  ggplot(aes(log2(`ALL`),log2(ID2)))+
  geom_point(pch=21,size=2.5,fill=id2color[2])+
  geom_smooth(method = 'lm')+
  scale_x_continuous(breaks = pretty_breaks(n = 8))+
  scale_y_continuous(breaks = pretty_breaks(n = 8))+
  labs(x="Total L1 insertions (log2)", y= "ID2 deletions (log2)" )+
  theme_ipsum_rc(base_size = 12,plot_title_size = 12,axis_title_size = 16,axis_title_just = "m",grid = FALSE,axis = "XY",ticks = TRUE)+
  panel_border(color = 'black')

ggsave(filename = 'L1_ID2.pdf',width = 5,height = 4,device = cairo_pdf)

tdata %>% 
  filter(ID2>0,`GERMLINE`>0) %>% 
  ggplot(aes(log2(`GERMLINE`),log2(ID2)))+
  geom_point(pch=21,size=2.5,fill=id2color[2])+
  geom_smooth(method = 'lm')+
  scale_x_continuous(breaks = pretty_breaks(n = 8))+
  scale_y_continuous(breaks = pretty_breaks(n = 8))+
  labs(x="L1 insertions due to germline L1 retrotranspostion (log2)", y= "ID2 deletions (log2)" )+
  theme_ipsum_rc(base_size = 12,plot_title_size = 12,axis_title_size = 16,axis_title_just = "m",grid = FALSE,axis = "XY",ticks = TRUE)+
  panel_border(color = 'black')

ggsave(filename = 'L1_ID2_germline.pdf',width = 5,height = 4,device = cairo_pdf)


tdata %>% 
  filter(ID2>0,`SOMATIC`>0) %>% 
  ggplot(aes(log2(`SOMATIC`),log2(ID2)))+
  geom_point(pch=21,size=2.5,fill=id2color[2])+
  geom_smooth(method = 'lm')+
  scale_x_continuous(breaks = pretty_breaks(n = 8))+
  scale_y_continuous(breaks = pretty_breaks(n = 8))+
  labs(x="L1 insertions due to germline L1 retrotranspostion (log2)", y= "ID2 deletions (log2)" )+
  theme_ipsum_rc(base_size = 12,plot_title_size = 12,axis_title_size = 16,axis_title_just = "m",grid = FALSE,axis = "XY",ticks = TRUE)+
  panel_border(color = 'black')

ggsave(filename = 'L1_ID2_somatic.pdf',width = 5,height = 4,device = cairo_pdf)









