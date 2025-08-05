# ------------------------------------------------------------------------------
# Script: Figure 4 - Tumor Evolution Analysis (LUAD cohort)
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
library(hablar)
library(circlize)
library(ggseqlogo)
library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)

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
load('./data/clinical_data.RData',verbose = T)
load('./data/sherlock_data_all.RData',verbose = T)
load('./data/sherlock_variable.RData',verbose = T)
load('./data/sp_group_data.RData',verbose = T)
load('./data/id2data.RData',verbose = T)
load('./data/tedata.RData',verbose = T)
load('./data/RNASeq_Exp.RData',verbose = T)
load('./data/Signature_Lumidl.RData',verbose = T)
load('./data/Signature_Lumidl_CN.RData',verbose = T)
load('./data/Signature_Lumidl_SV.RData',verbose = T)
load('./data/sherlock_profiles.RData',verbose = T)
load('./data/Mutation_Signature_Probability_ID.RData',verbose = T)
load('./data/Mutation_Signature_Probability_SBS.RData',verbose = T)
load('./data/Mutation_Signature_Probability_DBS.RData',verbose = T)
load('./data/Chronological_timing_short.RData',verbose = T)


# load function -----------------------------------------------------------
load('./data/ZTW_functions.RData')
source('./functions/Sherlock_functions.R')

# load analysis related data set
load('./data/trafic.RData')

## analysis limited to luad only
hq_samples2 <- covdata0 %>% filter(Histology == 'Adenocarcinoma', Tumor_Barcode %in% hq_samples) %>% pull(Tumor_Barcode)
rm(list=c('hq_samples'))


# Fig. 4a-b -----------------
#Fig. 4a: This panel illustrates a sample (NSLC-0622-T01) as an example of a tumor harboring L1 insertions from a germline source. In the Circos plot, an arrow indicates the direction from the location of L1 elements in the human germline genome to the position of L1 somatic insertions in the tumor genome. The gray line without an arrow in the Circos plot indicates L1 insertions with an unknown source. L1 retrotranspositions originating from the master L1 on chromosome 22q12.1 (highlighted in blue in the Circos plot) are zoomed in. Different partnered 3’ transductions are highlighted by triangles with various colors on chromosome 22, accompanied by a sequence depth plot (gray bars). A specific partnered 3’ transduction (including the L1 repetitive element and the adjacent unique sequence) between chromosome 22 and chromosome 2 serves as an example of one germline-source L1 retrotransposition.
#Fig. 4b: This section presents another example (tumor sample NSLC-0832-T01) predominantly characterized by somatic-source L1 insertions. In the Circos plot, a green arrow highlights multiple L1 retrotranspositions detected solely in the tumor genome. For improved visualization of these somatic-source L1 insertions, a zoomed-in view specifically focusing on the somatic retrotransposition between chromosome 3 and chromosome 11 is provided in the second Circos plot. Additionally, a specific partnered 3’ transduction serves to elucidate somatic-source L1 retrotransposition. 

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

pdf('./output/NSLC-0622-T01_Somatic_transductions.pdf',height = 5,width = 5,family = 'Roboto Condensed')
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

ggsave(filename = './output/tmp.pdf',width = 8,height = 5)

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

pdf('./output/NSLC-0832-T01_Somatic_transductions.pdf',height = 5,width = 5,family = 'Roboto Condensed')
circos.initializeWithIdeogram()
#circos.initializeWithIdeogram(chromosome.index = c("chr11", "chr3"))
circos.genomicLink(bed1, bed2, col = scols , border = NA,directional = -1,arr.length=0.3)
dev.off()

# for somatic only 
bedall <- bedall %>% filter(SRCTYPE == 'SOMATIC')
bed1 <- bedall %>% select(chr=CHROM,start=POS,value1=SRCTYPE) %>% mutate(end=start,start=start-100000,end=end+100000) %>% select(chr,start,end,value1) %>% as.data.frame()
bed2 <- bedall %>% select(chr=chrom,start,end,value2=SRCTYPE) %>% mutate(end=end,start=start-100000,end=end+100000) %>% select(chr,start,end,value2)%>% as.data.frame()
pdf('./output/NSLC-0832-T01_Somatic_transductions-somatic.pdf',height = 5,width = 5,family = 'Roboto Condensed')
#circos.initializeWithIdeogram()
circos.clear()
circos.initializeWithIdeogram(chromosome.index = c("chr3", "chr11"))
circos.genomicLink(bed1, bed2, col = '#01665e' , border = NA,directional = -1,arr.length=0.2,)
dev.off()


# Fig. 4c -----------------------------------------------------------------
#Distribution of retrotransposable sources of L1 insertions. The bottom pie chart displays the percentage of L1 insertions retrotransposed from germline, somatic, and unknown L1 elements. The top pie chart show the proportion of L1 insertions originating from specific germline L1 masters.
load('./data/trafic.RData')
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

source('./functions/ztw.R')
pdfhr2()
myggstyle()
PieDonut_ztw(trdatatmp %>% filter(SP_Group_New=='AS_N'),aes(pies=Status,count=n),mainCol = ncicolpal[c(3,1,8)],showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.05,pieLabelSize=5,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = './output/AS_N_TE_insertions_hq_luad.pdf',width = 4,height = 4)

PieDonut_ztw(trdatatmp %>% filter(SP_Group_New=='EU_N'),aes(pies=Status,count=n),mainCol = ncicolpal[c(3,1,8)],showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.01,pieLabelSize=5,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = './output/EU_N_TE_insertions_hq_luad.pdf',width = 4,height = 4)

PieDonut_ztw(trdatatmp %>% filter(SP_Group_New=='EU_S'),aes(pies=Status,count=n),mainCol = ncicolpal[c(3,1,8)],showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.01,pieLabelSize=5,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = './output/EU_S_TE_insertions_hq_luad.pdf',width = 4,height = 4)

PieDonut_ztw(trdatatmp %>% filter(SP_Group_New=='Others'),aes(pies=Status,count=n),mainCol = ncicolpal[c(3,1,8)],showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.01,pieLabelSize=5,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = './output/Others_TE_insertions_hq_luad.pdf',width = 4,height = 4)

PieDonut_ztw(trdatatmp %>% filter(SP_Group_New=='ALL'),aes(pies=Status,count=n),mainCol = ncicolpal[c(3,1,8)],showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.01,pieLabelSize=5,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = './output/ALL_TE_insertions_hq_luad.pdf',width = 4,height = 4)


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

ggsave('./output/Germline_L1_master_hq_luad_barplots.pdf',width = 7,height = 4,device = cairo_pdf)

trdatatmp <- trdata %>% 
  mutate(ID=if_else(SRCID %in% tmplevs,SRCID,"Others")) %>% 
  mutate(ID=factor(ID,levels=tmplevs)) %>% 
  count(SP_Group_New,ID,sort=T) %>% 
  mutate(SP_Group_New=factor(SP_Group_New,levels=rev(c('AS_N','EU_N','Others','EU_S','ALL'))))

tmpcolor <- c(pal_d3()(10),'#cccccc')
names(tmpcolor) <- rev(tmplevs)


trdatatmp0 <- trdatatmp %>% filter(SP_Group_New=='AS_N')
PieDonut_ztw(trdatatmp0,aes(pies=ID,count=n),mainCol = tmpcolor,color = 'black',showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.1,pieLabelSize=6,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = './output/tmp1.pdf',width = 4,height = 4)

trdatatmp0 <- trdatatmp %>% filter(SP_Group_New=='EU_N')
PieDonut_ztw(trdatatmp0,aes(pies=ID,count=n),mainCol = tmpcolor,color = 'black',showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.1,pieLabelSize=6,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = './output/tmp2.pdf',width = 4,height = 4)

trdatatmp0 <- trdatatmp %>% filter(SP_Group_New=='EU_S')
PieDonut_ztw(trdatatmp0,aes(pies=ID,count=n),mainCol = tmpcolor,color = 'black',showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.1,pieLabelSize=6,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = './output/tmp3.pdf',width = 4,height = 4)

trdatatmp0 <- trdatatmp %>% filter(SP_Group_New=='Others')
PieDonut_ztw(trdatatmp0,aes(pies=ID,count=n),mainCol = tmpcolor,color = 'black',showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.1,pieLabelSize=6,titlesize = 0,r0=0,title='', showPieName = FALSE,family = 'Roboto Condensed',labelpositionThreshold = 0,ThresholdToLable = T)
ggsave(filename = './output/tmp4.pdf',width = 4,height = 4)



# Fig. 4d -----------------------------------------------------------------
#Enrichment of mutational signatures in tumors with germline source L1 insertions. The horizontal lines indicate the significance threshold FDR < 0.05 in orange and FDR < 0.01 in red. 

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

ggsave(filename = './output/L1_regression_signatures_germline.pdf',width = 6.5,height = 5,device = cairo_pdf())
ggsave(filename = './output/L1_regression_signatures_somatic.pdf',width = 6.5,height = 5,device = cairo_pdf())


# Fig. 4e -----------------------------------------------------------------
#Pearson correlation between deletions attributed to signature ID2 and insertions attributed to signature ID1. Pearson correlation coefficients and corresponding p-values are displayed in the plot.

# ID2_ID1 correlation

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

ggsave(file='./output/ID2_ID1_correlation.pdf',width = 5,height = 4)

ludmil_activity_all %>% 
  filter(Tumor_Barcode %in% hq_samples2) %>% 
  filter(ID2>0,ID1>0) %>% 
  do(tidy(cor.test(log2(.$ID2),log2(.$ID1))))


# Fig. 4f-g -----------------------------------------------------------------
#Mutational signature profiles and motifs for mutational signature ID1 and ID2, respectively
# ID2 sequences motif 

tdata <- Mutation_Signature_Probability_ID %>%
  filter(ID2==1,MutationType=='1:Del:T:5') %>% #,str_detect(ID,"T:")
  select(Tumor_Barcode,ID,Gene_Name,Variant_Classification,AAChange,VAF,CCF,MutationType,ID2) 

## output for motif analysis
tdata %>% 
  filter(Tumor_Barcode %in% hq_samples) %>% 
  select(Tumor_Barcode,ID) %>% 
  separate(ID,into = c('chr','pos','ref','alt'),sep = ':') %>% 
  write_delim('./output/ID2_indels.txt',delim = '\t',col_names = F)


tdata <- Mutation_Signature_Probability_ID %>%
  filter(ID1==1,MutationType=='1:Ins:T:5') %>%
  select(Tumor_Barcode,ID,Gene_Name,Variant_Classification,AAChange,VAF,CCF,MutationType,ID1) 


tdata %>% 
  filter(Tumor_Barcode %in% hq_samples) %>% 
  select(Tumor_Barcode,ID) %>% 
  separate(ID,into = c('chr','pos','ref','alt'),sep = ':') %>% 
  write_delim('./output/ID1_indels.txt',delim = '\t',col_names = F)


tdata <- Mutation_Signature_Probability_DBS %>%
  filter(DBS2==1,MutationType=='CC>AA') %>%
  select(Tumor_Barcode,ID,MutationType,DBS2) 

tdata %>% 
  filter(Tumor_Barcode %in% hq_samples) %>% 
  select(Tumor_Barcode,ID) %>% 
  separate(ID,into = c('chr','pos'),sep = ':') %>% 
  write_delim('./output/DBS2_indels.txt',delim = '\t',col_names = F)

# seqlog by ggseqlogo

id1seq <- read_delim(file = './data/ID1_seq15.txt',delim = '\t',col_names = F) %>% pull(X1)

ggplot()+
  geom_logo(id1seq,font = 'roboto_bold',rev_stack_order = T)+
  theme_logo()+
  theme(axis.text.x = element_text(size = 8))

ggsave(file='./output/id1_motif.pdf',width = 6,height = 2.5,device = cairo_pdf())

id2seq <- read_delim(file = './data/ID2_seq15.txt',delim = '\t',col_names = F) %>% pull(X1)

ggplot()+
  geom_logo(id2seq,font = 'roboto_bold')+
  theme_logo()

ggsave(file='./output/id2_motif.pdf',width = 6,height = 2.5,device = cairo_pdf())


dbs2seq <- read_delim(file = './data/DBS2_seq15.txt',delim = '\t',col_names = F) %>% pull(X1)

ggplot()+
  geom_logo(dbs2seq,font = 'roboto_bold')+
  theme_logo()



# Fig. 4h -----------------------------------------------------------------
#Pearson correlation between deletions attributed to ID2 and total somatic L1 insertions. Pearson correlation coefficients and corresponding p-values are shown in the plot.
load('./data/trafic.RData')
# L1_ID2 
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

ggsave(filename = './output/L1_ID2.pdf',width = 5,height = 4,device = cairo_pdf)

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

ggsave(filename = './output/L1_ID2_germline.pdf',width = 5,height = 4,device = cairo_pdf)


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

ggsave(filename = './output/L1_ID2_somatic.pdf',width = 5,height = 4,device = cairo_pdf)









