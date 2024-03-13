set_wd()
libztw()
pdfhr2()
myggstyle()

source('/Users/zhangt8/NIH-Work/MutationSignature/mSigPortal/CBIIT/mSigPortal/Codes/Sigvisualfunc.R')
source('~/NIH-Work/R/ZTW_function/ztw.R')


# load common data files 
load('../../BBsolution_final3_short.RData')
load('../../sp_group_data.RData')
load('../../covdata0.RData')
load('../../Clinical/clinical_data.RData')
load('../../ZTW_functions.RData')
load('../../RDS/sherlock_data_all.RData')
load('../../RDS/sherlock_variable.RData')
load('../../MAFtools/sherlock_driver_mutations.RData')
load('../../sherlock_PGA.RData')

load('../Signature_Lumidl.RData',verbose = T)
load('../sherlock_profiles.RData',verbose = T)
load('../Signature_Lumidl_CN.RData',verbose = T)
load('../Signature_Lumidl_SV.RData',verbose = T)

load('../sampleset.RData')
load('../sigcol.RData')

drglist <- readRDS('../../../../Collaborators/Nuria/Update2/drivers_intogene.RDS') %>% pull(symbol)
drglist <- sherlock_driver_mutations %>% filter(Hugo_Symbol %in% drglist,Tumor_Barcode %in% sherlock_nonsmoker) %>% select(Hugo_Symbol,Tumor_Barcode) %>% unique() %>% count(Hugo_Symbol) %>% filter(n>2) %>% pull(Hugo_Symbol)


#The presence of a specific signature was considered for the enrichment analysis (or the dichotomization of values above and below the median of assigned mutations, if the signature was present in above 50% of the samples). 
ludmil_activity_all <- ludmil_activity_all %>% left_join(cn68_activity) %>% left_join(sv38_activity) %>% mutate_all(~replace_na(., 0L))

ludmil_activity_all <- ludmil_activity_all %>% filter(Tumor_Barcode %in% sherlock_nonsmoker)
ludmil_activity_all_obs1 <- ludmil_activity_all %>% mutate(across(where(is.numeric),~ . > 0))
tmpcolnames <- colnames(ludmil_activity_all_obs1)
excludesigs <- ludmil_activity_all_obs1 %>% pivot_longer(-Tumor_Barcode) %>% mutate(value=if_else(is.na(value),FALSE,value)) %>% count(name,value) %>% filter(value) %>% arrange(n) %>% filter(n==1) %>% pull(name)
tmpcolnames <- tmpcolnames[!tmpcolnames %in% excludesigs]

siglist1 <- ludmil_activity_all_obs1 %>% pivot_longer(-Tumor_Barcode) %>% mutate(value=if_else(is.na(value),FALSE,value)) %>% count(name,value) %>% group_by(name) %>% mutate(freq=n/sum(n)) %>% ungroup() %>% filter(value,freq<0.5) %>% pull(name)

ludmil_activity_all_obs2 <- ludmil_activity_all %>% pivot_longer(-Tumor_Barcode) %>% group_by(name) %>% mutate(mvalue=median(value,na.rm=T)) %>% ungroup() %>% mutate(value2=if_else(value>=mvalue,TRUE,FALSE)) %>% select(-value,-mvalue) %>% rename(value=value2) %>% filter(!(name %in% siglist1)) %>% pivot_wider()

ludmil_activity_all_obs <- left_join(
  ludmil_activity_all_obs1 %>% select(Tumor_Barcode,one_of(siglist1)),
  ludmil_activity_all_obs2
) %>% 
  select(one_of(tmpcolnames))

ludmil_activity_all_obs





# load pollution data
pollution_data <- read_csv('../../Exposure_data/Sherlock_Non-Smoker_PM2.5_groups_2023OCT6.csv',col_names = T)

pollution_data <- pollution_data %>% 
  select(Tumor_Barcode,Country, Country_pollution,PM25=Population.Weighted.PM2.5.ug.m3,Pollution_group2=Pollution_group) %>% 
  mutate(Pollution_group3=case_when(
    PM25<=10.01 ~ 'Low',
    PM25>10.01 & PM25<24.57  ~ 'Intermediate',
    PM25>=24.57 ~ 'High',
    TRUE ~ NA_character_
  )) %>% 
  mutate(Pollution_group2 = factor(Pollution_group2,levels=c('Low','High'))) %>% 
  mutate(Pollution_group3 = factor(Pollution_group3,levels=c('Low','Intermediate','High')))

pollution_colors <- c('#01665e','#947100','#BB0E3D')
names(pollution_colors) <- c('Low','Intermediate','High')




# ED_Figure 1 -------------------------------------------------------------

# 1a manual combined mulitple signature plot together


# 1b ----------------------------------------------------------------------

sbs288_denovo_mapping_detail_nonsmoker %>% 
  mutate(cosmic_signature = factor(cosmic_signature,levels=names(sigcol))) %>% 
  ggplot(aes(de_novo_extracted,contribution,fill=cosmic_signature))+
  geom_col( position = 'fill')+
  scale_fill_manual(values = sigcol)+
  geom_text(aes(label = paste0(cosmic_signature,"\n",percent_format(accuracy = 0.1)(contribution))), position = position_fill(vjust = .5),size=3)+
  theme_ipsum_rc(grid = FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank(),axis.text.x = element_text(angle = 45))+
  guides(fill = guide_legend(title = 'Cosmic Signatures',ncol = 1,title.position = 'top'))

ggsave(filename = 'Denovo_sig_mapping.pdf',width = 7.5,height = 6,device = cairo_pdf)


# 1c ----------------------------------------------------------------------
#genomesize <- genome2size('GRCh38')
tdata <- sbs288_denovo_activity_nonsmoker %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker) %>%
  pivot_longer(-Tumor_Barcode) %>% 
  left_join(covdata0 %>% select(Tumor_Barcode,Histology) %>% mutate(Histology = as.character(Histology))) %>% 
  arrange(name,Tumor_Barcode) %>% 
  select(Cancer_Type=name,Sample=Tumor_Barcode,Sample_Group=Histology,Burden=value) %>% 
  mutate(Burden=log10(Burden))

sgroupcol <- ncicolpal[1:4]
names(sgroupcol) <- c('Squamous cell carcinoma','Other','Carcinoid tumor','Adenocarcinoma')
p <- TMBplot2(tdata,sgroupcol = sgroupcol,ylab = 'Number of SBS mutations (log10)')+ theme(axis.text.x.bottom = element_text(size = 10,angle = 45,hjust = 1,vjust = 1))
ggsave(file='Denovo_sig_TMB.pdf',width = 5,height = 5,device = cairo_pdf)


# ED_Figure 2 -------------------------------------------------------------

# 1a ----------------------------------------------------------------------
id83_signature_denovo_nonsmoker 
tdata <- id83_signature_denovo_nonsmoker %>% select(MutationType,ID83A)
pa <- plot_id_83_profile(data = tdata)  
tdata <- id83_signature_denovo_nonsmoker %>% select(MutationType,ID83B)
pb <- plot_id_83_profile(data = tdata)  
tdata <- id83_signature_denovo_nonsmoker %>% select(MutationType,ID83C)
pc <- plot_id_83_profile(data = tdata)  
tdata <- id83_signature_denovo_nonsmoker %>% select(MutationType,ID83D)
pd <- plot_id_83_profile(data = tdata)  
tdata <- id83_signature_denovo_nonsmoker %>% select(MutationType,ID83E)
pe <- plot_id_83_profile(data = tdata)  
tdata <- id83_signature_denovo_nonsmoker %>% select(MutationType,ID83F)
pf <- plot_id_83_profile(data = tdata)  
tdata <- id83_signature_denovo_nonsmoker %>% select(MutationType,ID83G)
pg <- plot_id_83_profile(data = tdata)  
tdata <- id83_signature_denovo_nonsmoker %>% select(MutationType,ID83H)
ph <- plot_id_83_profile(data = tdata)  
tdata <- id83_signature_denovo_nonsmoker %>% select(MutationType,ID83I)
pi <- plot_id_83_profile(data = tdata)  
tdata <- id83_signature_denovo_nonsmoker %>% select(MutationType,ID83J)
pj <- plot_id_83_profile(data = tdata)  

p <- plot_grid(pa,pb,pc,pd,pe,pf,pg,ph,pi,pj,align = 'h',ncol = 2)

ggsave(filename = 'id83_signature_denovo_nonsmoker.pdf',p,width = 24,height = 18,device = cairo_pdf)

# 2b ----------------------------------------------------------------------
id83_denovo_mapping_detail_nonsmoker %>% 
  mutate(de_novo_extracted = toupper(de_novo_extracted)) %>% 
  mutate(cosmic_signature = factor(cosmic_signature,levels=names(sigcol))) %>% 
  ggplot(aes(de_novo_extracted,contribution,fill=cosmic_signature))+
  geom_col( position = 'fill')+
  scale_fill_manual(values = sigcol)+
  geom_text(aes(label = paste0(cosmic_signature,"\n",percent_format(accuracy = 0.1)(contribution))), position = position_fill(vjust = .5),size=3)+
  theme_ipsum_rc(grid = FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank(),axis.text.x = element_text(angle = 45))+
  guides(fill = guide_legend(title = 'Cosmic Signatures',ncol = 1,title.position = 'top'))

ggsave(filename = 'Denovo_sig_mapping_id83.pdf',width = 7,height = 6,device = cairo_pdf)


# 2c ----------------------------------------------------------------------
#genomesize <- genome2size('GRCh38')
tdata <- id83_denovo_activity_nonsmoker %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker) %>%
  pivot_longer(-Tumor_Barcode) %>% 
  left_join(covdata0 %>% select(Tumor_Barcode,Histology) %>% mutate(Histology = as.character(Histology))) %>% 
  arrange(name,Tumor_Barcode) %>% 
  select(Cancer_Type=name,Sample=Tumor_Barcode,Sample_Group=Histology,Burden=value) %>% 
  mutate(Burden=log10(Burden))

sgroupcol <- ncicolpal[1:4]
names(sgroupcol) <- c('Squamous cell carcinoma','Other','Carcinoid tumor','Adenocarcinoma')
p <- TMBplot2(tdata,sgroupcol = sgroupcol,ylab = 'Number of ID mutations (log10)')+ theme(axis.text.x.bottom = element_text(size = 10,angle = 45,hjust = 1,vjust = 1))
ggsave(file='Denovo_sig_TMB_ID.pdf',width = 5,height = 5,device = cairo_pdf)




# ED_Figure 3 -------------------------------------------------------------

# 3a ----------------------------------------------------------------------
dbs78_signature_denovo_nonsmoker 
tdata <- dbs78_signature_denovo_nonsmoker %>% select(MutationType,DBS78A)
pa <- plot_dbs_78_profile(data = tdata,percentage = TRUE)  
tdata <- dbs78_signature_denovo_nonsmoker %>% select(MutationType,DBS78B)
pb <- plot_dbs_78_profile(data = tdata,percentage = TRUE)  
tdata <- dbs78_signature_denovo_nonsmoker %>% select(MutationType,DBS78C)
pc <- plot_dbs_78_profile(data = tdata,percentage = TRUE)  

p <- plot_grid(pa+theme(plot.margin = margin(t = 10,b=10)),pb+theme(plot.margin = margin(t = 10,b=10)),pc+theme(plot.margin = margin(t = 10,b=10)),align = 'v',axis = 'lr',ncol = 1)

ggsave(filename = 'dbs78_signature_denovo_nonsmoker.pdf',p,width = 12,height = 10,device = cairo_pdf)

# 3b ----------------------------------------------------------------------
dbs78_denovo_mapping_detail_nonsmoker %>% 
  mutate(de_novo_extracted = toupper(de_novo_extracted)) %>% 
  mutate(cosmic_signature = factor(cosmic_signature,levels=names(sigcol))) %>% 
  ggplot(aes(de_novo_extracted,contribution,fill=cosmic_signature))+
  geom_col( position = 'fill')+
  scale_fill_manual(values = sigcol)+
  geom_text(aes(label = paste0(cosmic_signature,"\n",percent_format(accuracy = 0.1)(contribution))), position = position_fill(vjust = .5),size=3)+
  theme_ipsum_rc(grid = FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank(),axis.text.x = element_text(angle = 45))+
  guides(fill = guide_legend(title = 'Cosmic Signatures',ncol = 1,title.position = 'top'))

ggsave(filename = 'Denovo_sig_mapping_dbs78.pdf',width = 3.5,height = 6,device = cairo_pdf)


# 3c ----------------------------------------------------------------------
#genomesize <- genome2size('GRCh38')
tdata <- dbs78_denovo_activity_nonsmoker %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker) %>%
  pivot_longer(-Tumor_Barcode) %>% 
  left_join(covdata0 %>% select(Tumor_Barcode,Histology) %>% mutate(Histology = as.character(Histology))) %>% 
  arrange(name,Tumor_Barcode) %>% 
  select(Cancer_Type=name,Sample=Tumor_Barcode,Sample_Group=Histology,Burden=value) %>% 
  mutate(Burden=log10(Burden))

sgroupcol <- ncicolpal[1:4]
names(sgroupcol) <- c('Squamous cell carcinoma','Other','Carcinoid tumor','Adenocarcinoma')
p <- TMBplot2(tdata,sgroupcol = sgroupcol,ylab = 'Number of DBS (log10)')+ theme(axis.text.x.bottom = element_text(size = 10,angle = 45,hjust = 1,vjust = 1))
ggsave(file='Denovo_sig_TMB_DBS.pdf',width = 3.5,height = 5,device = cairo_pdf)




# ED_Figure 4 -------------------------------------------------------------

# 4a ----------------------------------------------------------------------
cn68_signature_denovo_nonsmoker
tdata <- cn68_signature_denovo_nonsmoker %>% select(MutationType,CN68A)
pa <- plot_cn_68_profile(data = tdata)  
tdata <- cn68_signature_denovo_nonsmoker %>% select(MutationType,CN68B)
pb <- plot_cn_68_profile(data = tdata)  
tdata <- cn68_signature_denovo_nonsmoker %>% select(MutationType,CN68C)
pc <- plot_cn_68_profile(data = tdata)  
tdata <- cn68_signature_denovo_nonsmoker %>% select(MutationType,CN68D)
pd <- plot_cn_68_profile(data = tdata)  

p <- plot_grid(pa,pb,pc,pd,align = 'h',ncol = 2)

ggsave(filename = 'cn68_signature_denovo_nonsmoker.pdf',p,width = 24,height = 8,device = cairo_pdf)

# 4b ----------------------------------------------------------------------
cn68_denovo_mapping_detail_nonsmoker %>% 
  mutate(de_novo_extracted = toupper(de_novo_extracted)) %>% 
  mutate(cosmic_signature = factor(cosmic_signature,levels=names(sigcol))) %>% 
  ggplot(aes(de_novo_extracted,contribution,fill=cosmic_signature))+
  geom_col( position = 'fill')+
  scale_fill_manual(values = sigcol)+
  geom_text(aes(label = paste0(cosmic_signature,"\n",percent_format(accuracy = 0.1)(contribution))), position = position_fill(vjust = .5),size=3)+
  theme_ipsum_rc(grid = FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank(),axis.text.x = element_text(angle = 45))+
  guides(fill = guide_legend(title = 'Cosmic Signatures',ncol = 1,title.position = 'top'))

ggsave(filename = 'Denovo_sig_mapping_id83.pdf',width = 7,height = 6,device = cairo_pdf)


# 4c ----------------------------------------------------------------------
#genomesize <- genome2size('GRCh38')
tdata <- cn68_denovo_activity_nonsmoker %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker) %>%
  pivot_longer(-Tumor_Barcode) %>% 
  left_join(covdata0 %>% select(Tumor_Barcode,Histology) %>% mutate(Histology = as.character(Histology))) %>% 
  arrange(name,Tumor_Barcode) %>% 
  select(Cancer_Type=name,Sample=Tumor_Barcode,Sample_Group=Histology,Burden=value) %>% 
  mutate(Burden=log10(Burden))

sgroupcol <- ncicolpal[1:4]
names(sgroupcol) <- c('Squamous cell carcinoma','Other','Carcinoid tumor','Adenocarcinoma')
p <- TMBplot2(tdata,sgroupcol = sgroupcol,ylab = 'Number of ID mutations (log10)')+ theme(axis.text.x.bottom = element_text(size = 10,angle = 45,hjust = 1,vjust = 1))
ggsave(file='Denovo_sig_TMB_CN68.pdf',width = 3,height = 5,device = cairo_pdf)




# ED_Figure 5 -------------------------------------------------------------

# 5a ----------------------------------------------------------------------
sv38_signature_denovo_nonsmoker
tdata <- sv38_signature_denovo_nonsmoker %>% select(MutationType,SV38A)
pa <- plot_sv_38_profile(data = tdata)  
tdata <- sv38_signature_denovo_nonsmoker %>% select(MutationType,SV38B)
pb <- plot_sv_38_profile(data = tdata)  
tdata <- sv38_signature_denovo_nonsmoker %>% select(MutationType,SV38C)
pc <- plot_sv_38_profile(data = tdata)  
tdata <- sv38_signature_denovo_nonsmoker %>% select(MutationType,SV38D)
pd <- plot_sv_38_profile(data = tdata)  

p <- plot_grid(pa,pb,pc,pd,align = 'h',ncol = 2)

ggsave(filename = 'sv38_signature_denovo_nonsmoker.pdf',p,width = 24,height = 8,device = cairo_pdf)

# 5b ----------------------------------------------------------------------
cn68_denovo_mapping_detail_nonsmoker %>% 
  mutate(de_novo_extracted = toupper(de_novo_extracted)) %>% 
  mutate(cosmic_signature = factor(cosmic_signature,levels=names(sigcol))) %>% 
  ggplot(aes(de_novo_extracted,contribution,fill=cosmic_signature))+
  geom_col( position = 'fill')+
  scale_fill_manual(values = sigcol)+
  geom_text(aes(label = paste0(cosmic_signature,"\n",percent_format(accuracy = 0.1)(contribution))), position = position_fill(vjust = .5),size=3)+
  theme_ipsum_rc(grid = FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank(),axis.text.x = element_text(angle = 45))+
  guides(fill = guide_legend(title = 'Cosmic Signatures',ncol = 1,title.position = 'top'))

ggsave(filename = 'Denovo_sig_mapping_id83.pdf',width = 7,height = 6,device = cairo_pdf)


# 5c ----------------------------------------------------------------------
#genomesize <- genome2size('GRCh38')
tdata <- sv38_denovo_activity_nonsmoker %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker) %>%
  pivot_longer(-Tumor_Barcode) %>% 
  left_join(covdata0 %>% select(Tumor_Barcode,Histology) %>% mutate(Histology = as.character(Histology))) %>% 
  arrange(name,Tumor_Barcode) %>% 
  select(Cancer_Type=name,Sample=Tumor_Barcode,Sample_Group=Histology,Burden=value) %>% 
  mutate(Burden=log10(Burden))

sgroupcol <- ncicolpal[1:4]
names(sgroupcol) <- c('Squamous cell carcinoma','Other','Carcinoid tumor','Adenocarcinoma')
p <- TMBplot2(tdata,sgroupcol = sgroupcol,ylab = 'Number of SV events (log10)')+ theme(axis.text.x.bottom = element_text(size = 10,angle = 45,hjust = 1,vjust = 1))
ggsave(file='Denovo_sig_TMB_SV38.pdf',width = 3,height = 5,device = cairo_pdf)



# ED_Figure 6 -------------------------------------------------------------


# 6a ----------------------------------------------------------------------
testdata <- ludmil_activity_all_obs %>%
  select(Tumor_Barcode,!starts_with('SBS')) %>% 
  pivot_longer(-Tumor_Barcode) %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% luad_nonsmoker,Assigned_Population %in% c('EAS','EUR')) %>% 
  mutate(Assigned_Population = factor(Assigned_Population, levels=c('EAS','EUR'))) %>% 
  group_by(name) %>% 
  #do(tidy(fisher.test(.$value,.$Assigned_Population)))
  mutate(value=as.factor(value)) %>% 
  #do(tresult = safely(stats::fisher.test)(.$value,.$Assigned_Population)) %>% 
  do(tresult = safely(stats::glm)(value ~ Assigned_Population + Gender + Age + Tumor_Purity,family = binomial(),data=.)) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(name,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='Assigned_PopulationEUR') %>% 
  filter(abs(log(estimate))<10) %>% 
  mutate(profile=str_remove_all(name,'[:digit:].*')) %>% 
  group_by(profile) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  ungroup() %>% 
  mutate(profile=factor(profile,levels=c('DBS','ID','CN','SV'),labels= c('DBS Signatures','ID Signatures','CN Signatures','SV Signatures'))) %>% 
  mutate(dir=if_else(estimate>1,'Enriched in EUR ancestry','Enriched in EAS ancestry'))

testdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR),fill=dir))+
  geom_hline(yintercept = -log10(0.05),col=ncicolpal[2],linetype=2)+
  geom_hline(yintercept = -log10(0.01),col=ncicolpal[1],linetype=2)+
  geom_vline(xintercept = 0,col="gray10",linetype=1,linewidth=0.5)+
  geom_point(pch=21,size=3.5,col='white',stroke=0.2)+
  ggrepel::geom_text_repel(data=testdata %>% filter(p.value<0.05),aes(label=name))+
  facet_wrap(~profile,nrow = 1)+
  scale_fill_nejm()+
  # scale_x_continuous(breaks = pretty_breaks(n = 7),limits = c(-4.05,4.05))+
  scale_y_continuous(breaks = pretty_breaks(n=5),limits = c(0,4))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 12,grid = F,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing.x = unit('0.25','cm'),strip.text.x = element_text(face = 'bold',hjust = 0.5),legend.position = 'bottom')+
  labs(x = 'Odd ratio (log2)', y = '-log10(FDR)',fill='Enrichment')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  coord_cartesian(clip = 'off')

ggsave(filename = 'LCINS_LUAD_otherSigs_Ancestry_association.pdf',width = 12,height = 3.5,device = cairo_pdf)



# 6b ----------------------------------------------------------------------
source('../../Sherlock_functions.R')
plotdata <- sherlock_data_full %>% 
  filter(Gene %in% c('TP53','KRAS','EGFR'), Type=='Mutation_Driver') %>%
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% luad_nonsmoker,Assigned_Population %in% c('EAS','EUR')) %>% 
  mutate(Assigned_Population = factor(Assigned_Population, levels=c('EAS','EUR'))) 

barplot_fisher(plotdata %>% filter(Gene=='EGFR'),'Assigned_Population','Alteration',var1lab =  expression("Ancestry"),var2lab = 'EGFR mutation status',filename = 'LCINS_LUAD_Ancestry_driver_EGFR.pdf')

barplot_fisher(plotdata %>% filter(Gene=='TP53'),'Assigned_Population','Alteration',var1lab =  expression("Ancestry"),var2lab = 'TP53 mutation status',filename = 'LCINS_LUAD_Ancestry_driver_TP53.pdf')

barplot_fisher(plotdata %>% filter(Gene=='KRAS'),'Assigned_Population','Alteration',var1lab =  expression("Ancestry"),var2lab = 'KRAS mutation status',filename = 'LCINS_LUAD_Ancestry_driver_KRAS.pdf')



# ED_Figure 7 -------------------------------------------------------------


# 7a ----------------------------------------------------------------------
testdata <- ludmil_activity_all_obs %>%
  #select(Tumor_Barcode,!starts_with('SBS')) %>% 
  pivot_longer(-Tumor_Barcode) %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% luad_nonsmoker) %>% 
  group_by(name) %>% 
  #do(tidy(fisher.test(.$value,.$Assigned_Population)))
  mutate(value=as.factor(value)) %>% 
  #do(tresult = safely(stats::fisher.test)(.$value,.$Assigned_Population)) %>% 
  do(tresult = safely(stats::glm)(value ~ Assigned_Population + Gender + Age + Tumor_Purity,family = binomial(),data=.)) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(name,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='GenderFemale') %>% 
  filter(abs(log(estimate))<10) %>% 
  mutate(profile=str_remove_all(name,'[:digit:].*')) %>% 
  group_by(profile) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  ungroup() %>% 
  mutate(profile=factor(profile,levels=c('SBS','DBS','ID','CN','SV'),labels= c('SBS Signatures','DBS Signatures','ID Signatures','CN Signatures','SV Signatures'))) %>% 
  mutate(dir=if_else(estimate>1,'Enriched in Female','Enriched in Male'))

testdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR),fill=dir))+
  geom_hline(yintercept = -log10(0.05),col=ncicolpal[2],linetype=2)+
  geom_hline(yintercept = -log10(0.01),col=ncicolpal[1],linetype=2)+
  geom_vline(xintercept = 0,col="gray10",linetype=1,linewidth=0.5)+
  geom_point(pch=21,size=3.5,col='black',stroke=0.2)+
  ggrepel::geom_text_repel(data=testdata %>% filter(p.value<0.05),aes(label=name))+
  facet_wrap(~profile,nrow = 1)+
  scale_fill_npg()+
  # scale_x_continuous(breaks = pretty_breaks(n = 7),limits = c(-4.05,4.05))+
  scale_y_continuous(breaks = pretty_breaks(n=5))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 12,grid = F,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing.x = unit('0.25','cm'),strip.text.x = element_text(face = 'bold',hjust = 0.5),legend.position = 'bottom')+
  labs(x = 'Odd ratio (log2)', y = '-log10(FDR)',fill='Enrichment')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  coord_cartesian(clip = 'off')

ggsave(filename = 'LUAD_LUAD_Sigs_Gender_association.pdf',width = 15,height = 3.5,device = cairo_pdf)


# 7b ----------------------------------------------------------------------
testdata <- sherlock_data_full %>% 
  filter(Gene %in% drglist, Type=='Mutation_Driver') %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% luad_nonsmoker)

plotdata <- testdata %>% 
  group_by(Gene) %>% 
  mutate(Alteration=as.factor(Alteration)) %>% 
  do(tresult = safely(glm)(Alteration ~ Assigned_Population + Gender + Age + Tumor_Purity, family='binomial',data=. )) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(Gene,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='GenderFemale') %>% 
  filter(abs(log(estimate))<10) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(dir=if_else(estimate>1,'Enriched in Female','Enriched in Male'))



plotdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR)))+
  geom_hline(yintercept = -log10(0.01),linetype=2,col=ncicolpal[1],size=0.5)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col=ncicolpal[2],size=0.5)+
  geom_vline(xintercept = 0,linetype=1,col='gray50',size=0.5)+
  geom_point(aes(fill=dir),pch=21,stroke=0.1,col='black',size=4)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_y_continuous(breaks = pretty_breaks())+
  scale_fill_npg()+
  ggrepel::geom_text_repel(data=plotdata %>% filter(p.value<0.05),aes(label=Gene),max.overlaps = 30,size=4.5)+
  labs(x='Odd ratio (log2)',y='-log10(FDR)',fill='Enrichment')+
  # guides(fill = "none")+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid=FALSE,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing.x = unit('0.25','cm'),strip.text.x = element_text(face = 'bold',hjust = 0.5),legend.position = 'top')+
  panel_border(color = 'black',size = 0.5)

ggsave(filename = 'LCINSLUAD_DriverGene_Gender_association.pdf',width = 5,height = 4,device = cairo_pdf)


# 7c ----------------------------------------------------------------------
source('../../Sherlock_functions.R')
plotdata <- sherlock_data_full %>% 
  filter(Gene %in% c('TP53','KRAS','EGFR'), Type=='Mutation_Driver') %>%
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% luad_nonsmoker)

barplot_fisher(plotdata %>% filter(Gene=='EGFR'),'Gender','Alteration',var1lab =  "Gender",var2lab = 'EGFR mutation status',filename = 'LCINS_LUAD_Gender_driver_EGFR.pdf')

barplot_fisher(plotdata %>% filter(Gene=='TP53'),'Gender','Alteration',var1lab =  "Gender",var2lab = 'TP53 mutation status',filename = 'LCINS_LUAD_Gender_driver_TP53.pdf')

barplot_fisher(plotdata %>% filter(Gene=='KRAS'),'Gender','Alteration',var1lab =  "Gender",var2lab = 'KRAS mutation status',filename = 'LCINS_LUAD_Gender_driver_KRAS.pdf')



# ED_Figure 8 -------------------------------------------------------------

# a-b ---------------------------------------------------------------------


library(rstatix)
library(ggbeeswarm)

sbs4negative <- ludmil_activity_all %>% filter(SBS4==0) %>% pull(Tumor_Barcode)
plotdata <- bind_rows(
  sherlock_sbs288_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='SBS'),
  sherlock_dbs78_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='DBS'),
  sherlock_id83_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='ID'),
  sherlock_cn48_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='CN segments'),
  sherlock_sv38_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='SV')
) %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,Tumor_Barcode %in% sbs4negative) %>% 
  mutate(Histology=factor(Histology,levels=c("Adenocarcinoma", "Carcinoid tumor", "Squamous cell carcinoma", "Other"),labels = c('LUAD','LUSC','Carcinoid tumors','Others'))) %>% 
  mutate(Value=log10(Value)) %>% 
  filter(is.finite(Value)) %>% 
  mutate(profile = factor(profile,levels=c('SBS','DBS','ID','CN segments','SV')))
  


my_comparisons <- list(c("LUAD", "LUSC"),c("LUAD", "Carcinoid tumors"),c("LUAD", "Others"),c("LUSC", "Carcinoid tumors"),c("LUSC", "Others"),c("Carcinoid tumors","Others"))

stat.test <- plotdata %>% 
  group_by(profile) %>% 
  wilcox_test(Value ~ Histology, comparisons = my_comparisons) %>%
  ungroup() %>% 
  add_xy_position(x = "Histology",scales = 'free_y',step.increase = 0.15) %>%
  mutate(myformatted.p = sprintf("P = %.2e",p)) %>% 
  mutate(Histology = "LUAD") %>% 
  mutate(y.position =y.position)

stat.test


sgroupcol <- ncicolpal[1:4]
names(sgroupcol) <- c('LUSC','Others','Carcinoid tumors','LUAD')

plotdata %>% 
  ggplot(aes(Histology,Value,fill=Histology))+
  geom_quasirandom(pch=21,size=2,width = 0.3,color="white",stroke=0.4)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.7)+
  facet_wrap(~profile,scales = 'free_y',nrow = 1)+
  scale_fill_manual(values = sgroupcol)+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 15,axis_title_just = 'm',axis_title_size = 16,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 18,ticks = T)+
  labs(x = NULL, y = 'Number of events (log10)')+
  theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=14,face = 'bold',hjust = 0.5),panel.spacing.x = unit(0.2,'cm'),axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_nonSBS4_burden.pdf',width = 16 ,height =8,device = cairo_pdf)  

# b
plotdata <- sherlock_variable %>% 
  filter(str_detect(name,'Tel')) %>% 
  pivot_wider()%>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,Tumor_Barcode %in% sbs4negative) %>% 
  mutate(Histology=factor(Histology,levels=c("Adenocarcinoma", "Carcinoid tumor", "Squamous cell carcinoma", "Other"),labels = c('LUAD','LUSC','Carcinoid tumors','Others'))) %>% 
  mutate(Value=Telseq_TL_Ratio) %>% 
  filter(is.finite(Value))

stat.test <- plotdata %>% 
  wilcox_test(Value ~ Histology, comparisons = my_comparisons) %>%
  ungroup() %>% 
  add_xy_position(x = "Histology",scales = 'free_y',step.increase = 0.2) %>%
  mutate(myformatted.p = sprintf("P = %.2e",p)) %>% 
  mutate(Histology = "LUAD")

plotdata %>% 
  ggplot(aes(Histology,Value,fill=Histology))+
  geom_quasirandom(pch=21,size=2,width = 0.3,color="white",stroke=0.4)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.7)+
  # facet_wrap(~profile,scales = 'free_y',nrow = 1)+
  scale_fill_manual(values = sgroupcol)+
  scale_y_continuous(breaks = pretty_breaks(n = 7))+
  theme_ipsum_rc(base_size = 15,axis_title_just = 'm',axis_title_size = 16,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 18,ticks = T)+
  labs(x = NULL, y = 'Telomere length ratio, log2(Tumor TL/Normal TL')+
  theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=14,face = 'bold',hjust = 0.5),panel.spacing.x = unit(0.2,'cm'),axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_nonSBS4_TL.pdf',width = 4 ,height =7,device = cairo_pdf)  



# ED_Figure 9 -------------------------------------------------------------
conflicts_prefer(cowplot::get_legend)
conflicts_prefer(dplyr::lag)
carcinoid_nonsmoker <- covdata0 %>% filter(Smoking == 'Non-Smoker',Histology == 'Carcinoid tumor') %>% pull(Tumor_Barcode)

# 9a ----------------------------------------------------------------------
#SBS 
sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% carcinoid_nonsmoker) %>% select(Samples=Tumor_Barcode,starts_with('SBS')) 
tmpx <- sigdata %>% pivot_longer(-Samples) %>% group_by(name) %>% summarise(value=sum(value)) %>% filter(value==0) %>% pull(name)
sigdata <- sigdata %>% select(-one_of(tmpx))
studydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Study=Assigned_Population))

puritydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Gender=Gender))
puritycol <- ncicolpal[5:6]
names(puritycol) <- levels(puritydata$Gender)

cosinedata <-  sigdata %>% select(Samples) %>% left_join(sbs288_decompsite %>% select(Samples=Sample_Names,Similarity=Cosine_similarity))


Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Population',studydata_cat = TRUE,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Gender",cosinedata = cosinedata,clustern=10,sampletext_size = 0,output_plot = 'carcinoid_nonsmoker_SBS.pdf' )

prevalence_plot(sigdata = sigdata,nmutation = 0,output_plot = paste0('Carcinoid_nonsmoker_SBS_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))

# DBS
sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% carcinoid_nonsmoker) %>% select(Samples=Tumor_Barcode,starts_with('DBS')) 
tmpx <- sigdata %>% pivot_longer(-Samples) %>% group_by(name) %>% summarise(value=sum(value)) %>% filter(value==0) %>% pull(name)
sigdata <- sigdata %>% select(-one_of(tmpx))
studydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Study=Assigned_Population))

puritydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Gender=Gender))
puritycol <- ncicolpal[5:6]
names(puritycol) <- levels(puritydata$Gender)

cosinedata <-  sigdata %>% select(Samples) %>% left_join(sbs288_decompsite %>% select(Samples=Sample_Names,Similarity=Cosine_similarity))


Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Population',studydata_cat = TRUE,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Gender",cosinedata = cosinedata,clustern=2,sampletext_size = 0,output_plot = 'carcinoid_nonsmoker_DBS.pdf' )

prevalence_plot(sigdata = sigdata,nmutation = 0,output_plot = paste0('Carcinoid_nonsmoker_DBS_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))


# 9b ----------------------------------------------------------------------


# ID
sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% carcinoid_nonsmoker) %>% select(Samples=Tumor_Barcode,starts_with('ID')) 
tmpx <- sigdata %>% pivot_longer(-Samples) %>% group_by(name) %>% summarise(value=sum(value)) %>% filter(value==0) %>% pull(name)
sigdata <- sigdata %>% select(-one_of(tmpx))
studydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Study=Assigned_Population))

puritydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Gender=Gender))
puritycol <- ncicolpal[5:6]
names(puritycol) <- levels(puritydata$Gender)

cosinedata <-  sigdata %>% select(Samples) %>% left_join(sbs288_decompsite %>% select(Samples=Sample_Names,Similarity=Cosine_similarity))


Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Population',studydata_cat = TRUE,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Gender",cosinedata = cosinedata,clustern=2,sampletext_size = 0,output_plot = 'carcinoid_nonsmoker_ID.pdf' )

prevalence_plot(sigdata = sigdata,nmutation = 0,output_plot = paste0('Carcinoid_nonsmoker_ID_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))



# 9c ----------------------------------------------------------------------
library(data.table)
source('../../RDS/Oncoplots_functions.R')
tmp <- read_csv('../../RDS/oncoplot_colors.csv')
landscape_colors <- tmp$Color
names(landscape_colors) <- tmp$Name

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


# Define the sample set
sample_level0 <-  carcinoid_nonsmoker

result_top <- oncoplot(data = data_top,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_scna <- oncoplot(data = data_focal,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_arm <- oncoplot(data = data_arm,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_feature <- oncoplot(data = data_feature,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)

sample_new_level <- result_top$sample_level %>% 
  left_join(
    result_fusion$sample_level
  ) %>% left_join(
    result_scna$sample_level
  ) %>% left_join(
    result_arm$sample_level
  ) %>% left_join(
    result_feature$sample_level
  ) %>% arrange(Mutation_Driver,SCNA_Focal_Gene,Genomic_Feature) %>% 
  pull(Tumor_Barcode)

result_top <- oncoplot(data = data_top,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar = -0.2,bmar = -0.05)
#result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05)
result_scna <- oncoplot(data = data_focal,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,)
#result_arm <- oncoplot(data = data_arm,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,)
result_feature <- oncoplot(data = data_feature,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE,removeAlterbyColor=TRUE)

result_purity <- oncoplot2(data = data_feature2 %>% filter(Gene == 'Tumor_Purity'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(),tmar=0,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_ratio <- oncoplot2(data = data_feature2 %>% filter(Gene=='Subclonal_Mutation_Ratio'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(option = "C"),tmar=-0.2,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_nrpcc <- oncoplot2(data = data_feature2 %>% filter(Gene =='NRPCC'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material(palette = 'green'),tmar=-0.2,bmar =0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_cnvratio <- oncoplot2(data = data_feature2 %>% filter(Gene =='PGA'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material(palette = 'blue'),tmar=-0.2,bmar =0,p2_axis_hidden = TRUE,p2_hidden = TRUE)

result_tmb <- oncoplot3(data = data_tmb,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar=0.2,bmar=-0.2,height = 6)

oncoplot_final <- oncoplot_combined(result_tmb,result_top,result_scna,result_feature,result_purity,result_ratio,result_nrpcc,result_cnvratio) #result_fusion

save_plot(filename = 'Genome_landscape_LCINS_carcinoid.pdf',plot = oncoplot_final,base_height = 5,base_width = 6,device=cairo_pdf)



# ED_Figure 10 ------------------------------------------------------------

# 10a ---------------------------------------------------------------------
testdata <- ludmil_activity_all_obs %>%
  #select(Tumor_Barcode,!starts_with('SBS')) %>% 
  pivot_longer(-Tumor_Barcode) %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker, Histology %in% c('Adenocarcinoma','Carcinoid tumor')) %>% 
  mutate(Histology = factor(Histology, levels = c('Adenocarcinoma','Carcinoid tumor'),labels=c('LUAD','Carcinoids'))) %>% 
  group_by(name) %>% 
  #do(tidy(fisher.test(.$value,.$Assigned_Population)))
  mutate(value=as.factor(value)) %>% 
  #do(tresult = safely(stats::fisher.test)(.$value,.$Assigned_Population)) %>% 
  do(tresult = safely(stats::glm)(value ~ Assigned_Population + Gender + Age + Tumor_Purity + Histology,family = binomial(),data=.)) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(name,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='HistologyCarcinoids') %>% 
  filter(abs(log(estimate))<10) %>% 
  mutate(profile=str_remove_all(name,'[:digit:].*')) %>% 
  group_by(profile) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  ungroup() %>% 
  mutate(profile=factor(profile,levels=c('SBS','DBS','ID','CN','SV'),labels= c('SBS Signatures','DBS Signatures','ID Signatures','CN Signatures','SV Signatures'))) %>% 
  mutate(dir=if_else(estimate>1,'Enriched in Carcinoid tumors','Enriched in LUADs'))

testdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR),fill=dir))+
  geom_hline(yintercept = -log10(0.05),col=ncicolpal[2],linetype=2)+
  geom_hline(yintercept = -log10(0.01),col=ncicolpal[1],linetype=2)+
  geom_vline(xintercept = 0,col="gray10",linetype=1,linewidth=0.5)+
  geom_point(pch=21,size=3.5,col='black',stroke=0.2)+
  ggrepel::geom_text_repel(data=testdata %>% filter(p.value<0.05),aes(label=name))+
  facet_wrap(~profile,nrow = 1)+
  scale_fill_manual(values = ncicolpal[c(3,4)])+
  # scale_x_continuous(breaks = pretty_breaks(n = 7),limits = c(-4.05,4.05))+
  scale_y_continuous(breaks = pretty_breaks(n=5))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 12,grid = F,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing.x = unit('0.25','cm'),strip.text.x = element_text(face = 'bold',hjust = 0.5),legend.position = 'bottom')+
  labs(x = 'Odd ratio (log2)', y = '-log10(FDR)',fill='Enrichment')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  coord_cartesian(clip = 'off')

ggsave(filename = 'LCINS_Sigs_LUAVCarcinode_association.pdf',width = 15,height = 3.5,device = cairo_pdf)


# 10b ----------------------------------------------------------------------
testdata <- sherlock_data_full %>% 
  filter(Gene %in% drglist, Type=='Mutation_Driver') %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,Histology %in% c('Adenocarcinoma','Carcinoid tumor')) %>% 
  mutate(Histology = factor(Histology, levels = c('Adenocarcinoma','Carcinoid tumor'),labels=c('LUAD','Carcinoids'))) 
  
  
plotdata <- testdata %>% 
  group_by(Gene) %>% 
  mutate(Alteration=as.factor(Alteration)) %>% 
  do(tresult = safely(glm)(Alteration ~ Assigned_Population + Gender + Age + Tumor_Purity + Histology, family='binomial',data=. )) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(Gene,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='HistologyCarcinoids') %>% 
  filter(abs(log(estimate))<10) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(dir=if_else(estimate>1,'Enriched in Carcinoid tumors','Enriched in LUADs'))



plotdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR)))+
  geom_hline(yintercept = -log10(0.01),linetype=2,col=ncicolpal[1],size=0.5)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col=ncicolpal[2],size=0.5)+
  geom_vline(xintercept = 0,linetype=1,col='gray50',size=0.5)+
  geom_point(aes(fill=dir),pch=21,stroke=0.1,col='black',size=4)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_y_continuous(breaks = pretty_breaks())+
  scale_fill_manual(values = ncicolpal[c(3,4)])+
  ggrepel::geom_text_repel(data=plotdata %>% filter(p.value<0.05),aes(label=Gene),max.overlaps = 30,size=4.5)+
  labs(x='Odd ratio (log2)',y='-log10(FDR)',fill='Enrichment')+
  # guides(fill = "none")+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid=FALSE,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing.x = unit('0.25','cm'),strip.text.x = element_text(face = 'bold',hjust = 0.5),legend.position = 'top')+
  panel_border(color = 'black',size = 0.5)

ggsave(filename = 'LCINS_LUADvsCarcinoide_DriverGene_Gender_association.pdf',width = 5,height = 4,device = cairo_pdf)


# 10c ----------------------------------------------------------------------
source('../../Sherlock_functions.R')
plotdata <- sherlock_data_full %>% 
  filter(Gene %in% c('TP53','ARID1A','ARID1B'), Type=='Mutation_Driver') %>%
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,Histology %in% c('Adenocarcinoma','Carcinoid tumor')) %>% 
  mutate(Histology = factor(Histology,levels=c('Adenocarcinoma','Carcinoid tumor'),labels=c('LUAD','Carcinoid tumor')))

barplot_fisher(plotdata %>% filter(Gene=='ARID1A'),'Histology','Alteration',var1lab =  "Histology",var2lab = 'ARID1A mutation status',filename = 'LCINS_LUADvsCarcinoide_ARID1A.pdf',textcol = 'black')

barplot_fisher(plotdata %>% filter(Gene=='TP53'),'Histology','Alteration',var1lab =  "Histology",var2lab = 'TP53 mutation status',filename = 'LCINS_LUADvsCarcinoide_TP53.pdf',textcol = 'black')

barplot_fisher(plotdata %>% filter(Gene=='ARID1B'),'Histology','Alteration',var1lab =  "Histology",var2lab = 'ARID1B mutation status',filename = 'LCINS_LUADvsCarcinoide_ARID1B.pdf',textcol = 'black')




# ED_Figure 11 -------------------------------------------------------------
conflicts_prefer(cowplot::get_legend)
conflicts_prefer(dplyr::lag)
lusc_nonsmoker <- covdata0 %>% filter(Smoking == 'Non-Smoker',Histology == 'Squamous cell carcinoma') %>% pull(Tumor_Barcode)

# 11a ----------------------------------------------------------------------
#SBS 
sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% lusc_nonsmoker) %>% select(Samples=Tumor_Barcode,starts_with('SBS')) 
tmpx <- sigdata %>% pivot_longer(-Samples) %>% group_by(name) %>% summarise(value=sum(value)) %>% filter(value==0) %>% pull(name)
sigdata <- sigdata %>% select(-one_of(tmpx))
studydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Study=Assigned_Population))

puritydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Gender=Gender))
puritycol <- ncicolpal[5:6]
names(puritycol) <- levels(puritydata$Gender)

cosinedata <-  sigdata %>% select(Samples) %>% left_join(sbs288_decompsite %>% select(Samples=Sample_Names,Similarity=Cosine_similarity))


Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Population',studydata_cat = TRUE,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Gender",cosinedata = cosinedata,clustern=10,sampletext_size = 0,output_plot = 'LUSC_nonsmoker_SBS.pdf' )

prevalence_plot(sigdata = sigdata,nmutation = 0,output_plot = paste0('LUSC_nonsmoker_SBS_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))


# 11b ----------------------------------------------------------------------

# DBS
sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% lusc_nonsmoker) %>% select(Samples=Tumor_Barcode,starts_with('DBS')) 
tmpx <- sigdata %>% pivot_longer(-Samples) %>% group_by(name) %>% summarise(value=sum(value)) %>% filter(value==0) %>% pull(name)
sigdata <- sigdata %>% select(-one_of(tmpx))
studydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Study=Assigned_Population))

puritydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Gender=Gender))
puritycol <- ncicolpal[5:6]
names(puritycol) <- levels(puritydata$Gender)

cosinedata <-  sigdata %>% select(Samples) %>% left_join(sbs288_decompsite %>% select(Samples=Sample_Names,Similarity=Cosine_similarity))


Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Population',studydata_cat = TRUE,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Gender",cosinedata = cosinedata,clustern=2,sampletext_size = 0,output_plot = 'LUSC_nonsmoker_DBS.pdf' )

prevalence_plot(sigdata = sigdata,nmutation = 0,output_plot = paste0('LUSC_nonsmoker_DBS_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))


# 11c ----------------------------------------------------------------------
# ID
sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% lusc_nonsmoker) %>% select(Samples=Tumor_Barcode,starts_with('ID')) 
tmpx <- sigdata %>% pivot_longer(-Samples) %>% group_by(name) %>% summarise(value=sum(value)) %>% filter(value==0) %>% pull(name)
sigdata <- sigdata %>% select(-one_of(tmpx))
studydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Study=Assigned_Population))

puritydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Gender=Gender))
puritycol <- ncicolpal[5:6]
names(puritycol) <- levels(puritydata$Gender)

cosinedata <-  sigdata %>% select(Samples) %>% left_join(sbs288_decompsite %>% select(Samples=Sample_Names,Similarity=Cosine_similarity))


Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Population',studydata_cat = TRUE,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Gender",cosinedata = cosinedata,clustern=2,sampletext_size = 0,output_plot = 'LUSC_nonsmoker_ID.pdf' )

prevalence_plot(sigdata = sigdata,nmutation = 0,output_plot = paste0('LUSC_nonsmoker_ID_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))



# 11d ----------------------------------------------------------------------
library(data.table)
source('../../RDS/Oncoplots_functions.R')
tmp <- read_csv('../../RDS/oncoplot_colors.csv')
landscape_colors <- tmp$Color
names(landscape_colors) <- tmp$Name

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


# Define the sample set
sample_level0 <-  lusc_nonsmoker

result_top <- oncoplot(data = data_top,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_scna <- oncoplot(data = data_focal,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_arm <- oncoplot(data = data_arm,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_feature <- oncoplot(data = data_feature,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)

sample_new_level <- result_top$sample_level %>% 
   left_join(
    result_scna$sample_level
  ) %>% 
  left_join(
    result_feature$sample_level
  ) %>% arrange(Mutation_Driver,SCNA_Focal_Gene,Genomic_Feature) %>% 
  pull(Tumor_Barcode)

result_top <- oncoplot(data = data_top,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar = -0.2,bmar = -0.05)
#result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05)
result_scna <- oncoplot(data = data_focal,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,)
#result_arm <- oncoplot(data = data_arm,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,)
result_feature <- oncoplot(data = data_feature,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE,removeAlterbyColor=TRUE)

result_purity <- oncoplot2(data = data_feature2 %>% filter(Gene == 'Tumor_Purity'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(),tmar=0,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_ratio <- oncoplot2(data = data_feature2 %>% filter(Gene=='Subclonal_Mutation_Ratio'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(option = "C"),tmar=-0.2,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_nrpcc <- oncoplot2(data = data_feature2 %>% filter(Gene =='NRPCC'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material(palette = 'green'),tmar=-0.2,bmar =0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_cnvratio <- oncoplot2(data = data_feature2 %>% filter(Gene =='PGA'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material(palette = 'blue'),tmar=-0.2,bmar =0,p2_axis_hidden = TRUE,p2_hidden = TRUE)

result_tmb <- oncoplot3(data = data_tmb,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar=0.2,bmar=-0.2,height = 6)

oncoplot_final <- oncoplot_combined(result_tmb,result_top,result_scna,result_feature,result_purity,result_ratio,result_nrpcc,result_cnvratio) #

save_plot(filename = 'Genome_landscape_LCINS_LUSC.pdf',plot = oncoplot_final,base_height = 8,base_width = 6,device=cairo_pdf)




# ED_Figure 12 ------------------------------------------------------------

# 12a ---------------------------------------------------------------------
testdata <- ludmil_activity_all_obs %>%
  #select(Tumor_Barcode,!starts_with('SBS')) %>% 
  pivot_longer(-Tumor_Barcode) %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker, Histology %in% c('Adenocarcinoma','Squamous cell carcinoma')) %>% 
  mutate(Histology = factor(Histology, levels = c('Adenocarcinoma','Squamous cell carcinoma'),labels=c('LUAD','LUSC'))) %>% 
  group_by(name) %>% 
  #do(tidy(fisher.test(.$value,.$Assigned_Population)))
  mutate(value=as.factor(value)) %>% 
  #do(tresult = safely(stats::fisher.test)(.$value,.$Assigned_Population)) %>% 
  do(tresult = safely(stats::glm)(value ~ Assigned_Population + Gender + Age + Tumor_Purity + Histology,family = binomial(),data=.)) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(name,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='HistologyLUSC') %>% 
  filter(abs(log(estimate))<10) %>% 
  mutate(profile=str_remove_all(name,'[:digit:].*')) %>% 
  group_by(profile) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  ungroup() %>% 
  mutate(profile=factor(profile,levels=c('SBS','DBS','ID','CN','SV'),labels= c('SBS Signatures','DBS Signatures','ID Signatures','CN Signatures','SV Signatures'))) %>% 
  mutate(dir=if_else(estimate>1,'Enriched in LUSC','Enriched in LUAD'))

testdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR),fill=dir))+
  geom_hline(yintercept = -log10(0.05),col=ncicolpal[2],linetype=2)+
  geom_hline(yintercept = -log10(0.01),col=ncicolpal[1],linetype=2)+
  geom_vline(xintercept = 0,col="gray10",linetype=1,linewidth=0.5)+
  geom_point(pch=21,size=3.5,col='black',stroke=0.2)+
  ggrepel::geom_text_repel(data=testdata %>% filter(p.value<0.05),aes(label=name))+
  facet_wrap(~profile,nrow = 1)+
  scale_fill_manual(values = ncicolpal[c(4,1)])+
  # scale_x_continuous(breaks = pretty_breaks(n = 7),limits = c(-4.05,4.05))+
  scale_y_continuous(breaks = pretty_breaks(n=5))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 12,grid = F,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing.x = unit('0.25','cm'),strip.text.x = element_text(face = 'bold',hjust = 0.5),legend.position = 'bottom')+
  labs(x = 'Odd ratio (log2)', y = '-log10(FDR)',fill='Enrichment')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  coord_cartesian(clip = 'off')

ggsave(filename = 'LCINS_Sigs_LUADvsLUSC_association.pdf',width = 15,height = 3.5,device = cairo_pdf)


# 12b ----------------------------------------------------------------------
testdata <- sherlock_data_full %>% 
  filter(Gene %in% drglist, Type=='Mutation_Driver') %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker, Histology %in% c('Adenocarcinoma','Squamous cell carcinoma')) %>% 
  mutate(Histology = factor(Histology, levels = c('Adenocarcinoma','Squamous cell carcinoma'),labels=c('LUAD','LUSC'))) 


plotdata <- testdata %>% 
  group_by(Gene) %>% 
  mutate(Alteration=as.factor(Alteration)) %>% 
  do(tresult = safely(glm)(Alteration ~ Assigned_Population + Gender + Age + Tumor_Purity + Histology, family='binomial',data=. )) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(Gene,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='HistologyLUSC') %>% 
  filter(abs(log(estimate))<10) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(dir=if_else(estimate>1,'Enriched in LUSC','Enriched in LUAD'))



plotdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR)))+
  geom_hline(yintercept = -log10(0.01),linetype=2,col=ncicolpal[1],size=0.5)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col=ncicolpal[2],size=0.5)+
  geom_vline(xintercept = 0,linetype=1,col='gray50',size=0.5)+
  geom_point(aes(fill=dir),pch=21,stroke=0.1,col='black',size=4)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_y_continuous(breaks = pretty_breaks())+
  scale_fill_manual(values = ncicolpal[c(4,1)])+
  ggrepel::geom_text_repel(data=plotdata %>% filter(p.value<0.05),aes(label=Gene),max.overlaps = 30,size=4.5)+
  labs(x='Odd ratio (log2)',y='-log10(FDR)',fill='Enrichment')+
  # guides(fill = "none")+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid=FALSE,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing.x = unit('0.25','cm'),strip.text.x = element_text(face = 'bold',hjust = 0.5),legend.position = 'top')+
  panel_border(color = 'black',size = 0.5)

ggsave(filename = 'LCINS_LUADvsLUSC_DriverGene_Gender_association.pdf',width = 5,height = 4,device = cairo_pdf)


# 12c ----------------------------------------------------------------------
source('../../Sherlock_functions.R')
plotdata <- sherlock_data_full %>% 
  filter(Gene %in% c('TP53','EGFR','LRP1B','PTEN'), Type=='Mutation_Driver') %>%
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,Histology %in% c('Adenocarcinoma','Squamous cell carcinoma')) %>% 
  mutate(Histology = factor(Histology,levels=c('Adenocarcinoma','Squamous cell carcinoma'),labels=c('LUAD','LUSC')))


barplot_fisher(plotdata %>% filter(Gene=='LRP1B'),'Histology','Alteration',var1lab =  "Histology",var2lab = 'LRP1B mutation status',filename = 'LCINS_LUADvsLUSC_LRP1B.pdf',textcol = 'black')

barplot_fisher(plotdata %>% filter(Gene=='TP53'),'Histology','Alteration',var1lab =  "Histology",var2lab = 'TP53 mutation status',filename = 'LCINS_LUADvsLUSC_TP53.pdf')

barplot_fisher(plotdata %>% filter(Gene=='EGFR'),'Histology','Alteration',var1lab =  "Histology",var2lab = 'EGFR mutation status',filename = 'LCINS_LUADvsLUSC_EGFR.pdf')

barplot_fisher(plotdata %>% filter(Gene=='PTEN'),'Histology','Alteration',var1lab =  "Histology",var2lab = 'PTEN mutation status',filename = 'LCINS_LUADvsLUSC_PTEN.pdf',textcol = 'black')






# ED_Figure 13 -------------------------------------------------------------
conflicts_prefer(cowplot::get_legend)
conflicts_prefer(dplyr::lag)

covdata_rare <- covdata0 %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,Histology %in% 'Other') %>% 
  select(Tumor_Barcode) %>% 
  left_join(clinical_data %>% select(Tumor_Barcode,Histology,Histology_Detail)) %>% 
  mutate(Histology_Detail = toupper(Histology_Detail),Histology=if_else(str_detect(Histology_Detail,'LARGE'),'Large Cell Carcinomas',Histology)) %>%
  mutate(Histology = if_else(is.na(Histology),"Others",Histology)) %>% 
  select(Tumor_Barcode,Histology) %>% 
  left_join(covdata0 %>% select(-Histology))

rare_nonsmoker <- covdata_rare %>% pull(Tumor_Barcode)

# 13a ----------------------------------------------------------------------
#SBS 
sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% rare_nonsmoker) %>% select(Samples=Tumor_Barcode,starts_with('SBS')) 
tmpx <- sigdata %>% pivot_longer(-Samples) %>% group_by(name) %>% summarise(value=sum(value)) %>% filter(value==0) %>% pull(name)
sigdata <- sigdata %>% select(-one_of(tmpx))
studydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Study=Assigned_Population))

puritydata <- sigdata %>% select(Samples) %>% left_join(covdata_rare %>% select(Samples=Tumor_Barcode,Histology))
puritycol <- ncicolpal[c(20,18,2)]
names(puritycol) <- sort(unique(puritydata$Histology))

cosinedata <-  sigdata %>% select(Samples) %>% left_join(sbs288_decompsite %>% select(Samples=Sample_Names,Similarity=Cosine_similarity))


Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Ancestry',studydata_cat = TRUE,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Histology",cosinedata = cosinedata,clustern=10,sampletext_size = 0,output_plot = 'Rare_nonsmoker_SBS.pdf' )

prevalence_plot(sigdata = sigdata,nmutation = 0,output_plot = paste0('Rare_nonsmoker_SBS_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))


# 13b ----------------------------------------------------------------------

# DBS
sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% rare_nonsmoker) %>% select(Samples=Tumor_Barcode,starts_with('DBS')) 
tmpx <- sigdata %>% pivot_longer(-Samples) %>% group_by(name) %>% summarise(value=sum(value)) %>% filter(value==0) %>% pull(name)
sigdata <- sigdata %>% select(-one_of(tmpx))
studydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Study=Assigned_Population))

puritydata <- sigdata %>% select(Samples) %>% left_join(covdata_rare %>% select(Samples=Tumor_Barcode,Histology))
puritycol <- ncicolpal[c(20,18,2)]
names(puritycol) <- sort(unique(puritydata$Histology))

cosinedata <-  sigdata %>% select(Samples) %>% left_join(sbs288_decompsite %>% select(Samples=Sample_Names,Similarity=Cosine_similarity))


Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Ancestry',studydata_cat = TRUE,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Histology",cosinedata = cosinedata,clustern=2,sampletext_size = 0,output_plot = 'Rare_nonsmoker_DBS.pdf' )

prevalence_plot(sigdata = sigdata,nmutation = 0,output_plot = paste0('Rare_nonsmoker_DBS_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))


# 13c ----------------------------------------------------------------------

# ID
sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% rare_nonsmoker) %>% select(Samples=Tumor_Barcode,starts_with('ID')) 
tmpx <- sigdata %>% pivot_longer(-Samples) %>% group_by(name) %>% summarise(value=sum(value)) %>% filter(value==0) %>% pull(name)
sigdata <- sigdata %>% select(-one_of(tmpx))
studydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Study=Assigned_Population))

puritydata <- sigdata %>% select(Samples) %>% left_join(covdata_rare %>% select(Samples=Tumor_Barcode,Histology))
puritycol <- ncicolpal[c(20,18,2)]
names(puritycol) <- sort(unique(puritydata$Histology))

cosinedata <-  sigdata %>% select(Samples) %>% left_join(sbs288_decompsite %>% select(Samples=Sample_Names,Similarity=Cosine_similarity))


Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Ancestry',studydata_cat = TRUE,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Histology",cosinedata = cosinedata,clustern=2,sampletext_size = 0,output_plot = 'Rare_nonsmoker_ID.pdf' )

prevalence_plot(sigdata = sigdata,nmutation = 0,output_plot = paste0('Rare_nonsmoker_ID_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))


# 13d ----------------------------------------------------------------------

library(data.table)
source('../../RDS/Oncoplots_functions.R')
tmp <- read_csv('../../RDS/oncoplot_colors.csv')
landscape_colors <- tmp$Color
names(landscape_colors) <- tmp$Name

landscape_colors <-  c(landscape_colors,puritycol)

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
  select(Subject,Tumor_Barcode,Assigned_Population,Gender,Smoking,WGD_Status,SCNA_Group,MMR_Status,HR_Status,HRD_Type,HLA_LOH,HRDetect_Status,Kataegis,EBV) %>% 
  left_join(covdata_rare %>% select(Tumor_Barcode,Histology)) %>% 
  pivot_longer(cols = -c(Subject,Tumor_Barcode)) %>% 
  filter(!value %in% c('No','None','N')) %>% 
  mutate(Type='Genomic_Feature') %>% 
  select(Subject,Tumor_Barcode,Gene=name,Alteration = value,Type) 

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


# Define the sample set
sample_level0 <-  rare_nonsmoker

result_top <- oncoplot(data = data_top,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_scna <- oncoplot(data = data_focal,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
#result_arm <- oncoplot(data = data_arm,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
#result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)
result_feature <- oncoplot(data = data_feature,landscape_colors = landscape_colors,sample_level0 = sample_level0,GeneSortOnly = TRUE)

sample_new_level <- result_top$sample_level %>% 
  left_join(
    result_scna$sample_level
  ) %>% 
  left_join(
    result_feature$sample_level
  ) %>% arrange(Mutation_Driver,SCNA_Focal_Gene,Genomic_Feature) %>% 
  pull(Tumor_Barcode)

result_top <- oncoplot(data = data_top,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar = -0.2,bmar = -0.05)
#result_fusion <- oncoplot(data = data_fusion,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05)
result_scna <- oncoplot(data = data_focal,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,)
#result_arm <- oncoplot(data = data_arm,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = -0.05,)
result_feature <- oncoplot(data = data_feature,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,bmar = 0,p2_axis_hidden = TRUE,p2_hidden = TRUE,removeAlterbyColor=TRUE)

result_purity <- oncoplot2(data = data_feature2 %>% filter(Gene == 'Tumor_Purity'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(),tmar=0,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_ratio <- oncoplot2(data = data_feature2 %>% filter(Gene=='Subclonal_Mutation_Ratio'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_viridis_c(option = "C"),tmar=-0.2,bmar = -0.1,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_nrpcc <- oncoplot2(data = data_feature2 %>% filter(Gene =='NRPCC'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material(palette = 'green'),tmar=-0.2,bmar =0,p2_axis_hidden = TRUE,p2_hidden = TRUE)
result_cnvratio <- oncoplot2(data = data_feature2 %>% filter(Gene =='PGA'),sample_level0 = sample_level0,sample_level = sample_new_level,scale_fill_ztw = scale_fill_material(palette = 'blue'),tmar=-0.2,bmar =0,p2_axis_hidden = TRUE,p2_hidden = TRUE)

result_tmb <- oncoplot3(data = data_tmb,landscape_colors = landscape_colors,sample_level0 = sample_level0,sample_level = sample_new_level,tmar=0.2,bmar=-0.2,height = 6)

oncoplot_final <- oncoplot_combined(result_tmb,result_top,result_scna,result_feature,result_purity,result_ratio,result_nrpcc,result_cnvratio) #

save_plot(filename = 'Genome_landscape_LCINS_Rare.pdf',plot = oncoplot_final,base_height = 9,base_width = 6,device=cairo_pdf)






# ED_Figure 14 ------------------------------------------------------------
covdata_rare <- covdata0 %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,Histology %in% 'Other') %>% 
  select(Tumor_Barcode) %>% 
  left_join(clinical_data %>% select(Tumor_Barcode,Histology,Histology_Detail)) %>% 
  mutate(Histology_Detail = toupper(Histology_Detail),Histology=if_else(str_detect(Histology_Detail,'LARGE'),'Large Cell Carcinomas',Histology)) %>%
  mutate(Histology = if_else(is.na(Histology),"Others",Histology)) %>% 
  select(Tumor_Barcode,Histology) %>% 
  left_join(covdata0 %>% select(-Histology))

rare_nonsmoker <- covdata_rare %>% pull(Tumor_Barcode)


# 14a ---------------------------------------------------------------------
testdata <- ludmil_activity_all_obs %>%
  #select(Tumor_Barcode,!starts_with('SBS')) %>% 
  pivot_longer(-Tumor_Barcode) %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker, Histology %in% c('Adenocarcinoma','Other')) %>% 
  mutate(Histology = factor(Histology, levels = c('Adenocarcinoma','Other'),labels=c('LUAD','Others'))) %>% 
  group_by(name) %>% 
  #do(tidy(fisher.test(.$value,.$Assigned_Population)))
  mutate(value=as.factor(value)) %>% 
  #do(tresult = safely(stats::fisher.test)(.$value,.$Assigned_Population)) %>% 
  do(tresult = safely(stats::glm)(value ~ Assigned_Population + Gender + Age + Tumor_Purity + Histology,family = binomial(),data=.)) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(name,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='HistologyOthers') %>% 
  filter(abs(log(estimate))<10) %>% 
  mutate(profile=str_remove_all(name,'[:digit:].*')) %>% 
  group_by(profile) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  ungroup() %>% 
  mutate(profile=factor(profile,levels=c('SBS','DBS','ID','CN','SV'),labels= c('SBS Signatures','DBS Signatures','ID Signatures','CN Signatures','SV Signatures'))) %>% 
  mutate(dir=if_else(estimate>1,'Enriched in Others','Enriched in LUAD'))

testdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR),fill=dir))+
  geom_hline(yintercept = -log10(0.05),col=ncicolpal[2],linetype=2)+
  geom_hline(yintercept = -log10(0.01),col=ncicolpal[1],linetype=2)+
  geom_vline(xintercept = 0,col="gray10",linetype=1,linewidth=0.5)+
  geom_point(pch=21,size=3.5,col='black',stroke=0.2)+
  ggrepel::geom_text_repel(data=testdata %>% filter(p.value<0.05),aes(label=name))+
  facet_wrap(~profile,nrow = 1)+
  scale_fill_manual(values = ncicolpal[c(4,18)])+
  # scale_x_continuous(breaks = pretty_breaks(n = 7),limits = c(-4.05,4.05))+
  scale_y_continuous(breaks = pretty_breaks(n=5))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 12,grid = F,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing.x = unit('0.25','cm'),strip.text.x = element_text(face = 'bold',hjust = 0.5),legend.position = 'bottom')+
  labs(x = 'Odd ratio (log2)', y = '-log10(FDR)',fill='Enrichment')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  coord_cartesian(clip = 'off')

ggsave(filename = 'LCINS_Sigs_LUADvsOthers_association.pdf',width = 15,height = 3.5,device = cairo_pdf)


# 14b ----------------------------------------------------------------------
testdata <- sherlock_data_full %>% 
  filter(Gene %in% drglist, Type=='Mutation_Driver') %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker, Histology %in% c('Adenocarcinoma','Other')) %>% 
  mutate(Histology = factor(Histology, levels = c('Adenocarcinoma','Other'),labels=c('LUAD','Others')))


plotdata <- testdata %>% 
  group_by(Gene) %>% 
  mutate(Alteration=as.factor(Alteration)) %>% 
  do(tresult = safely(glm)(Alteration ~ Assigned_Population + Gender + Age + Tumor_Purity + Histology, family='binomial',data=. )) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(Gene,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='HistologyOthers') %>% 
  filter(abs(log(estimate))<10) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(dir=if_else(estimate>1,'Enriched in Others','Enriched in LUAD'))



plotdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR)))+
  geom_hline(yintercept = -log10(0.01),linetype=2,col=ncicolpal[1],size=0.5)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col=ncicolpal[2],size=0.5)+
  geom_vline(xintercept = 0,linetype=1,col='gray50',size=0.5)+
  geom_point(aes(fill=dir),pch=21,stroke=0.1,col='black',size=4)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_y_continuous(breaks = pretty_breaks())+
  scale_fill_manual(values = ncicolpal[c(4,18)])+
  ggrepel::geom_text_repel(data=plotdata %>% filter(p.value<0.05),aes(label=Gene),max.overlaps = 30,size=4.5)+
  labs(x='Odd ratio (log2)',y='-log10(FDR)',fill='Enrichment')+
  # guides(fill = "none")+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid=FALSE,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing.x = unit('0.25','cm'),strip.text.x = element_text(face = 'bold',hjust = 0.5),legend.position = 'top')+
  panel_border(color = 'black',size = 0.5)

ggsave(filename = 'LCINS_LUADvsOthers_DriverGene_Gender_association.pdf',width = 5,height = 4,device = cairo_pdf)


# 14c ----------------------------------------------------------------------
source('../../Sherlock_functions.R')
plotdata <- sherlock_data_full %>% 
  filter(Gene %in% c('TP53','EGFR'), Type=='Mutation_Driver') %>%
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,Histology %in% c('Adenocarcinoma','Other')) %>% 
  mutate(Histology = factor(Histology,levels=c('Adenocarcinoma','Other'),labels=c('LUAD','Others')))


barplot_fisher(plotdata %>% filter(Gene=='TP53'),'Histology','Alteration',var1lab =  "Histology",var2lab = 'TP53 mutation status',filename = 'LCINS_LUADvsOthers_TP53.pdf')

barplot_fisher(plotdata %>% filter(Gene=='EGFR'),'Histology','Alteration',var1lab =  "Histology",var2lab = 'EGFR mutation status',filename = 'LCINS_LUADvsOthers_EGFR.pdf')


# ED_Figure 15 ------------------------------------------------------------

# 15a ---------------------------------------------------------------------
library(rstatix)
library(ggbeeswarm)

plotdata <- bind_rows(
  sherlock_sbs288_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='SBS'),
  sherlock_dbs78_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='DBS'),
  sherlock_id83_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='ID'),
  sherlock_cn48_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='CN segments'),
  sherlock_sv38_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='SV')
) %>% 
 left_join(ludmil_activity_all_obs %>% select(Tumor_Barcode,SBS4) %>% mutate(SBS4=if_else(SBS4,'SBS4','non-SBS4'))) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker) %>% 
  mutate(Value=log10(Value)) %>% 
  filter(is.finite(Value)) %>% 
  mutate(profile = factor(profile,levels=c('SBS','DBS','ID','CN segments','SV')))

my_comparisons <- list(c("SBS4", "non-SBS4"))

stat.test <- plotdata %>% 
  group_by(profile) %>% 
  wilcox_test(Value ~ SBS4, comparisons = my_comparisons) %>%
  ungroup() %>% 
  add_xy_position(x = "SBS4",scales = 'free_y',step.increase = 0.5) %>%
  mutate(myformatted.p = sprintf("P = %.2e",p)) %>% 
  mutate(SBS4 = "SBS4")

stat.test

plotdata %>% 
  ggplot(aes(SBS4,Value,fill=SBS4))+
  geom_quasirandom(pch=21,size=2,width = 0.3,color="black",stroke=0.4)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.7)+
  facet_wrap(~profile,scales = 'free_y',nrow = 1)+
  scale_fill_manual(values = c('#78c679','#2EBAED'))+
  scale_y_continuous(breaks = pretty_breaks(n=7))+
  theme_ipsum_rc(base_size = 15,axis_title_just = 'm',axis_title_size = 16,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 18,ticks = T)+
  labs(x = NULL, y = 'Number of events (log10)')+
  theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=14,face = 'bold',hjust = 0.5),panel.spacing.x = unit(0.2,'cm'),axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_SBS4vsnonSBS4_burden.pdf',width = 9 ,height =5,device = cairo_pdf)  


# 15b ----------------------------------------------------------------------
conflicts_prefer(cowplot::get_legend)
conflicts_prefer(dplyr::lag)
#SBS 
sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% sherlock_nonsmoker_SBS4) %>% select(Samples=Tumor_Barcode,starts_with('SBS')) 
tmpx <- sigdata %>% pivot_longer(-Samples) %>% group_by(name) %>% summarise(value=sum(value)) %>% filter(value==0) %>% pull(name)
sigdata <- sigdata %>% select(-one_of(tmpx))
studydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Study=Assigned_Population))

puritydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Gender=Gender))
puritycol <- ncicolpal[5:6]
names(puritycol) <- levels(puritydata$Gender)

cosinedata <-  sigdata %>% select(Samples) %>% left_join(sbs288_decompsite %>% select(Samples=Sample_Names,Similarity=Cosine_similarity))

Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Ancestry',studydata_cat = TRUE,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Gender",cosinedata = cosinedata,clustern=10,sampletext_size = 0,output_plot = 'SBS4_nonsmoker_SBS.pdf' )

studydata <- sigdata %>% select(Samples) %>% left_join(clinical_data %>% select(Samples=Tumor_Barcode,Passive_Smoking=Passive_Smoking))
studycol <- ncicolpal[10:11]
names(studycol) <- c('Yes','No')

Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Passive Smoking',studydata_cat = TRUE,studycolor = studycol,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Gender",cosinedata = cosinedata,clustern=10,sampletext_size = 0,output_plot = 'tmp1.pdf' )

prevalence_plot(sigdata = sigdata,nmutation = 0,output_plot = paste0('SBS4_nonsmoker_SBS_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))


# 13b ----------------------------------------------------------------------

# DBS
sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% sherlock_nonsmoker_SBS4) %>% select(Samples=Tumor_Barcode,starts_with('DBS')) 
tmpx <- sigdata %>% pivot_longer(-Samples) %>% group_by(name) %>% summarise(value=sum(value)) %>% filter(value==0) %>% pull(name)
sigdata <- sigdata %>% select(-one_of(tmpx))
studydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Study=Assigned_Population))

puritydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Gender=Gender))
puritycol <- ncicolpal[5:6]
names(puritycol) <- levels(puritydata$Gender)

cosinedata <-  sigdata %>% select(Samples) %>% left_join(sbs288_decompsite %>% select(Samples=Sample_Names,Similarity=Cosine_similarity))

Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Ancestry',studydata_cat = TRUE,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Gender",cosinedata = cosinedata,clustern=10,sampletext_size = 0,output_plot = 'SBS4_nonsmoker_DBS.pdf' )

studydata <- sigdata %>% select(Samples) %>% left_join(clinical_data %>% select(Samples=Tumor_Barcode,Passive_Smoking=Passive_Smoking))
studycol <- ncicolpal[10:11]
names(studycol) <- c('Yes','No')

Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Passive Smoking',studydata_cat = TRUE,studycolor = studycol,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Gender",cosinedata = cosinedata,clustern=10,sampletext_size = 0,output_plot = 'tmp1.pdf' )

prevalence_plot(sigdata = sigdata,nmutation = 0,output_plot = paste0('SBS4_nonsmoker_DBS_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))



# 13c ----------------------------------------------------------------------

# ID
sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% sherlock_nonsmoker_SBS4) %>% select(Samples=Tumor_Barcode,starts_with('ID')) 
tmpx <- sigdata %>% pivot_longer(-Samples) %>% group_by(name) %>% summarise(value=sum(value)) %>% filter(value==0) %>% pull(name)
sigdata <- sigdata %>% select(-one_of(tmpx))
studydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Study=Assigned_Population))

puritydata <- sigdata %>% select(Samples) %>% left_join(covdata0 %>% select(Samples=Tumor_Barcode,Gender=Gender))
puritycol <- ncicolpal[5:6]
names(puritycol) <- levels(puritydata$Gender)

cosinedata <-  sigdata %>% select(Samples) %>% left_join(sbs288_decompsite %>% select(Samples=Sample_Names,Similarity=Cosine_similarity))

Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Ancestry',studydata_cat = TRUE,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Gender",cosinedata = cosinedata,clustern=10,sampletext_size = 0,output_plot = 'SBS4_nonsmoker_ID.pdf' )

studydata <- sigdata %>% select(Samples) %>% left_join(clinical_data %>% select(Samples=Tumor_Barcode,Passive_Smoking=Passive_Smoking))
studycol <- ncicolpal[10:11]
names(studycol) <- c('Yes','No')

Exposure_Clustering(sigdata = sigdata,studydata = studydata,sigcolor = sigcol[colnames(sigdata)],studyname = 'Passive Smoking',studydata_cat = TRUE,studycolor = studycol,puritydata = puritydata,puritydata_cat = T,puritycol = puritycol,purityname = "Gender",cosinedata = cosinedata,clustern=10,sampletext_size = 0,output_plot = 'tmp1.pdf' )

prevalence_plot(sigdata = sigdata,nmutation = 0,output_plot = paste0('SBS4_nonsmoker_ID_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))




# ED_Figure 17 ---------------------------------------------------------------
library(rstatix)
library(ggbeeswarm)
tdata <- tibble(Tumor_Barcode=sherlock_nonsmoker) %>% 
  left_join(
    sherlock_sbs288_profile %>% group_by(Tumor_Barcode) %>% summarise(Total_SBS=log10(sum(Contribution)))
  ) %>% 
  left_join(
    sherlock_variable %>% filter(name=='Telseq_TL_Ratio') %>% pivot_wider()
  ) %>% 
  left_join(covdata0) %>% 
  left_join(clinical_data %>% select(Tumor_Barcode, Passive_Smoking)) %>% 
  filter(!is.na(Passive_Smoking)) %>% 
  filter(Tumor_Barcode %in% luad_nonsmoker)

my_comparisons <- list(c("No", "Yes"))


# 17a ---------------------------------------------------------------------
stat.test <- tdata %>% 
  wilcox_test(Total_SBS ~ Passive_Smoking, comparisons = my_comparisons) %>%
  add_xy_position(x = "Passive_Smoking") %>%
  mutate(myformatted.p = sprintf("P = %.4f",p)) %>% 
  mutate(Passive_Smoking = NA)

stat.test

tdata %>% 
  ggplot(aes(Passive_Smoking,Total_SBS,fill=Passive_Smoking))+
  geom_quasirandom(pch=21,size=2,width = 0.3,color="white",stroke=0.5)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.7)+
  scale_fill_jama()+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 13,axis_title_just = 'm',axis_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 15,ticks = T)+
  labs(x = "Passive smoking", y = 'Total number of SBS mutations (log10)',title = 'SBS')+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_LUAD_passive_smoking_SBS.pdf',width = 3 ,height =6,device = cairo_pdf)  

# 17b ---------------------------------------------------------------------

stat.test <- tdata %>% 
  wilcox_test(Telseq_TL_Ratio ~ Passive_Smoking, comparisons = my_comparisons) %>%
  add_xy_position(x = "Passive_Smoking") %>%
  mutate(myformatted.p = sprintf("P = %.2e",p)) %>% 
  mutate(Passive_Smoking = NA)

stat.test

tdata %>% 
  ggplot(aes(Passive_Smoking,Telseq_TL_Ratio,fill=Passive_Smoking))+
  geom_quasirandom(pch=21,size=2,width = 0.3,color="white",stroke=0.5)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.7)+
  scale_fill_jama()+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 13,axis_title_just = 'm',axis_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 15,ticks = T)+
  labs(x = "Passive smoking", y = 'Telomere length ratio, log2(Tumor TL/Normal TL)',title = 'Telomere length')+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_LUAD_passive_smoking_telomere.pdf',width = 3 ,height =6,device = cairo_pdf)  


# 17c ---------------------------------------------------------------------

plotdata <- tdata %>% 
  do(tidy(lm(Total_SBS ~ Passive_Smoking + Age + Gender + Assigned_Population + Tumor_Purity, data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term != '(Intercept)') %>% 
  #mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(label=paste0('β = ',round(estimate,2), ', q-value = ',scientific_format(digits = 3)(p.value))) %>% 
  arrange(term)

plotdata$term
plotdata$term <-  c('Age','Ethnicity: EAS\nRef=EUR', 'Ethnicity: Others\nRef=EUR','Sex: Female\nRef=Male','Passive Smoking: Yes\nRef=No','Tumor Purity')

p <- plotdata %>% 
  mutate(term=fct_reorder(term,estimate)) %>% 
  #mutate(p.value = if_else(p.value>0.05,NA,p.value)) %>% 
  mutate(p.value = if_else(str_detect(term, 'Passive'),p.value,NA)) %>% 
  ggplot(aes(estimate, term, xmin = conf.low, xmax = conf.high, height = 0,color=-log10(p.value))) +
  geom_vline(xintercept = 0, lty = 5) +
  geom_errorbarh(height = .3)+
  geom_point(pch=19,size=4) +
  scale_color_gradient(low = "#984ea3" ,high = "#984ea3",na.value = ncicolpal[8])+
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  ggrepel::geom_text_repel(aes(label=label),size=3.5,nudge_y = 0.4)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid_col = 'gray85',plot_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing = unit(0.2,"cm"),plot.title = element_text(hjust = 0.5))+
  labs(x = "Linear regression coefficient", y = NULL,title='SBS')+
  guides(color="none")+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = 'LCINS_LUAD_passive_smoking_SBS_multivariable.pdf',plot = p,width = 7,height = 4,device = cairo_pdf)

# 17d ---------------------------------------------------------------------

plotdata <- tdata %>% 
  do(tidy(lm(Telseq_TL_Ratio ~ Passive_Smoking + Age + Gender + Assigned_Population  + Tumor_Purity, data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term != '(Intercept)') %>% 
  #mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(label=paste0('β = ',round(estimate,2), ', q-value = ',scientific_format(digits = 3)(p.value))) %>% 
  arrange(term)

plotdata$term
plotdata$term <-  c('Age','Ethnicity: EAS\nRef=EUR', 'Ethnicity: Others\nRef=EUR','Sex: Female\nRef=Male','Passive Smoking: Yes\nRef=No','Tumor Purity')

p <- plotdata %>% 
  mutate(term=fct_reorder(term,estimate)) %>% 
  #mutate(p.value = if_else(p.value>0.05,NA,p.value)) %>% 
  mutate(p.value = if_else(str_detect(term, 'Passive'),p.value,NA)) %>%   ggplot(aes(estimate, term, xmin = conf.low, xmax = conf.high, height = 0,color=-log10(p.value))) +
  geom_vline(xintercept = 0, lty = 5) +
  geom_errorbarh(height = .3)+
  geom_point(pch=19,size=4) +
  scale_color_gradient(low = "#984ea3" ,high = "#984ea3",na.value = ncicolpal[8])+
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  ggrepel::geom_text_repel(aes(label=label),size=3.5,nudge_y = 0.4)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid_col = 'gray85',plot_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing = unit(0.2,"cm"),plot.title = element_text(hjust = 0.5))+
  labs(x = "Linear regression coefficient", y = NULL,title='Telomere length')+
  guides(color="none")+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = 'LCINS_LUAD_passive_smoking_telomere_multivariable.pdf',plot = p,width = 7,height = 4,device = cairo_pdf)



# ED_Figure 18 ------------------------------------------------------------
library(rstatix)
library(ggbeeswarm)

tdata <- bind_rows(
  #sherlock_sbs288_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='SBS'),
  sherlock_dbs78_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='DBS'),
  sherlock_id83_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='ID'),
  sherlock_cn48_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='CN segments'),
  sherlock_sv38_profile %>% group_by(Tumor_Barcode) %>% summarise(Value=sum(Contribution)) %>% mutate(profile='SV')
) %>% 
  left_join(covdata0) %>% 
  left_join(clinical_data %>% select(Tumor_Barcode, Passive_Smoking)) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker) %>% 
  mutate(Value=log10(Value)) %>% 
  filter(is.finite(Value),!is.na(Passive_Smoking)) %>% 
  mutate(profile = factor(profile,levels=c('DBS','ID','CN segments','SV')))

# 18a ---------------------------------------------------------------------
my_comparisons <- list(c("No", "Yes"))

stat.test <- tdata %>% 
  group_by(profile) %>% 
  wilcox_test(Value ~ Passive_Smoking, comparisons = my_comparisons) %>%
  ungroup() %>% 
  add_xy_position(x = "Passive_Smoking",scales = 'free_y',step.increase = 1) %>%
  mutate(myformatted.p = sprintf("P = %.3f",p)) %>% 
  mutate(Passive_Smoking = "Yes")

stat.test

tdata %>% 
  ggplot(aes(Passive_Smoking,Value,fill=Passive_Smoking))+
  geom_quasirandom(pch=21,size=2,width = 0.3,color="gray90",stroke=0.4)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.7)+
  facet_wrap(~profile,scales = 'free_y',nrow = 1)+
  scale_fill_jama()+
  scale_y_continuous(breaks = pretty_breaks(n=7))+
  theme_ipsum_rc(base_size = 15,axis_title_just = 'm',axis_title_size = 16,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 18,ticks = T)+
  labs(x = NULL, y = 'Number of events (log10)')+
  theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=14,face = 'bold',hjust = 0.5),panel.spacing.x = unit(0.2,'cm'),axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_passive_smoking_burden.pdf',width = 7 ,height =5,device = cairo_pdf)  


# 18b ---------------------------------------------------------------------
plotdata <- tdata %>% filter(profile == 'DBS') %>% 
  do(tidy(lm(Value ~ Passive_Smoking + Age + Gender + Assigned_Population + Tumor_Purity + Histology, data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term != '(Intercept)') %>% 
  #mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(label=paste0('β = ',round(estimate,2), ', q-value = ',scientific_format(digits = 3)(p.value))) %>% 
  arrange(term)

plotdata$term
plotdata$term <-  c('Age','Ethnicity: EAS\nRef=EUR', 'Ethnicity: Others\nRef=EUR','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Passive Smoking: Yes\nRef=No','Tumor Purity')


p <- plotdata %>% 
  mutate(term=fct_reorder(term,estimate)) %>% 
  #mutate(p.value = if_else(p.value>0.05,NA,p.value)) %>% 
  mutate(p.value = if_else(str_detect(term, 'Passive'),p.value,NA)) %>% 
  ggplot(aes(estimate, term, xmin = conf.low, xmax = conf.high, height = 0,color=-log10(p.value))) +
  geom_vline(xintercept = 0, lty = 5) +
  geom_errorbarh(height = .3)+
  geom_point(pch=19,size=4) +
  scale_color_gradient(low = "#984ea3" ,high = "#984ea3",na.value = ncicolpal[8])+
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  ggrepel::geom_text_repel(aes(label=label),size=3.5,nudge_y = 0.4)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid_col = 'gray85',plot_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing = unit(0.2,"cm"),plot.title = element_text(hjust = 0.5))+
  labs(x = "Linear regression coefficient", y = NULL,title='DBS')+
  guides(color="none")+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = 'LCINS_passive_smoking_DBS_multivariable.pdf',plot = p,width = 7,height = 4.6,device = cairo_pdf)



# 18c ---------------------------------------------------------------------

plotdata <- tdata %>% filter(profile == 'ID') %>% 
  do(tidy(lm(Value ~ Passive_Smoking + Age + Gender + Assigned_Population + Tumor_Purity + Histology, data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term != '(Intercept)') %>% 
  #mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(label=paste0('β = ',round(estimate,2), ', q-value = ',scientific_format(digits = 3)(p.value))) %>% 
  arrange(term)

plotdata$term
plotdata$term <-  c('Age','Ethnicity: EAS\nRef=EUR', 'Ethnicity: Others\nRef=EUR','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Passive Smoking: Yes\nRef=No','Tumor Purity')


p <- plotdata %>% 
  mutate(term=fct_reorder(term,estimate)) %>% 
  #mutate(p.value = if_else(p.value>0.05,NA,p.value)) %>% 
  mutate(p.value = if_else(str_detect(term, 'Passive'),p.value,NA)) %>% 
  ggplot(aes(estimate, term, xmin = conf.low, xmax = conf.high, height = 0,color=-log10(p.value))) +
  geom_vline(xintercept = 0, lty = 5) +
  geom_errorbarh(height = .3)+
  geom_point(pch=19,size=4) +
  scale_color_gradient(low = "#984ea3" ,high = "#984ea3",na.value = ncicolpal[8])+
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  ggrepel::geom_text_repel(aes(label=label),size=3.5,nudge_y = 0.4)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid_col = 'gray85',plot_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing = unit(0.2,"cm"),plot.title = element_text(hjust = 0.5))+
  labs(x = "Linear regression coefficient", y = NULL,title='ID')+
  guides(color="none")+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = 'LCINS_passive_smoking_ID_multivariable.pdf',plot = p,width = 7,height = 4.6,device = cairo_pdf)



# 18d ---------------------------------------------------------------------

plotdata <- tdata %>% filter(profile == 'CN segments') %>% 
  do(tidy(lm(Value ~ Passive_Smoking + Age + Gender + Assigned_Population + Tumor_Purity + Histology, data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term != '(Intercept)') %>% 
  #mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(label=paste0('β = ',round(estimate,2), ', q-value = ',scientific_format(digits = 3)(p.value))) %>% 
  arrange(term)

plotdata$term
plotdata$term <-  c('Age','Ethnicity: EAS\nRef=EUR', 'Ethnicity: Others\nRef=EUR','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Passive Smoking: Yes\nRef=No','Tumor Purity')


p <- plotdata %>% 
  mutate(term=fct_reorder(term,estimate)) %>% 
  #mutate(p.value = if_else(p.value>0.05,NA,p.value)) %>% 
  mutate(p.value = if_else(str_detect(term, 'Passive'),p.value,NA)) %>% 
  ggplot(aes(estimate, term, xmin = conf.low, xmax = conf.high, height = 0,color=-log10(p.value))) +
  geom_vline(xintercept = 0, lty = 5) +
  geom_errorbarh(height = .3)+
  geom_point(pch=19,size=4) +
  scale_color_gradient(low = "#984ea3" ,high = "#984ea3",na.value = ncicolpal[8])+
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  ggrepel::geom_text_repel(aes(label=label),size=3.5,nudge_y = 0.4)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid_col = 'gray85',plot_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing = unit(0.2,"cm"),plot.title = element_text(hjust = 0.5))+
  labs(x = "Linear regression coefficient", y = NULL,title='CN segements')+
  guides(color="none")+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = 'LCINS_passive_smoking_CN_multivariable.pdf',plot = p,width = 7,height = 4.6,device = cairo_pdf)


# 18e ---------------------------------------------------------------------

plotdata <- tdata %>% filter(profile == 'SV') %>% 
  do(tidy(lm(Value ~ Passive_Smoking + Age + Gender + Assigned_Population + Tumor_Purity + Histology, data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term != '(Intercept)') %>% 
  #mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(label=paste0('β = ',round(estimate,2), ', q-value = ',scientific_format(digits = 3)(p.value))) %>% 
  arrange(term)

plotdata$term
plotdata$term <-  c('Age','Ethnicity: EAS\nRef=EUR', 'Ethnicity: Others\nRef=EUR','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Passive Smoking: Yes\nRef=No','Tumor Purity')


p <- plotdata %>% 
  mutate(term=fct_reorder(term,estimate)) %>% 
  #mutate(p.value = if_else(p.value>0.05,NA,p.value)) %>% 
  mutate(p.value = if_else(str_detect(term, 'Passive'),p.value,NA)) %>% 
  ggplot(aes(estimate, term, xmin = conf.low, xmax = conf.high, height = 0,color=-log10(p.value))) +
  geom_vline(xintercept = 0, lty = 5) +
  geom_errorbarh(height = .3)+
  geom_point(pch=19,size=4) +
  scale_color_gradient(low = "#984ea3" ,high = "#984ea3",na.value = ncicolpal[8])+
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  ggrepel::geom_text_repel(aes(label=label),size=3.5,nudge_y = 0.4)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid_col = 'gray85',plot_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing = unit(0.2,"cm"),plot.title = element_text(hjust = 0.5))+
  labs(x = "Linear regression coefficient", y = NULL,title='SV')+
  guides(color="none")+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = 'LCINS_passive_smoking_SV_multivariable.pdf',plot = p,width = 7,height = 4.6,device = cairo_pdf)


# 18f ---------------------------------------------------------------------
testdata <- ludmil_activity_all_obs %>%
  select(Tumor_Barcode,!starts_with('SBS')) %>% 
  pivot_longer(-Tumor_Barcode) %>% 
  left_join(clinical_data %>% select(Tumor_Barcode,Passive_Smoking)) %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,!is.na(Passive_Smoking)) %>% 
  group_by(name) %>% 
  # do(tresult = safely(stats::fisher.test)(.$value,.$Passive_Smoking)) %>% 
  do(tresult = safely(stats::glm)(value ~ Passive_Smoking + Assigned_Population + Histology + Gender + Age + Tumor_Purity,family = binomial(),data=.)) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(name,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='Passive_SmokingYes') %>% 
  filter(abs(log(estimate))<10) %>% 
  mutate(profile=str_remove_all(name,'[:digit:].*')) %>% 
  group_by(profile) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  ungroup() %>% 
  mutate(profile=factor(profile,levels=c('DBS','ID','CN','SV'),labels= c('DBS Signatures','ID Signatures','CN Signatures','SV Signatures'))) %>% 
  mutate(dir=if_else(estimate>1,'Enriched in passive smokers','Enriched in non-passive smokers'))

testdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR),fill=dir))+
  geom_hline(yintercept = -log10(0.05),col=ncicolpal[2],linetype=2)+
  geom_hline(yintercept = -log10(0.01),col=ncicolpal[1],linetype=2)+
  geom_vline(xintercept = 0,col="gray10",linetype=1,linewidth=0.5)+
  geom_point(pch=21,size=3.5,col='black',stroke=0.2)+
  ggrepel::geom_text_repel(data=testdata %>% filter(p.value<0.05),aes(label=name),force=20)+
  facet_wrap(~profile,nrow = 1)+
  scale_fill_jama()+
  scale_color_manual(values = sigcol)+
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  scale_y_continuous(breaks = pretty_breaks(n=7))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 12,grid = F,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing.x = unit('0.25','cm'),strip.text.x = element_text(face = 'bold',hjust = 0.5),legend.position = 'bottom')+
  labs(x = 'Odd ratio (log2)', y = '-log10(FDR)',fill='Enrichment')+
  #guides(fill="none",color='none')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  coord_cartesian(clip = 'off')

ggsave(filename = 'LCINS_passive_smoking_signature_enrichment_others.pdf',width = 15,height = 4,device = cairo_pdf)


# ED_Figure 19 ------------------------------------------------------------
tdata <- bind_rows(
  sherlock_cn68_profile %>% mutate(Profile='CN segements'),
  sherlock_sv38_profile %>% mutate(Profile='SV')
) %>% 
  group_by(Tumor_Barcode,Profile) %>% 
  summarise(Burden=log10(sum(Contribution))) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Profile,values_from = Burden) %>% 
  left_join(pollution_data) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker) %>% 
  filter(!is.na(Pollution_group2))

my_comparisons <- list(c("Low", "High"))

# 19a ---------------------------------------------------------------
plotdata <- tdata %>% select(Tumor_Barcode,`CN segements`, SV,Pollution_group2) %>% pivot_longer(cols = -c(Tumor_Barcode,Pollution_group2)) %>% filter(is.finite(value)) %>% 
  mutate(name=factor(name,levels=c('CN segements','SV')))

stat.test <- plotdata %>% 
  group_by(name) %>% 
  wilcox_test(value ~ Pollution_group2, comparisons = my_comparisons) %>%
  ungroup() %>% 
  add_xy_position(x = "Pollution_group2") %>%
  mutate(myformatted.p = sprintf("P = %.2e",p)) %>% 
  mutate(Pollution_group2 = "High")

stat.test

plotdata %>% 
  ggplot(aes(Pollution_group2,value,fill=Pollution_group2))+
  geom_quasirandom(pch=21,size=2,width = 0.3,color="white",stroke=0.4)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.7)+
  facet_wrap(~name,scales = 'free')+
  scale_fill_manual(values = c('#01665e','#BB0E3D'))+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 13,axis_title_just = 'm',axis_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 15,ticks = T)+
  labs(x = expression("Pollution group (Population weighted PM"[2.5]*")"), y = 'Total number of events (log)')+
  theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=14,face = 'bold',hjust = 0.5),panel.spacing = unit(0.4,'cm'))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_Pollution_mutations_others.pdf',width = 4.6 ,height =6,device = cairo_pdf)  


# 19b ---------------------------------------------------------------------
plotdata <- tdata %>% select(Tumor_Barcode,`CN segements`,SV,Pollution_group2) %>% pivot_longer(cols = -c(Tumor_Barcode,Pollution_group2)) %>% filter(is.finite(value)) %>% 
  mutate(name=factor(name,levels=c('CN segements','SV'))) %>% 
  left_join(covdata0)

plotdata <- plotdata %>% 
  group_by(name) %>% 
  do(tidy(lm(value ~ Pollution_group2 + Age + Gender + Assigned_Population + Histology + Tumor_Purity, data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term != '(Intercept)') %>% 
  #mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(label=paste0('β = ',round(estimate,2), ', q-value = ',scientific_format(digits = 3)(p.value))) %>% 
  arrange(term) %>% 
  ungroup()

plotdata$term
plotdata$term <-  rep(c('Age','Ethnicity: EAS\nRef=EUR', 'Ethnicity: Others\nRef=EUR','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Pollution Group: High\nRef=Low','Tumor Purity'),each=2)

p <- plotdata %>% 
  mutate(term=fct_reorder(term,estimate)) %>% 
  #mutate(p.value = if_else(p.value>0.05,NA,p.value)) %>% 
  mutate(p.value = if_else(str_detect(term,'Pollution'),p.value,NA)) %>% 
  ggplot(aes(estimate, term, xmin = conf.low, xmax = conf.high, height = 0,color=-log10(p.value))) +
  geom_vline(xintercept = 0, lty = 5) +
  geom_errorbarh(height = .3)+
  geom_point(pch=19,size=3) +
  facet_wrap(~name,nrow = 1)+
  scale_color_gradient(low = "#B71B1BFF" ,high = "#B71B1BFF",na.value = ncicolpal[8])+ #"#EE5250FF"
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  ggrepel::geom_text_repel(aes(label=label),size=3.5,nudge_y = 0.4)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid_col = 'gray85',plot_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing = unit(0.3,"cm"),plot.title = element_text(hjust = 0.5),strip.text.x = element_text(face = 'bold',hjust = 0.5))+
  labs(x = "Linear regression coefficient", y = NULL)+
  guides(color="none")+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = 'LCINS_Pollution_mutations_multivariable_others.pdf',plot = p,width = 12,height = 5,device = cairo_pdf)

# Figure 19c ---------------------------------------------------------------
testdata <- ludmil_activity_all_obs %>%
  pivot_longer(-Tumor_Barcode) %>% 
  mutate(Profile=str_remove_all(name,"[0-9a-z]")) %>% 
  filter(Profile %in% c('CN','SV')) %>% 
  left_join(pollution_data) %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,!is.na(PM25)) %>% 
  mutate(PM25 = PM25/10) %>% 
  #mutate(value=as.factor(value)) %>% 
  group_by(name,Profile) %>% 
  do(tresult = safely(stats::glm)(value ~ PM25+Histology+Assigned_Population + Gender + Age + Tumor_Purity,family = binomial(),data=.)) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate = TRUE))) %>% 
  select(name,Profile,fit) %>% 
  unnest(cols = c(fit)) %>%
  filter(term=='PM25') %>% 
  group_by(Profile) %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  mutate(Profile=factor(Profile,levels=c('CN','SV'),labels=c('CN signatures','SV signatures'))) %>% 
  mutate(dir=if_else(estimate>1,'Enrich in high PM2.5','Enriched in low PM2.5')) 
  #mutate(dir=if_else(estimate>1,expression("Enriched in high PM"[2.5]*),expression("Enriched in low PM"[2.5]*)))

testdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR),fill=dir))+
  geom_hline(yintercept = -log10(0.05),col=ncicolpal[2],linetype=2)+
  geom_hline(yintercept = -log10(0.01),col=ncicolpal[1],linetype=2)+
  geom_vline(xintercept = 0,col="gray60",linetype=1,linewidth=0.5)+
  geom_point(pch=21,size=3.5,col='black',stroke=0.2)+
  ggrepel::geom_text_repel(data=testdata %>% filter(p.value<0.05),aes(label=name),force=20)+
  facet_wrap(~Profile,nrow = 1)+
  scale_fill_manual(values = as.character(pollution_colors[c(3,1)]))+
  scale_color_manual(values = sigcol)+
  scale_x_continuous(breaks = pretty_breaks(n = 5))+
  scale_y_continuous(breaks = pretty_breaks(n=7))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid = F,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_text(hjust = 0.5,face = 'bold',size = 14),panel.spacing = unit(0.3,'cm'),legend.position = 'bottom')+
  labs(x = expression("Odd ratio (" * 10 ~ mu*g/m^3 * " of PM" [2.5] * ", log2)"), y = '-log10(FDR)', fill='Enrichement')+
 # guides(fill="none",color='none')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  coord_cartesian(clip = 'off')

ggsave(filename = 'LCINS_pollution_signature_enrichment_others.pdf',width = 7,height = 4,device = cairo_pdf)


# ED_Figure 21 ------------------------------------------------------------

tdata <- sherlock_sbs96_profile %>%
  group_by(Tumor_Barcode) %>% 
  summarise(Burden=(sum(Contribution))) %>% 
  ungroup() %>% 
  left_join(ludmil_activity_all %>% select(Tumor_Barcode,SBS4)) %>%
  mutate(SBS4=replace_na(SBS4,0)) %>% 
  mutate(value=log10(Burden - SBS4)) %>% 
  select(Tumor_Barcode,value) %>% 
  left_join(pollution_data) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker) %>% 
  filter(!is.na(Pollution_group2),!is.na(value))

my_comparisons <- list(c("Low", "High"))

# 21a ---------------------------------------------------------------
plotdata <- tdata %>% select(Tumor_Barcode,value,Pollution_group2) 
  #mutate(name=factor(name,levels=c('SBS','DBS','ID')))
stat.test <- plotdata %>% 
  wilcox_test(value ~ Pollution_group2, comparisons = my_comparisons) %>%
  ungroup() %>% 
  add_xy_position(x = "Pollution_group2") %>%
  mutate(myformatted.p = sprintf("P = %.4f",p)) %>% 
  mutate(Pollution_group2 = "High")

stat.test

plotdata %>% 
  ggplot(aes(Pollution_group2,value,fill=Pollution_group2))+
  geom_quasirandom(pch=21,size=2,width = 0.3,color="white",stroke=0.4)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.7)+
  scale_fill_manual(values = c('#01665e','#BB0E3D'))+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 13,axis_title_just = 'm',axis_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 15,ticks = T)+
  labs(x = expression("Pollution group (Population weighted PM"[2.5]*")"), y = 'Total number of mutations (log10)',title = 'SBS (excluded SBS4 mutations) ')+
  theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=14,face = 'bold',hjust = 0.5),panel.spacing = unit(0.4,'cm'))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_Pollution_mutations_nonSBS4_mutations.pdf',width = 4.2 ,height =6,device = cairo_pdf)  


# 21b ---------------------------------------------------------------------
plotdata <- tdata %>% select(Tumor_Barcode,value,Pollution_group2) %>% left_join(covdata0)

plotdata <- plotdata %>% 
  do(tidy(lm(value ~ Pollution_group2 + Age + Gender + Assigned_Population + Histology + Tumor_Purity, data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term != '(Intercept)') %>% 
  #mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(label=paste0('β = ',round(estimate,2), ', q-value = ',scientific_format(digits = 3)(p.value))) %>% 
  arrange(term) %>% 
  ungroup()

plotdata$term
plotdata$term <-  rep(c('Age','Ethnicity: EAS\nRef=EUR', 'Ethnicity: Others\nRef=EUR','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Pollution Group: High\nRef=Low','Tumor Purity'),each=1)

p <- plotdata %>% 
  mutate(term=fct_reorder(term,estimate)) %>% 
  #mutate(p.value = if_else(p.value>0.05,NA,p.value)) %>% 
  mutate(p.value = if_else(str_detect(term,'Pollution'),p.value,NA)) %>% 
  ggplot(aes(estimate, term, xmin = conf.low, xmax = conf.high, height = 0,color=-log10(p.value))) +
  geom_vline(xintercept = 0, lty = 5) +
  geom_errorbarh(height = .3)+
  geom_point(pch=19,size=3) +
  scale_color_gradient(low = "#B71B1BFF" ,high = "#B71B1BFF",na.value = ncicolpal[8])+ #"#EE5250FF"
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  ggrepel::geom_text_repel(aes(label=label),size=3.5,nudge_y = 0.4)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid_col = 'gray85',plot_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing = unit(0.3,"cm"),plot.title = element_text(hjust = 0.5),strip.text.x = element_text(face = 'bold',hjust = 0.5))+
  labs(x = "Linear regression coefficient", y = NULL,title='SBS (exclude SBS4 mutations)')+
  guides(color="none")+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = 'LCINS_Pollution_mutations_multivariable_nonSBS4_mutations.pdf',plot = p,width = 6.5,height = 5,device = cairo_pdf)

# 21c ---------------------------------------------------------------
my_comparisons <- list(c("Low", "High"),c("Low", "Intermediate"),c("High", "Intermediate"))

plotdata <- tdata %>% select(Tumor_Barcode,value,Pollution_group3) 
#mutate(name=factor(name,levels=c('SBS','DBS','ID')))
stat.test <- plotdata %>% 
  wilcox_test(value ~ Pollution_group3, comparisons = my_comparisons) %>%
  ungroup() %>% 
  add_xy_position(x = "Pollution_group3",step.increase = 0.3) %>%
  mutate(myformatted.p = sprintf("P = %.2e",p)) %>% 
  mutate(Pollution_group3 = "High")

stat.test

plotdata %>% 
  ggplot(aes(Pollution_group3,value,fill=Pollution_group3))+
  geom_quasirandom(pch=21,size=2,width = 0.3,color="white",stroke=0.4)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.7)+
  scale_fill_manual(values = pollution_colors)+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 13,axis_title_just = 'm',axis_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 15,ticks = T)+
  labs(x = expression("Pollution group (Population weighted PM"[2.5]*")"), y = 'Total number of mutations (log10)',title = 'SBS (excluded SBS4 mutations)')+
  theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=14,face = 'bold',hjust = 0.5),panel.spacing = unit(0.4,'cm'))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_Pollution_mutations_nonSBS4_mutations2.pdf',width = 4.5 ,height =6,device = cairo_pdf)  


# 21d ---------------------------------------------------------------------
testdata <- tdata %>% left_join(covdata0)

plabel <- testdata %>% 
  do(tidy(lm(value ~ PM25 + Age + Gender + Assigned_Population + Histology + Tumor_Purity, data = .))) %>% 
  filter(term=='PM25') %>% 
  mutate(label=paste0('Multiple linear regression\nβ = ',round(estimate,2),sprintf(", q-value = %.2e",p.value))) %>% 
  pull(label)

tdata_group <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(value,na.rm = T),nsample=n_distinct(Tumor_Barcode)) %>% 
  mutate(Burden= if_else(Country_pollution == 'Azerbaijan',4.2,Burden))

tmpbreak <- pretty_breaks()(tdata_group$Burden)
tmplabel <- tmpbreak
tmplabel[length(tmplabel)] <- paste0('>',tmplabel[length(tmplabel)])

tdata_group %>% 
  ggplot(aes(mean_pollution,(Burden)))+
  stat_smooth(method="lm",fullrange=TRUE)+
  geom_point(aes(size=nsample),pch=21,fill='gray20',col='white',stroke=0.2)+
  ggrepel::geom_text_repel(aes(label=Country_pollution),size=3.2,col='gray30',force = 30,segment.color='gray35',max.overlaps = 30)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_size_binned(breaks = c(0,40,80,120,160,200))+
  scale_y_continuous(breaks = tmpbreak,labels = tmplabel)+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = FALSE,axis = 'XY',ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  labs(x = expression("Average population weighted PM"[2.5] ~ (mu*g/m^3)), y = 'Average number of SBS mutations (log10)',size='Number of samples',title = 'SBS (excluded SBS4 mutations)')+
  theme(legend.position = 'top',legend.key.width = unit(1,'cm'),legend.margin=margin(0,0,-12,0),legend.box.margin=margin(0,0,0,0),plot.title = element_text(hjust = 0.5,size = 15))+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  annotate("text",x=28,y=3,label=plabel,family= 'Roboto Condensed')

ggsave(filename = 'Pollution_assocaition_SBS_noSBS4_mutations.pdf',width = 5,height = 5.2,device = cairo_pdf)



