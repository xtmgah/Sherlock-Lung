set_wd()
libztw()
pdfhr2()
myggstyle()

library(gtsummary)
sci_notation <- function(x) ifelse(is.na(x), NA, format(x, digits = 3, scientific = TRUE))


source('/Users/zhangt8/NIH-Work/MutationSignature/mSigPortal/CBIIT/mSigPortal/Codes/Sigvisualfunc.R')
source('~/NIH-Work/R/ZTW_function/ztw.R')


# load common data files  -------------------------------------------------

load('../../BBsolution_final3_short.RData')
load('../../sp_group_data.RData')
load('../../covdata0.RData')
load('../../Clinical/clinical_data.RData')
load('../../ZTW_functions.RData')
load('../../sherlock_PGA.RData')
load('../../RDS/sherlock_data_all.RData')
load('../../RDS/sherlock_variable.RData')
load('../../MAFtools/sherlock_driver_mutations.RData')

load('../Signature_Lumidl.RData',verbose = T)
load('../sherlock_profiles.RData',verbose = T)
load('../Signature_Lumidl_CN.RData',verbose = T)
load('../Signature_Lumidl_SV.RData',verbose = T)

load('../sampleset.RData')
load('../sigcol.RData')

load('../../verifybamID.RData')

verifybamID <- wgs_groups_info %>% select(Tumor_Barcode,Barcode=Normal_Barcode) %>% left_join(verifybamID) %>% select(Tumor_Barcode,PC1,PC2)


drglist <- readRDS('../../../../Collaborators/Nuria/Update2/drivers_intogene.RDS') %>% pull(symbol)
drglist <- sherlock_driver_mutations %>% filter(Hugo_Symbol %in% drglist,Tumor_Barcode %in% sherlock_nonsmoker) %>% select(Hugo_Symbol,Tumor_Barcode) %>% unique() %>% count(Hugo_Symbol) %>% mutate(freq=n/length(sherlock_nonsmoker)) %>% filter(freq>0.02) %>% pull(Hugo_Symbol)


#The presence of a specific signature was considered for the enrichment analysis (or the dichotomization of values above and below the median of assigned mutations, if the signature was present in above 50% of the samples). 
ludmil_activity_all <- ludmil_activity_all %>% left_join(cn68_activity) %>% left_join(sv38_activity) %>% mutate_all(~replace_na(., 0L)) %>% filter(Tumor_Barcode %in% sherlock_nonsmoker)

# remove signature not present in non-smokers
exclude_sigs <- ludmil_activity_all %>% pivot_longer(-Tumor_Barcode) %>% group_by(name) %>% summarise(value=sum(value)) %>% filter(value==0) %>% pull(name)
ludmil_activity_all <- ludmil_activity_all %>% select(-one_of(exclude_sigs))


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


# define the geographical regions
pollution_data %>% select(Tumor_Barcode, Country) %>% left_join(wgs_groups_info %>% select(Tumor_Barcode,Assigned_Population)) %>% count(Country,Assigned_Population)

East <- c("Hong Kong","Korea","Philippines","Taiwan")
Unknown <- c('Azerbaijan','Tunisia','Uzbekistan')

country_group <- pollution_data %>% 
  select(Tumor_Barcode,Country) %>% 
  mutate(Country_Group=case_when(
    Country %in% East ~ 'AS',
    Country %in% Unknown ~ 'Other',
    is.na(Country) ~ 'Other',
    TRUE ~ 'NA/EU'
  ))

# remove the EAS Canada 
country_group2 <- country_group %>% 
  left_join(
    wgs_groups_info %>% select(Tumor_Barcode,Assigned_Population)
  ) %>% 
  filter(!(Country == 'Canada_Ontario' & Assigned_Population == 'AS')) 



# Figure 1 ----------------------------------------------------------------

# Figure 1a ---------------------------------------------------------------
library(maps)
library(sf)

df <- readRDS('../../Signature_ludmil2/Barcode2map.RDS')
tmpx <- df %>% 
  left_join(covdata0) %>% 
  filter(Smoking=='Non-Smoker') %>% 
  count(Country,CURRENT_LAST_RES,City) %>% 
  filter(!is.na(Country))


tmpx.1 <- tmpx %>% filter(City=="Yes") %>% filter(n<5) %>% group_by(Country) %>% summarise(n=sum(n)) %>% mutate(CURRENT_LAST_RES=NA_character_, City="No")

tmpx.2 <- tmpx %>% filter(City!="Yes") %>% bind_rows(tmpx.1) %>% group_by(Country,CURRENT_LAST_RES,City) %>% summarise(n=sum(n))

tmpx <- tmpx %>% filter(City=="Yes") %>% filter(n>=5) %>% bind_rows(tmpx.2) %>% arrange(Country)


tmpx %>% 
  mutate(name=if_else(City=='Yes',CURRENT_LAST_RES,Country)) %>% 
  filter(!is.na(name)) %>% 
  arrange(Country) %>% 
  mutate(name=fct_inorder(name)) %>% 
  ggplot(aes("x",name,fill=n))+
  geom_point(pch=21,size=7)+
  scale_fill_viridis_c(direction = -1,limits = c(5, 200))+
  geom_blank()+
  theme(legend.position = 'top',legend.key.width = unit(2.5,'cm'),legend.key.height = unit(0.3,'cm'))

ggsave(filename = 'tmp.pdf',width = 8,height = 12,device = cairo_pdf)




df <- df %>% count(Country,SP_Group_New) %>% pivot_wider(names_from = SP_Group_New,values_from = n,values_fill = 0) %>% filter(!is.na(Country)) %>% mutate(sample_size= Others + AS_N + EU_N + EU_S)


sdf <- rnaturalearthdata::countries50 %>% 
  st_as_sf() %>% 
  # st_make_valid() %>% 
  #st_crop(xmin = -24, xmax = 31, ymin = 33, ymax = 73) %>% 
  #filter(admin %in% df$Country) %>% 
  left_join(df, by = c("admin" = "Country"))

ranking <- st_geometry(sdf) %>% 
  st_point_on_surface() %>% 
  st_coordinates() %>% 
  as_tibble() %>% 
  bind_cols(tibble(country = sdf$admin, AS_N=sdf$AS_N, EU_N=sdf$EU_N, EU_S=sdf$EU_S, Others=sdf$Others,sample_size=sdf$sample_size )) %>% #SP_Group_New=sdf$SP_Group_New
  filter(country %in% df$Country)



# map
p <- ggplot() + 
  geom_sf(data = sdf,fill='gray70', size = .2, color = 'white') +
  geom_point(data = ranking, aes(x = X, y = Y,fill=sample_size),size = 4.5,pch=21,stroke=0.2,col='black')+
  #ggrepel::geom_text_repel(data = ranking, aes(x = X, y = Y, label = country),size=3.5,fontface=2) +
  #geom_text(data = ranking, aes(x = X-.5, y = Y, label = country, color = SP_Group_New), hjust = 1, size = 2.5, nudge_y = .5) +
  scale_fill_viridis_c(direction = -1)+
  #coord_sf(xlim=c(-150,150),ylim=c(5, 60),clip = "off") +
  coord_sf(xlim=c(-120,120),ylim=c(30, 60),clip = "off") +
  theme_void(base_family = "Roboto Condensed") +
  guides(fill = "none")+
  theme(plot.margin = margin(.5, 1, .5, .5, "cm"),
        legend.position = "top", legend.title = element_text(face = 2),
        #plot.background = element_rect(fill = "black"),
        plot.caption = element_text(color = "gray40"),
        plot.title = element_text(color = "gray40", size = 16, family = "Roboto Condensed", face = "bold"),
        plot.subtitle = element_text(color = "gray40", size = 8))

ggsave("WorldMap_Samples0.pdf", plot=p, width=10, height=6, device=cairo_pdf())



p <- ggplot() + 
  geom_sf(data = sdf, size = .2, fill = "transparent", color = 'gray60') +
  #geom_point(data = ranking, aes(x = X, y = Y), fill=ncicolpal[2],size = 4,pch=21,stroke=0.2)+
  geom_rect(data = ranking, aes(xmin = X-4, xmax=X-2, ymin = Y, ymax=Y+AS_N*0.25), fill= sp_group_color_new['AS_N'])+
  geom_rect(data = ranking, aes(xmin = X-2, xmax=X, ymin = Y, ymax=Y+EU_N*0.25), fill= sp_group_color_new['EU_N'])+
  geom_rect(data = ranking, aes(xmin = X, xmax=X+2, ymin = Y, ymax=Y+EU_S*0.25), fill= sp_group_color_new['EU_S'])+
  geom_rect(data = ranking, aes(xmin = X+2, xmax=X+4, ymin = Y, ymax=Y+Others*0.25), fill= sp_group_color_new['Others'])+
  geom_rect(data = ranking, aes(xmin = X-4, xmax=X+4, ymin = Y, ymax=Y+0.1), fill= ncicolpal[7])+
  ggrepel::geom_text_repel(data = ranking, aes(x = X, y = Y, label = country),size=3,force = 2,nudge_y = -1) +
  #geom_text(data = ranking, aes(x = X-.5, y = Y, label = country), hjust = 0.5, size = 2.5, nudge_y = -1) +
  scale_fill_manual(values = sp_group_color_new)+
  coord_sf(xlim=c(-120,120),ylim=c(30, 60),clip = "off") +
  theme_void() +
  theme(plot.margin = margin(.5, 1, .5, .5, "cm"),
        #legend.position = "none",
        #plot.background = element_rect(fill = "black"),
        plot.caption = element_text(color = "gray40"),
        plot.title = element_text(color = "gray40", size = 16, family = "Roboto Condensed", face = "bold"),
        plot.subtitle = element_text(color = "gray40", size = 8))




ggsave("WorldMap_Samples.pdf", plot=p, width=12, height=8, device=cairo_pdf())



# Figure 1b  -------------------------------------------------------------
# Sanity plot for mutational signature paper N=871
library(ggsankey)
tdata0 <- covdata0 %>% 
  mutate(Histology=if_else(Histology == 'Other','Others',Histology)) %>% 
  left_join(
    ludmil_activity_all_obs %>% select(Tumor_Barcode,SBS4) %>% mutate(SBS4=if_else(SBS4,'SBS4','Non-SBS4'))
  ) %>% 
  left_join(
    clinical_data %>% select(Tumor_Barcode,Passive_Smoking) %>% mutate(Passive_Smoking = if_else(is.na(Passive_Smoking),'Unknown',Passive_Smoking))
  ) %>% 
  left_join(country_group) %>% 
  mutate(Country_Group=if_else(Country_Group == 'Other','Other_region',Country_Group)) %>% 
  select(Tumor_Barcode,Histology,Assigned_Population,Region=Country_Group,Gender,Smoking,Passive_Smoking) %>% 
  #mutate(Smoking_Signature=SBS4) %>% 
  filter(Smoking=='Non-Smoker') 
#%>%  mutate(Histology=if_else(Histology == 'Other','Others',Histology))

tdata0 %>% count(Histology,Assigned_Population,Region,Gender,Passive_Smoking) %>% View()

tdata <- tdata0 %>% make_long(Histology,Assigned_Population,Region,Gender,Passive_Smoking)
tdata_nr <- 
  tdata %>% 
  filter(!is.na(node)) %>% 
  group_by(x, node)%>% 
  summarise(count = n())

tdata <- tdata %>% left_join(tdata_nr)

tdata$node <- factor(tdata$node,levels = c('Adenocarcinoma','Squamous cell carcinoma','Carcinoid tumor','Others', 'EUR','EAS','Other',"NA/EU","AS","Other_region",'Female','Male','Yes','No','Unknown'))

ggplot(tdata, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6, node.color = "gray30",space = 35,width = 0.13) +
  geom_sankey_text(size = 4, color = "black") +
  geom_sankey_text(aes(label = count), size = 3.5, check_overlap = TRUE,) +
  scale_fill_manual(values = ncicolpal[c(4,1,3,2,5,6,11,13,15,12,18,10,9,20,8)])+
  #scale_fill_d3(palette = 'category20') +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme_ipsum_rc(axis = FALSE,axis_text_size = 16)+
  theme(axis.text.y = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

ggsave(filename = 'Sherlock_871_WGS_dataset.pdf',width = 8,height = 7,device = cairo_pdf)





# Figure 1c ----------------------------------------------------------------------
readRDS(file = '../../Signature_ludmil2/Barcode2map.RDS')-> b2map

b2map <- b2map %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker) %>% 
  mutate(Country = if_else(Country == 'Hong Kong S.A.R.','Hong Kong',Country)) %>% 
  mutate(Country = if_else(Country == 'United States of America','USA',Country)) %>% 
  mutate(Country = if_else(Country == 'Czech Republic','Czech',Country)) %>% 
  mutate(Country = if_else(Country == 'Republic of Serbia','Serbia',Country)) %>% 
  mutate(Country = if_else(Country %in% c('Hong Kong','Taiwan'),Country,str_to_upper(Country))) %>% 
  mutate(Country = if_else(Country == 'Taiwan','Taipei', Country)) %>% 
  mutate(CURRENT_LAST_RES = if_else(CURRENT_LAST_RES == 'CT','Connecticut',CURRENT_LAST_RES)) %>% 
  mutate(CURRENT_LAST_RES = if_else(CURRENT_LAST_RES == 'MA','Massachusetts',CURRENT_LAST_RES)) %>% 
  mutate(CURRENT_LAST_RES = if_else(CURRENT_LAST_RES == 'NEW YORK','New York',CURRENT_LAST_RES)) %>% 
  mutate(CURRENT_LAST_RES = if_else(CURRENT_LAST_RES == 'MINNESOTA','Minnesota',CURRENT_LAST_RES)) %>% 
  mutate(CURRENT_LAST_RES = if_else(CURRENT_LAST_RES == 'FLORIDA','Florida',CURRENT_LAST_RES))


tmp <- b2map %>% count(Country,sort=T) %>% filter(!is.na(Country),n>=5) %>% pull(Country)

b2map <- b2map %>% 
  mutate(Country = if_else(!is.na(Country) & Country %in% tmp, Country, 'Others')) %>% 
  mutate(CURRENT_LAST_RES = if_else(Country == 'USA',CURRENT_LAST_RES,NA_character_)) %>% 
  mutate(CURRENT_LAST_RES = if_else(is.na(CURRENT_LAST_RES) | CURRENT_LAST_RES %in% c('Connecticut','Massachusetts','New York','Florida','Minnesota'),CURRENT_LAST_RES,'Others')) %>% 
  mutate(CURRENT_LAST_RES = if_else(Country == 'USA' & is.na(CURRENT_LAST_RES),"Others",CURRENT_LAST_RES)) %>% 
  select(Tumor_Barcode,Country,State=CURRENT_LAST_RES) 


genomesize <-  genome2size('GRCh38')
b2map <- b2map %>% 
  left_join(
    BBsolution4 %>% select(Tumor_Barcode,Mutations) %>% mutate(Burden=log10((Mutations)/genomesize))
  ) %>% 
  left_join(
    covdata0 %>% select(Tumor_Barcode,Histology)
  ) %>% 
  left_join(sherlock_PGA) %>% 
  left_join(
    sherlock_variable %>% filter(name=='SV') %>% pivot_wider() %>% mutate(SV=log2(SV))
  )

sgroupcol <- ncicolpal[1:4]
names(sgroupcol) <- c('Squamous cell carcinoma','Other','Carcinoid tumor','Adenocarcinoma')

tdata <-  b2map %>% select(Cancer_Type=Country,Sample=Tumor_Barcode,Sample_Group=Histology,Burden)
p1 <- TMBplot2(tdata,sgroupcol = sgroupcol )+ theme(axis.text.x.bottom = element_text(size = 8,angle = 45,hjust = 1,vjust = 1))
global_levels <- tdata %>% group_by(Cancer_Type) %>% summarise(Burden=median(Burden)) %>% arrange((Burden)) %>% pull(Cancer_Type)
tdata <-  b2map %>% filter(!is.na(State))%>% select(Cancer_Type=State,Sample=Tumor_Barcode,Sample_Group=Histology,Burden)
p2 <- TMBplot2(tdata,sgroupcol = sgroupcol )+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+theme(axis.text.x.bottom = element_text(size = 8,angle = 45,hjust = 1,vjust = 1))
usa_levels <- tdata %>% group_by(Cancer_Type) %>% summarise(Burden=median(Burden)) %>% arrange((Burden)) %>% pull(Cancer_Type)
#p <- plot_grid(p1,NULL,p2,align = 'h',nrow = 1,axis = 'tb',rel_widths = c(3,-0.36,1.5))
#ggsave(filename = 'tmp.pdf',plot = p,width = 8,height = 4.5,device = cairo_pdf)

# PGA
tdata <-  b2map %>% select(Cancer_Type=Country,Sample=Tumor_Barcode,Sample_Group=Histology,Burden=PGA_WGD)
p3 <- TMBplot2(tdata,sgroupcol = sgroupcol,cancer_type_levels = global_levels )+ theme(axis.text.x.bottom = element_text(size = 8,angle = 45,hjust = 1,vjust = 1))+scale_y_continuous(limits = c(0,1),expand = c(0,0),breaks = seq(0,1,0.2))+ylab('Percentage of the Genome Aberrated\n')
tdata <-  b2map %>% filter(!is.na(State))%>% select(Cancer_Type=State,Sample=Tumor_Barcode,Sample_Group=Histology,Burden=PGA_WGD)
p4 <- TMBplot2(tdata,sgroupcol = sgroupcol,cancer_type_levels = usa_levels )+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+theme(axis.text.x.bottom = element_text(size = 8,angle = 45,hjust = 1,vjust = 1))+scale_y_continuous(limits = c(0,1),expand = c(0,0),breaks = seq(0,1,0.2))
#p <- plot_grid(p1,NULL,p2,align = 'h',nrow = 1,axis = 'tb',rel_widths = c(3,-0.36,1.5))
#ggsave(filename = 'tmp.pdf',plot = p,width = 8,height = 4.5,device = cairo_pdf)

# SV
tdata <-  b2map %>% select(Cancer_Type=Country,Sample=Tumor_Barcode,Sample_Group=Histology,Burden=SV)
p5 <- TMBplot2(tdata,sgroupcol = sgroupcol,cancer_type_levels = global_levels )+ theme(axis.text.x.bottom = element_text(size = 8,angle = 45,hjust = 1,vjust = 1))+scale_y_continuous(limits = c(0,10.5),expand = c(0,0),breaks = seq(0,10,2))+ylab('Number of SVs (log2)\n')
tdata <-  b2map %>% filter(!is.na(State))%>% select(Cancer_Type=State,Sample=Tumor_Barcode,Sample_Group=Histology,Burden=SV)
p6 <- TMBplot2(tdata,sgroupcol = sgroupcol,cancer_type_levels = usa_levels )+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+theme(axis.text.x.bottom = element_text(size = 8,angle = 45,hjust = 1,vjust = 1))+scale_y_continuous(limits = c(0,10.5),expand = c(0,0),breaks = seq(0,10,2))


pall1 <- plot_grid(p1,NULL,p2,align = 'h',nrow = 1,axis = 'tb',rel_widths = c(3,-0.08,1.5))
pall2 <- plot_grid(p3,NULL,p4,align = 'h',nrow = 1,axis = 'tb',rel_widths = c(3,-0.08,1.5))
pall3 <- plot_grid(p5,NULL,p6,align = 'h',nrow = 1,axis = 'tb',rel_widths = c(3,-0.08,1.5))
p <- plot_grid(pall1,pall2, pall3, align = 'v',ncol = 1,axis = 'lr')

ggsave(filename = 'TMB_PGA_by_Country.pdf',plot = p,width = 8,height = 11,device = cairo_pdf)
#ggsave(filename = 'TMB_PGA_by_Country.pdf',plot = p,width = 6.5,height = 9,device = cairo_pdf)

# LUAD only
tdata <-  b2map %>% select(Cancer_Type=Country,Sample=Tumor_Barcode,Sample_Group=Histology,Burden) %>% filter(Sample_Group == 'Adenocarcinoma')
p1 <- TMBplot2(tdata,sgroupcol = sgroupcol )+ theme(axis.text.x.bottom = element_text(size = 8,angle = 45,hjust = 1,vjust = 1))
global_levels <- tdata %>% group_by(Cancer_Type) %>% summarise(Burden=median(Burden)) %>% arrange((Burden)) %>% pull(Cancer_Type)
tdata <-  b2map %>% filter(!is.na(State))%>% select(Cancer_Type=State,Sample=Tumor_Barcode,Sample_Group=Histology,Burden) %>% filter(Sample_Group == 'Adenocarcinoma')
p2 <- TMBplot2(tdata,sgroupcol = sgroupcol )+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+theme(axis.text.x.bottom = element_text(size = 8,angle = 45,hjust = 1,vjust = 1))
usa_levels <- tdata %>% group_by(Cancer_Type) %>% summarise(Burden=median(Burden)) %>% arrange((Burden)) %>% pull(Cancer_Type)
#p <- plot_grid(p1,NULL,p2,align = 'h',nrow = 1,axis = 'tb',rel_widths = c(3,-0.36,1.5))
#ggsave(filename = 'tmp.pdf',plot = p,width = 8,height = 4.5,device = cairo_pdf)

# PGA
tdata <-  b2map %>% select(Cancer_Type=Country,Sample=Tumor_Barcode,Sample_Group=Histology,Burden=PGA_WGD)%>% filter(Sample_Group == 'Adenocarcinoma')
p3 <- TMBplot2(tdata,sgroupcol = sgroupcol,cancer_type_levels = global_levels )+ theme(axis.text.x.bottom = element_text(size = 8,angle = 45,hjust = 1,vjust = 1))+scale_y_continuous(limits = c(0,1),expand = c(0,0),breaks = seq(0,1,0.2))+ylab('Percentage of the Genome Aberrated\n')
tdata <-  b2map %>% filter(!is.na(State))%>% select(Cancer_Type=State,Sample=Tumor_Barcode,Sample_Group=Histology,Burden=PGA_WGD)%>% filter(Sample_Group == 'Adenocarcinoma')
p4 <- TMBplot2(tdata,sgroupcol = sgroupcol,cancer_type_levels = usa_levels )+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+theme(axis.text.x.bottom = element_text(size = 8,angle = 45,hjust = 1,vjust = 1))+scale_y_continuous(limits = c(0,1),expand = c(0,0),breaks = seq(0,1,0.2))
#p <- plot_grid(p1,NULL,p2,align = 'h',nrow = 1,axis = 'tb',rel_widths = c(3,-0.36,1.5))
#ggsave(filename = 'tmp.pdf',plot = p,width = 8,height = 4.5,device = cairo_pdf)

# SV
tdata <-  b2map %>% select(Cancer_Type=Country,Sample=Tumor_Barcode,Sample_Group=Histology,Burden=SV)%>% filter(Sample_Group == 'Adenocarcinoma')
p5 <- TMBplot2(tdata,sgroupcol = sgroupcol,cancer_type_levels = global_levels )+ theme(axis.text.x.bottom = element_text(size = 8,angle = 45,hjust = 1,vjust = 1))+scale_y_continuous(limits = c(0,10.5),expand = c(0,0),breaks = seq(0,10,2))+ylab('Number of SVs (log2)\n')
tdata <-  b2map %>% filter(!is.na(State))%>% select(Cancer_Type=State,Sample=Tumor_Barcode,Sample_Group=Histology,Burden=SV)%>% filter(Sample_Group == 'Adenocarcinoma')
p6 <- TMBplot2(tdata,sgroupcol = sgroupcol,cancer_type_levels = usa_levels )+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+theme(axis.text.x.bottom = element_text(size = 8,angle = 45,hjust = 1,vjust = 1))+scale_y_continuous(limits = c(0,10.5),expand = c(0,0),breaks = seq(0,10,2))

pall1 <- plot_grid(p1,NULL,p2,align = 'h',nrow = 1,axis = 'tb',rel_widths = c(3,-0.08,1.5))
pall2 <- plot_grid(p3,NULL,p4,align = 'h',nrow = 1,axis = 'tb',rel_widths = c(3,-0.08,1.5))
pall3 <- plot_grid(p5,NULL,p6,align = 'h',nrow = 1,axis = 'tb',rel_widths = c(3,-0.08,1.5))
p <- plot_grid(pall1,pall2, pall3, align = 'v',ncol = 1,axis = 'lr')

ggsave(filename = 'TMB_PGA_by_Country_LUAD.pdf',plot = p,width = 8,height = 11,device = cairo_pdf)

# polar plot
ggdata <-  b2map %>% select(Cancer_Type=Country,Sample=Tumor_Barcode,Sample_Group=Histology,TMB=Burden) %>% 
  group_by(Cancer_Type,Sample_Group) %>% summarise(TMB=median(TMB)) %>% ungroup()
tmp <- ggdata %>% filter(Sample_Group=='Adenocarcinoma') %>% arrange(TMB) %>% pull(Cancer_Type)
ggdata <- ggdata %>% mutate(Cancer_Type = factor(Cancer_Type,levels=tmp)) %>% arrange(Cancer_Type)
summary(ggdata$TMB)

p1 <- ggplot(data = ggdata, aes(x = factor(Cancer_Type), y = TMB, group = Sample_Group,color=Sample_Group)) + 
  #ylim(0, NA)+
  geom_point(stat = 'identity') + 
  geom_polygon(fill=NA,linewidth=0.5) + 
  coord_polar(start = - pi * 1/24)+
  theme_ipsum_rc(ticks = FALSE,grid = "Xx")+
  scale_color_manual(values = sgroupcol)+
  theme(axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid = element_line(linetype = 2))+
  geom_hline(yintercept = seq(-0.8,1.8,by=0.4), color = "gray90",size=0.15) +
  annotate('text', x = 1, y =seq(-0.8,1.8,by=0.4),label=seq(-0.8,1.8,by=0.4),size=3,col='#542788')+
  guides(color="none")


ggdata <-  b2map %>% select(Cancer_Type=Country,Sample=Tumor_Barcode,Sample_Group=Histology,TMB=PGA_WGD) %>% 
  group_by(Cancer_Type,Sample_Group) %>% summarise(TMB=median(TMB)) %>% ungroup()
#tmp <- ggdata %>% filter(Sample_Group=='Adenocarcinoma') %>% arrange(TMB) %>% pull(Cancer_Type)
ggdata <- ggdata %>% mutate(Cancer_Type = factor(Cancer_Type,levels=tmp)) %>% arrange(Cancer_Type)
summary(ggdata$TMB)

p2 <- ggplot(data = ggdata, aes(x = factor(Cancer_Type), y = TMB, group = Sample_Group,color=Sample_Group)) + 
  #ylim(0, NA)+
  geom_point(stat = 'identity') + 
  geom_polygon(fill=NA,linewidth=0.5) + 
  coord_polar(start = - pi * 1/24)+
  theme_ipsum_rc(ticks = FALSE,grid = "Xx")+
  scale_color_manual(values = sgroupcol)+
  theme(axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid = element_line(linetype = 2))+
  geom_hline(yintercept = seq(0,1,by=0.2), color = "gray90",size=0.15) +
  annotate('text', x = 1, y = seq(0,1,by=0.2),label=seq(0,1,by=0.2),size=3,col='#542788')+
  guides(color="none")


ggdata <-  b2map %>% select(Cancer_Type=Country,Sample=Tumor_Barcode,Sample_Group=Histology,TMB=SV) %>% 
  group_by(Cancer_Type,Sample_Group) %>% summarise(TMB=median(TMB)) %>% ungroup()
#tmp <- ggdata %>% filter(Sample_Group=='Adenocarcinoma') %>% arrange(TMB) %>% pull(Cancer_Type)
ggdata <- ggdata %>% mutate(Cancer_Type = factor(Cancer_Type,levels=tmp)) %>% arrange(Cancer_Type)
summary(ggdata$TMB)

p3 <- ggplot(data = ggdata, aes(x = factor(Cancer_Type), y = TMB, group = Sample_Group,color=Sample_Group)) + 
  ylim(0,10)+
  geom_point(stat = 'identity') + 
  geom_polygon(fill=NA,linewidth=0.5) + 
  coord_polar(start = - pi * 1/24)+
  theme_ipsum_rc(ticks = FALSE,grid = "Xx")+
  scale_color_manual(values = sgroupcol)+
  theme(axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid = element_line(linetype = 2),axis.line.y  = element_blank())+
  geom_hline(yintercept = seq(0,10,by=2), color = "gray90",size=0.15) +
  annotate('text', x = 1, y = seq(0,10,by=2),label=seq(0,10,by=2),size=3,col='#542788')+
  guides(color="none")

p <- plot_grid(p1,p2, p3, align = 'v',ncol = 1,axis = 'lr')
ggsave(filename = 'TMB_PGA_by_Country_polar.pdf',plot = p,width = 8,height = 12,device = cairo_pdf)


# Figure 1d ---------------------------------------------------------------
conflicts_prefer(dplyr::lag)
nmutation_input <- 0
freqdata_all <-  NULL
# SBS
freqdata_SBS <- tibble(Type=NULL,Signature=NULL,Value=NULL, Percentage=NULL)
#histmp <- 'Adenocarcinoma'
for(histmp in unique(covdata0$Histology)){
  samtmp <- covdata0 %>% filter(Histology == histmp) %>% pull(Tumor_Barcode)
  sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% sherlock_nonsmoker,Tumor_Barcode %in% samtmp) %>% select(Samples=Tumor_Barcode,starts_with('SBS')) 
  tmp <- prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,colset = sigcol)
  freqdata_SBS <- freqdata_SBS %>% bind_rows(tmp$freq_data %>% mutate(Histology=histmp))
}
freqdata_all <- freqdata_SBS %>% filter(Value>0) %>% mutate(Profile='SBS') %>% bind_rows(freqdata_all)

# DBS
freqdata_DBS <- tibble(Type=NULL,Signature=NULL,Value=NULL, Percentage=NULL)
#histmp <- 'Adenocarcinoma'
for(histmp in unique(covdata0$Histology)){
  samtmp <- covdata0 %>% filter(Histology == histmp) %>% pull(Tumor_Barcode)
  sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% sherlock_nonsmoker,Tumor_Barcode %in% samtmp) %>% select(Samples=Tumor_Barcode,starts_with('DBS')) 
  tmp <- prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,colset = sigcol)
  freqdata_DBS <- freqdata_DBS %>% bind_rows(tmp$freq_data %>% mutate(Histology=histmp))
}
freqdata_all <- freqdata_DBS %>% filter(Value>0) %>% mutate(Profile='DBS') %>% bind_rows(freqdata_all)

# ID
freqdata_ID <- tibble(Type=NULL,Signature=NULL,Value=NULL, Percentage=NULL)
#histmp <- 'Adenocarcinoma'
for(histmp in unique(covdata0$Histology)){
  samtmp <- covdata0 %>% filter(Histology == histmp) %>% pull(Tumor_Barcode)
  sigdata <- ludmil_activity_all %>% filter(Tumor_Barcode %in% sherlock_nonsmoker,Tumor_Barcode %in% samtmp) %>% select(Samples=Tumor_Barcode,starts_with('ID')) 
  tmp <- prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,colset = sigcol)
  freqdata_ID <- freqdata_ID %>% bind_rows(tmp$freq_data %>% mutate(Histology=histmp))
}
freqdata_all <- freqdata_ID %>% filter(Value>0) %>% mutate(Profile='ID') %>% bind_rows(freqdata_all)

# CN
freqdata_CN <- tibble(Type=NULL,Signature=NULL,Value=NULL, Percentage=NULL)
#histmp <- 'Adenocarcinoma'
for(histmp in unique(covdata0$Histology)){
  samtmp <- covdata0 %>% filter(Histology == histmp) %>% pull(Tumor_Barcode)
  sigdata <- cn68_activity %>% filter(Tumor_Barcode %in% sherlock_nonsmoker,Tumor_Barcode %in% samtmp) %>% select(Samples=Tumor_Barcode,starts_with('CN')) 
  tmp <- prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,colset = sigcol)
  freqdata_CN <- freqdata_CN %>% bind_rows(tmp$freq_data %>% mutate(Histology=histmp))
}
freqdata_all <- freqdata_CN %>% filter(Value>0) %>% mutate(Profile='CN') %>% bind_rows(freqdata_all)

# SV
freqdata_SV <- tibble(Type=NULL,Signature=NULL,Value=NULL, Percentage=NULL)
#histmp <- 'Adenocarcinoma'
for(histmp in unique(covdata0$Histology)){
  samtmp <- covdata0 %>% filter(Histology == histmp) %>% pull(Tumor_Barcode)
  sigdata <- sv38_activity %>% filter(Tumor_Barcode %in% sherlock_nonsmoker,Tumor_Barcode %in% samtmp) %>% select(Samples=Tumor_Barcode,starts_with('SV')) 
  tmp <- prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,colset = sigcol)
  freqdata_SV <- freqdata_SV %>% bind_rows(tmp$freq_data %>% mutate(Histology=histmp))
}
freqdata_all <- freqdata_SV %>% filter(Value>0) %>% mutate(Profile='SV') %>% bind_rows(freqdata_all)

# visualization
freqdata_all <- freqdata_all %>% 
  mutate(Profile=factor(Profile,levels=c('SBS','DBS','ID','CN','SV'))) %>% 
  mutate(Histology = factor(Histology, levels=c('Adenocarcinoma','Squamous cell carcinoma','Carcinoid tumor','Other'),labels=c('LUAD','LUSC','Carcinoid tumor','Others')))

siglevs <- freqdata_all %>% select(Profile,Signature) %>% mutate(Channel=extract_numeric(Signature)) %>% arrange(Profile,Channel) %>% unique() %>% pull(Signature)

freqdata_all %>% 
  mutate(Signature = factor(Signature,levels=siglevs)) %>% 
  mutate(Type = str_replace(Type,"Prevalence by","By")) %>% 
  mutate(lab=if_else(Value>0.05,as.character(percent(Value,accuracy=1)), NA_character_)) %>% 
  mutate(Value=Value*100) %>% 
  mutate(Type=str_remove(Type,'s$')) %>% 
  ggplot(aes(Signature,Histology,fill=Value))+
  geom_tile(col='white',linewidth=0.2)+
  geom_text(aes(label=lab),col='gray50',size=3)+
  facet_grid(Type~Profile,scales = 'free',space = 'free')+
  scale_fill_viridis_c(limit=c(0,100))+
  labs(x=NULL,y=NULL,fill='Prevalence')+
  theme_ipsum_rc(grid = FALSE,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),panel.spacing = unit(0.2,'cm'),strip.text.x = element_text(hjust = 0.5,face = 'bold'),strip.text.y = element_text(hjust = 0.5,face = 'bold'),legend.key.height = unit(1,'cm'))+
  panel_border(color = 'black',size = 0.4)


freqdata_all %>% 
  mutate(Signature = factor(Signature,levels=siglevs)) %>% 
  mutate(Type = str_replace(Type,"Prevalence by","By")) %>% 
  mutate(lab=if_else(Value>0.05,as.character(percent(Value,accuracy=1)), NA_character_)) %>% 
  mutate(Value=Value*100) %>% 
  mutate(Type=str_remove(Type,'s$')) %>% 
  ggplot(aes(Signature,Histology,fill=Value))+
  #geom_tile(col='white',linewidth=0.2)+
  geom_point(aes(size=Value),pch=21,stroke=0.5)+
  #geom_text(aes(label=lab),col='gray50',size=3)+
  facet_grid(Type~Profile,scales = 'free',space = 'free')+
  scale_fill_viridis_c(limit=c(0,100),breaks=pretty_breaks())+
  scale_size(limit=c(0,100),breaks = pretty_breaks())+
  labs(x=NULL,y=NULL,fill='Prevalence')+
  theme_ipsum_rc(grid = FALSE,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),panel.spacing = unit(0.2,'cm'),strip.text.x = element_text(hjust = 0.5,face = 'bold'),strip.text.y = element_text(hjust = 0.5,face = 'bold'),legend.key.height = unit(0.5,'cm'),legend.position = 'top',legend.key.width = unit(2,'cm'))+
  guides(size = guide_legend(nrow = 1))+
  panel_border(color = 'black',size = 0.4)

ggsave(filename = 'Prevalence_overall2.pdf',width = 16,height = 3.6,device = cairo_pdf)



# Figure 2 ----------------------------------------------------------------

# Figure 2a ---------------------------------------------------------------

ggdata <- sherlock_variable %>% 
  filter(name=='Mutations') %>% 
  pivot_wider() %>% 
  left_join(sherlock_PGA) %>% 
  filter(Tumor_Barcode %in% luad_nonsmoker) %>% 
  left_join(
    ludmil_activity_all %>% dplyr::select(Tumor_Barcode,starts_with('SBS')) %>% pivot_longer(-Tumor_Barcode) %>% group_by(Tumor_Barcode) %>% arrange(desc(value)) %>% slice(1)
  ) 

tmp <- names(sigcol[names(sigcol) %in% unique(ggdata$name)])

ggdata %>% 
  mutate(Mutations = log10(Mutations)) %>% 
  ggplot(aes(PGA_WGD,(Mutations)))+
  geom_point(aes(fill=name),pch=21,size=3,col='white',stroke=0.2)+
  geom_smooth(method = 'lm')+
  scale_fill_manual(values = sigcol,breaks = tmp)+
  scale_x_continuous(breaks = pretty_breaks(),labels = percent_format())+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = FALSE,ticks = T)+
  labs(x = 'Percentage of Genome Aberrated', y = 'Number of mutations (log10)',fill='Most prevalent\nsignature')+
  theme(legend.position = 'top',legend.margin=margin(0,0,-12,0),legend.box.margin=margin(0,0,0,0))+
  panel_border(color = 'black',linetype = 1,size=0.5)

ggsave(filename = 'TMB_PGA.pdf',width = 6,height = 5,device = cairo_pdf)

ggdata %>% do(tidy(cor.test(log10(.$Mutations), .$PGA)))





# Figure 2b ---------------------------------------------------------------
genomesize <-  genome2size('GRCh38')

data_input <- ludmil_activity_all %>% filter(Tumor_Barcode %in% luad_nonsmoker) %>% pivot_longer(cols = -Tumor_Barcode) %>% filter(str_starts(name,'SBS')) %>% select(Cancer_Type=name,Sample=Tumor_Barcode,Burden=value) %>% mutate(Burden=log10((Burden)/genomesize)) %>% mutate(Sample_Group=Cancer_Type)

#data_input <- data_input %>% filter(Cancer_Type %in% tmpx)
TMBplot2(data_input,sgroupcol = sigcol,output_plot = 'LUAD_SBS_TMplot.pdf',plot_width = 6,plot_height = 3.6)

# Figure 2c-d ---------------------------------------------------------------
nmutation_input <-  0
conflicts_prefer(dplyr::lag)

# for(tmpg in c('N_A')){
#   tmp <- wgs_groups_info %>% left_join(covdata0 %>% select(Tumor_Barcode, Histology)) %>% filter(Assigned_Population == 'EAS',Histology=='Adenocarcinoma',Smoking=='Non-Smoker') %>% pull(Tumor_Barcode)
#   sigdata <- ludmil_activity_all %>% select(Samples=Tumor_Barcode,starts_with('SBS')) %>% filter(Samples %in% tmp)
#   data1 <- prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,output_plot = paste0('SBS_N_LUAD_EAS_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))
# }
# 
# for(tmpg in c('N_U')){
#   tmp <- wgs_groups_info %>% left_join(covdata0 %>% select(Tumor_Barcode, Histology)) %>% filter(Assigned_Population == 'EUR',Histology=='Adenocarcinoma',Smoking=='Non-Smoker') %>% pull(Tumor_Barcode)
#   sigdata <- ludmil_activity_all %>% select(Samples=Tumor_Barcode,starts_with('SBS')) %>% filter(Samples %in% tmp)
#   data2 <- prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,output_plot = paste0('SBS_N_LUAD_EUR_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))
# }

tmp <- country_group %>% filter(Tumor_Barcode %in% luad_nonsmoker,Country_Group == 'East') %>% pull(Tumor_Barcode)
sigdata <- ludmil_activity_all %>% select(Samples=Tumor_Barcode,starts_with('SBS')) %>% filter(Samples %in% tmp)
data1 <- prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,output_plot = paste0('SBS_N_LUAD_East_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))

tmp <- country_group %>% filter(Tumor_Barcode %in% luad_nonsmoker,Country_Group == 'NA/EU') %>% pull(Tumor_Barcode)

sigdata <- ludmil_activity_all %>% select(Samples=Tumor_Barcode,starts_with('SBS')) %>% filter(Samples %in% tmp)
data2 <- prevalence_plot(sigdata = sigdata,nmutation = nmutation_input,output_plot = paste0('SBS_N_LUAD_West_prevalence.pdf'),colset = sigcol,rel_widths = c(1,4))

tdata <- bind_rows(
  data1$freq_data %>% mutate(Group='AS'),
  data2$freq_data %>% mutate(Group='NA/EU')
)

myggstyle()

plotdata <- tdata %>% 
  filter(Type == 'Prevalence by samples') %>% 
  arrange(desc(Value)) %>% 
  mutate(Signature = fct_inorder(Signature)) %>% 
  mutate(Value2=if_else(Value<0.1,0.1,Value)) %>% 
  mutate(Value = if_else(Group=='East',-1*Value,Value)) %>% 
  mutate(Value2 = if_else(Group=='East',-1*Value2,Value2)) 

p1 <- plotdata %>% 
  ggplot(aes(Value, y=Signature,fill=Signature))+
  geom_col(width = 0.75,col='gray10',linewidth=0.3)+
  geom_vline(xintercept = 0,col='gray10')+
  geom_text(data = plotdata %>% filter(Group=='AS'),aes(label=Percentage),vjust=0.5,hjust=1.2)+
  geom_text(data = plotdata %>% filter(Group=='NA/EU'),aes(label=Percentage),vjust=0.5,hjust=-0.2)+
  scale_fill_manual(values = sigcol)+
  scale_x_continuous(limits = c(-1.2,1.2),breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),labels = c("100%","75%","50%","25%","0%","25%","50%","75%","100%"))+
  labs(x='Prevalence by samples (%)',y=NULL,fill='Group')+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 12,grid = 'Y',ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  panel_border(color = 'gray10',size = 1)+
  theme(legend.position = 'none',panel.grid.major.y = element_line(linetype = 2,colour = 'gray90'))


testdata <- ludmil_activity_all_obs %>%
  select(Tumor_Barcode,starts_with('SBS')) %>% 
  pivot_longer(-Tumor_Barcode) %>% 
  left_join(country_group) %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% luad_nonsmoker,Country_Group %in% c('AS','NA/EU')) %>% 
  mutate(Country_Group = factor(Country_Group, levels=c('AS','NA/EU'))) %>% 
  group_by(name) %>% 
  #do(tidy(fisher.test(.$value,.$Assigned_Population)))
  mutate(value=as.factor(value)) %>% 
  #do(tresult = safely(stats::fisher.test)(.$value,.$Assigned_Population)) %>% 
  do(tresult = safely(stats::glm)(value ~ Country_Group + Gender + Age + Tumor_Purity,family = binomial(),data=.)) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(name,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='Country_GroupNA/EU') %>% 
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

plot_grid(p2,p1,align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2,3))

ggsave(filename = 'LUAD_SBS_Region_prevalence_tmp.pdf',width = 6,height = 6,device = cairo_pdf)

# piechart
library(webr)

plotdata <- tdata %>% filter(Type == 'Prevalence by mutations') %>% mutate(Value=round(Value,3))

PieDonut_ztw(plotdata %>% filter(Group=='AS') %>% arrange(Group,desc(Signature)) %>% mutate(Signature = fct_inorder(Signature)),aes(pies=Signature,count=Value),mainCol = sigcol,start=3*pi/2,showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.02,showDonutName = T,pieLabelSize=3.5,titlesize = 5,r0=0.3,title='AS_N', showPieName = FALSE,family = 'Roboto Condensed', labelpositionThreshold =0.02)
ggsave(filename = 'SBS_N_LUAD_East_prevalence2.pdf',width = 5,height = 5,device = cairo_pdf)

PieDonut_ztw(plotdata %>% filter(Group=='NA/EU') %>% arrange(Group,desc(Signature)) %>% mutate(Signature = fct_inorder(Signature)),aes(pies=Signature,count=Value),mainCol = sigcol,start=3*pi/2,showNum = FALSE,showRatioPie = TRUE,showRatioThreshold = 0.02,showDonutName = T,pieLabelSize=3.5,titlesize = 5,r0=0.3,title='AS_N', showPieName = FALSE,family = 'Roboto Condensed', labelpositionThreshold =0.02)
ggsave(filename = 'SBS_N_LUAD_West_prevalence2.pdf',width = 5,height = 5,device = cairo_pdf)







# Figure 2e-f ---------------------------------------------------------------
testdata <- sherlock_data_full %>% 
  filter(Gene %in% drglist, Type=='Mutation_Driver') %>% 
  left_join(covdata0) %>% 
  left_join(country_group) %>% 
  filter(Tumor_Barcode %in% luad_nonsmoker, Country_Group %in% c('AS','NA/EU')) %>% 
  mutate(Country_Group = factor(Country_Group,levels=c('AS','NA/EU')))

plotdata <- testdata %>% 
  group_by(Gene) %>% 
  mutate(Alteration=as.factor(Alteration)) %>% 
  do(tresult = safely(glm)(Alteration ~ Country_Group + Gender + Age + Tumor_Purity, family='binomial',data=. )) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(Gene,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='Country_GroupNA/EU') %>% 
  filter(abs(log(estimate))<10) %>% 
  mutate(FDR=p.adjust(p.value,method='BH'))

plotdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR)))+
  geom_hline(yintercept = -log10(0.01),linetype=2,col=ncicolpal[1],size=0.5)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col=ncicolpal[2],size=0.5)+
  geom_vline(xintercept = 0,linetype=1,col='gray50',size=0.5)+
  geom_point(aes(fill=estimate>1),pch=21,stroke=0.1,col='black',size=4)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_y_continuous(breaks = pretty_breaks())+
  scale_fill_nejm()+
  ggrepel::geom_text_repel(data=plotdata %>% filter(FDR<0.05),aes(label=Gene),max.overlaps = 30,size=4.5)+
  labs(x='Odd ratio (log2)',y='-log10(FDR)')+
  guides(fill = "none")+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid=FALSE,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  panel_border(color = 'black',size = 0.5)

ggsave(filename = 'LUAD_N_driver_region.pdf',width = 5,height = 4,device = cairo_pdf)

plotdata <- testdata %>% 
  filter(Gene %in% c('EGFR','KRAS','TP53')) %>% 
  count(Gene,Country_Group,Alteration) %>% 
  group_by(Gene,Country_Group) %>% 
  mutate(freq=n/(sum(n))) %>% 
  filter(Alteration == 'Yes') %>% 
  ungroup()

plotdata %>% 
  mutate(Gene= fct_rev(fct_reorder(Gene,freq))) %>% 
  ggplot(aes(Gene,freq,group=Country_Group,fill=Country_Group))+
  geom_bar (stat="identity", position = position_dodge(width = 1),width = 0.8)+
  scale_y_continuous(breaks = pretty_breaks(n=7),labels =percent_format())+
  scale_fill_nejm()+
  labs(x=NULL,y="Mutation frequency",fill='Group')+
  theme_ipsum_rc(grid = "Y",ticks = TRUE,axis_title_just = 'm',axis_title_size = 14,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(legend.position = c(0.7,0.9),axis.text.x = element_text(face = 'italic'))+
  guides(fill=guide_legend(title.position = 'left',nrow = 1))+
  panel_border(color = 'black',size = 0.5)
ggsave(filename = 'LUAD_N_driver_region2.pdf',width = 3,height = 4,device = cairo_pdf)


# Figure 2g ---------------------------------------------------------------
load('../Mutation_Signature_Probability_SBS.RData',verbose = T)
load('../../MAFtools/sherlock_driver_mutations.RData',verbose = T)

tmp <- sherlock_driver_mutations %>% filter(Hugo_Symbol %in% c('EGFR','KRAS','TP53')) %>% mutate(ID=str_remove(ID,"^chr")) %>% pull(ID)
plotdata <- Mutation_Signature_Probability_SBS %>% filter(Gene_Name %in% c('EGFR','KRAS','TP53'),ID %in% tmp)

plotdata <- plotdata %>% select(Tumor_Barcode,ID,Gene_Name,starts_with('SBS')) %>% unique() %>% pivot_longer(-c(Tumor_Barcode,ID,Gene_Name))

plotdata <- plotdata %>% filter(Tumor_Barcode %in% luad_nonsmoker)
tmp <- plotdata %>% select(Tumor_Barcode,ID,Gene_Name) %>% unique() %>% count(Gene_Name) %>% mutate(Label=paste0(Gene_Name,"\n(N=",n,")")) %>% select(-n)

plotdata <- plotdata  %>% group_by(Gene_Name,name) %>% summarise(value=sum(value)) %>% ungroup() %>% left_join(tmp) %>% filter(value>0)

plotdata %>% 
  mutate(name=factor(name,levels=names(sigcol))) %>% 
  mutate(Label=factor(Label,levels=tmp$Label[c(1,3,2)])) %>% 
  ggplot(aes(Label,value,fill=name))+
  geom_col(position = "fill",width = 0.75)+
  scale_fill_manual(values = sigcol,breaks = names(sigcol))+
  scale_y_continuous(breaks = pretty_breaks(n=7),labels = percent_format(),expand = c(0,0))+
  labs(x=NULL,y='Proportion of driver mutations (SBS)',title = 'SBS signatures assigned to driver mutations',fill='Signature')+
  theme_ipsum_rc(plot_title_size = 16,base_size = 14,grid = "Y",ticks = TRUE,axis_title_just = 'm',axis_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(axis.text.x = element_text(face = 'italic'),plot.title = element_text(hjust = 0.5))

ggsave(filename = 'LUAD_N_driver_signature_assignment.pdf',width = 5,height = 7,device = cairo_pdf)

# Figure 2h ---------------------------------------------------------------
load('../Mutation_Signature_Probability_ID.RData',verbose = T)
load('../../MAFtools/sherlock_driver_mutations.RData',verbose = T)

tmp <- sherlock_driver_mutations %>% filter(Hugo_Symbol %in% c('EGFR','KRAS','TP53')) %>% mutate(ID=str_remove(ID,"^chr")) %>% pull(ID)
plotdata <- Mutation_Signature_Probability_ID %>% filter(Gene_Name %in% c('EGFR','KRAS','TP53'),ID %in% tmp)

plotdata <- plotdata %>% select(Tumor_Barcode,ID,Gene_Name,starts_with('ID')) %>% unique() %>% pivot_longer(-c(Tumor_Barcode,ID,Gene_Name))

plotdata <- plotdata %>% filter(Tumor_Barcode %in% luad_nonsmoker)
tmp <- plotdata %>% select(Tumor_Barcode,ID,Gene_Name) %>% unique() %>% count(Gene_Name) %>% mutate(Label=paste0(Gene_Name,"\n(N=",n,")")) %>% select(-n)

plotdata <- plotdata  %>% group_by(Gene_Name,name) %>% summarise(value=sum(value)) %>% ungroup() %>% left_join(tmp) %>% filter(value>0)

plotdata %>% 
  mutate(name=factor(name,levels=names(sigcol))) %>% 
  mutate(Label=factor(Label,levels=tmp$Label[c(1,3,2)])) %>% 
  ggplot(aes(Label,value,fill=name))+
  geom_col(position = "fill",width = 0.75)+
  scale_fill_manual(values = sigcol,breaks = names(sigcol))+
  scale_y_continuous(breaks = pretty_breaks(n=7),labels = percent_format(),expand = c(0,0))+
  labs(x=NULL,y='Proportion of driver mutations (ID)',title = 'ID signatures assigned to driver mutations',fill='Signature')+
  theme_ipsum_rc(plot_title_size = 16,base_size = 14,grid = "Y",ticks = TRUE,axis_title_just = 'm',axis_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(axis.text.x = element_text(face = 'italic'),plot.title = element_text(hjust = 0.5))

ggsave(filename = 'LUAD_N_driver_signature_assignment_ID.pdf',width = 3.8,height = 7,device = cairo_pdf)



# Figure 3a-d ---------------------------------------------------------------
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
  left_join(verifybamID)

my_comparisons <- list(c("No", "Yes"))

# 3a
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

ggsave(filename = 'LCINS_passive_smoking_SBS.pdf',width = 3 ,height =6,device = cairo_pdf)  

# 3b
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

ggsave(filename = 'LCINS_passive_smoking_telomere.pdf',width = 3 ,height =6,device = cairo_pdf)  


#3c
#3c
tdata %>% 
  lm(Total_SBS ~ Passive_Smoking + Age + Gender + PC1 + PC2 + Histology + Tumor_Purity, data=.) %>% tbl_regression(pvalue_fun = sci_notation) %>% as_gt() %>% gt::gtsave(filename = 'tmp.docx')

plotdata <- tdata %>% 
  do(tidy(lm(Total_SBS ~ Passive_Smoking + Age + Gender + Assigned_Population + Histology + Tumor_Purity, data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term != '(Intercept)') %>% 
  #mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(label=paste0('β = ',round(estimate,2), ', q-value = ',scientific_format(digits = 3)(p.value))) %>% 
  arrange(term)

plotdata$term
#plotdata$term <-  c('Age','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Ancestry: PC1', 'Ancestry: PC2','Passive Smoking: Yes\nRef=No','Tumor Purity')

plotdata$term <-  c('Age','Ancestry: EAS\nRef=EUR', 'Ancestry: Others\nRef=EUR','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Passive Smoking: Yes\nRef=No','Tumor Purity')


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

ggsave(filename = 'LCINS_passive_smoking_SBS_multivariable.pdf',plot = p,width = 7,height = 5,device = cairo_pdf)

# 3d 
tdata %>% 
  lm(Telseq_TL_Ratio ~ Passive_Smoking + Age + Gender + PC1 + PC2 + Histology + Tumor_Purity, data=.) %>% tbl_regression(pvalue_fun = sci_notation) %>% as_gt() %>% gt::gtsave(filename = 'tmp.docx')


plotdata <- tdata %>% 
  do(tidy(lm(Telseq_TL_Ratio ~ Passive_Smoking + Age + Gender + Assigned_Population + Histology + Tumor_Purity, data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term != '(Intercept)') %>% 
  #mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(label=paste0('β = ',round(estimate,2), ', q-value = ',scientific_format(digits = 3)(p.value))) %>% 
  arrange(term)

plotdata$term
#plotdata$term <-  c('Age','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Ancestry: PC1', 'Ancestry: PC2','Passive Smoking: Yes\nRef=No','Tumor Purity')
plotdata$term <-  c('Age','Ancestry: EAS\nRef=EUR', 'Ancestry: Others\nRef=EUR','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Passive Smoking: Yes\nRef=No','Tumor Purity')

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

ggsave(filename = 'LCINS_passive_smoking_telomere_multivariable.pdf',plot = p,width = 7,height = 5,device = cairo_pdf)


# Figure 3e ---------------------------------------------------------------
testdata <- ludmil_activity_all_obs %>%
  select(Tumor_Barcode,starts_with('SBS')) %>% 
  pivot_longer(-Tumor_Barcode) %>% 
  left_join(clinical_data %>% select(Tumor_Barcode,Passive_Smoking)) %>% 
  left_join(covdata0) %>% 
  left_join(verifybamID) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,!is.na(Passive_Smoking)) %>% 
  group_by(name) %>% 
  # do(tresult = safely(stats::fisher.test)(.$value,.$Passive_Smoking)) %>% 
  do(tresult = safely(stats::glm)(value ~ Passive_Smoking + PC1 + PC2 + Histology + Gender + Age + Tumor_Purity,family = binomial(),data=.)) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(name,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='Passive_SmokingYes') %>% 
  filter(abs(log(estimate))<10) %>% 
  mutate(FDR=p.adjust(p.value,method='BH'))

testdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR),fill=name))+
  geom_hline(yintercept = -log10(0.05),col=ncicolpal[2],linetype=2)+
  geom_hline(yintercept = -log10(0.01),col=ncicolpal[1],linetype=2)+
  geom_vline(xintercept = 0,col="gray10",linetype=1,linewidth=0.5)+
  geom_point(pch=21,size=3.5,col='black',stroke=0.2)+
  ggrepel::geom_text_repel(data=testdata %>% filter(FDR<1),aes(label=name,col=name),force=20)+
  scale_fill_manual(values = sigcol)+
  scale_color_manual(values = sigcol)+
  scale_x_continuous(breaks = pretty_breaks(n = 7),limits = c(-2,2))+
  scale_y_continuous(breaks = pretty_breaks(n=7))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 12,grid = F,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x = 'Odd ratio (log2)', y = '-log10(FDR)',title='SBS mutational signatures')+
  guides(fill="none",color='none')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  coord_cartesian(clip = 'off')

ggsave(filename = 'LCINS_passive_smoking_signature_enrichment.pdf',width = 6,height = 5,device = cairo_pdf)


# Figure 3f ---------------------------------------------------------------
plotdata <- sherlock_sbs96_profile %>% 
  mutate(Type=str_sub(MutationType,3,5)) %>% 
  group_by(Tumor_Barcode,Type) %>% 
  summarise(Contribution = log10(sum(Contribution))) %>% 
  ungroup() %>% 
  left_join(clinical_data %>% select(Tumor_Barcode,Passive_Smoking)) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,!is.na(Passive_Smoking)) %>% 
  mutate(Passive_Smoking = fct_drop(Passive_Smoking))

my_comparisons <- list(c("No", "Yes"))
stat.test <- plotdata %>% 
  group_by(Type) %>% 
  wilcox_test(Contribution ~ Passive_Smoking, comparisons = my_comparisons) %>%
  mutate(FDR=p.adjust(p,method='BH')) %>% 
  add_xy_position(x = "Type") %>%
  mutate(myformatted.p = sprintf("FDR = %.3f",FDR)) %>% 
  mutate(Passive_Smoking = "Yes")

stat.test

stat.test <- plotdata %>%
  left_join(covdata0) %>% 
  left_join(verifybamID) %>% 
  group_by(Type) %>% 
  do(tidy(lm(Contribution ~ Passive_Smoking + Age + Gender + PC1 + PC2 + Histology + Tumor_Purity, data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term=='Passive_SmokingYes') %>% 
  mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  #mutate(label=paste0('β = ',round(estimate,2), ', FDR = ',scientific_format(digits = 3)(FDR))) %>%
  arrange(term) %>% 
  ungroup() %>% 
  mutate(group1='No',group2='Yes',x=as.integer(as.factor(Type)),xmin=x-0.2,xmax=x+0.2,y.position=4.5) %>% 
  mutate(myformatted.p = sprintf("FDR=%.2f",FDR)) %>% 
  mutate(Passive_Smoking = "Yes")

stat.test

plotdata %>% 
  ggplot(aes(Type,Contribution,fill=Passive_Smoking))+
  geom_quasirandom(pch=21,size=1.6,width = 0.12,color="white",stroke=0.25,dodge.width=0.7)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.6,position = position_dodge(width = 0.7))+
  scale_fill_jama()+
  scale_y_continuous(breaks = pretty_breaks(n = 7),limits = c(1.2,4.6))+
  theme_ipsum_rc(base_size = 13,axis_title_just = 'm',axis_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 15,ticks = T,grid='Y')+
  labs(x = NULL, y = 'Total number of SBS mutations (log10)',title = 'SBS mutation types',fill='Passive smoking')+
  theme(plot.title = element_text(hjust = 0.5))+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_passive_smoking_mutationtypes.pdf',width = 8 ,height =5,device = cairo_pdf)  

# Figure 3g ---------------------------------------------------------------
testdata <- sherlock_data_full %>% 
  filter(Gene %in% drglist, Type=='Mutation_Driver') %>%
  left_join(clinical_data %>% select(Tumor_Barcode,Passive_Smoking)) %>% 
  left_join(covdata0) %>% 
  left_join(verifybamID) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker, !is.na(Passive_Smoking)) 

plotdata <- testdata %>% 
  group_by(Gene) %>% 
  mutate(Alteration=as.factor(Alteration)) %>% 
  do(tresult = safely(glm)(Alteration ~ Passive_Smoking + Histology + PC1 + PC2 + Gender + Age + Tumor_Purity, family='binomial',data=. )) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(Gene,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='Passive_SmokingYes') %>% 
  filter(abs(log(estimate))<10) %>%
  mutate(FDR=p.adjust(p.value,method='BH'))


plotdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR)))+
  geom_hline(yintercept = -log10(0.01),linetype=2,col=ncicolpal[1],size=0.5)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col=ncicolpal[2],size=0.5)+
  geom_vline(xintercept = 0,linetype=1,col='gray50',size=0.5)+
  geom_point(aes(fill=log2(estimate)>0),pch=21,stroke=0.1,col='black',size=4)+
  scale_x_continuous(breaks = pretty_breaks(n = 7),limits = c(-2.7,2.7))+
  scale_y_continuous(breaks = pretty_breaks())+
  scale_fill_jama()+
  ggrepel::geom_text_repel(data=plotdata %>% filter(FDR<1),aes(label=Gene),force = 10,max.overlaps = 30,size=4.5)+
  labs(x='Odd ratio (log2)',y='-log10(FDR)',title = 'Driver genes')+
  guides(fill = "none")+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid=FALSE,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 15)+
  theme(plot.title = element_text(hjust = 0.5))+
  panel_border(color = 'black',size = 0.5)

ggsave(filename = 'LCINS_passive_smoking_driver_mutations.pdf',width = 6,height = 5,device = cairo_pdf)


# Figure 4 ----------------------------------------------------------------
tdata <- bind_rows(
  sherlock_sbs96_profile %>% mutate(Profile='SBS'),
  sherlock_dbs78_profile %>% mutate(Profile='DBS'),
  sherlock_id83_profile %>% mutate(Profile='ID')
) %>% 
  group_by(Tumor_Barcode,Profile) %>% 
  summarise(Burden=log10(sum(Contribution))) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Profile,values_from = Burden) %>% 
  left_join(
    sherlock_variable %>% filter(name=='Telseq_TL_Ratio') %>% pivot_wider()
  ) %>% 
  left_join(pollution_data) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker) %>% 
  filter(!is.na(Pollution_group2))

my_comparisons <- list(c("Low", "High"))


# Figure 4a ---------------------------------------------------------------
plotdata <- tdata %>% select(Tumor_Barcode,SBS,DBS,ID,Pollution_group2) %>% pivot_longer(cols = -c(Tumor_Barcode,Pollution_group2)) %>% filter(is.finite(value)) %>% 
  mutate(name=factor(name,levels=c('SBS','DBS','ID')))


stat.test <- plotdata %>% 
  group_by(name) %>% 
  wilcox_test(value ~ Pollution_group2, comparisons = my_comparisons) %>%
  ungroup() %>% 
  add_xy_position(x = "Pollution_group2") %>%
  mutate(myformatted.p = sprintf("P = %.3f",p)) %>% 
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
  labs(x = expression("Pollution group (Population weighted PM"[2.5]*")"), y = 'Total number of mutations (log10)')+
  theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=14,face = 'bold',hjust = 0.5),panel.spacing = unit(0.4,'cm'))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_Pollution_mutations.pdf',width = 6 ,height =6,device = cairo_pdf)  


# Figure 4b ---------------------------------------------------------------
plotdata <- tdata %>% select(Tumor_Barcode,Telseq_TL_Ratio,Pollution_group2) %>%  filter(!is.na(Telseq_TL_Ratio))


stat.test <- plotdata %>% 
  wilcox_test(Telseq_TL_Ratio ~ Pollution_group2, comparisons = my_comparisons) %>%
  ungroup() %>% 
  add_xy_position(x = "Pollution_group2") %>%
  mutate(myformatted.p = sprintf("P = %.2e",p)) %>% 
  mutate(Pollution_group2 = "High")

stat.test

plotdata %>% 
  ggplot(aes(Pollution_group2,Telseq_TL_Ratio,fill=Pollution_group2))+
  geom_quasirandom(pch=21,size=2,width = 0.3,color="white",stroke=0.4)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.7)+
  scale_fill_manual(values = c('#01665e','#BB0E3D'))+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 13,axis_title_just = 'm',axis_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 15,ticks = T)+
  labs(x = expression("Pollution group (Population weighted PM"[2.5]*")"), y = 'Telomere length ratio, log2(Tumor TL/Normal TL)',title='Telomere length')+
  theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=14,face = 'bold',hjust = 0.5),panel.spacing = unit(0.4,'cm'))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_Pollution_telomere.pdf',width = 2.5 ,height =6,device = cairo_pdf)  


# Figure 4c ---------------------------------------------------------------

plotdata <- tdata %>% select(Tumor_Barcode,SBS,DBS,ID,Pollution_group2) %>% pivot_longer(cols = -c(Tumor_Barcode,Pollution_group2)) %>% filter(is.finite(value)) %>% 
  mutate(name=factor(name,levels=c('SBS','DBS','ID'))) %>% 
  left_join(covdata0) %>% 
  left_join(verifybamID)

plotdata %>% filter(name == 'SBS') %>% lm(value ~ Pollution_group2 + Age + Gender + PC1 + PC2 + Histology + Tumor_Purity, data=.) %>% tbl_regression(pvalue_fun = sci_notation) %>% as_gt() %>% gt::gtsave(filename = 'tmp.docx')
plotdata %>% filter(name == 'DBS') %>% lm(value ~ Pollution_group2 + Age + Gender + PC1 + PC2 + Histology + Tumor_Purity, data=.) %>% tbl_regression(pvalue_fun = sci_notation) %>% as_gt() %>% gt::gtsave(filename = 'tmp.docx')
plotdata %>% filter(name == 'ID') %>% lm(value ~ Pollution_group2 + Age + Gender + PC1 + PC2 + Histology + Tumor_Purity, data=.) %>% tbl_regression(pvalue_fun = sci_notation) %>% as_gt() %>% gt::gtsave(filename = 'tmp.docx')


plotdata <- plotdata %>% 
  group_by(name) %>% 
  do(tidy(lm(value ~ Pollution_group2 + Age + Gender + Assigned_Population + Histology + Tumor_Purity, data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term != '(Intercept)') %>% 
  #mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(label=paste0('β = ',round(estimate,2), ', q-value = ',scientific_format(digits = 3)(p.value))) %>% 
  arrange(term) %>% 
  ungroup()

plotdata$term
#plotdata$term <-  rep(c('Age','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Ancestry: PC1', 'Ancestry: PC2','Pollution Group: High\nRef=Low','Tumor Purity'),each=3)
plotdata$term <-  rep(c('Age','Ancestry: EAS\nRef=EUR', 'Ancestry: Others\nRef=EUR','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Pollution Group: High\nRef=Low','Tumor Purity'),each=3)

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

ggsave(filename = 'LCINS_Pollution_mutations_multivariable.pdf',plot = p,width = 14,height = 5,device = cairo_pdf)


# Figure 4d ---------------------------------------------------------------

plotdata <- tdata %>% select(Tumor_Barcode,Telseq_TL_Ratio,Pollution_group2) %>%  filter(!is.na(Telseq_TL_Ratio)) %>%  left_join(covdata0) %>% left_join(verifybamID)

plotdata %>% lm(Telseq_TL_Ratio ~ Pollution_group2 + Age + Gender + PC1 + PC2 + Histology + Tumor_Purity, data=.) %>% tbl_regression(pvalue_fun = sci_notation) %>% as_gt() %>% gt::gtsave(filename = 'tmp.docx')

plotdata <- plotdata %>% 
  do(tidy(lm(Telseq_TL_Ratio ~ Pollution_group2 + Age + Gender + PC1 + PC2 + Histology + Tumor_Purity, data=.),conf.int = TRUE, conf.level = 0.95)) %>% 
  filter(term != '(Intercept)') %>% 
  #mutate(FDR=p.adjust(p.value,method='BH')) %>% 
  mutate(label=paste0('β = ',round(estimate,2), ', q-value = ',scientific_format(digits = 3)(p.value))) %>% 
  arrange(term) %>% 
  ungroup()

plotdata$term
plotdata$term <-  rep(c('Age','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Ancestry: PC1', 'Ancestry: PC2','Pollution Group: High\nRef=Low','Tumor Purity'),each=1)

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
  labs(x = "Linear regression coefficient", y = NULL,title = 'Telomere length')+
  guides(color="none")+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = 'LCINS_Pollution_telomere_multivariable.pdf',plot = p,width = 6,height = 5,device = cairo_pdf)


# Figure 4e ---------------------------------------------------------------
my_comparisons <- list(c("Low", "Intermediate"),c("Intermediate", "High"),c("Low", "High"))

plotdata <- tdata %>% select(Tumor_Barcode,SBS,DBS,ID,Pollution_group3) %>% pivot_longer(cols = -c(Tumor_Barcode,Pollution_group3)) %>% filter(is.finite(value)) %>% 
  mutate(name=factor(name,levels=c('SBS','DBS','ID')))


stat.test <- plotdata %>% 
  group_by(name) %>% 
  wilcox_test(value ~ Pollution_group3, comparisons = my_comparisons) %>%
  ungroup() %>% 
  add_xy_position(x = "Pollution_group3") %>%
  mutate(myformatted.p = sprintf("P = %.2e",p)) %>% 
  mutate(Pollution_group3 = "High")

stat.test

plotdata %>% 
  ggplot(aes(Pollution_group3,value,fill=Pollution_group3))+
  geom_quasirandom(pch=21,size=2,width = 0.3,color="white",stroke=0.4)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.7)+
  facet_wrap(~name,scales = 'free')+
  scale_fill_manual(values = pollution_colors)+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 13,axis_title_just = 'm',axis_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 15,ticks = T)+
  labs(x = expression("Pollution group (Population weighted PM"[2.5]*")"), y = 'Total number of mutations (log10)')+
  theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=14,face = 'bold',hjust = 0.5),panel.spacing = unit(0.4,'cm'))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_Pollution_mutations2.pdf',width = 8 ,height =6,device = cairo_pdf)  


# Figure 4f ---------------------------------------------------------------
plotdata <- tdata %>% select(Tumor_Barcode,Telseq_TL_Ratio,Pollution_group3) %>%  filter(!is.na(Telseq_TL_Ratio))

stat.test <- plotdata %>% 
  wilcox_test(Telseq_TL_Ratio ~ Pollution_group3, comparisons = my_comparisons) %>%
  ungroup() %>% 
  add_xy_position(x = "Pollution_group3") %>%
  mutate(myformatted.p = sprintf("P = %.2e",p)) %>% 
  mutate(Pollution_group3 = "High")

stat.test

plotdata %>% 
  ggplot(aes(Pollution_group3,Telseq_TL_Ratio,fill=Pollution_group3))+
  geom_quasirandom(pch=21,size=2,width = 0.3,color="white",stroke=0.4)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.7)+
  scale_fill_manual(values = pollution_colors)+
  scale_y_continuous(breaks = pretty_breaks())+
  theme_ipsum_rc(base_size = 13,axis_title_just = 'm',axis_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 15,ticks = T)+
  labs(x = expression("Pollution group (Population weighted PM"[2.5]*")"), y = 'Telomere length ratio, log2(Tumor TL/Normal TL)',title='Telomere length')+
  theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_text(size=14,face = 'bold',hjust = 0.5),panel.spacing = unit(0.4,'cm'))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_Pollution_telomere2.pdf',width = 3 ,height =6,device = cairo_pdf)  


# Figure 4g ---------------------------------------------------------------
#library(ggbreak)
testdata <- tdata %>% left_join(covdata0) %>% filter(is.finite(SBS))

plabel <- testdata %>% 
  do(tidy(lm(SBS ~ PM25 + Age + Gender + Assigned_Population + Histology + Tumor_Purity, data = .))) %>% 
  filter(term=='PM25') %>% 
  mutate(label=paste0('Multiple linear regression\nβ = ',round(estimate,2),sprintf(", q-value = %.2e",p.value))) %>% 
  pull(label)

tdata_group <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(SBS,na.rm = T),nsample=n_distinct(Tumor_Barcode)) %>%
  mutate(Burden= if_else(Country_pollution == 'Azerbaijan',4.2,Burden))

tdata_group2 <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(SBS,na.rm = T),nsample=n_distinct(Tumor_Barcode)) 


tmpbreak <- pretty_breaks()(tdata_group$Burden)
tmplabel <- tmpbreak
tmplabel[length(tmplabel)] <- paste0('>',tmplabel[length(tmplabel)])

tdata_group %>% 
  ggplot(aes(mean_pollution,(Burden)))+
  stat_smooth(data = tdata_group2,method="lm",fullrange=TRUE)+
  geom_point(aes(size=nsample),pch=21,fill='gray20',col='white',stroke=0.2)+
  ggrepel::geom_text_repel(aes(label=Country_pollution),size=3.2,col='gray30',force = 30,segment.color='gray35',max.overlaps = 30)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_size_binned(breaks = c(0,40,80,120,160,200))+
  scale_y_continuous(breaks = tmpbreak,labels = tmplabel,limits = c(NA,4.21))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = FALSE,axis = 'XY',ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  labs(x = expression("Average population weighted PM"[2.5] ~ (mu*g/m^3)), y = 'Average number of SBS mutations (log10)',size='Number of samples')+
  theme(legend.position = 'top',legend.key.width = unit(1,'cm'),legend.margin=margin(0,0,-12,0),legend.box.margin=margin(0,0,0,0))+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  annotate("text",x=28,y=3,label=plabel,family= 'Roboto Condensed')

ggsave(filename = 'Pollution_assocaition_SBS.pdf',width = 5,height = 5,device = cairo_pdf)


# Figure 4h ---------------------------------------------------------------
#library(ggbreak)
testdata <- tdata %>% left_join(covdata0) %>% filter(is.finite(DBS))

plabel <- testdata %>% 
  do(tidy(lm(DBS ~ PM25 + Age + Gender + Assigned_Population + Histology + Tumor_Purity, data = .))) %>% 
  filter(term=='PM25') %>% 
  mutate(label=paste0('Multiple linear regression\nβ = ',round(estimate,2),sprintf(", q-value = %.2e",p.value))) %>% 
  pull(label)

tdata_group <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(DBS,na.rm = T),nsample=n_distinct(Tumor_Barcode)) %>%
  mutate(Burden= if_else(Country_pollution == 'Azerbaijan',2,Burden))

tdata_group2 <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(DBS,na.rm = T),nsample=n_distinct(Tumor_Barcode)) 

tmpbreak <- pretty_breaks()(tdata_group$Burden)
tmplabel <- tmpbreak
tmplabel[length(tmplabel)] <- paste0('>',tmplabel[length(tmplabel)])

tdata_group %>% 
  ggplot(aes(mean_pollution,(Burden)))+
  stat_smooth(data = tdata_group2,method="lm",fullrange=TRUE)+
  geom_point(aes(size=nsample),pch=21,fill='gray20',col='white',stroke=0.2)+
  ggrepel::geom_text_repel(aes(label=Country_pollution),size=3.2,col='gray30',force = 30,segment.color='gray35',max.overlaps = 30)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_size_binned(breaks = c(0,40,80,120,160,200))+
  scale_y_continuous(breaks = tmpbreak,labels = tmplabel,limits = c(NA,2.01))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = FALSE,axis = 'XY',ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  labs(x = expression("Average population weighted PM"[2.5] ~ (mu*g/m^3)), y = 'Average number of DBS mutations (log10)',size='Number of samples')+
  theme(legend.position = 'top',legend.key.width = unit(1,'cm'),legend.margin=margin(0,0,-12,0),legend.box.margin=margin(0,0,0,0))+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  annotate("text",x=28,y=0.5,label=plabel,family= 'Roboto Condensed')

ggsave(filename = 'Pollution_assocaition_DBS.pdf',width = 5,height = 5,device = cairo_pdf)


# Figure 4i ---------------------------------------------------------------
#library(ggbreak)
testdata <- tdata %>% left_join(covdata0) %>% filter(is.finite(ID))

plabel <- testdata %>% 
  do(tidy(lm(ID ~ PM25 + Age + Gender + Assigned_Population + Histology + Tumor_Purity, data = .))) %>% 
  filter(term=='PM25') %>% 
  mutate(label=paste0('Multiple linear regression\nβ = ',round(estimate,2),sprintf(", q-value = %.2e",p.value))) %>% 
  pull(label)

tdata_group <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(ID,na.rm = T),nsample=n_distinct(Tumor_Barcode)) %>%
  mutate(Burden= if_else(Country_pollution == 'Azerbaijan',2.8,Burden))

tdata_group2 <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(ID,na.rm = T),nsample=n_distinct(Tumor_Barcode))

tdata_group %>% arrange(desc(Burden))

tmpbreak <- pretty_breaks()(tdata_group$Burden)
tmplabel <- tmpbreak
tmplabel[length(tmplabel)] <- paste0('>',tmplabel[length(tmplabel)])

tdata_group %>% 
  ggplot(aes(mean_pollution,(Burden)))+
  stat_smooth(data = tdata_group2,method="lm",fullrange=TRUE)+
  geom_point(aes(size=nsample),pch=21,fill='gray20',col='white',stroke=0.2)+
  ggrepel::geom_text_repel(aes(label=Country_pollution),size=3.2,col='gray30',force = 30,segment.color='gray35',max.overlaps = 30)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_size_binned(breaks = c(0,40,80,120,160,200))+
  scale_y_continuous(breaks = tmpbreak,labels = tmplabel,limits = c(NA,2.81))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = FALSE,axis = 'XY',ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  labs(x = expression("Average population weighted PM"[2.5] ~ (mu*g/m^3)), y = 'Average number of ID mutations (log10)',size='Number of samples')+
  theme(legend.position = 'top',legend.key.width = unit(1,'cm'),legend.margin=margin(0,0,-12,0),legend.box.margin=margin(0,0,0,0))+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  annotate("text",x=28,y=1.55,label=plabel,family= 'Roboto Condensed')

ggsave(filename = 'Pollution_assocaition_ID.pdf',width = 5,height = 5,device = cairo_pdf)


# Figure 4j ---------------------------------------------------------------
#library(ggbreak)
testdata <- tdata %>% left_join(covdata0) %>% filter(is.finite(ID))

plabel <- testdata %>% 
  do(tidy(lm(Telseq_TL_Ratio ~ PM25 + Age + Gender + Assigned_Population + Histology + Tumor_Purity, data = .))) %>% 
  filter(term=='PM25') %>% 
  mutate(label=paste0('Multiple linear regression\nβ = ',round(estimate,2),sprintf(", q-value = %.2e",p.value))) %>% 
  pull(label)

tdata_group <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(Telseq_TL_Ratio,na.rm = T),nsample=n_distinct(Tumor_Barcode)) %>%
  mutate(Burden= if_else(Country_pollution == 'Portugal',0.5,Burden))

tdata_group2 <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(Telseq_TL_Ratio,na.rm = T),nsample=n_distinct(Tumor_Barcode))

tdata_group %>% arrange(desc(Burden))

tmpbreak <- pretty_breaks()(tdata_group$Burden)
tmplabel <- tmpbreak
tmplabel[length(tmplabel)] <- paste0('>',tmplabel[length(tmplabel)])

tdata_group %>% 
  ggplot(aes(mean_pollution,(Burden)))+
  stat_smooth(data = tdata_group2,method="lm",fullrange=TRUE)+
  geom_point(aes(size=nsample),pch=21,fill='gray20',col='white',stroke=0.2)+
  ggrepel::geom_text_repel(aes(label=Country_pollution),size=3.2,col='gray30',force = 30,segment.color='gray35',max.overlaps = 30)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_size_binned(breaks = c(0,40,80,120,160,200))+
  scale_y_continuous(breaks = tmpbreak,labels = tmplabel,limits = c(NA,0.51))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = FALSE,axis = 'XY',ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  labs(x = expression("Average population weighted PM"[2.5] ~ (mu*g/m^3)), y = 'Telomere length ratio, log2(Tumor TL/Normal TL)',size='Number of samples')+
  theme(legend.position = 'top',legend.key.width = unit(1,'cm'),legend.margin=margin(0,0,-12,0),legend.box.margin=margin(0,0,0,0))+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  annotate("text",x=28.2,y=0.42,label=plabel,family= 'Roboto Condensed')

ggsave(filename = 'Pollution_assocaition_telomere.pdf',width = 5,height = 5,device = cairo_pdf)



# Figure 5a ---------------------------------------------------------------
testdata <- ludmil_activity_all_obs %>%
  pivot_longer(-Tumor_Barcode) %>% 
  mutate(Profile=str_remove_all(name,"[0-9a-z]")) %>% 
  filter(Profile %in% c('SBS','DBS','ID')) %>% 
  left_join(pollution_data) %>% 
  left_join(covdata0) %>% 
  left_join(verifybamID) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,!is.na(PM25)) %>% 
  mutate(PM25 = PM25/10) %>% 
  #mutate(value=as.factor(value)) %>% 
  group_by(name,Profile) %>% 
  do(tresult = safely(stats::glm)(value ~ PM25 + Histology + PC1 + PC2 + Gender + Age + Tumor_Purity,family = binomial(),data=.)) %>% 
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
  mutate(Profile=factor(Profile,levels=c('SBS','DBS','ID'),labels=c('SBS signatures','DBS signatures','ID signatures')))


testdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR),fill=name))+
  geom_hline(yintercept = -log10(0.05),col=ncicolpal[2],linetype=2)+
  geom_hline(yintercept = -log10(0.01),col=ncicolpal[1],linetype=2)+
  geom_vline(xintercept = 0,col="gray60",linetype=1,linewidth=0.5)+
  geom_point(pch=21,size=3.5,col='black',stroke=0.2)+
  ggrepel::geom_text_repel(data=testdata %>% filter(FDR<1),aes(label=name,col=name),force=20)+
  facet_wrap(~Profile,nrow = 1)+
  scale_fill_manual(values = sigcol)+
  scale_color_manual(values = sigcol)+
  scale_x_continuous(breaks = pretty_breaks(n = 5),limits = c(-2,2))+
  scale_y_continuous(breaks = pretty_breaks(n=7))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid = F,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_text(hjust = 0.5,face = 'bold',size = 14),panel.spacing = unit(0.3,'cm'))+
  labs(x = expression("Odd ratio (" * 10 ~ mu*g/m^3 * " of PM" [2.5] * ", log2)"), y = '-log10(FDR)')+
  guides(fill="none",color='none')+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  coord_cartesian(clip = 'off')

ggsave(filename = 'LCINS_pollution_signature_enrichment.pdf',width = 9,height = 3,device = cairo_pdf)


# Figure 5b ---------------------------------------------------------------
plotdata <- ludmil_activity_all_obs %>%
  pivot_longer(-Tumor_Barcode) %>% 
  mutate(Profile=str_remove_all(name,"[0-9a-z]")) %>% 
  left_join(pollution_data) %>% 
  left_join(covdata0) %>% 
  left_join(verifybamID) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,!is.na(PM25)) %>% 
  mutate(PM25 = PM25/10) %>% 
  filter(name %in% c('SBS4','SBS5','ID3')) 

plotdata %>% filter(name =='SBS4') %>% glm(value ~ PM25+Histology+ PC1 + PC2 + Gender + Age + Tumor_Purity,family = binomial(),data=.) %>% tbl_regression(pvalue_fun = sci_notation) %>% as_gt() %>% gt::gtsave(filename = 'tmp.docx')

plotdata %>% filter(name =='SBS5') %>% glm(value ~ PM25+Histology+ PC1 + PC2 + Gender + Age + Tumor_Purity,family = binomial(),data=.) %>% tbl_regression(pvalue_fun = sci_notation) %>% as_gt() %>% gt::gtsave(filename = 'tmp.docx')

plotdata %>% filter(name =='ID3') %>% glm(value ~ PM25+Histology+ PC1 + PC2 + Gender + Age + Tumor_Purity,family = binomial(),data=.) %>% tbl_regression(pvalue_fun = sci_notation) %>% as_gt() %>% gt::gtsave(filename = 'tmp.docx')

plotdata %>% 
  #mutate(value=as.factor(value)) %>% 
  group_by(name,Profile) %>% 
  do(tresult = safely(stats::glm)(value ~ PM25+Histology+ Assigned_Population + Gender + Age + Tumor_Purity,family = binomial(),data=.)) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],conf.int = TRUE, conf.level = 0.95,exponentiate = TRUE))) %>% 
  select(name,Profile,fit) %>% 
  unnest(cols = c(fit)) %>%
  filter(term!='(Intercept)') %>% 
  mutate(label=paste0('OR = ',round(estimate,2), ', q-value = ',scientific_format(digits = 3)(p.value))) %>% 
  arrange(term) %>% 
  ungroup()

plotdata$term
# plotdata$term <-  rep(c('Age','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Ancestry: PC1', 'Ancestry: PC2','PM[2.5] ~ mu*g/m^3','Tumor Purity'),each=3)
# 
# plotdata$term <-  rep(c('Age','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','Ancestry: PC1', 'Ancestry: PC2',"PM~2~.~5~ 10 μg/m^3^",'Tumor Purity'),each=3)

plotdata$term <-  rep(c('Age','Ancestry: EAS\nRef=EUR', 'Ancestry: Others\nRef=EUR','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD','PM[2.5] ~ mu*g/m^3','Tumor Purity'),each=3)

plotdata$term <-  rep(c('Age','Ancestry: EAS\nRef=EUR', 'Ancestry: Others\nRef=EUR','Sex: Female\nRef=Male','Histology: Carcinoid Tumor\nRef=LUAD','Histology: Others\nRef=LUAD','Histology: LUSC\nRef=LUAD',"PM~2~.~5~ 10 μg/m^3^",'Tumor Purity'),each=3)

termlevel <- plotdata %>% filter(name=='SBS4') %>% arrange(desc(term)) %>% pull(term)
p <- plotdata %>% 
  mutate(term=factor(term,levels=termlevel)) %>% 
  mutate(name=factor(name,levels=c('SBS4','SBS5','ID3'))) %>% 
  #mutate(p.value = if_else(p.value>0.05,NA,p.value)) %>% 
  mutate(p.value = if_else(str_detect(term,'PM'),p.value,NA)) %>% 
  ggplot(aes(log2(estimate), term, xmin = log2(conf.low), xmax = log2(conf.high), height = 0,color=-log10(p.value))) +
  geom_vline(xintercept = 0, lty = 5) +
  geom_errorbarh(height = .3)+
  geom_point(pch=19,size=3) +
  facet_wrap(~name,nrow = 1)+
  scale_color_gradient( low = "#B71B1BFF",high = "#B71B1BFF",na.value = ncicolpal[8])+ #low = "#EE5250FF"
  scale_x_continuous(breaks = pretty_breaks(n = 7))+
  #scale_y_discrete(labels = parse_format())+
  ggrepel::geom_text_repel(aes(label=label),size=3.5,nudge_y = 0.4)+
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid_col = 'gray85',plot_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5))+
  theme(panel.spacing = unit(0.3,"cm"),plot.title = element_text(hjust = 0.5),strip.text.x = element_text(face = 'bold',hjust = 0.5),axis.text.y = element_markdown())+
  #labs(x =  expression("Odd ratio (log2) PM" [2.5] * 10 ~ mu*g/m^3), y = NULL)+
  labs(x =  "Odd ratio (log2)", y = NULL)+
  guides(color="none")+
  panel_border(color = 'black',linetype = 1)

ggsave(filename = 'LCINS_pollution_signature_enrichment_examples_tmp.pdf',plot = p,width = 14,height = 4,device = cairo_pdf)



# Figure 5c ----------------------------------------------------------------
testdata <- ludmil_activity_all %>% 
  select(Tumor_Barcode, SBS4) %>% 
  left_join(pollution_data) %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,!is.na(PM25)) %>% 
  filter(SBS4>0) %>% 
  mutate(Burden=log10(SBS4))

plabel <- testdata %>% 
  do(tidy(lm(Burden ~ PM25 + Age + Gender + Assigned_Population + Histology + Tumor_Purity, data = .))) %>% 
  filter(term=='PM25') %>% 
  mutate(label=paste0('Multiple linear regression\nβ = ',round(estimate,2),sprintf(", q-value = %.2e",p.value))) %>% 
  pull(label)

tdata_group <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(Burden,na.rm = T),nsample=n_distinct(Tumor_Barcode)) 
# mutate(Burden= if_else(Country_pollution == 'Azerbaijan',2,Burden))

# tmpbreak <- pretty_breaks()(tdata_group$Burden)
# tmplabel <- tmpbreak
# tmplabel[length(tmplabel)] <- paste0('>',tmplabel[length(tmplabel)])

tdata_group %>% 
  ggplot(aes(mean_pollution,(Burden)))+
  stat_smooth(method="lm",fullrange=TRUE)+
  geom_point(aes(size=nsample),pch=21,fill='gray20',col='white',stroke=0.2)+
  ggrepel::geom_text_repel(aes(label=Country_pollution),size=3.2,col='gray30',force = 30,segment.color='gray35',max.overlaps = 30)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_size_binned(breaks = c(0,3,6,9,12))+
  #scale_y_continuous(breaks = tmpbreak,labels = tmplabel)+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = FALSE,axis = 'XY',ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  labs(x = expression("Average population weighted PM"[2.5] ~ (mu*g/m^3)), y = 'Average number of SBS4 mutations (log10)',size='Number of samples')+
  theme(legend.position = 'top',legend.key.width = unit(1,'cm'),legend.margin=margin(0,0,-12,0),legend.box.margin=margin(0,0,0,0))+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  annotate("text",x=25,y=2.8,label=plabel,family= 'Roboto Condensed')

ggsave(filename = 'Pollution_assocaition_SBS4.pdf',width = 5,height = 5,device = cairo_pdf)


# Figure 5d ----------------------------------------------------------------
testdata <- ludmil_activity_all %>% 
  select(Tumor_Barcode, SBS5) %>% 
  left_join(pollution_data) %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,!is.na(PM25)) %>% 
  filter(SBS5>0) %>% 
  mutate(Burden=log10(SBS5))

plabel <- testdata %>% 
  do(tidy(lm(Burden ~ PM25 + Age + Gender + Assigned_Population + Histology + Tumor_Purity, data = .))) %>% 
  filter(term=='PM25') %>% 
  mutate(label=paste0('Multiple linear regression\nβ = ',round(estimate,2),sprintf(", q-value = %.2e",p.value))) %>% 
  pull(label)

tdata_group <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(Burden,na.rm = T),nsample=n_distinct(Tumor_Barcode)) %>% 
  mutate(Burden= if_else(Country_pollution == 'Azerbaijan',3.6,Burden))

tdata_group2 <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(Burden,na.rm = T),nsample=n_distinct(Tumor_Barcode))  


tmpbreak <- pretty_breaks()(tdata_group$Burden)
tmplabel <- tmpbreak
tmplabel[length(tmplabel)] <- paste0('>',tmplabel[length(tmplabel)])

tdata_group %>% 
  ggplot(aes(mean_pollution,(Burden)))+
  stat_smooth(data = tdata_group2,method="lm",fullrange=TRUE)+
  geom_point(aes(size=nsample),pch=21,fill='gray20',col='white',stroke=0.2)+
  ggrepel::geom_text_repel(aes(label=Country_pollution),size=3.2,col='gray30',force = 30,segment.color='gray35',max.overlaps = 30)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_size_binned(breaks = c(0,40,80,120,160,200))+
  scale_y_continuous(breaks = tmpbreak,labels = tmplabel,limits = c(NA,3.61))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = FALSE,axis = 'XY',ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  labs(x = expression("Average population weighted PM"[2.5] ~ (mu*g/m^3)), y = 'Average number of SBS5 mutations (log10)',size='Number of samples')+
  theme(legend.position = 'top',legend.key.width = unit(1,'cm'),legend.margin=margin(0,0,-12,0),legend.box.margin=margin(0,0,0,0))+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  annotate("text",x=28,y=2.7,label=plabel,family= 'Roboto Condensed')

ggsave(filename = 'Pollution_assocaition_SBS5.pdf',width = 5,height = 5,device = cairo_pdf)




# Figure 5e ----------------------------------------------------------------
testdata <- ludmil_activity_all %>% 
  select(Tumor_Barcode, ID3) %>% 
  left_join(pollution_data) %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,!is.na(PM25)) %>% 
  filter(ID3>0) %>% 
  mutate(Burden=log10(ID3))

plabel <- testdata %>% 
  do(tidy(lm(Burden ~ PM25 + Age + Gender + Assigned_Population + Histology + Tumor_Purity, data = .))) %>% 
  filter(term=='PM25') %>% 
  mutate(label=paste0('Multiple linear regression\nβ = ',round(estimate,2),sprintf(", q-value = %.2e",p.value))) %>% 
  pull(label)

tdata_group <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(Burden,na.rm = T),nsample=n_distinct(Tumor_Barcode)) %>% 
  mutate(Burden= if_else(Country_pollution == 'Azerbaijan',2.6,Burden))

tdata_group2 <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(Burden,na.rm = T),nsample=n_distinct(Tumor_Barcode)) 

tmpbreak <- pretty_breaks()(tdata_group$Burden)
tmplabel <- tmpbreak
tmplabel[length(tmplabel)] <- paste0('>',tmplabel[length(tmplabel)])

tdata_group %>% 
  ggplot(aes(mean_pollution,(Burden)))+
  stat_smooth(data = tdata_group2,method="lm",fullrange=TRUE)+
  geom_point(aes(size=nsample),pch=21,fill='gray20',col='white',stroke=0.2)+
  ggrepel::geom_text_repel(aes(label=Country_pollution),size=3.2,col='gray30',force = 30,segment.color='gray35',max.overlaps = 30)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_size_binned(breaks = c(0,20,40,60,80,100))+
  scale_y_continuous(breaks = tmpbreak, labels = tmplabel,limits = c(NA,2.61))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = FALSE,axis = 'XY',ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  labs(x = expression("Average population weighted PM"[2.5] ~ (mu*g/m^3)), y = 'Average number of ID3 mutations (log10)',size='Number of samples')+
  theme(legend.position = 'top',legend.key.width = unit(1,'cm'),legend.margin=margin(0,0,-12,0),legend.box.margin=margin(0,0,0,0))+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  annotate("text",x=13,y=2.5,label=plabel,family= 'Roboto Condensed')

ggsave(filename = 'Pollution_assocaition_ID3.pdf',width = 5,height = 5,device = cairo_pdf)



# Figure 5 (ID8) ----------------------------------------------------------------
testdata <- ludmil_activity_all %>% 
  select(Tumor_Barcode, ID8) %>% 
  left_join(pollution_data) %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker,!is.na(PM25)) %>% 
  filter(ID8>0) %>% 
  mutate(Burden=log10(ID8))

plabel <- testdata %>% 
  do(tidy(lm(Burden ~ PM25 + Age + Gender + Assigned_Population + Histology + Tumor_Purity, data = .))) %>% 
  filter(term=='PM25') %>% 
  mutate(label=paste0('Multiple linear regression\nβ = ',round(estimate,2),sprintf(", q-value = %.2e",p.value))) %>% 
  pull(label)

tdata_group <- testdata %>% 
  group_by(Country_pollution) %>% 
  summarise(mean_pollution=mean(PM25,na.rm = T),Burden = mean(Burden,na.rm = T),nsample=n_distinct(Tumor_Barcode)) #%>% 
#mutate(Burden= if_else(Country_pollution == 'Azerbaijan',3.6,Burden))

# tmpbreak <- pretty_breaks()(tdata_group$Burden)
# tmplabel <- tmpbreak
# tmplabel[length(tmplabel)] <- paste0('>',tmplabel[length(tmplabel)])

tdata_group %>% 
  ggplot(aes(mean_pollution,(Burden)))+
  stat_smooth(method="lm",fullrange=TRUE)+
  geom_point(aes(size=nsample),pch=21,fill='gray20',col='white',stroke=0.2)+
  ggrepel::geom_text_repel(aes(label=Country_pollution),size=3.2,col='gray30',force = 30,segment.color='gray35',max.overlaps = 30)+
  scale_x_continuous(breaks = pretty_breaks())+
  scale_size_binned(breaks = c(0,40,80,120,160,200))+
  scale_y_continuous(breaks = pretty_breaks())+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14,grid = FALSE,axis = 'XY',ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  labs(x = expression("Average population weighted PM"[2.5] ~ (mu*g/m^3)), y = 'Average number of SBS5 mutations',size='Number of samples')+
  theme(legend.position = 'top',legend.key.width = unit(1,'cm'),legend.margin=margin(0,0,-12,0),legend.box.margin=margin(0,0,0,0))+
  panel_border(color = 'black',linetype = 1,size=0.5)+
  annotate("text",x=28,y=1.2,label=plabel,family= 'Roboto Condensed')

ggsave(filename = 'Pollution_assocaition_ID8.pdf',width = 5,height = 5,device = cairo_pdf)


# Figure 5f --------------------------------------------------------------
testdata <- sherlock_data_full %>% 
  filter(Gene %in% drglist, Type=='Mutation_Driver') %>%
  left_join(pollution_data) %>% 
  left_join(covdata0) %>% 
  left_join(verifybamID) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker, !is.na(PM25)) %>% 
  mutate(PM25=PM25/10)

plotdata <- testdata %>% 
  group_by(Gene) %>% 
  mutate(Alteration=as.factor(Alteration)) %>% 
  do(tresult = safely(glm)(Alteration ~ PM25 + Histology + PC1 + PC2 + Gender + Age + Tumor_Purity, family='binomial',data=. )) %>% 
  mutate(tresult_null = map_lgl(tresult['result'], is.null)) %>% 
  filter(!tresult_null) %>% 
  mutate(fit = list(tidy(tresult[['result']],exponentiate=TRUE))) %>% 
  select(Gene,fit) %>% 
  unnest(cols = c(fit)) %>% 
  ungroup() %>% 
  arrange(p.value) %>% 
  filter(term=='PM25') %>% 
  filter(abs(log(estimate))<10) %>%
  mutate(FDR=p.adjust(p.value,method='BH'))


plotdata %>% 
  ggplot(aes(log2(estimate),-log10(FDR)))+
  geom_hline(yintercept = -log10(0.01),linetype=2,col=ncicolpal[1],size=0.5)+
  geom_hline(yintercept = -log10(0.05),linetype=2,col=ncicolpal[2],size=0.5)+
  geom_vline(xintercept = 0,linetype=1,col='gray50',size=0.5)+
  geom_point(aes(fill=log2(estimate)>0),pch=21,stroke=0.1,col='black',size=4)+
  scale_x_continuous(breaks = pretty_breaks(n = 7),limits = c(-4,4))+
  scale_y_continuous(breaks = pretty_breaks())+
  scale_fill_manual(values = as.character(pollution_colors[c(1,3)]))+
  ggrepel::geom_text_repel(data=plotdata %>% filter(FDR<0.5),aes(label=Gene),force = 10,max.overlaps = 30,size=4.5)+
  labs(x = expression("Odd ratio (" * 10 ~ mu*g/m^3 * " of PM" [2.5] * ", log2)"), y = '-log10(FDR)')+  guides(fill = "none")+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid=FALSE,ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 15)+
  theme(plot.title = element_text(hjust = 0.5))+
  panel_border(color = 'black',size = 0.5)

ggsave(filename = 'LCINS_Pollution_driver_mutations.pdf',width = 6,height = 5,device = cairo_pdf)


# figure 5g-h -------------------------------------------------------------
source('../../Sherlock_functions.R')

plotdata <- sherlock_data_full %>% 
  filter(Gene %in% drglist, Type=='Mutation_Driver') %>%
  left_join(pollution_data) %>% 
  left_join(covdata0) %>% 
  filter(Tumor_Barcode %in% sherlock_nonsmoker, !is.na(Pollution_group2)) 

barplot_fisher(plotdata %>% filter(Gene=='TP53'),'Pollution_group2','Alteration',var1lab =  expression("Pollution group (Population weighted PM"[2.5]*")"),var2lab = 'TP53 mutation status',filename = 'LCINS_Pollution_driver_TP53.pdf')

barplot_fisher(plotdata %>% filter(Gene=='CTNNB1'),'Pollution_group2','Alteration',var1lab =  expression("Pollution group (Population weighted PM"[2.5]*")"),var2lab = 'CTNNB1 mutation status',filename = 'LCINS_Pollution_driver_CTNNB1.pdf')




# Save R object for Macros ------------------------------------------------
wgs_groups_info <- wgs_groups_info %>% filter(Tumor_Barcode %in% sherlock_nonsmoker)
BBsolution4 <- BBsolution4 %>% filter(Tumor_Barcode %in% sherlock_nonsmoker)
covdata0 <- covdata0 %>% filter(Tumor_Barcode %in% sherlock_nonsmoker)
sherlock_PGA <-  sherlock_PGA %>% filter(Tumor_Barcode %in% sherlock_nonsmoker)
clinical_data <-  clinical_data %>% filter(Tumor_Barcode %in% sherlock_nonsmoker)
sherlock_variable <- sherlock_variable %>% filter(Tumor_Barcode %in% sherlock_nonsmoker)
sherlock_driver_mutations <- sherlock_driver_mutations %>% filter(Tumor_Barcode %in% sherlock_nonsmoker)
Mutation_Signature_Probability_SBS <- Mutation_Signature_Probability_SBS %>% filter(Tumor_Barcode %in% sherlock_nonsmoker,Gene_Name %in% drglist)
sherlock_data_full <- sherlock_data_full %>% filter(Gene %in% drglist, Type=='Mutation_Driver')


wgs_groups_info 
BBsolution4
covdata0
sherlock_PGA
clinical_data
sherlock_variable
sherlock_driver_mutations
Mutation_Signature_Probability_SBS
sherlock_data_full
drglist
ludmil_activity_all_obs
ludmil_activity_all
pollution_data
pollution_colors
Barcode2map
covdata0
sp_group_data2
sp_group_color_new
b2map
genome2size
ncicolpal
TMBplot2
sherlock_sbs288_profile
sherlock_sbs96_profile
sherlock_dbs78_profile
sherlock_id83_profile
prevalence_plot
sv38_activity
cn68_activity
PieDonut_ztw
barplot_fisher


save(b2map,Barcode2map,barplot_fisher,BBsolution4,clinical_data,cn68_activity,covdata0,drglist,genome2size,ludmil_activity_all,ludmil_activity_all_obs,Mutation_Signature_Probability_SBS,ncicolpal,PieDonut_ztw,pollution_colors,pollution_data,prevalence_plot,sherlock_data_full,sherlock_dbs78_profile,sherlock_driver_mutations,sherlock_id83_profile,sherlock_PGA,sherlock_sbs288_profile,sherlock_sbs96_profile,sherlock_variable,sp_group_color_new,sp_group_data2,sv38_activity,TMBplot2,wgs_groups_info,file='Main_Figures.RData')


