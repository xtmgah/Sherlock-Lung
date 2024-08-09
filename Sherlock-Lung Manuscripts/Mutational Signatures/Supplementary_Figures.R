set_wd()
libztw()
pdfhr2()
myggstyle()


load('../sampleset.RData')

# SFigure 1 ---------------------------------------------------------------

# CHORD vs HRDetect
load('/Volumes/data/NSLC2/CHORD/Sherlock_1217/chord_output.RData')
load('/Volumes/data/NSLC2/Signature_tools_lib/Sherlock_1217/Signature_tools_lib.RData')

hrdetect_result <-
  as.data.frame(res$hrdetect_output) %>% rownames_to_column(var = 'Tumor_Barcode') %>%
  arrange(desc(Probability)) %>%
  mutate(HRDetect_Status = if_else(Probability >0.7, 'High',if_else(Probability > 0.005, 'Moderate','Low')))

left_join(
  hrdetect_result %>% select(Tumor_Barcode, HRDetect=Probability),
  chord_output %>% select(Tumor_Barcode=sample,CHORD=p_hrd,hrd_type)
) %>%
  left_join(ludmil_activity_all_obs %>% select(Tumor_Barcode, SBS3) %>% mutate(SBS3=if_else(SBS3,'Present','Absent')))%>%
  left_join(covdata0) %>%
  filter(Tumor_Barcode %in% sherlock_nonsmoker) %>% 
  arrange(SBS3) %>% 
  ggplot(aes(HRDetect, CHORD, fill=Histology,shape=SBS3))+geom_point(size=2.5,stroke=0.25)+
  scale_x_continuous(breaks = pretty_breaks(),limits = c(0,1))+
  scale_y_continuous(breaks = pretty_breaks(),limits = c(0,1))+
  scale_fill_manual(name='Histology',values = ncicolpal[1:4])+
  scale_shape_manual(name='Mutational signature SBS3',values = c(21,25))+
  geom_hline(yintercept = 0.5,linetype=2)+
  geom_vline(xintercept = 0.7,linetype=2)+
  theme_ipsum_rc(base_size = 14,axis_title_just = 'm',axis_title_size = 16,grid = FALSE,ticks = TRUE)+
  guides(fill = guide_legend(override.aes = list(pch=21,size=4)),shape = guide_legend(override.aes = list(size=4,fill='black'))) +
  panel_border(color = 'black')+
  labs(x='HRDetect probability Score', y = 'CHORD probability Score')

ggsave(filename = 'HRDetect_vs_chord3.pdf',width = 8,height = 5,device = cairo_pdf)


# SFigure 2 ---------------------------------------------------------------




# # Sfigure 3 -------------------------------------------------------------
# pollution_data <- read_csv('../../Exposure_data/Sherlock_Non-Smoker_PM2.5_groups_2023OCT6.csv',col_names = T)
# 
# pollution_data <- pollution_data %>% 
#   select(Tumor_Barcode,Country, Country_pollution,PM25=Population.Weighted.PM2.5.ug.m3,Pollution_group2=Pollution_group) %>% 
#   mutate(Pollution_group3=case_when(
#     PM25<=10.01 ~ 'Low',
#     PM25>10.01 & PM25<24.57  ~ 'Intermediate',
#     PM25>=24.57 ~ 'High',
#     TRUE ~ NA_character_
#   )) %>% 
#   mutate(Pollution_group2 = factor(Pollution_group2,levels=c('Low','High'))) %>% 
#   mutate(Pollution_group3 = factor(Pollution_group3,levels=c('Low','Intermediate','High')))
# 
# pollution_colors <- c('#01665e','#947100','#BB0E3D')
# names(pollution_colors) <- c('Low','Intermediate','High')
# 
# pollution_data <- pollution_data %>% filter(Tumor_Barcode %in% sherlock_nonsmoker) %>% filter(!is.na(PM25))

#dens <- density(pollution_data$PM25)

mdata <- read_delim('~/NIH-Work/EAGLE_NSLC/SecondPaper/Mutational_Signature_Manuscript/Code_and_Data_Sherlock_mutational_signatures_manuscript_2024FEB8/Met_Pol_TMB_PGA_Sigs_SBS_ID_DBS_CN_SV_2024FEB7.tsv',delim = '\t',col_names = T)

dens <- density(mdata$Population.Weighted.PM2.5.ug.m3,na.rm = T)
df <- data.frame(x=dens$x, y=dens$y)
probs <- c(0,0.2,0.4,0.5,0.6,0.8,1)
quantiles <- quantile(mdata$Population.Weighted.PM2.5.ug.m3, prob=probs,na.rm = T)
df$quant <- factor(findInterval(df$x,quantiles[c(2,4,6)]))

tmp <- df %>% group_by(quant) %>% slice(1) %>% ungroup() %>% slice(2:4)
tmp

ggplot(df, aes(x,y)) + 
  geom_line(col='black') + 
  geom_area(aes(fill=quant)) +
  geom_segment(data = tmp,aes(x=x,xend=x,y=0,yend=y),col='white')+
  scale_x_continuous(breaks = pretty_breaks(n=7),expand = c(0,0))+
  scale_y_continuous(breaks = pretty_breaks(n=7),expand = expansion(add=c(0,0.005)))+
  scale_fill_manual(values = c(c('#01665e','#947100','#947100','#BB0E3D')))+
  annotate("text",x=tmp$x[1],y=0,label='Split in 3 groups (low-intermediate)',angle=90,hjust=-0.05,vjust=1.5,col='white')+
  annotate("text",x=tmp$x[2],y=0,label='Split in 2 groups (low-high)',angle=90,hjust=-0.05,vjust=1.5,col='white')+
annotate("text",x=tmp$x[3],y=0,label='Split in 3 groups (intermediate-high)',angle=90,hjust=-0.05,vjust=1.5,col='white')+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid = 'Yy',ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  guides(fill="none")+
  labs(x=expression("Average population weighted PM"[2.5] ~ (mu*g/m^3)),y="Sample density")+
  panel_border(color = 'black',linetype = 1,size=0.5)

ggsave(filename = 'Pollution_groups.pdf',width = 6,height = 5,device = cairo_pdf)


df$quant <- factor(findInterval(df$x,quantiles[c(4)]))
tmp <- df %>% group_by(quant) %>% slice(1) %>% ungroup() %>% slice(2)
tmp

ggplot(df, aes(x,y)) + 
  geom_line(col='black') + 
  geom_area(aes(fill=quant)) +
  geom_segment(data = tmp,aes(x=x,xend=x,y=0,yend=y),col='white')+
  scale_x_continuous(breaks = pretty_breaks(n=7),expand = c(0,0))+
  scale_y_continuous(breaks = pretty_breaks(n=7),expand = expansion(add=c(0,0.005)))+
  scale_fill_manual(values = c(c('#01665e','#BB0E3D')))+
  annotate("text",x=tmp$x[1],y=0,label='Split in 2 groups (low-high)',angle=90,hjust=-0.05,vjust=1.5,col='white')+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid = 'Yy',ticks = T,plot_margin=margin(5.5,5.5,5.5,5.5))+
  guides(fill="none")+
  labs(x=expression("Average population weighted PM"[2.5] ~ (mu*g/m^3)),y="Sample density")+
  panel_border(color = 'black',linetype = 1,size=0.5)

ggsave(filename = 'Pollution_groups2.pdf',width = 6,height = 5,device = cairo_pdf)



# TMB vs EGFR in the context of TP53 mutations ----------------------------
library(ggbeeswarm)
library(rstatix)

load('../../RDS/sherlock_variable.RData')
load('../../RDS/sherlock_data_all.RData')

my_comparisons <- list(c("No", "Yes"))

tdata <- sherlock_data_full %>%
  filter(Gene %in% c('EGFR','TP53'), Type=='Mutation_Driver') %>% 
  mutate(Alteration = if_else(Alteration == 'Yes','Mutant','Wild-type')) %>% 
  mutate(Alteration = factor(Alteration,levels=c('Wild-type','Mutant'))) %>% 
  pivot_wider(names_from = Gene,values_from = Alteration) %>% 
  mutate(TP53=factor(TP53,labels =c('TP53 Wild-type','TP53 Mutant'))) %>% 
  left_join(
    sherlock_variable %>% filter(name == 'TMB') %>% pivot_wider()
  ) %>% 
  mutate(TMB=log2(TMB)) %>% 
  filter(Tumor_Barcode %in% luad_nonsmoker ) %>% 
  left_join(covdata0)

tdata %>% group_by(TP53) %>% do(tidy(lm(TMB ~ EGFR + Tumor_Purity + Assigned_Population + Gender + Age, data=.)))
  
stat.test <- tdata %>% 
  group_by(TP53) %>% 
  wilcox_test(TMB ~ EGFR, comparisons = my_comparisons) %>%
  add_xy_position(x = "EGFR") %>%
  mutate(myformatted.p = sprintf("P = %.2e",p)) %>% 
  mutate(EGFR = NA)

stat.test$myformatted.p <- c("q-value = 1.15e-04", "q-value = 0.15")

tdata %>% 
  ggplot(aes(EGFR,TMB,fill=EGFR))+
  geom_quasirandom(pch=21,size=2,width = 0.3,color="white",stroke=0.5)+
  geom_boxplot(width=0.5,outlier.shape = NA,alpha =0.6,size=0.7)+
  scale_fill_jama()+
  scale_y_continuous(breaks = pretty_breaks())+
  facet_wrap(~TP53)+
  theme_ipsum_rc(base_size = 13,axis_title_just = 'm',axis_title_size = 15,plot_margin=margin(5.5,5.5,5.5,5.5),plot_title_size = 15,ticks = T)+
  labs(x = "EGFR mutation status", y = 'Tumor mutational burden (log2)')+
  theme(plot.title = element_text(hjust = 0.5),strip.text.x = element_text(face = 'bold',hjust = 0.5))+
  guides(fill="none")+
  panel_border(color = 'black',linetype = 1)+
  stat_pvalue_manual(stat.test, label = "myformatted.p")

ggsave(filename = 'LCINS_EGFR_TMB.pdf',width = 4 ,height =6,device = cairo_pdf)  



# Save object for Macros --------------------------------------------------

hrdetect_result <- hrdetect_result %>% filter(Tumor_Barcode %in% sherlock_nonsmoker)
chord_output <- chord_output %>% filter(sample %in% sherlock_nonsmoker)
mdata

save(mdata,hrdetect_result, chord_output,sherlock_nonsmoker, file='Supplementary_Figures.RData')
