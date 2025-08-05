require(tidyverse)
require(ggtext)
require(ggforce)
library(janitor)
library(ggpubr)
library(cowplot)
library(hrbrthemes)
library(ggsci)
library(ggrepel)




# Conflicted package ------------------------------------------------------
resolve_conflicts <- function(){
  require(conflicted)
  conflict_prefer("select", "dplyr")
  conflict_prefer("filter", "dplyr")
  #filter <- dplyr::filter
  #select <- dplyr::select
}


### reformat the signature name ###
signames <- 
  tibble(
    SBSname=c(
      "SBS1",
      "SBS2",
      "SBS3",
      "SBS4",
      "SBS5",
      "SBS6",
      "SBS7a",
      "SBS7b",
      "SBS7c",
      "SBS7d",
      "SBS8",
      "SBS9",
      "SBS10a",
      "SBS10b",
      "SBS10c",
      "SBS11",
      "SBS12",
      "SBS13",
      "SBS14",
      "SBS15",
      "SBS16",
      "SBS17a",
      "SBS17b",
      "SBS18",
      "SBS19",
      "SBS20",
      "SBS21",
      "SBS22",
      "SBS23",
      "SBS24",
      "SBS25",
      "SBS26",
      "SBS27",
      "SBS28",
      "SBS29",
      "SBS30",
      "SBS31",
      "SBS32",
      "SBS33",
      "SBS34",
      "SBS35",
      "SBS36",
      "SBS37",
      "SBS38",
      "SBS39",
      "SBS40",
      "SBS41",
      "SBS42",
      "SBS43",
      "SBS44",
      "SBS45",
      "SBS46",
      "SBS47",
      "SBS48",
      "SBS49",
      "SBS50",
      "SBS51",
      "SBS52",
      "SBS53",
      "SBS54",
      "SBS55",
      "SBS56",
      "SBS57",
      "SBS58",
      "SBS59",
      "SBS60",
      "SBS84",
      "SBS85",
      "SBS92",
      "SBS-others"
    ),
    Subsname=c(
      "Signature Subs-01",
      "Signature Subs-02",
      "Signature Subs-03",
      "Signature Subs-04",
      "Signature Subs-05",
      "Signature Subs-06",
      "Signature Subs-07a",
      "Signature Subs-07b",
      "Signature Subs-07c",
      "Signature Subs-07d",
      "Signature Subs-08",
      "Signature Subs-09",
      "Signature Subs-10a",
      "Signature Subs-10b",
      "Signature Subs-10c",
      "Signature Subs-11",
      "Signature Subs-12",
      "Signature Subs-13",
      "Signature Subs-14",
      "Signature Subs-15",
      "Signature Subs-16",
      "Signature Subs-17a",
      "Signature Subs-17b",
      "Signature Subs-18",
      "Signature Subs-19",
      "Signature Subs-20",
      "Signature Subs-21",
      "Signature Subs-22",
      "Signature Subs-23",
      "Signature Subs-24",
      "Signature Subs-25",
      "Signature Subs-26",
      "Signature Subs-27",
      "Signature Subs-28",
      "Signature Subs-29",
      "Signature Subs-30",
      "Signature Subs-31",
      "Signature Subs-32",
      "Signature Subs-33",
      "Signature Subs-34",
      "Signature Subs-35",
      "Signature Subs-36",
      "Signature Subs-37",
      "Signature Subs-38",
      "Signature Subs-39",
      "Signature Subs-40",
      "Signature Subs-41",
      "Signature Subs-42",
      "Signature Subs-43",
      "Signature Subs-44",
      "Signature Subs-45",
      "Signature Subs-46",
      "Signature Subs-47",
      "Signature Subs-48",
      "Signature Subs-49",
      "Signature Subs-50",
      "Signature Subs-51",
      "Signature Subs-52",
      "Signature Subs-53",
      "Signature Subs-54",
      "Signature Subs-55",
      "Signature Subs-56",
      "Signature Subs-57",
      "Signature Subs-58",
      "Signature Subs-59",
      "Signature Subs-60",
      "Signature Subs-84",
      "Signature Subs-85",
      "Signature Subs-92",
      "Signature Subs-others"
      
    )
  )



SBS2Subs <- signames$Subsname
names(SBS2Subs) <- signames$SBSname
Subs2SBS <- signames$SBSname
names(Subs2SBS) <- signames$Subsname


### define the signature funciton ####
Subscolor <- c(
  'Signature Subs-01'='#4a9855',
  'Signature Subs-02'='#e2a8ab',
  'Signature Subs-03'='#40004b',
  'Signature Subs-04'='#5aa1ca',
  'Signature Subs-05'='#305d39',
  'Signature Subs-06'='#785940',
  "Signature Subs-07a"='#6e70b7',
  "Signature Subs-07b"='#ff7f00',
  "Signature Subs-07c"='#fec44f',
  "Signature Subs-07d"='#846a2a',
  "Signature Subs-08"='#cab2d6',
  "Signature Subs-09"='#f4a582',
  "Signature Subs-10a"='#8dd3c7',
  "Signature Subs-10b"='#5e4fa2',
  "Signature Subs-10c"='#761429',
  'Signature Subs-12'='#ffed6f',
  'Signature Subs-13'='#e41a1c',
  'Signature Subs-14'='#ffffbf',
  'Signature Subs-15'='#4d4d4d',
  'Signature Subs-16'='#513276',
  'Signature Subs-17a'='#df4c7d',
  'Signature Subs-17b'='#08519c',
  'Signature Subs-18'='#b3de69',
  'Signature Subs-19'='#dfc27d',
  'Signature Subs-20'='#b2182b',
  'Signature Subs-22'='#01665e',
  'Signature Subs-21'='#9ecae1',
  'Signature Subs-24'='#1c9099',
  'Signature Subs-25'='#35978f',
  'Signature Subs-26'='#ec7014',
  'Signature Subs-28'='#de77ae',
  'Signature Subs-30'='#d9d9d9',
  'Signature Subs-31'='#f781bf',
  'Signature Subs-32'='#dd1c77',
  'Signature Subs-33'='#b25d7e',
  'Signature Subs-35'='#fc8d59',
  'Signature Subs-36'='yellow',
  'Signature Subs-39'='#636363',
  'Signature Subs-40'='#b15928',
  'Signature Subs-41'='#fccde5',
  'Signature Subs-44'='#8c6bb1',
  'Signature Subs-46'='#e6f598',
  'Signature Subs-47'='#bababa',
  'Signature Subs-42'='#ae017e',
  'Signature Subs-54'='#fcc5c0',
  'Signature Subs-56'='#8c510a',
  'Signature Subs-84'='#063C3C',
  'Signature Subs-85'='#AA9139',
  'Signature Subs-92'='#0E1844',
  'Signature Subs-others'='#cececa'
)


tmpcolor <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2','#40004b','#762a83','#9970ab','#c2a5cf','#e7d4e8','#d9f0d3','#a6dba0','#5aae61','#1b7837','#00441b')
names(tmpcolor) <- c('Signature Subs-11','Signature Subs-23','Signature Subs-27','Signature Subs-29','Signature Subs-34','Signature Subs-37','Signature Subs-38','Signature Subs-43','Signature Subs-45','Signature Subs-48','Signature Subs-49','Signature Subs-50','Signature Subs-51','Signature Subs-52','Signature Subs-53','Signature Subs-55','Signature Subs-57','Signature Subs-58','Signature Subs-59','Signature Subs-60')

Subscolor <- c(Subscolor,tmpcolor)

SBScolor <- Subscolor
names(SBScolor) <- Subs2SBS[names(SBScolor)]


### color 12 for the clustering
color12 <- rev(c(
  "#a6cee3",
  "#1f78b4",
  "#b2df8a",
  "#33a02c",
  "#fb9a99",
  "#e31a1c",
  "#fdbf6f",
  "#ff7f00",
  "#cab2d6",
  "#6a3d9a",
  "#b15928"
))


### SBS cluster funciton ####
# p4 proportion bar

SBS96_Clustering <- function(sigdata,studydata=NULL,cosinedata,cosine_cutoff=0,sigcolor=NULL,studycolor=NULL,puritydata=NULL,puritydata_cat=FALSE, puritycol=NULL,purity_cutoff=0,clustern,highlight=NULL,legendnrow=3,sampletext_size=6, filename='SBS96_clustering.pdf',height=10,width=20,hc_func='hclust',hc_metric = 'euclidean',hc_method = 'ward.D2',stand=TRUE){
  
  require(tidyverse)
  require(scales)
  require(janitor)
  require(ggsci)
  require(ggpubr)
  require(factoextra)
  require(cowplot)
  
  blankcol <- NA
  names(blankcol) <- "        "
  
  if(is.null(sigcolor)){
    sigcolorindex <- as.character(Subscolor[colnames(sigdata)[-1]])
    names(sigcolorindex) <- colnames(sigdata)[-1]
  }else{
    sigcolorindex <- sigcolor
  }
  
  # colorall <- c(sigcolorindex,"white",studycolor)
  # names(colorall) <- c(names(sigcolorindex),"        ",names(studycolor))
  # 
  colorall <- c(sigcolorindex,blankcol,studycolor)
  
  tmp=sigdata %>% adorn_percentages('row') 
  mdata=as.matrix(tmp[,-1])
  rownames(mdata) <- tmp$Samples
  
  #fviz_nbclust(mdata, kmeans, method = "gap_stat")
  kcolors <- pal_d3("category20")(clustern)
  
  res <- hcut(mdata,k = clustern,hc_func = hc_func,hc_metric = hc_metric,hc_method = hc_method,stand=stand)
  p1 <- fviz_dend(res, rect = TRUE, cex = 0.5,k_colors = kcolors,lwd = 0.5,show_labels = F)+scale_x_discrete(expand = c(0,0))+ theme(plot.margin=margin(b=-0.7,unit="cm"),title = element_blank())
  #plot.margin=margin(b=-1,unit="cm")
  
  if(!is.null(studydata)){
    studydata <- studydata %>% filter(Samples %in% res$labels) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order]))
    p2 <- studydata %>% ggplot(aes(Samples,1,fill=factor(Study,levels = names(studycolor))))+geom_tile(col="black")+scale_fill_manual(values =studycolor,drop=FALSE)+theme_minimal()+theme(legend.position = "none",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t=0.4,unit="cm"),title = element_blank())+ylim(c(0,2))
  }
  p2.1 <- cosinedata %>%  mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% mutate(Similarity=if_else(Similarity<cosine_cutoff,NA_real_,Similarity)) %>% ggplot(aes(Samples,1,fill=Similarity))+geom_tile(col="black")+scale_fill_viridis_c(na.value = "#cccccc",option = 'C',limits = c(0.6, 1), oob = scales::squish)+theme_minimal()+theme(legend.position = "bottom",legend.key.width =unit(3, "cm"),legend.key.height = unit(0.3,"cm"), panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(title = element_blank())+ylim(c(0,2))
  p2.2 <- as_ggplot(get_legend(p2.1))
  p2.2 <- p2.2+theme(plot.margin = margin(b = 30))
  p2.1 <- cosinedata %>%  mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% mutate(Similarity=if_else(Similarity<cosine_cutoff,NA_real_,Similarity)) %>% ggplot(aes(Samples,1,fill=Similarity))+geom_tile(col="black")+scale_fill_viridis_c(na.value = "#cccccc",option = 'C',limits = c(0.6, 1), oob = scales::squish)+theme_minimal()+theme(legend.position = "none", panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t = 0.5,unit="cm"),title = element_blank())+ylim(c(0,2))
  
  p2.5 <- p2.2
  
  if(!is.null(puritydata)){
    puritydata <- puritydata %>% filter(Samples %in% res$labels) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) 
    
    if(!puritydata_cat) {
      puritydata <- puritydata %>% mutate(Purity=if_else(Purity<purity_cutoff,NA_real_,Purity))
      p2.3 <- puritydata  %>% ggplot(aes(Samples,1,fill=Purity))+geom_tile(col="black")+scale_fill_viridis_c(na.value = "#cccccc",option = 'D',breaks=pretty_breaks(5))+theme_minimal()+theme(legend.position = "bottom",legend.key.width =unit(3, "cm"),legend.key.height = unit(0.3,"cm"), panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(title = element_blank())+ylim(c(0,2))
      p2.4 <- as_ggplot(get_legend(p2.3))
      p2.4 <- p2.4+theme(plot.margin = margin(b = 30))
      p2.3 <- puritydata %>% ggplot(aes(Samples,1,fill=Purity))+geom_tile(col="black")+scale_fill_viridis_c(na.value = "#cccccc",option = 'D',breaks=pretty_breaks(5))+theme_minimal()+theme(legend.position = "none", panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t = 0.5,unit="cm"),title = element_blank())+ylim(c(0,2))
      #
      p2.5 <- plot_grid(p2.2, p2.4, align = "h", axis = "b", rel_widths = c(1, 1))
      
    }else { 
      names(blankcol) <- "    "
      colorall <- c(colorall,blankcol,puritycol)
      p2.3 <- puritydata %>% ggplot(aes(Samples,1,fill=factor(Purity,levels = names(puritycol))))+geom_tile(col="black")+scale_fill_manual(values =puritycol)+theme_minimal()+theme(legend.position = "none",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t=0.5,unit="cm"),title = element_blank())+ylim(c(0,2))
      
    }
  }
  
  
  p3 <- sigdata %>% gather(Signature,Weight,-Samples) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,(Weight),fill=factor(Signature,levels = names(colorall))))+geom_bar(stat="identity",col="gray95",width=1,size=0.1)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = "xy")+theme(legend.title = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_blank(),panel.grid.major.x=element_blank(),legend.position = "none",legend.box.spacing = unit(0,"cm"),plot.margin=margin(b=-1,t = 1,unit="cm"))+scale_y_continuous(expand = c(0, 0),breaks = pretty_breaks(),labels = comma)+scale_fill_manual(values = colorall,drop=FALSE)+guides(fill=guide_legend(nrow=2,byrow=TRUE))+xlab("")+ylab("Number of mutations \n")
  #+theme(plot.margin=margin(b=4,unit="pt"))
  
  p4 <- sigdata %>% gather(Signature,Weight,-Samples) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,Weight,fill=factor(Signature,levels = names(colorall))))+geom_bar(stat="identity",position="fill",col="gray95",width = 1,size=0.1)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12)+theme(legend.title = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_text(size = sampletext_size, angle = 90, vjust = 0.5, hjust = 1),panel.grid.major=element_line(),legend.position = "bottom",legend.box.background = element_blank(),legend.box.spacing = unit(-0.5,"cm"),legend.key = element_rect(size = 0),axis.ticks.y = element_line(colour = "black"),legend.key.size = unit(0.25, "cm"),legend.key.width =unit(1, "cm"))+scale_y_continuous(expand = c(0, 0),breaks = pretty_breaks())+xlab("")+scale_fill_manual(values = colorall,drop=FALSE)+guides(fill=guide_legend(nrow=legendnrow,byrow=TRUE,label.position = "bottom"))+ylab("Signature contribution\n")
  #+theme(plot.margin=margin(t=4,unit="pt"))
  
  
  if(!is.null(highlight)){
    samhigh <- sigdata %>% mutate(Samples_high=if_else(Samples %in% highlight,paste0("*",Samples),Samples))
    p4 <-  p4+scale_x_discrete(breaks=samhigh$Samples,labels=samhigh$Samples_high)
  }
  
  #p3 <- flush_ticks(p3,flush = "Y",plot = FALSE)
  #p4 <- flush_ticks(p4,flush = "Y",plot = FALSE)
  
  
  if(!is.null(studydata)){
    if(!is.null(puritydata)){
      pall <- ggarrange(p1,p2,p2.3,p3,p2.1,p4,p2.5,ncol=1,nrow = 7,align = 'v',heights = c(2,0.1,0.1,4,0.1,7,0.1))
    }else {
      pall <- ggarrange(p1,p2,p3,p2.1,p4,p2.5,ncol=1,nrow = 6,align = 'v',heights = c(2,0.1,4,0.1,7,0.1))
    }
  }else{
    if(!is.null(puritydata)){
      pall <- ggarrange(p1,p2.3,p3,p2.1,p4,p2.5,ncol=1,nrow = 6,align = 'v',heights = c(2,0.1,4,0.1,7,0.1))
    }
  }
  
  ggsave(filename = filename,plot = pall,height = height,width = width,device = cairo_pdf)
}




### calculate_similarities ####

# cos_sim <- function(a, b) {
#   if (sum(a) == 0 | sum(b) == 0) {
#     value = 0
#   }
#   dot_product = a %*% b
#   norm_a <- norm(as.matrix(a), "2")
#   norm_b <- norm(as.matrix(b), "2")
#   value = dot_product / (norm_a * norm_b)
#   return(value)
# }

## from MutationalPattern
cos_sim <- function (x, y) 
{
  res = x %*% y/(sqrt(x %*% x) * sqrt(y %*% y))
  res = as.numeric(res)
  return(res)
}

# add the Correlation
calculate_similarities <- function(orignal_genomes, signature, signature_activaties) {
  require(entropy)
  data2 <- as.data.frame(orignal_genomes)
  #data2 <- read.delim(orignal_genomes_file,header = T,check.names = F,stringsAsFactors = F)
  #data2 <- data.frame(data2)
  data3 <-  as.data.frame(signature)
  data4 <-  as.data.frame(signature_activaties)
  #colnames(data2[,2:ncol(data2)])
  
  if(dim(data2)[1]==0 | dim(data3)[1]==0 | dim(data4)[1]==0){
    return(NA)
  }else {
    data2 <- data2[,!is.na(colSums(data2 != 0)) & colSums(data2 != 0) > 0]
    ## filter out the samples and signatures
    sample_name_total <- data4$Sample
    sample_name_tmp<- colnames(data2)[-1]
    sample_name <- sample_name_total[sample_name_total %in% sample_name_tmp]
    
    signature_name_total <- colnames(data4)[-1]
    signature_name_tmp <- colnames(data3)[-1]
    signature_name <- signature_name_total[signature_name_total %in% signature_name_tmp]
    
    data2 <- data2 %>% select(MutationType,one_of(sample_name))
    data3 <- data3 %>% select(MutationType,one_of(signature_name))
    data4 <- data4 %>% select(Sample,one_of(signature_name)) %>% filter(Sample %in% sample_name)
    
    genomes <- data2[, 2:length(data2)]
    
    #data3 <- read.delim(signature_file,header = T,check.names = F,stringsAsFactors = F)
    
    data3 <- data3[,2:ncol(data3)]
    #data4 <- read.delim(signature_activaties_file,header = T,check.names = F,stringsAsFactors = F)
    
    data4 <- data4[,2:ncol(data4)]
    est_genomes <- as.data.frame(as.matrix(data3) %*% as.matrix(t(data4)))
    
    cosine_sim_list = c()
    correlation_list = c()
    kl_divergence_list = c()
    l1_norm_list = c()
    l2_norm_list = c()
    total_mutations = c()
    relative_l1_list = c()
    relative_l2_list = c()
    for (i in 1:ncol(genomes)) {
      p_i <- as.numeric(genomes[, i])
      q_i = (est_genomes[, i])
      if(sum(q_i)>0){
        cosine_sim_list = append(cosine_sim_list, round(cos_sim(p_i, q_i), digits=3))
        correlation_list = append(correlation_list, round(cor(p_i, q_i), digits=3))
        kl_divergence_list = append(kl_divergence_list, round(KL.empirical(p_i, q_i), digits=4))
        l1_norm_list = append(l1_norm_list, round(norm(as.matrix(p_i-q_i), "1"), digits=3))
        relative_l1_list = append(relative_l1_list, round((dplyr::last(l1_norm_list)/norm(as.matrix(p_i), "1"))*100, digits=3))
        l2_norm_list = append(l2_norm_list, round(norm(as.matrix(p_i-q_i), "2"), digits=3))
        relative_l2_list = append(relative_l2_list, round((dplyr::last(l2_norm_list)/norm(as.matrix(p_i), "2"))*100, digits=3))
        total_mutations = append(total_mutations, sum(p_i))
      }else{
        cosine_sim_list = append(cosine_sim_list, NA_real_)
        correlation_list = append(correlation_list, NA_real_)
        kl_divergence_list = append(kl_divergence_list,NA_real_)
        l1_norm_list = append(l1_norm_list, NA_real_)
        relative_l1_list = append(relative_l1_list, NA_real_)
        l2_norm_list = append(l2_norm_list, NA_real_)
        relative_l2_list = append(relative_l2_list, NA_real_)
        total_mutations = append(total_mutations, sum(p_i))
      }
      
      
    }
    kl_divergence_list[!is.na(kl_divergence_list) & !is.finite(kl_divergence_list)] <- 1000
    similarities_dataframe = data.frame("Sample_Names"=sample_name,
                                        "Total_Mutations"=total_mutations,
                                        "Cosine_similarity"=cosine_sim_list,
                                        "L1_Norm"=l1_norm_list,
                                        `L1_Norm_%`=relative_l1_list,
                                        `100-L1_Norm_%`=100-relative_l1_list,
                                        "L2_Norm"=l2_norm_list,
                                        `L2_Norm_%`=relative_l2_list,
                                        `100-L2_Norm_%`=100-relative_l2_list,
                                        "KL_Divergence"= kl_divergence_list,
                                        "Correlation" = correlation_list,
                                        check.names = F)
    #write.csv(similarities_dataframe, file="MyData.csv")
    #write.table(similarities_dataframe, file="MyData.txt", sep="\t", row.names = FALSE)
    
    similarities_dataframe <- tibble("Sample_Names"=sample_name_total) %>% left_join(similarities_dataframe)
    return(similarities_dataframe)
  }
}

#calculate_similarities("Test/orignal_genomes.txt", "Test/Wsignature_sigs.txt", "Test/Wsignature_activaties.txt")
L1_normal_displot <- function(similarities_dataframe,...){
  require(ggplot2)
  require(hrbrthemes)
  similarities_dataframe %>% ggplot(aes(100-`L1_Norm_%`))+geom_histogram(color="white",binwidth = 2)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 16,axis = "xy")+scale_x_continuous(name ="\n% mutations explained by extracted signatures (100-L1_Norm %)",expand = c(0,0),limits = c(0,100))+scale_y_continuous(name = "Number of sample",expand = c(0,0))
}




#Generation of probabilities for each processes given to A mutation type 
#probabilities('Decomposed_Solution_Signatures.txt','Decomposed_Solution_Activities.txt')
probabilities <-  function(W, H){ 
  require(dplyr)
  #W=read.delim(Wfile,header = T,check.names = F,stringsAsFactors = F)
  #H=read.delim(Hfile,header = T,check.names = F,stringsAsFactors = F)
  
  MutationType <- W$MutationType
  allcolnames <- H$Samples
  
  W <- as.matrix(W[,2:ncol(W)])
  H <- t(as.matrix(H[,2:ncol(H)]))
  genomes <- W %*% H
  
  probs_all <- NULL
  
  for(i in 1:dim(H)[2]){
    probs <-(W*(H[,i])[col(W)])/(genomes[,i])[row(W)]
    rownames(probs) <- MutationType
    probs <- as.data.frame(probs) %>% rownames_to_column(var = "MutationType") %>% mutate(Sample_Names=allcolnames[i]) %>% select(Sample_Names,MutationType,everything())
    probs_all <- bind_rows(probs_all,probs)
  }
  return(probs_all)
  
}



# SBS96_plot --------------------------------------------------------------

plot_sbs_96_profile <- function(data,samplename=NULL,totalmut=NULL,samplename_plot=TRUE, totalmut_plot=TRUE,filename=NULL,percentage=TRUE,ytitle=NULL){
  
  require(tidyverse)
  require(hrbrthemes)
  require(scales)
  require(ggpubr)  
  
  #### data is a frame with the following columns: Type, SubType, MutationType, Value ###
  
  if(dim(data)[2]==2){
    data <- data %>% mutate(Type=str_sub(MutationType,3,5),SubType=paste0(str_sub(MutationType,1,1),str_sub(MutationType,3,3),str_sub(MutationType,7,7))) %>% select(Type,SubType,MutationType,everything()) %>% arrange(Type,SubType)
  }
  
  if(is.null(samplename)){
    samplename <- colnames(data)[4]
  }
  
  colnames(data)[4] <- "Value"
  
  data <- data %>% mutate(Seq=seq_along(MutationType)) %>% mutate(Type=fct_inorder(Type),SubType=fct_inorder(SubType),MutationType=fct_inorder(MutationType))
  
  stype <- unique(data$Type)
  
  if(!is.null(totalmut)){
    totalmut <- comma_format()(totalmut)
  }
  
  COLORS6 = c(
    "#2EBAED", "#000000", "#DE1C14",
    "#D4D2D2", "#ADCC54", "#F0D0CE")
  names(COLORS6) <- stype
  
  data0 <- bind_cols(
    data %>% group_by(Type) %>% dplyr::slice(1) %>% ungroup() %>% mutate(Value1=Seq+0.25) %>% select(Value1,Type),
    data %>% group_by(Type) %>% dplyr::slice(16) %>% ungroup() %>% mutate(Value2=Seq+0.25) %>% select(Value2)
  )
  
  
  if(is.null(totalmut)){
    totalmut <- sum(data$Value)
  }  
  
  
  if(percentage){
    if(sum(data$Value)>10){
      data$Value <- data$Value/totalmut
    }
    ymax <- max(data$Value)/0.9
    ymax <- 0.04*ceiling(ymax/0.04)
    labelx <- percent
    #ylabp <- percent(pretty_breaks(n = 5)(ymax))
    ylabp <- percent(seq(0,ymax,length.out = 5))
    if(is.null(ytitle)) {
      ytext <- "Percentage of Single Base Subsitutions"
    }else{
      ytext <- ytitle
    }
    
    
  }else{
    ymax <- max(data$Value)/0.9
    ymax <- ceiling(ymax)
    labelx <- comma
    #ylabp <- comma(pretty_breaks(n = 5)(ymax))
    ylabp <- comma(seq(0,ymax,length.out = 5))
    if(is.null(ytitle)) {
      ytext <- "Number of Single Base Subsitutions"
    }else{
      ytext <- ytitle
    }
  }
  
  
  
  #ylabp <- percent(pretty_breaks(n = 5)(data$Value))
  ylabp <- ylabp[length(ylabp)]
  
  data00 <- data %>%  mutate(B1=str_sub(SubType,1,1),B2=str_sub(SubType,2,2),B3=str_sub(SubType,3,3))
  
  p1 <- data0 %>% ggplot(aes(xmin=Value1,xmax=Value2,ymin=0,ymax=1,fill=Type))+
    geom_rect()+
    scale_fill_manual(values=COLORS6)+
    scale_x_continuous(expand = c(0,0),sec.axis = dup_axis(labels = NULL,name = NULL),limits=c(1,96.5))+
    scale_y_continuous(expand = c(0,0),labels =ylabp,breaks = 1,limits = c(0,2.2),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y=" ")+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,color = "white"),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank(),
          #axis.line = element_blank(),
          axis.line.x = element_line(colour = 'white',size=0.6),
          axis.line.y = element_line(colour = 'white',size=0.6),
          plot.margin=margin(b=-1.2,l = 0.2,unit="cm")
    )+
    annotate("text", x = 8.75+16*0, y = 1.5, label = stype[1],size=5,fontface =2)+
    annotate("text", x = 8.75+16*1, y = 1.5, label = stype[2],size=5,fontface =2)+
    annotate("text", x = 8.75+16*2, y = 1.5, label = stype[3],size=5,fontface =2)+
    annotate("text", x = 8.75+16*3, y = 1.5, label = stype[4],size=5,fontface =2)+
    annotate("text", x = 8.75+16*4, y = 1.5, label = stype[5],size=5,fontface =2)+
    annotate("text", x = 8.75+16*5, y = 1.5, label = stype[6],size=5,fontface =2)
  
  
  
  p2 <- data %>% 
    ggplot(aes(Seq,Value,fill=Type))+
    geom_col(color="white",width = 0.5,size=0)+
    scale_fill_manual(values=COLORS6)+
    scale_x_continuous(expand = expansion(add = c(0.25,0.3)),labels = data$SubType,breaks = data$Seq,sec.axis = dup_axis(labels = NULL,name = NULL))+
    scale_y_continuous(expand = c(0,0),labels = labelx,breaks = seq(0,ymax,length.out = 5),limits = c(0,ymax),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y=ytext)+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10),
          #axis.text.x=element_text(angle = 90,hjust = 1,vjust = 0.5),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "gray90",size = 0.3,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "gray96",size = 0.3,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour = 'gray80',size=0.6),
          axis.line.y = element_line(colour = 'gray80',size=0.6),
          plot.margin=margin(t=-1.2,b=-1,l=0.2,unit="cm")
          
    )
  
  if(samplename_plot & !is.null(samplename)){
    p2 <- p2+annotate("text", x = 1.5, y = ymax*0.9, label = samplename,size=7,fontface =2,hjust = 0)
  }
  
  if(totalmut_plot & !is.null(totalmut)){
    p2 <- p2+annotate("text", x = 78, y = ymax*0.88, label = paste0("Total Mutations: ",totalmut),size=6.5,fontface =2,hjust = 0)
  }
  
  #invisible(capture.output(p2 = flush_ticks(p2)+theme(axis.text.x = element_blank())))
  p2 = flush_ticks(p2)+theme(axis.text.x = element_blank())
  
  p3 <- data00 %>% 
    ggplot(aes(Seq))+
    geom_text(aes(y=0.4,label=B1),col="gray60",angle=90,hjust = 0,vjust = 0.5,size=3.5)+
    geom_text(aes(y=1,label=B2,col=Type),angle=90,hjust = 0,vjust = 0.5,size=3.5,family = "Arial")+
    geom_text(aes(y=1.6,label=B3),col="gray60",angle=90,hjust = 0,vjust = 0.5,size=3.5)+
    scale_color_manual(values=COLORS6)+
    scale_x_continuous(expand = c(0,0),sec.axis = dup_axis(labels = NULL,name = NULL),limits = c(0.5,96.5))+
    scale_y_continuous(expand = c(0,0),labels =ylabp,breaks = 1,limits = c(0,2.2),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y=" ")+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,colour = "white"),
          #axis.text.x=element_text(angle = 90,hjust = 1,vjust = 0.5),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "white",size = 0.4,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "white",size = 0.4,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour = 'white',size=0.6),
          axis.line.y = element_line(colour = 'white',size=0.6),
          plot.margin=margin(t=-1,l=0.2,unit="cm")
    )
  #,text = element_text(family = "Arial")
  
  
  #cairo_pdf(file = 'tmp.pdf',width = 16,height = 4)
  pcomb <- ggarrange(p1,p2,p3,nrow = 3,align = "h",heights = c(1.8,10,1.5))
  #pcomb <- plot_grid(p1,p2,p3,nrow = 3,align = "v",axis = c("tb"),rel_heights = c(1.8,10,1.5))
  
  if(is.null(filename)){
    return(pcomb)
  } else {
    ggsave(filename,width = 18,height = 4.5)
  }
  
}



# PLOT_ID_83 -------------------------------------------------------------

plot_id_83_profile <- function(data,samplename=NULL,totalmut=NULL,samplename_plot=TRUE, totalmut_plot=TRUE,filename=NULL,percentage=TRUE,ytitle=NULL){
  
  require(tidyverse)
  require(hrbrthemes)
  require(scales)
  require(ggpubr)  
  
  if(is.null(samplename)){
    samplename <- colnames(data)[2]
  }
  
  
  colnames(data) <- c('MutationType','Value')
  
  data <- profile_format_df(data,factortype = T) 
  
  COLORS16=c("#FCBD6F", "#FE8002", "#AFDC8A", "#36A02E", "#FCC9B4", "#FB896A", "#F04432", "#BB191A", "#CFE0F1", "#93C3DE", "#4A97C8", "#1764AA", "#E1E1EE", "#B5B5D7", "#8582BC", "#62409A")
  names(COLORS16) <- as.character(unique(data$Type))
  
  ###
  if(percentage){
    if(sum(data$Value)>2){
      data$Value <- data$Value/totalmut
    }
    ymax <- max(data$Value)/0.9
    ymax <- ceiling(ymax*100/4)*0.04
    ymax <- if_else(ymax>1,1,ymax)
    
    if(is.null(ytitle)) { ytitle='Percentage of Indels' }
    
    totalmut_plot <-  FALSE
    labely <- percent_format()
    
  }else{
    ymax <- max(data$Value)/0.9
    ymax <- ceiling(ymax)
    
    if(is.null(ytitle)) { ytitle='Percentage of Indels' }
    labely <- comma_format()
  }
  
  
  #### 
  bdata <- data %>% mutate(seq=seq_along(MutationType)) %>% select(Type,SubType,seq) %>% group_by(Type) %>% summarise(xmin=min(seq)-0.2,xmax=max(seq)+0.2,rcol=NA_character_) 
  bdata$rcol = COLORS16[bdata$Type]
  
  tdata <- data %>% 
    mutate(Group=if_else(str_starts(Type,'1'),str_remove(Type,".*:"),str_remove(Type,":.*"))) %>% 
    mutate(seq=seq_along(MutationType)) %>% 
    group_by(Type) %>%
    summarise(Group=unique(Group),seq=(min(seq)+max(seq))/2) %>% 
    mutate(tcol=if_else(Group %in% c('T','5+'),'white','black'))
  
  if(is.null(totalmut)){
    totalmut <- sum(data$Value)
  }  
  
  if(!is.null(totalmut)){
    totalmut <- comma_format()(totalmut)
  }
  
  
  
  p <- data %>% 
    ggplot(aes(MutationType,Value,fill=Type))+
    geom_col(color="white",width = 0.4,size=0)+
    scale_fill_manual(values=COLORS16)+
    #scale_x_discrete(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0),breaks = seq(0,ymax,length.out = 5),labels = labely)+
    labs(x="",y=ytitle)+
    theme_ipsum_rc(axis = FALSE,grid = 'Y',axis_title_just = 'm',plot_margin=margin(60,5.5,40,5.5))+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10),
          axis.text.x=element_blank(),
          # axis.text.x=element_text(size=11,angle = 90,hjust = 1,vjust = 0.5,colour = COLORS10[data$Type]),
          #axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "gray90",size = 0.3,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "gray96",size = 0.3,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank()
    )+
    panel_border(color = 'gray20',size = 0.3)+
    annotate("rect", xmin = bdata$xmin,xmax = bdata$xmax, ymin = -0.05*ymax,ymax = -0.01*ymax,fill=bdata$rcol)+
    annotate("text", x=1:length(data$SubType),y=-0.05*ymax,label=data$SubType,hjust=0.5,vjust=1.2,size=3,fontface =2)+
    annotate("text", x=c(6.5,18.5,36.5,60.5,78),y = -0.1*ymax,label=c('Homopolymer Length','Homopolymer Length','Number of Repeat Units','Number of Repeat Units','Microhomology Length'),hjust=0.5,vjust=1.2,size=4.5,fontface =2)+
    #annotate("segment", x=c(6.5,18.5,36.5,60.5,78),xend=c(6.5,18.5,36.5,60.5,78),y = -0.1, yend=1)+
    annotate("rect", xmin = bdata$xmin,xmax = bdata$xmax, ymin = ymax+0.01*ymax,ymax = ymax+0.07*ymax,fill=bdata$rcol)+
    annotate("text", x=tdata$seq,y=ymax+0.04*ymax,label=tdata$Group,color=tdata$tcol,hjust=0.5,vjust=0.5,size=4,fontface =2)+
    annotate("text", x=c(6.5,18.5,36.5,60.5,78),y = ymax+0.1*ymax,label=c('1bp Deletion','1bp Insertion','>1bp Deletion at Repeats\n(Deletion Length)','>1bp Insertion at Repeats\n(Insertion Length)','Microhomology\n(Deletion Length)'),hjust=0.5,vjust=-0.2,size=4.5,fontface =2)+
    coord_cartesian(ylim = c(0,ymax),clip = 'off')
  
  
  if(samplename_plot & !is.null(samplename)){
    p <- p+annotate("text", x = 1.5, y = ymax*0.95, label = samplename,size=6,fontface =2,hjust = 0)
  }
  
  if(totalmut_plot & !is.null(totalmut)){
    p <- p+annotate("text", x =66, y = ymax*0.95, label = paste0("Total Indels: ",totalmut),size=6,fontface =2,hjust = 0)
  }
  
  if(is.null(filename)){
    return(p)
  } else {
    ggsave(filename,width = 14,height = 4,device = cairo_pdf)
  }
  
}

# PLOT_CN_68 -------------------------------------------------------------

plot_cn_68_profile <- function(data,samplename=NULL,totalmut=NULL,samplename_plot=TRUE, totalmut_plot=TRUE,filename=NULL,percentage=TRUE,ytitle=NULL){
  
  require(tidyverse)
  require(hrbrthemes)
  require(scales)
  require(ggpubr)  
  
  if(is.null(samplename)){
    samplename <- colnames(data)[2]
  }
  
  colnames(data) <- c('MutationType','Value')
  data <- profile_format_df(data,factortype = T) 
  
  COLORS3 <- c('#0000F2','#bababa','#3D4551')
  names(COLORS3) <-  c('homdel','LOH','het')
  
  COLORS10 <- c('#0000F2','gray20','#4daf4a','#984ea3','#947100','#BB0E3D','#4daf4a','#984ea3','#947100','#BB0E3D')
  names(COLORS10) <-  c("0:homdel","1:LOH","2:LOH","3-4:LOH","5-8:LOH","9+:LOH","2:het","3-4:het","5-8:het","9+:het")
  
  COLORS68 <- c(alpha('#0000F2',seq(0.1,1,length=5)),alpha('gray20',seq(0.1,1,length=7)),alpha('#4daf4a',seq(0.1,1,length=7)),alpha('#984ea3',seq(0.1,1,length=7)),alpha('#947100',seq(0.1,1,length=7)),alpha('#BB0E3D',seq(0.1,1,length=7)),alpha('#4daf4a',seq(0.1,1,length=7)),alpha('#984ea3',seq(0.1,1,length=7)),alpha('#947100',seq(0.1,1,length=7)),alpha('#BB0E3D',seq(0.1,1,length=7)))
  names(COLORS68) <- levels(data$MutationType) 
  
  ###
  if(percentage){
    if(sum(data$Value)>2){
      data$Value <- data$Value/totalmut
    }
    ymax <- max(data$Value)/0.9
    ymax <- ceiling(ymax*100/4)*0.04
    ymax <- if_else(ymax>1,1,ymax)
    
    if(is.null(ytitle)) { ytitle='Proportion' }
    
    totalmut_plot <-  FALSE
    labely <- percent_format()
    
  }else{
    ymax <- max(data$Value)/0.9
    ymax <- ceiling(ymax)
    
    if(is.null(ytitle)) { ytitle='Number of CNV Segements' }
    labely <- comma_format()
  }
  
  
  #### 
  tdata <- data %>% mutate(seq=seq_along(MutationType)) %>% select(Type,SubType,seq) %>% group_by(Type) %>% summarise(xmin=min(seq)-0.5,xmax=max(seq)+0.5,rcol=NA_character_) %>% mutate(label=str_remove(Type,":.*"))
  tdata$rcol = COLORS10[tdata$Type]
  
  if(is.null(totalmut)){
    totalmut <- sum(data$Value)
  }  
  
  if(!is.null(totalmut)){
    totalmut <- comma_format()(totalmut)
  }
  
  
  
  p <- data %>% 
    ggplot(aes(MutationType,Value,fill=MutationType))+
    geom_col(color="black",width = 0.6,size=0.1)+
    scale_fill_manual(values=COLORS68)+
    #scale_x_discrete(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0),breaks = seq(0,ymax,length.out = 5),labels = labely)+
    labs(x="",y=ytitle)+
    theme_ipsum_rc(axis = FALSE,grid = 'Y',axis_title_just = 'm',plot_margin=margin(40,5.5,50,5.5))+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10),
          axis.text.x=element_blank(),
          # axis.text.x=element_text(size=11,angle = 90,hjust = 1,vjust = 0.5,colour = COLORS10[data$Type]),
          #axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "gray90",size = 0.3,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "gray96",size = 0.3,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank()
    )+
    panel_border(color = 'gray20',size = 0.3)+
    annotate("text", x=1:length(data$SubType),y=-0.015*ymax,label=data$SubType,angle=90,size=3,fontface =2,hjust=1)+
    annotate("rect", xmin = tdata$xmin,xmax = tdata$xmax, ymin = ymax+0.01*ymax,ymax = ymax+0.09*ymax,fill=tdata$rcol,col='gray20',linewidth=0.2)+
    annotate("text", x=(tdata$xmin+tdata$xmax)/2,y = ymax+0.05*ymax,label=tdata$label,hjust=0.5,vjust=0.5,size=4.5,fontface =2,col='white')+
    annotate("rect", xmin = tdata$xmin[c(1,2,7)],xmax = tdata$xmax[c(1,6,10)], ymin = ymax+0.09*ymax,ymax = ymax+0.17*ymax,fill=COLORS3,col='gray20',linewidth=0.2)+
    annotate("text", x=c(3,23,54.5),y = ymax+0.13*ymax,label=c('HD','LOH','Het'),hjust=0.5,vjust=0.5,size=4.5,fontface =2,col=c('white','black','white'))+
    coord_cartesian(ylim = c(0,ymax),clip = 'off')
  
  
  if(samplename_plot & !is.null(samplename)){
    p <- p+annotate("text", x = 1.5, y = ymax*0.95, label = samplename,size=6,fontface =2,hjust = 0)
  }
  
  if(totalmut_plot & !is.null(totalmut)){
    p <- p+annotate("text", x =66, y = ymax*0.95, label = paste0("Total Indels: ",totalmut),size=6,fontface =2,hjust = 0)
  }
  
  if(is.null(filename)){
    return(p)
  } else {
    ggsave(filename,width = 12,height = 4,device = cairo_pdf)
  }
  
}


# PLOT_DBS_78 -------------------------------------------------------------

plot_dbs_78_profile <- function(data,samplename=NULL,samplename_plot=TRUE, filename=NULL,totalmut=NULL,totalmut_plot=TRUE,percentage=TRUE,ytitle=NULL){
  require(tidyverse)
  require(hrbrthemes)
  require(scales)
  require(cowplot)  
  
  data <- data %>% mutate(Type=paste0(str_sub(MutationType,1,3),"NN"),SubType=str_sub(MutationType,4,5)) %>% select(Type,SubType,MutationType,everything()) %>% arrange(Type,SubType)
  
  if(is.null(samplename)){
    samplename <- colnames(data)[4]
  }
  
  colnames(data)[4] <- "Value"
  data <- data %>% mutate(Seq=seq_along(MutationType)) %>% mutate(Type=fct_inorder(Type),SubType=fct_inorder(SubType),MutationType=fct_inorder(MutationType))
  stype <- unique(data$Type)
  
  if(!is.null(totalmut)){
    totalmut <- comma_format()(totalmut)
  }
  
  if(percentage){
    if(sum(data$Value)>10){
      totalmut <- sum(data$Value)
      data$Value <- data$Value/totalmut
    }
    
    ymax <- max(data$Value)/0.9
    ymax <- ceiling(ymax*100/4)*0.04
    ymax <- if_else(ymax>1,1,ymax)
    
    labely <- percent_format()
    
    if(is.null(ytitle)) {
      ytext <- "Percentage of Double Base Subsitutions"
    }else{
      ytext <- ytitle
    }
  }else{
    if(is.null(ytitle)) {
      ytext <- "Number of Double Base Subsitutions"
    }else{
      ytext <- ytitle
    }  
    labely <- comma_format()
    
  }
  
  print(ymax)
  
  COLORS10 = c("#03BCEE", "#0366CB", "#A1CE63", "#016601", "#FE9898", "#E32926", "#FEB166", "#FE8001", "#CB98FE", "#4C0198")
  names(COLORS10) <- stype
  
  # ymax <- max(data$Value)/0.9
  # ymax <- 4*ceiling(ymax/4)
  
  #ylabp <- pretty_breaks(n = 5)(data$Value)
  ylabp <- comma(seq(0,ymax,length.out = 5))
  ylabp <- ylabp[length(ylabp)]
  
  data0 <- bind_cols(
    data %>% group_by(Type) %>% filter(row_number() == 1) %>% ungroup() %>% mutate(Value1=Seq+0.25) %>% select(Value1,Type),
    data %>% group_by(Type) %>% "["(.,c(which(data$Type != lag(data$Type))-1, 78),) %>% ungroup() %>% mutate(Value2=Seq+0.25) %>% select(Value2)
  )
  
  #ylabp <- pretty_breaks(n = 5)(data$Value)
  #ylabp <- ylabp[length(ylabp)]
  
  anno.x=(data0$Value1+data0$Value2)/2
  anno.lab=data0$Type
  
  p1 <- 
    data0 %>% ggplot(aes(xmin=Value1,xmax=Value2,ymin=0,ymax=0.9,fill=Type))+
    geom_rect()+
    scale_fill_manual(values=COLORS10)+
    scale_x_continuous(expand = c(0,0),limits=c(1,78.5))+
    scale_y_continuous(expand = c(0,0),breaks = 1,labels = ylabp,limits = c(0,2.2),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(y=" ",x=NULL)+annotate("text",x=anno.x,y=1.5,label=anno.lab,size=5,fontface =2,hjust=0.5)+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,color = "white"),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank(),
          #axis.line.x = element_line(colour = 'white',size=0.6),
          #axis.line.y = element_line(colour = 'white',size=0.6),
          axis.line = element_blank(),
          plot.margin=margin(b=0,l = 0.2,unit="cm"),
    )
  #+theme(plot.background = element_rect(fill = "darkblue"))
  
  p2 <- data %>% 
    ggplot(aes(Seq,Value,fill=Type))+
    geom_col(color="white",width = 0.4,size=0)+
    scale_fill_manual(values=COLORS10)+
    scale_x_continuous(expand = c(0,0),labels = data$SubType,breaks = data$Seq,sec.axis = dup_axis(labels = NULL,name = NULL))+
    scale_y_continuous(expand = c(0,0),breaks = seq(0,ymax,length.out = 5),labels = labely,limits = c(0,ymax),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y=ytext)+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=11,angle = 90,hjust = 1,vjust = 0.5,colour = COLORS10[data$Type]),
          #axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "gray90",size = 0.3,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "gray96",size = 0.3,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour = 'gray80',size=0.6),
          axis.line.y = element_line(colour = 'gray80',size=0.6),
          plot.margin=margin(t=0,l=0.2,unit="cm")
    )
  
  # if(samplename_plot & !is.null(samplename)){
  #   p2 <- p2+annotate("text", x = 1, y = ymax*0.9, label = paste0(samplename,": ",totalmut," double subs"),size=7,fontface =2,hjust = 0)
  # }
  # 
  
  if(samplename_plot & !is.null(samplename)){
    p2 <- p2+annotate("text", x = 1.5, y = ymax*0.9, label = samplename,size=7,fontface =2,hjust = 0)
  }
  
  if(totalmut_plot & !is.null(totalmut)){
    p2 <- p2+annotate("text", x = 78, y = ymax*0.88, label = paste0("Total Mutations: ",totalmut),size=6.5,fontface =2,hjust = 0)
  }
  
  
  #invisible(capture.output(p2 <- flush_ticks(p2,flush = "X")))
  p2 <- p2+ theme(axis.text.y=element_text(vjust=c(0, rep(0.5, 3), 1)))
  
  
  #cairo_pdf(file = 'tmp.pdf',width = 16,height = 4)
  #pcomb <- ggarrange(p1,p2,nrow = 2,align = "h",heights = c(3,15))
  #pcomb <- grid.arrange(p1,p2,nrow=2)
  pcomb <- plot_grid(p1, p2, align = "v", nrow = 2, rel_heights = c(1/6,5/6))
  
  if(is.null(filename)){
    return(pcomb)
  } else {
    ggsave(filename,pcomb,width = 18,height = 4,device = cairo_pdf)
  }
  
}
plot_dbs_78_profile_old <- function(data,samplename=NULL,samplename_plot=TRUE, filename=NULL){
  require(tidyverse)
  require(hrbrthemes)
  require(scales)
  require(cowplot)  
  
  data <- data %>% mutate(Type=paste0(str_sub(MutationType,1,3),"NN"),SubType=str_sub(MutationType,4,5)) %>% select(Type,SubType,MutationType,everything()) %>% arrange(Type,SubType)
  
  if(is.null(samplename)){
    samplename <- colnames(data)[4]
  }
  
  colnames(data)[4] <- "Value"
  data <- data %>% mutate(Seq=seq_along(MutationType)) %>% mutate(Type=fct_inorder(Type),SubType=fct_inorder(SubType),MutationType=fct_inorder(MutationType))
  stype <- unique(data$Type)
  
  totalmut <- sum(data$Value)
  totalmut <- comma_format()(totalmut)
  
  COLORS10 = c("#03BCEE", "#0366CB", "#A1CE63", "#016601", "#FE9898", "#E32926", "#FEB166", "#FE8001", "#CB98FE", "#4C0198")
  names(COLORS10) <- stype
  
  ymax <- max(data$Value)/0.9
  ymax <- 4*ceiling(ymax/4)
  
  #ylabp <- pretty_breaks(n = 5)(data$Value)
  ylabp <- comma(seq(0,ymax,length.out = 5))
  ylabp <- ylabp[length(ylabp)]
  
  data0 <- bind_cols(
    data %>% group_by(Type) %>% filter(row_number() == 1) %>% ungroup() %>% mutate(Value1=Seq+0.25) %>% select(Value1,Type),
    data %>% group_by(Type) %>% "["(.,c(which(data$Type != lag(data$Type))-1, 78),) %>% ungroup() %>% mutate(Value2=Seq+0.25) %>% select(Value2)
  )
  
  #ylabp <- pretty_breaks(n = 5)(data$Value)
  #ylabp <- ylabp[length(ylabp)]
  
  anno.x=(data0$Value1+data0$Value2)/2
  anno.lab=data0$Type
  
  p1 <- 
    data0 %>% ggplot(aes(xmin=Value1,xmax=Value2,ymin=0,ymax=0.9,fill=Type))+
    geom_rect()+
    scale_fill_manual(values=COLORS10)+
    scale_x_continuous(expand = c(0,0),limits=c(1,78.5))+
    scale_y_continuous(expand = c(0,0),breaks = 1,labels = ylabp,limits = c(0,2.2),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(y=" ",x=NULL)+annotate("text",x=anno.x,y=1.5,label=anno.lab,size=5,fontface =2,hjust=0.5)+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,color = "white"),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank(),
          #axis.line.x = element_line(colour = 'white',size=0.6),
          #axis.line.y = element_line(colour = 'white',size=0.6),
          axis.line = element_blank(),
          plot.margin=margin(b=0,l = 0.2,unit="cm"),
    )
  #+theme(plot.background = element_rect(fill = "darkblue"))
  
  p2 <- data %>% 
    ggplot(aes(Seq,Value,fill=Type))+
    geom_col(color="white",width = 0.4,size=0)+
    scale_fill_manual(values=COLORS10)+
    scale_x_continuous(expand = c(0,0),labels = data$SubType,breaks = data$Seq,sec.axis = dup_axis(labels = NULL,name = NULL))+
    scale_y_continuous(expand = c(0,0),breaks = seq(0,ymax,length.out = 5),limits = c(0,ymax),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y="Number of Double Base Subsitutions")+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=11,angle = 90,hjust = 1,vjust = 0.5,colour = COLORS10[data$Type]),
          #axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "gray90",size = 0.3,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "gray96",size = 0.3,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour = 'gray80',size=0.6),
          axis.line.y = element_line(colour = 'gray80',size=0.6),
          plot.margin=margin(t=0,l=0.2,unit="cm")
    )
  
  if(samplename_plot & !is.null(samplename)){
    p2 <- p2+annotate("text", x = 1, y = ymax*0.9, label = paste0(samplename,": ",totalmut," double subs"),size=7,fontface =2,hjust = 0)
  }
  
  #invisible(capture.output(p2 <- flush_ticks(p2,flush = "X")))
  p2 <- p2+ theme(axis.text.y=element_text(vjust=c(0, rep(0.5, 3), 1)))
  
  
  #cairo_pdf(file = 'tmp.pdf',width = 16,height = 4)
  #pcomb <- ggarrange(p1,p2,nrow = 2,align = "h",heights = c(3,15))
  #pcomb <- grid.arrange(p1,p2,nrow=2)
  pcomb <- plot_grid(p1, p2, align = "v", nrow = 2, rel_heights = c(1/6,5/6))
  
  if(is.null(filename)){
    return(pcomb)
  } else {
    ggsave(filename,pcomb,width = 18,height = 4,device = cairo_pdf)
  }
  
}
#plot_dbs_78_profile(data = data_dbs_exm1,filename = "tmp.pdf")


#save(data_sbs_exm1,data_dbs_exm1,file='../Example.RData')



# PLOT_SV_38 -------------------------------------------------------------

plot_sv_38_profile <- function(data,samplename=NULL,totalmut=NULL,samplename_plot=TRUE, totalmut_plot=TRUE,filename=NULL,percentage=TRUE,ytitle=NULL){
  
  require(tidyverse)
  require(hrbrthemes)
  require(scales)
  require(ggpubr)  
  
  if(is.null(samplename)){
    samplename <- colnames(data)[2]
  }
  
  colnames(data) <- c('MutationType','Value')
  data <- profile_format_df(data,factortype = T) 
  
  COLORS2 <- c('#cccccc','#3D4551')
  names(COLORS2) <-  c('Clustered','Non-Clustered')
  
  
  COLORS8 <- c('#BB0E3D','#ff7f00','#984ea3','#71767A','#BB0E3D','#ff7f00','#984ea3','#71767A')
  names(COLORS8) <-  c("DEL:clustered",    "DUP:clustered",    "INV:clustered",    "TLOC:clustered",   "DEL:unclustered",  "DUP:unclustered",  "INV:unclustered",  "TLOC:unclustered")
  
  COLORS38 <- c(alpha('#BB0E3D',seq(0.1,1,length=6)),alpha('#ff7f00',seq(0.1,1,length=6)),alpha('#984ea3',seq(0.1,1,length=6)),'#71767A',alpha('#BB0E3D',seq(0.1,1,length=6)),alpha('#ff7f00',seq(0.1,1,length=6)),alpha('#984ea3',seq(0.1,1,length=6)),'#71767A')
  names(COLORS38) <- levels(data$MutationType) 
  
  ###
  if(percentage){
    if(sum(data$Value)>2){
      data$Value <- data$Value/totalmut
    }
    ymax <- max(data$Value)/0.9
    ymax <- ceiling(ymax*100/4)*0.04
    ymax <- if_else(ymax>1,1,ymax)
    
    if(is.null(ytitle)) { ytitle='Percentage (%)' }
    
    totalmut_plot <-  FALSE
    labely <- percent_format()
    
  }else{
    ymax <- max(data$Value)/0.9
    ymax <- ceiling(ymax)
    
    if(is.null(ytitle)) { ytitle='Number of SV events' }
    labely <- comma_format()
  }
  
  
  #### 
  tdata <- data %>% mutate(seq=seq_along(MutationType)) %>% select(Type,SubType,seq) %>% group_by(Type) %>% summarise(xmin=min(seq)-0.5,xmax=max(seq)+0.5,rcol=NA_character_) %>% mutate(label=str_remove(Type,":.*")) %>% mutate(label=str_replace(label,'TLOC','T'))
  tdata$rcol = COLORS8[tdata$Type]
  
  if(is.null(totalmut)){
    totalmut <- sum(data$Value)
  }  
  
  if(!is.null(totalmut)){
    totalmut <- comma_format()(totalmut)
  }
  
  
  p <- data %>% 
    ggplot(aes(MutationType,Value,fill=MutationType))+
    geom_col(color="black",width = 0.6,size=0.1)+
    scale_fill_manual(values=COLORS38)+
    #scale_x_discrete(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0),breaks = seq(0,ymax,length.out = 5),labels = labely)+
    labs(x="",y=ytitle)+
    theme_ipsum_rc(axis = FALSE,grid = 'Y',axis_title_just = 'm',plot_margin=margin(40,5.5,50,5.5))+
    theme(axis.title.y=element_text(size=12,face = "bold"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10),
          axis.text.x=element_blank(),
          # axis.text.x=element_text(size=11,angle = 90,hjust = 1,vjust = 0.5,colour = COLORS10[data$Type]),
          #axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "gray90",size = 0.3,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "gray96",size = 0.3,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank()
    )+
    panel_border(color = 'gray20',size = 0.3)+
    annotate("text", x=1:length(data$SubType),y=-0.015*ymax,label=data$SubType,angle=90,size=3,fontface =2,hjust=1,vjust=0.5)+
    annotate("rect", xmin = tdata$xmin,xmax = tdata$xmax, ymin = ymax+0.01*ymax,ymax = ymax+0.09*ymax,fill=tdata$rcol,col='gray20',linewidth=0.2)+
    annotate("text", x=(tdata$xmin+tdata$xmax)/2,y = ymax+0.05*ymax,label=tdata$label,hjust=0.5,vjust=0.5,size=4.5,fontface =2,col='white')+
    annotate("rect", xmin = tdata$xmin[c(1,5)],xmax = tdata$xmax[c(4,8)], ymin = ymax+0.09*ymax,ymax = ymax+0.17*ymax,fill=COLORS2,col='gray20',linewidth=0.2)+
    annotate("text", x=c(10,29),y = ymax+0.13*ymax,label=c('Clustered','Non-Clustered'),hjust=0.5,vjust=0.5,size=4.5,fontface =2,col=c('black','white'))+
    coord_cartesian(ylim = c(0,ymax),clip = 'off')
  
  
  if(samplename_plot & !is.null(samplename)){
    p <- p+annotate("text", x = 1.5, y = ymax*0.95, label = samplename,size=6,fontface =2,hjust = 0)
  }
  
  if(totalmut_plot & !is.null(totalmut)){
    p <- p+annotate("text", x =66, y = ymax*0.95, label = paste0("Total Indels: ",totalmut),size=6,fontface =2,hjust = 0)
  }
  
  if(is.null(filename)){
    return(p)
  } else {
    ggsave(filename,width = 8,height = 4,device = cairo_pdf)
  }
  
}







# Plot_Profile_Logo -------------------------------------------------------
plot_profile_logo <- function (profile, colors = NULL, condensed = FALSE,output_plot = NULL,plot_width=NULL, plot_height=NULL) 
{
  ## profile 1 and profile 2 will be the dataframe with two columns: MutationType and value
  colnames(profile) <- c('MutationType','Value')
  
  COLORS2 = c('#1B4564','#D03D32')
  COLORS6 = c("#03BCEE", "#010101", "#E32926", "#CAC9C9", "#A1CE63", "#EBC6C4")
  COLORS10 = c("#03BCEE", "#0366CB", "#A1CE63", "#016601", "#FE9898", "#E32926", "#FEB166", "#FE8001", "#CB98FE", "#4C0198")
  COLORS16=c("#FCBD6F", "#FE8002", "#AFDC8A", "#36A02E", "#FCC9B4", "#FB896A", "#F04432", "#BB191A", "#CFE0F1", "#93C3DE", "#4A97C8", "#1764AA", "#E1E1EE", "#B5B5D7", "#8582BC", "#62409A")
  
  typelength = dim(profile)[1]
  if (is.null(colors)) {
    if(typelength == 192){colors = COLORS2; }
    if(typelength == 96){colors = COLORS6; }
    if(typelength == 78){colors = COLORS10;}
    if(typelength == 83){colors = COLORS16;}
  }
  
  ## make sure the order profile 1 and profile 2 are the same order
  
  indel_short <-  FALSE
  if(typelength == 83) {indel_short = TRUE}
  profile <- profile_format_df(profile,factortype = TRUE,indel_short = indel_short) 
  
  if(typelength == 192){
    names(colors) <- levels(profile$Strand) 
    profile[,5] <- profile[,5]/sum(profile[,5])  
    plot <- ggplot(data = profile, aes(x = SubType, y = Value,  fill = Strand, width = 0.5)) + 
      geom_bar(stat = "identity", position = "dodge", colour = NA, size = 0) + 
      scale_fill_manual(values = colors) + 
      facet_grid(. ~ Type, scales = "free",space = "free_x") +
      guides(fill = "none") + 
      scale_x_discrete(expand = expansion(add = 0))+
      scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
      theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),strip.text.x = element_blank(),panel.spacing.x = unit(0, "lines"),axis.line.x = element_line(colour = 'black',size = 0.05))
    
    
  }else{
    names(colors) <- levels(profile$Type)
    profile[,4] <- profile[,4]/sum(profile[,4])  
    plot <- ggplot(data = profile, aes(x = SubType, y = Value,  fill = Type, width = 0.5)) + 
      geom_bar(stat = "identity", position = "identity", colour = NA, size = 0) + 
      scale_fill_manual(values = colors) + 
      facet_grid(. ~ Type, scales = "free",space = "free_x") +
      guides(fill = "none") + 
      scale_x_discrete(expand = expansion(add = 0))+
      scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
      theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),strip.text.x = element_blank(),panel.spacing.x = unit(0, "lines"),axis.line.x = element_line(colour = 'black',size = 0.05))
  }
  
  
  #print(colors)
  
  
  
  #element_text(size = 0,margin = margin(1.5,0,1.5,0))
  
  # ## add background color for strip
  # require(grid)
  # g <- ggplot_gtable(ggplot_build(plot))
  # strip_top <- which(grepl('strip-t', g$layout$name))
  # fills <- colors
  # k <- 1
  # for (i in strip_top) {
  #   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  #   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  #   k <- k+1
  # }
  # plot <- as_ggplot(g)
  # 
  if(is.null(output_plot)){
    return(plot)
  }else{
    xleng <- 2
    yleng <- 0.8
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = plot,width = plot_width,height = plot_height)
  }
  
}



# ID29_Plot ---------------------------------------------------------------
id29_col1 <- c(
  "[+C]A",
  "[+C]G",
  "[+C]T",
  "[+C]C",
  "[+C]CC",
  "[+C]LR",
  "[+T]A",
  "[+T]C",
  "[+T]G",
  "[+T]T",
  "[+T]TT",
  "[+T]LR",
  "[+>1]NonR",
  "[+>1]Rep",
  "[-C]A",
  "[-C]G",
  "[-C]T",
  "[-C]C",
  "[-C]CC",
  "[-C]LR",
  "[-T]A",
  "[-T]C",
  "[-T]G",
  "[-T]T",
  "[-T]TT",
  "[-T]LR",
  "[->1]NonR",
  "[->1]Rep",
  "[-]MH")

id29_col2 <- c(
  "[+C]A",
  "[+C]G",
  "[+C]T",
  "[+C]C",
  "[+C]CC",
  "[+C]LR",
  "[+T]A",
  "[+T]C",
  "[+T]G",
  "[+T]T",
  "[+T]TT",
  "[+T]LR",
  "[+>1]NonR",
  "[+>1]Rep",
  "[-C]A",
  "[-C]G",
  "[-C]T",
  "[-C]C",
  "[-C]CC",
  "[-C]LR",
  "[-T]A",
  "[-T]C",
  "[-T]G",
  "[-T]T",
  "[-T]TT",
  "[-T]LR",
  "[->1]NonR",
  "[->1]Rep",
  "[-]MH")

id29_colors <- c("[+C]"="#0C5CA2","[+T]"="#CA4904","[+>1]"="#BF6196","[-C]"="#47A4E3","[-T]"="#DF8E05","[->1]"="#139060","[-]"="#7C00A5")
id29 <- tibble(MutationType1=id29_col1,MutationType2=id29_col2) %>% 
  mutate(Type=str_replace(MutationType1,"].*","]"),SubType=str_remove(MutationType1,".*]")) %>% 
  mutate(Type=factor(Type,levels = names(id29_colors))) %>% 
  mutate(MutationType2=factor(MutationType2,levels = id29_col2))


plot_profile_id29 <- function (profile, colors = id29_colors,output_plot = NULL,plot_width=NULL, plot_height=NULL) 
{
  ## make sure the order profile 1 and profile 2 are the same order
  
  ## normalize 
  colnames(profile) <- c('MutationType1','Value')
  profile <- profile %>% left_join(id29) %>% select(MutationType2,Type,SubType,Value)
  
  profile[,4] <- profile[,4]/sum(profile[,4])  
  
  plot <- ggplot(data = profile, aes(x = MutationType2, y = Value,  fill = Type, width = 0.7)) + 
    geom_bar(stat = "identity", position = "identity", colour = NA, size = 0) + 
    scale_fill_manual(values = colors) + 
    facet_grid(. ~ Type, scales = "free",space = "free_x") +
    guides(fill = FALSE) + 
    scale_x_discrete(expand = expansion(add = 0))+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = pretty_breaks(5))+
    labs(x="",y="")+
    theme_ipsum_rc(grid = 'yY',axis_text_size = 10)+
    theme(axis.text.x = element_text(angle = 330,hjust = 0,vjust = 1),strip.text.x = element_text(hjust = 0.5,face = "bold",colour = "white"),panel.spacing.x = unit(0.2, "lines"),panel.spacing.y = unit(0.2, "lines"),axis.line.x.bottom = element_line(colour = 'black',size = 0.2),strip.background = element_rect(size = 0))
  
  ## add background color for strip
  require(grid)
  g <- ggplot_gtable(ggplot_build(plot))
  strip_top <- which(grepl('strip-t', g$layout$name))
  fills <- colors
  k <- 1
  for (i in strip_top) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  plot <- as_ggplot(g)
  
  if(is.null(output_plot)){
    return(plot)
  }else{
    xleng <- 14
    yleng <- 5
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = plot,width = plot_width,height = plot_height)
  }
  
}

plot_profile_id29_logo <- function (profile, colors = id29_colors,output_plot = NULL,plot_width=NULL, plot_height=NULL) 
{
  ## make sure the order profile 1 and profile 2 are the same order
  
  ## normalize 
  colnames(profile) <- c('MutationType1','Value')
  profile <- profile %>% left_join(id29) %>% select(MutationType2,Type,SubType,Value)
  
  profile[,4] <- profile[,4]/sum(profile[,4])  
  
  plot <- ggplot(data = profile, aes(x = MutationType2, y = Value,  fill = Type, width = 0.7)) + 
    geom_bar(stat = "identity", position = "identity", colour = NA, size = 0) + 
    scale_fill_manual(values = colors) + 
    facet_grid(. ~ Type, scales = "free",space = "free_x") +
    guides(fill = FALSE) + 
    scale_x_discrete(expand = expansion(add = 0))+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = pretty_breaks(5))+
    labs(x="",y="")+
    theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),strip.text.x = element_blank(),panel.spacing.x = unit(0, "lines"),axis.line.x = element_line(colour = 'black',size = 0.05))
  
  ## add background color for strip
  # require(grid)
  # g <- ggplot_gtable(ggplot_build(plot))
  # strip_top <- which(grepl('strip-t', g$layout$name))
  # fills <- colors
  # k <- 1
  # for (i in strip_top) {
  #   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  #   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  #   k <- k+1
  # }
  # plot <- as_ggplot(g)
  
  if(is.null(output_plot)){
    return(plot)
  }else{
    xleng <- 2
    yleng <- 0.8
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = plot,width = plot_width,height = plot_height)
  }
  
}





# CN48 logo ---------------------------------------------------------------

cn48_mutationtype <- c(
  '0:homdel:0-100kb',
  '0:homdel:100kb-1Mb',
  '0:homdel:>1Mb',
  '1:LOH:0-100kb',
  '1:LOH:100kb-1Mb',
  '1:LOH:1Mb-10Mb',
  '1:LOH:10Mb-40Mb',
  '1:LOH:>40Mb',
  '2:LOH:0-100kb',
  '2:LOH:100kb-1Mb',
  '2:LOH:1Mb-10Mb',
  '2:LOH:10Mb-40Mb',
  '2:LOH:>40Mb',
  '3-4:LOH:0-100kb',
  '3-4:LOH:100kb-1Mb',
  '3-4:LOH:1Mb-10Mb',
  '3-4:LOH:10Mb-40Mb',
  '3-4:LOH:>40Mb',
  '5-8:LOH:0-100kb',
  '5-8:LOH:100kb-1Mb',
  '5-8:LOH:1Mb-10Mb',
  '5-8:LOH:10Mb-40Mb',
  '5-8:LOH:>40Mb',
  '9+:LOH:0-100kb',
  '9+:LOH:100kb-1Mb',
  '9+:LOH:1Mb-10Mb',
  '9+:LOH:10Mb-40Mb',
  '9+:LOH:>40Mb',
  '2:het:0-100kb',
  '2:het:100kb-1Mb',
  '2:het:1Mb-10Mb',
  '2:het:10Mb-40Mb',
  '2:het:>40Mb',
  '3-4:het:0-100kb',
  '3-4:het:100kb-1Mb',
  '3-4:het:1Mb-10Mb',
  '3-4:het:10Mb-40Mb',
  '3-4:het:>40Mb',
  '5-8:het:0-100kb',
  '5-8:het:100kb-1Mb',
  '5-8:het:1Mb-10Mb',
  '5-8:het:10Mb-40Mb',
  '5-8:het:>40Mb',
  '9+:het:0-100kb',
  '9+:het:100kb-1Mb',
  '9+:het:1Mb-10Mb',
  '9+:het:10Mb-40Mb',
  '9+:het:>40Mb'
)


cn48_color_mapping = c('0:0-100kb'='#F0F8FF', '0:100kb-1Mb'='#787CE6', '0:>1Mb'='#0000CD', 
                       '1:0-100kb'='#EBEBEB', '1:100kb-1Mb'='#C5C5C5', '1:1Mb-10Mb'='#9F9F9F', '1:10Mb-40Mb'='#797979', 
                       '1:>40Mb'='#545454', '2:0-100kb'='#F5FFFA', '2:100kb-1Mb'='#C0E2C3', 
                       '2:1Mb-10Mb'='#8BC48E', '2:10Mb-40Mb'='#56A858', '2:>40Mb'='#228B22', 
                       '3-4:0-100kb'='#FFF0F5', '3-4:100kb-1Mb'='#DEBDEB', '3-4:1Mb-10Mb'='#BE8BE1', 
                       '3-4:10Mb-40Mb'='#9D58D7', '3-4:>40Mb'='#7D26CD', '5-8:0-100kb'='#FFFAF0', 
                       '5-8:100kb-1Mb'='#F2DCB3', '5-8:1Mb-10Mb'='#E6BF78', '5-8:10Mb-40Mb'='#D9A23C', 
                       '5-8:>40Mb'='#CD8500', '9+:0-100kb'='#FFE4E1', '9+:100kb-1Mb'='#E2ADBC', 
                       '9+:1Mb-10Mb'='#C47798', '9+:10Mb-40Mb'='#A84074', '9+:>40Mb'='#8B0A50')
#colors = ['#0000CD', '#545454', '#228B22', '#7D26CD','#CD8500', '#8B0A50']

cn48_mutationtype
cn48_signature <- tibble(MutationType=cn48_mutationtype) %>% mutate(tmp=MutationType) %>% separate(col = tmp,into = c('a','b','c'),sep = ':') %>% mutate(key=paste0(a,":",c)) %>% select(MutationType,Class=b,Type=a,SubType=c,key) %>% mutate(color=NA_character_)
cn48_signature$color <- cn48_color_mapping[cn48_signature$key]
cn48_signature <- cn48_signature %>% select(-key) 

cn48_signature_color <- cn48_signature$color
names(cn48_signature_color) <- cn48_signature$MutationType

plot_profile_logo_cn48 <- function (profile, output_plot = NULL,plot_width=NULL, plot_height=NULL) 
{
  ## profile 1 and profile 2 will be the dataframe with two columns: MutationType and value
  colnames(profile) <- c('MutationType','Value')
  
  profile <- cn48_signature %>% left_join(profile) %>% mutate(MutationType = fct_inorder(MutationType))
  
  #profile[,2] <- profile[,2]/sum(profile[,2])  
  plot <- 
    ggplot(data = profile, aes(x = MutationType, y = Value,  fill = MutationType, width = 0.7)) + 
    geom_bar(stat = "identity", position = "identity", colour = 'black', size = 0.05) + 
    scale_fill_manual(values = cn48_signature_color) + 
    #facet_grid(. ~ Type, scales = "free",space = "free_x") +
    guides(fill = "none") + 
    scale_x_discrete(expand = expansion(add = 0))+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),strip.text.x = element_blank(),panel.spacing.x = unit(0, "lines"),axis.line.x = element_line(colour = 'black',size = 0.05))
  
  # 
  if(is.null(output_plot)){
    return(plot)
  }else{
    xleng <- 2
    yleng <- 0.8
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = plot,width = plot_width,height = plot_height)
  }
  
}




# DBS_signal_plot ---------------------------------------------------------
# dbs78_signal_colors <- c('AA>NN'='#E4446E','AC>NN'='#63C26F','AG>NN'='#FEE74C','AT>NN'='#4195CF','CA>NN'='#F3935C','CC>NN'='#9F47BE','CG>NN'='#6FF3F3','GA>NN'='#EF50E8','GC>NN'='#D9F964','TA>NN'='#F8C6C6')

dbs78_signal_colors <- c('AC>NN'="#03BCEE",'AT>NN'="#0366CB",'CC>NN'="#A1CE63",'CG>NN'="#016601",'CT>NN'="#FE9898",'GC>NN'="#E32926",'TA>NN'="#FEB166",'TC>NN'="#FE8001",'TG>NN'="#CB98FE",'TT>NN'="#4C0198")

plot_profile_dbs78_signal <- function (profile, colors = dbs78_signal_colors,output_plot = NULL,plot_width=NULL, plot_height=NULL) 
{
  ## make sure the order profile 1 and profile 2 are the same order
  
  ## normalize 
  colnames(profile) <- c('MutationType','Value')
  profile <- profile_format_df(profile,factortype = TRUE,indel_short = indel_short) 
  names(colors) <- levels(profile$Type)
  
  
  profile[,4] <- profile[,4]/sum(profile[,4])  
  
  plot <- ggplot(data = profile, aes(x = MutationType, y = Value,  fill = Type, width = 0.7)) + 
    geom_bar(stat = "identity", position = "identity", colour = NA, size = 0) + 
    scale_fill_manual(values = colors) + 
    facet_grid(. ~ Type, scales = "free",space = "free_x") +
    guides(fill = FALSE) + 
    scale_x_discrete(expand = expansion(add = 0))+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = pretty_breaks(5))+
    labs(x="",y="")+
    theme_ipsum_rc(grid = 'yY',axis_text_size = 10)+
    theme(axis.text.x = element_text(angle = 270,hjust = 0,vjust = 0.5),strip.text.x = element_text(hjust = 0.5,face = "bold"),panel.spacing.x = unit(0.2, "lines"),panel.spacing.y = unit(0.2, "lines"),axis.line.x.bottom = element_line(colour = 'black',size = 0.2),strip.background = element_rect(size = 0))
  
  ## add background color for strip
  require(grid)
  g <- ggplot_gtable(ggplot_build(plot))
  strip_top <- which(grepl('strip-t', g$layout$name))
  fills <- colors
  k <- 1
  for (i in strip_top) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  plot <- as_ggplot(g)
  
  if(is.null(output_plot)){
    return(plot)
  }else{
    xleng <- 14
    yleng <- 5
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = plot,width = plot_width,height = plot_height)
  }
  
}




plot_profile_dbs78_signal_logo <- function (profile, colors = dbs78_signal_colors,output_plot = NULL,plot_width=NULL, plot_height=NULL) 
{
  ## make sure the order profile 1 and profile 2 are the same order
  
  ## normalize 
  colnames(profile) <- c('MutationType','Value')
  profile <- profile_format_df(profile,factortype = TRUE,indel_short = indel_short) 
  names(colors) <- levels(profile$Type)
  
  
  profile[,4] <- profile[,4]/sum(profile[,4])  
  
  plot <- ggplot(data = profile, aes(x = MutationType, y = Value,  fill = Type, width = 0.7)) + 
    geom_bar(stat = "identity", position = "identity", colour = NA, size = 0) + 
    scale_fill_manual(values = colors) + 
    facet_grid(. ~ Type, scales = "free",space = "free_x") +
    guides(fill = FALSE) + 
    scale_x_discrete(expand = expansion(add = 0))+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)),breaks = pretty_breaks(5))+
    labs(x="",y="")+
    theme(axis.ticks = element_blank(),axis.title = element_blank(), axis.text = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),strip.text.x = element_blank(),panel.spacing.x = unit(0, "lines"),axis.line.x = element_line(colour = 'black',size = 0.05))
  
  ## add background color for strip
  # require(grid)
  # g <- ggplot_gtable(ggplot_build(plot))
  # strip_top <- which(grepl('strip-t', g$layout$name))
  # fills <- colors
  # k <- 1
  # for (i in strip_top) {
  #   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  #   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  #   k <- k+1
  # }
  # plot <- as_ggplot(g)
  
  if(is.null(output_plot)){
    return(plot)
  }else{
    xleng <- 2
    yleng <- 0.8
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = plot,width = plot_width,height = plot_height)
  }
  
}



# Heatmap of profiles -----------------------------------------------------

profile_heatmap_plot <- function(data,output_plot = NULL,plot_width=NULL, plot_height=NULL){
  
  # data format: Sample  Profile Mutations
  
  # merge different profile 
  #if(length(unique(data_input$Sample)))
  data <- data %>% mutate(Key=str_remove(Profile,"\\d+")) %>% group_by(Sample,Mutations,Key) %>%mutate(Matrix=parse_number(Profile)) %>%arrange(Key,Matrix) %>% mutate(Profile=paste0(Matrix,collapse = '/')) %>% mutate(Profile=paste0(Key,": ",Profile)) %>% select(-Matrix) %>% ungroup() %>% unique() %>% select(-Key) 
  
  profile_tmp <- data %>% group_by(Profile) %>% summarise(mean=mean(Mutations)) %>% arrange(desc(mean))
  samleve_tmp <- data %>% filter(Profile==profile_tmp$Profile[1]) %>% arrange(Mutations) %>% pull(Sample)
  
  # typedata_tmp <- seqmatrix_refdata_public %>% group_by(Profile,MutationType) %>% summarise(mean=mean(Mutations,na.rm = TRUE)) %>% ungroup() %>% arrange(desc(mean)) %>% group_by(Profile) %>% slice(1:5) %>% mutate(Seq2=6-seq_along(MutationType)) %>% ungroup() %>% left_join(profile_tmp %>% select(-mean)) %>% mutate(Seq3=Seq+0.1666667*Seq2)
  # 
  # data2 <- seqmatrix_refdata_public %>% 
  #   left_join(typedata_tmp %>% select(Profile,MutationType,Seq,Seq2)) %>% 
  #   filter(!is.na(Seq)) %>% 
  #   left_join(
  #     seqmatrix_refdata_public %>% group_by(Sample,Profile) %>% summarise(total=sum(Mutations))
  #   ) %>% 
  #   mutate(Ratio=Mutations/total)
  # 
  # max_tmp <- data2 %>% group_by(Profile) %>% summarise(max=max(Ratio))
  # 
  # data2 <- data2 %>% left_join(max_tmp) %>% mutate(value=Ratio/max+Seq)
  name_max <- max(str_length(unique(data$Sample)))
  name_len <- length(unique(data$Sample))
  angle_max <- if_else(name_max < 15, 90, 30)
  vjust_max <- if_else(name_max < 15, 0.5, 1)
  
  if(name_len>0){
    
    plot_final <- data %>% left_join(profile_tmp) %>% 
      mutate(Sample=factor(Sample,levels = samleve_tmp),Profile=factor(Profile,levels = profile_tmp$Profile)) %>% 
      ggplot(aes(Sample,log10(Mutations),group=Profile,col=Profile))+
      geom_line()+
      geom_point()+
      scale_color_d3()+
      scale_y_continuous(breaks = pretty_breaks())+
      labs(x="Sample index",y="log10(Mutations)",color="Profile")+
      #title = 'Number of Mutations Per Sample with Regard to Mutational Profile',
      theme_ipsum_rc(grid = "XYy",ticks = TRUE,axis = FALSE,axis_title_just = 'm',axis_title_size = 14)+
      theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),plot.title = element_text(hjust = 0.5))+
      panel_border(color = 'black')
    #,legend.direction = "vertical",legend.position = "top"
    
  }else{
    
    plot_final <- data %>% left_join(profile_tmp) %>% 
      mutate(Sample=factor(Sample,levels = samleve_tmp),Profile=factor(Profile,levels = profile_tmp$Profile)) %>% 
      ggplot(aes(Sample,Profile,fill=(Mutations)))+
      geom_tile(col="white")+
      scale_fill_viridis_c(trans = "log10",label=comma_format(),na.value = 'black')+
      labs(x="",y="",fill="Number of mutations\n")+
      #title = 'Number of Mutations Per Sample with Regard to Mutational Profile'
      theme_ipsum_rc(grid = FALSE,ticks = FALSE,axis = FALSE)+
      theme(legend.key.width =unit(2, "cm"),legend.position = "top",plot.title = element_text(hjust = 0.5))
    
    if(name_len<=80){
      plot_final <- plot_final + theme(axis.text.x = element_text(angle = angle_max,hjust = 1,vjust = vjust_max,size = 10))
    } else{
      plot_final <- plot_final + theme(axis.text.x = element_blank())
    }
  }
  
  ## define the length of x and y
  leng0 <- 2.5
  leng_ratio <-  0.2
  
  xleng <- leng_ratio*(name_len)+leng0+4
  
  xleng <- if_else(xleng>15,15,xleng)
  #yleng <- leng_ratio*name_max+2.5+leng0
  yleng <- 6
  
  if(is.null(output_plot)){
    return(plot_final)
  }else{
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = plot_final,width = plot_width,height = plot_height)
  }
  
}



# Mutational Patten funcitons ---------------------------------------------

profile_format_df <- function(data,factortype=FALSE,indel_short=FALSE){
  # format, sorting and factor the signature dataframe for SBS96,DBS78 and ID83
  
  if(dim(data)[1]==96){
    data <- data %>% mutate(Type=str_sub(MutationType,3,5),SubType=paste0(str_sub(MutationType,1,1),str_sub(MutationType,3,3),str_sub(MutationType,7,7))) %>% select(Type,SubType,MutationType,everything()) %>% arrange(Type,SubType)
  }
  
  if(dim(data)[1]==78){
    data <- data %>% mutate(Type=paste0(str_sub(MutationType,1,3),"NN"),SubType=str_sub(MutationType,4,5)) %>% select(Type,SubType,MutationType,everything()) %>% arrange(Type,SubType)
  }
  
  if(dim(data)[1]==83){
    idtypeorder <- c("1:Del:C","1:Del:T","1:Ins:C","1:Ins:T","2:Del:R","3:Del:R","4:Del:R","5:Del:R","2:Ins:R","3:Ins:R","4:Ins:R","5:Ins:R","2:Del:M","3:Del:M","4:Del:M","5:Del:M")
    data <- data %>% mutate(Type=str_sub(MutationType,1,7),SubType=str_sub(MutationType,9,9)) %>% select(Type,SubType,MutationType,everything()) %>% mutate(Type=factor(Type,levels = idtypeorder)) %>% arrange(Type,SubType) %>% mutate(Type=as.character(Type))
  }
  
  if(dim(data)[1]==68){
    idtypeorder <- c("0:homdel","1:LOH","2:LOH","3-4:LOH","5-8:LOH","9+:LOH","2:het","3-4:het","5-8:het","9+:het")
    idsubtypeorder <- c('0-1kb','1-10kb','10-100kb','100kb-1Mb','>1Mb','1-10Mb','10-40Mb','>40Mb')
    idsubtypelabel <- c('0-1kb','1kb-10kb','10kb-100kb','100kb-1Mb','>1Mb','1Mb-10Mb','10Mb-40Mb','>40Mb')
    data <- data %>% mutate(tmp=MutationType) %>% separate(tmp,into = c('a','b','c'),sep = ':') %>% mutate(Type=paste0(a,":",b)) %>% select(-a,-b) %>% select(Type,SubType=c,MutationType,everything()) %>% mutate(Type=factor(Type,levels = idtypeorder),SubType=factor(SubType,levels = idsubtypeorder,labels=idsubtypelabel)) %>% arrange(Type,SubType) %>% mutate(Type=as.character(Type))
  }
  
  if(dim(data)[1]==38){
    idtypeorder <- c("DEL:clustered",'DUP:clustered','INV:clustered','TLOC:clustered',"DEL:unclustered",'DUP:unclustered','INV:unclustered','TLOC:unclustered')
    idsubtypeorder <-  c('NA','<1kb','1-10kb','10-100kb','100kb-1Mb','1Mb-10Mb','>10Mb')
    idsubtypelabel <-  c('','<1kb','1kb-10kb','10kb-100kb','100kb-1Mb','1Mb-10Mb','>10Mb')
    data <- data %>% mutate(Type=str_replace(MutationType,":.*:",":"),SubType=str_remove(str_remove(MutationType,":[^:]*$"),'^.*:')) %>% select(Type,SubType,MutationType,everything()) %>% mutate(Type=factor(Type,levels = idtypeorder),SubType=factor(SubType,levels = idsubtypeorder,labels=idsubtypelabel)) %>% arrange(Type,SubType) %>% mutate(Type=as.character(Type))
  }
  
  if(dim(data)[1]==192){
    data <- data %>% separate(MutationType,into = c('Strand','MutationType'),sep = ':')%>% mutate(Type=str_sub(MutationType,3,5),SubType=paste0(str_sub(MutationType,1,1),str_sub(MutationType,3,3),str_sub(MutationType,7,7))) %>% select(Type,SubType,MutationType,everything()) %>% arrange(Type,SubType)
  }
  
  if(factortype){
    
    if(dim(data)[1]==83){
      tmplev <- idtypeorder
      tmplab <- str_replace(tmplev,"5","5+")
      
      if(indel_short){
        tmplab = if_else(tmplab %in% c('2:Del:M','3:Del:M','4:Del:M'),str_remove(tmplab,":.*"),tmplab)
      }
      
      data <- data %>% 
        mutate(Type=factor(Type,levels = tmplev,labels = tmplab)) %>% 
        mutate(SubType = if_else(str_detect(Type,'Del') & !str_detect(Type,"M"),as.character(as.integer(SubType)+1L),SubType)) %>% 
        mutate(SubType = if_else(str_detect(Type,'Del') & !str_detect(Type,"M"), str_replace(SubType,"6","6+"),str_replace(SubType,"5","5+")))
    }else {
      tmplev <- unique(data$Type)
      data <- data %>% mutate(Type=factor(Type,levels = tmplev))
    }
    
    if(dim(data)[1]==192){ data <- data %>% mutate(Strand=factor(Strand, levels = c('T','U')))}
    
    #tmplev <- unique(data$SubType)
    #data <- data %>% mutate(SubType=factor(SubType,levels = tmplev))
    tmplev <- unique(data$MutationType)
    data <- data %>% mutate(MutationType=factor(MutationType,levels = tmplev))
  }
  
  return(data)
  
}


## need to define profile_format_df2 ## 






# Calculate cosine similarity between two signature in dataframe format --------
cos_sim_df_old <- function (mut_df1, mut_df2, output_matrix=FALSE) 
{
  colnames(mut_df1)[1] <- "MutationType"
  colnames(mut_df2)[1] <- "MutationType"
  
  mut_matrix1 <- mut_df1 %>% 
    arrange(MutationType) %>% 
    select(-MutationType) %>% 
    as.matrix()
  
  mut_matrix2 <- mut_df2 %>% 
    arrange(MutationType) %>% 
    select(-MutationType) %>% 
    as.matrix()
  
  n_samples1 = ncol(mut_matrix1)
  n_samples2 = ncol(mut_matrix2)
  res_matrix = matrix(nrow = n_samples1, ncol = n_samples2)
  for (s in 1:n_samples1) {
    signal1 = mut_matrix1[, s]
    cos_sim_vector = c()
    for (i in 1:n_samples2) {
      signal2 = mut_matrix2[, i]
      cos_sim_vector[i] = cos_sim(signal1, signal2)
    }
    res_matrix[s, ] = cos_sim_vector
  }
  rownames(res_matrix) = colnames(mut_matrix1)
  colnames(res_matrix) = colnames(mut_matrix2)
  
  if(output_matrix){
    return(res_matrix)
  }else {
    res_matrix %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()
  }
}


# Calculate cosine similarity between two signature in dataframe format using coop package --------
cos_sim_df <- function (mut_df1, mut_df2, output_matrix=FALSE) 
{
  require(coop)
  
  colnames(mut_df1)[1] <- "MutationType"
  colnames(mut_df2)[1] <- "MutationType"
  
  ### fix bug if only 1 profile 
  if(dim(mut_df1)[2]==2) {
    mut_df1 <- mut_df1[,c(1,2,2)]
    colnames(mut_df1)[3] <- 'XXXXX'
  }
  
  if(dim(mut_df2)[2]==2) {
    mut_df2 <- mut_df2[,c(1,2,2)]
    colnames(mut_df2)[3] <- 'YYYYY'
  }
  
  
  if(identical(mut_df1,mut_df2)){
    # two identical matrix/df
    mut_matrix1 <- mut_df1 %>% 
      arrange(MutationType) %>% 
      select(-MutationType) %>% 
      as.matrix()
    res_matrix <- coop::cosine(mut_matrix1) 
    
  }else{
    # two differnt matrix/df
    df_name1 <- colnames(mut_df1)[-1]
    df_name2 <- colnames(mut_df2)[-1]
    new_df_name1 <- paste0("N1_",seq_along(df_name1))
    new_df_name2 <- paste0("N2_",seq_along(df_name2))
    colnames(mut_df1)[-1] <- new_df_name1
    colnames(mut_df2)[-1] <- new_df_name2  
    mut_df <- left_join(mut_df1,mut_df2)
    mut_matrix <- mut_df %>% 
      arrange(MutationType) %>% 
      select(-MutationType) %>% 
      as.matrix()
    
    res_matrix <- coop::cosine(mut_matrix) 
    res_matrix <- res_matrix[new_df_name1,new_df_name2]
    rownames(res_matrix) <- df_name1
    colnames(res_matrix) <- df_name2
  }
  
  if(output_matrix){
    return(res_matrix)
  }else {
    res_matrix %>% as.data.frame() %>% rownames_to_column() %>% as_tibble() %>% select(-contains('XXXXX'),-contains('YYYYY')) %>% filter(!rowname %in% c('XXXXX','YYYYY'))
  }
  
}



# Heatmap of cosine similairty  -------------------------------------------
plot_cosine_heatmap_df <- function (cos_sim_df, col_order, cluster_rows = TRUE, method = "complete", nmax = 200L, plot_values = FALSE,output_plot = NULL,plot_width=NULL, plot_height=NULL) 
{
  colnames(cos_sim_df)[1] <- "Sample"
  # if (class(cos_sim_matrix) != "matrix") {
  #   stop("cos_sim_matrix must be a matrix")
  # }
  if (length(colnames(cos_sim_df)) == 0) {
    stop("cos_sim_df is missing colnames")
  }
  if (length(rownames(cos_sim_df)) == 0) {
    stop("cos_sim_df is missing rownames")
  }
  
  # limited to max row/col 
  
  if(dim(cos_sim_df)[1] > nmax | dim(cos_sim_df)[2] > (nmax+1)){
    nrow <- if_else(dim(cos_sim_df)[1]>nmax, nmax, dim(cos_sim_df)[1])
    ncol <- if_else(dim(cos_sim_df)[2]>nmax+1L, nmax+1L, dim(cos_sim_df)[2])
    cos_sim_df <- cos_sim_df[1:nrow,1:ncol]
  }
  
  
  # covert to matrix
  cos_sim_matrix <- as.matrix(cos_sim_df[,-1])
  rownames(cos_sim_matrix) <- cos_sim_df[[1]]
  
  if (missing(col_order)) {
    #col_order = colnames(cos_sim_df)[-1]
    hc.sample = hclust(dist(t(cos_sim_matrix)), method = method)
    col_order = rownames(t(cos_sim_matrix))[hc.sample$order]
  }
  
  if (class(col_order) != "character") {
    stop("col_order must be a character vector")
  }
  if (length(col_order) != ncol(cos_sim_df)-1) {
    stop("col_order must have the same length as the number of signatures in the explained df")
  }
  
  if (cluster_rows == TRUE) {
    hc.sample = hclust(dist(cos_sim_matrix), method = method)
    sample_order = rownames(cos_sim_matrix)[hc.sample$order]
  }
  else {
    sample_order = rownames(cos_sim_matrix)
  }
  Cosine.sim = NULL
  Signature = NULL
  Sample = NULL
  x = NULL
  y = NULL
  xend = NULL
  yend = NULL
  #cos_sim_matrix.m = melt(cos_sim_matrix)
  #colnames(cos_sim_matrix.m) = c("Sample", "Signature", "Cosine.sim")
  cos_sim_matrix.m <- cos_sim_df %>% pivot_longer(-1, names_to="Signature",values_to="Cosine.sim") %>% select(Sample,Signature,Cosine.sim)
  cos_sim_matrix.m$Signature = factor(cos_sim_matrix.m$Signature,  levels = col_order)
  cos_sim_matrix.m$Sample = factor(cos_sim_matrix.m$Sample,  levels = sample_order)
  
  cos_sim_matrix.m <- cos_sim_matrix.m %>% mutate(coslab=round(Cosine.sim, 2))
  cos_sim_matrix.m$Cosine.sim <- round(cos_sim_matrix.m$Cosine.sim,digits = 2)
  mincosine <- min(cos_sim_matrix.m$Cosine.sim)
  maxcosine <- max(cos_sim_matrix.m$Cosine.sim)
  
  if(mincosine< -0.1){
    plot_values <- TRUE
    cos_sim_matrix.m$Cosine.sim <- abs(cos_sim_matrix.m$Cosine.sim)
    mincosine <- min(cos_sim_matrix.m$Cosine.sim)
    cos_sim_matrix.m$coslab <- if_else(cos_sim_matrix.m$coslab>0,"","-")
  }
  mincosine <- if_else(mincosine<0.1,0,mincosine)
  maxcosine <- if_else(maxcosine>0.9,1,maxcosine)
  
  
  ## define the length of x and y
  name_max_x <- max(str_length(unique(cos_sim_matrix.m$Signature)))
  name_max_y <- max(str_length(unique(cos_sim_matrix.m$Sample)))
  name_len_x <- length(unique(cos_sim_matrix.m$Signature))
  name_len_y <- length(unique(cos_sim_matrix.m$Sample))
  
  angle_max_x <- if_else(name_max_x < 15, 90, 30)
  vjust_max_x <- if_else(name_max_x < 15, 0.5, 1)
  
  leng0 <- 2.5
  leng_ratio <-  0.2
  lenx <- name_max_x*0.4
  leny <- name_max_y*0.4
  xleng <- leng_ratio*name_len_x+leng0+3+leny
  yleng <- leng_ratio*name_len_y+leng0+1+lenx
  
  xleng <- if_else(xleng>30,30,xleng)
  yleng <- if_else(yleng>25,25,yleng)
  
  heatmap = ggplot(cos_sim_matrix.m, aes(x = Signature, y = Sample,  fill = Cosine.sim, order = Sample)) + 
    geom_tile(color = "white") + 
    scale_fill_viridis_c(name = "Cosine similarity\n", limits = c(mincosine, maxcosine),breaks=scales::pretty_breaks())+
    # scale_fill_distiller(palette = "YlGnBu", direction = 1,name = "Cosine similarity", limits = c(0, 1)) + 
    theme_ipsum_rc(grid = FALSE,ticks = T)+
    #theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) +
    labs(x = NULL, y = NULL)+
    scale_y_discrete(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    theme(legend.position = "top",legend.key.width = unit(2,"cm"))
  #legend.key.width = unit(1.5,"cm"),legend.key.height = unit(0.3,"cm")
  
  
  if(name_len_x<=80){
    heatmap <- heatmap + theme(axis.text.x = element_text(angle = angle_max_x,hjust = 1,vjust = vjust_max_x,size = 10))
  } else{
    heatmap <- heatmap + theme(axis.text.x = element_blank())
  }
  
  if(name_len_y>80){
    heatmap <- heatmap + theme(axis.text.y = element_blank())
  }
  
  
  if (plot_values) {
    heatmap = heatmap + geom_text(aes(label = coslab), size = 3,col="red")
  }
  if (cluster_rows == TRUE) {
    dhc = as.dendrogram(hc.sample)
    ddata = ggdendro::dendro_data(dhc, type = "rectangle")
    dendrogram = ggplot(ggdendro::segment(ddata)) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() + 
      scale_y_reverse(expand = c(0.1, 0)) + scale_x_continuous(expand = expansion(add=0.5)) + ggdendro::theme_dendro()
    plot_final = plot_grid(dendrogram+theme(plot.margin = margin(r = -0.2,unit = "cm")), heatmap+theme(plot.margin = margin(l = -0.2,r=0.5,unit = "cm")), align = "h",axis = "tb",rel_widths = c(0.15, 1))
  } else {
    plot_final = heatmap + ylim(rev(levels(factor(cos_sim_matrix.m$Sample))))
  }
  
  
  if(is.null(output_plot)){
    return(plot_final)
  }else{
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = plot_final,width = plot_width,height = plot_height)
  }
  
}




# Plot two profile difference for SBS96, ID83 and DBS78 ---------------------------------------------
plot_compare_profiles_diff <- function (profile1, profile2, profile_names = NULL, profile_ymax = NULL, diff_ylim = NULL, colors = NULL, condensed = FALSE,output_plot = NULL,plot_width=NULL, plot_height=NULL, output_data=NULL, normalize=TRUE) 
{
  ## profile 1 and profile 2 will be the dataframe with two columns: MutationType and value
  
  COLORS2 = c("#1F77B4FF","#D62728FF")
  COLORS6 = c("#03BCEE", "#010101", "#E32926", "#CAC9C9", "#A1CE63", "#EBC6C4")
  COLORS10 = c("#03BCEE", "#0366CB", "#A1CE63", "#016601", "#FE9898", "#E32926", "#FEB166", "#FE8001", "#CB98FE", "#4C0198")
  COLORS16=c("#FCBD6F", "#FE8002", "#AFDC8A", "#36A02E", "#FCC9B4", "#FB896A", "#F04432", "#BB191A", "#CFE0F1", "#93C3DE", "#4A97C8", "#1764AA", "#E1E1EE", "#B5B5D7", "#8582BC", "#62409A")
  
  typelength = dim(profile1)[1]
  typelength2 = dim(profile2)[1]
  
  if(!(typelength %in% c(78,83,96,192))){
    if(typelength == 0){
      stop('Error: Sample 1 have 0 mutation of selected profile')
    }else{
      stop('Error: Sample 1 detected unsupported profile types. Currently, Profile Comparision only supports SBS78, SBS192, DBS78, and ID83.')
    }
  }
  
  if(!(typelength2 %in% c(78,83,96,192))){
    if(typelength2 == 0){
      stop('Error: Sample 1 have 0 mutation of selected profile')
    }else{
      stop('Error: Sample 1 detected unsupported profile types. Currently, Profile Comparision only supports SBS78, SBS192, DBS78, and ID83.')
    }
  }
  
  
  
  
  if (is.null(colors)) {
    if(typelength == 192){colors = COLORS6;}
    if(typelength == 96){colors = COLORS6; }
    if(typelength == 78){colors = COLORS10;}
    if(typelength == 83){colors = COLORS16;}
  }
  
  ## make sure the order profile 1 and profile 2 are the same order
  
  indel_short <-  FALSE
  if(typelength == 83) {indel_short = TRUE}
  profile1 <- profile_format_df(profile1,factortype = TRUE,indel_short = indel_short) 
  profile2 <- profile_format_df(profile2,factortype = TRUE,indel_short = indel_short)
  
  names(colors) <- levels(profile1$Type)
  
  dim2 <- dim(profile1)[2]
  
  #print(colors)
  if(normalize){
    profile1[,dim2] <- profile1[,dim2]/sum(profile1[,dim2])  
    profile2[,dim2] <- profile2[,dim2]/sum(profile2[,dim2]) 
  }
  
  diff = profile1[[dim2]] - profile2[[dim2]]
  RSS = sum(diff^2)
  RSS = format(RSS, scientific = TRUE, digits = 3)
  cosine_sim = cos_sim(profile1[[dim2]], profile2[[dim2]])
  cosine_sim = round(cosine_sim, 3)
  if(colnames(profile1)[dim2] == colnames(profile2)[dim2]){
    colnames(profile1)[dim2] <-  "Profile1"
    colnames(profile2)[dim2] <-  "Profile2"
  }
  
  df <-  profile1 %>% left_join(profile2) %>% mutate(Difference=diff)
  if(is.null(profile_names)){
    profile_names <- colnames(df)[dim2:(dim2+1)]
  }
  colnames(df)[dim2:(dim2+1)] <- profile_names
  df <- df %>% pivot_longer(cols = -c(1:(dim2-1))) %>% mutate(name=factor(name,levels = c(profile_names,"Difference"))) %>% mutate(Type=factor(Type,levels = names(colors)))
  
  if(is.null(profile_ymax)){
    profile_ymax <- max(c(profile1[[dim2]],profile2[[dim2]]))*1.1
  }
  
  if(is.null(diff_ylim)){
    #diff_ylim <- range(diff)*1.5
    diff_ylim <- range(diff)
    diff_ylim[1] <- if_else(diff_ylim[1] > -0.02, -0.02,diff_ylim[1])
    diff_ylim[2] <- if_else(diff_ylim[2] < 0.02, 0.02,diff_ylim[2])
  }
  dftmp = tibble(Type = rep(levels(profile1$Type)[1], 4), SubType = rep((profile1$SubType)[1], 4), Strand = 'T', name = c(profile_names, "Difference", "Difference"), value = c(profile_ymax, profile_ymax, diff_ylim[1], diff_ylim[2])) %>% mutate(name=factor(name,levels = c(profile_names,"Difference"))) %>% mutate(Type=factor(Type,levels = names(colors)))
  
  
  ## if stranded
  if('Strand' %in% colnames(df)){
    
    StrandColor <- c("#1F77B4FF", "#D62728FF","#2CA02CFF")
    names(StrandColor) <- c('Transcribed Strand','Untranscribed Strand','Nontranscribed Strand')
    df <- df %>% mutate(Strand = factor(Strand,levels = c('T','U','N'),labels = c('Transcribed Strand','Untranscribed Strand','Nontranscribed Strand')))
    dftmp <- dftmp %>% mutate(Strand = factor(Strand,levels = c('T','U','N'),labels = c('Transcribed Strand','Untranscribed Strand','Nontranscribed Strand')))
    df$Strand <- fct_drop(df$Strand)
    dftmp$Strand <- fct_drop(dftmp$Strand)
    StrandColor <- StrandColor[levels(df$Strand)]
  }
  
  if (condensed) {
    
    if('Strand' %in% colnames(df)){
      plot <- ggplot(data = df, aes(x = SubType, y = value,  fill = Strand, width = 1)) + 
        geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), colour = "black", size = 0.2) + 
        geom_point(data = dftmp , aes(x = SubType, y = value), alpha = 0,size=0) + 
        scale_fill_manual(values = StrandColor) + 
        facet_grid(name ~ Type, scales = "free",space = "free_x") +
        ylab("Relative contribution") +
        labs(x="",fill="")+
        scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
        theme_ipsum_rc(axis_title_just = "m",grid = "Y",axis = TRUE) + 
        ggtitle(paste("RSS = ", RSS, "; Cosine Similarity = ", cosine_sim, sep = ""))+
        theme(axis.title.y = element_text(size = 14, vjust = 1), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8,face = 'bold',angle = 90, vjust = 0.5), strip.text.x = element_text(size = 10,hjust = 0.5,face = 'bold',colour = "white", margin = margin(0.1,0.1,0.1,0.1,unit = 'cm')), strip.text.y = element_text(size = 11,hjust = 0.5,margin = margin(0.1,0.15,0.1,0.15,unit = 'cm')),strip.background = element_rect(fill = "#f0f0f0"), panel.grid.major.x = element_blank(), panel.spacing.x = unit(0, "lines"),panel.spacing.y = unit(0.2, "lines"),plot.title = element_text(size=15,face = 'bold',hjust = 0.5),axis.line.y = element_line(colour = 'black',size = 0.25),legend.position = c(0.9,0.94),legend.text = element_text(size = 12),legend.key.size = unit(0.7,"line"))+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,colour = 'black',size = 0.5)+geom_vline(xintercept = Inf,colour = 'black',size = 0.5)
      
    }else{
      plot = ggplot(data = df, aes(x = SubType, y = value,  fill = Type, width = 1)) + 
        geom_bar(stat = "identity", position = "identity", colour = "black", size = 0.2) + 
        geom_point(data = dftmp, aes(x = SubType, y = value), alpha = 0,size=0) + 
        scale_fill_manual(values = colors) + 
        facet_grid(name ~ Type, scales = "free",space = "free_x") +
        ylab("Relative contribution") +
        guides(fill = FALSE) + 
        labs(x="")+
        scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
        #scale_y_continuous(expand = c(0,0))+
        theme_ipsum_rc(axis_title_just = "m",grid = "Y",axis = TRUE) + 
        ggtitle(paste("RSS = ", RSS, "; Cosine Similarity = ", cosine_sim, sep = ""))+
        theme(axis.title.y = element_text(size = 14, vjust = 1), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8,face = 'bold',angle = 90, vjust = 0.5), strip.text.x = element_text(size = 10,hjust = 0.5,face = 'bold',colour = "white", margin = margin(0.1,0.1,0.1,0.1,unit = 'cm')), strip.text.y = element_text(size = 11,hjust = 0.5,margin = margin(0.1,0.15,0.1,0.15,unit = 'cm')),strip.background = element_rect(fill = "#f0f0f0"), panel.grid.major.x = element_blank(), panel.spacing.x = unit(0, "lines"),panel.spacing.y = unit(0.2, "lines"),plot.title = element_text(size=15,face = 'bold',hjust = 0.5),axis.line.y = element_line(colour = 'black',size = 0.25))+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,colour = 'black',size = 0.5)+geom_vline(xintercept = Inf,colour = 'black',size = 0.5)#axis.line.x = element_line(colour = 'black',size = 0.25),
      #panel_border(color = gray(0.5),size = 0.3)
    }
  }
  else {
    
    if('Strand' %in% colnames(df)){
      plot = ggplot(data = df, aes(x = SubType, y = value,  fill = Strand, width = 0.7)) + 
        geom_bar(stat = "identity", position = position_dodge2(preserve = "single"), colour = "black", size = 0) + 
        geom_point(data = dftmp, aes(x = SubType, y = value), alpha = 0,size=0) + 
        scale_fill_manual(values = StrandColor) + 
        facet_grid(name ~ Type, scales = "free",space = "free_x") +
        ylab("Relative contribution") +
        labs(x="",fill="")+
        scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
        #scale_y_continuous(expand = c(0,0))+
        theme_ipsum_rc(axis_title_just = "m",grid = "Y",axis = TRUE) + 
        ggtitle(paste("RSS = ", RSS, "; Cosine Similarity = ", cosine_sim, sep = ""))+
        theme(axis.title.y = element_text(size = 14, vjust = 1), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8,face = 'bold', angle = 90, vjust = 0.5), strip.text.x = element_text(size = 10,hjust = 0.5,face = 'bold',colour = "white", margin = margin(0.1,0.1,0.1,0.1,unit = 'cm')), strip.text.y = element_text(size = 11,hjust = 0.5,margin = margin(0.1,0.15,0.1,0.15,unit = 'cm')),strip.background = element_rect(fill = "#f0f0f0",), panel.grid.major.x = element_blank(), panel.spacing.x = unit(0, "lines"),panel.spacing.y = unit(0.2, "lines"),plot.title = element_text(size=15,face = 'bold',hjust = 0.5),axis.line.y = element_line(colour = 'black',size = 0.25),legend.position = c(0.9,0.94),legend.text = element_text(size = 12),legend.key.size = unit(0.7,"line"))+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,colour = 'black',size = 0.5) +geom_vline(xintercept = Inf,colour = 'black',size = 0.5) #,axis.line.x = element_line(colour = 'black',size = 0.25)
      #     panel_border(color = gray(0.5),size = 0.3)
      
    }else{
      plot = ggplot(data = df, aes(x = SubType, y = value,  fill = Type, width = 0.7)) + 
        geom_bar(stat = "identity", position = "identity", colour = "black", size = 0) + 
        geom_point(data = dftmp, aes(x = SubType, y = value), alpha = 0,size=0) + 
        scale_fill_manual(values = colors) + 
        facet_grid(name ~ Type, scales = "free",space = "free_x") +
        ylab("Relative contribution") +
        guides(fill = 'none') + 
        labs(x="")+
        scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
        #scale_y_continuous(expand = c(0,0))+
        theme_ipsum_rc(axis_title_just = "m",grid = "Y",axis = TRUE) + 
        ggtitle(paste("RSS = ", RSS, "; Cosine Similarity = ", cosine_sim, sep = ""))+
        theme(axis.title.y = element_text(size = 14, vjust = 1), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8,face = 'bold', angle = 90, vjust = 0.5), strip.text.x = element_text(size = 10,hjust = 0.5,face = 'bold',colour = "white", margin = margin(0.1,0.1,0.1,0.1,unit = 'cm')), strip.text.y = element_text(size = 11,hjust = 0.5,margin = margin(0.1,0.15,0.1,0.15,unit = 'cm')),strip.background = element_rect(fill = "#f0f0f0",), panel.grid.major.x = element_blank(), panel.spacing.x = unit(0, "lines"),panel.spacing.y = unit(0.2, "lines"),plot.title = element_text(size=15,face = 'bold',hjust = 0.5),axis.line.y = element_line(colour = 'black',size = 0.25))+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,colour = 'black',size = 0.5) +geom_vline(xintercept = Inf,colour = 'black',size = 0.5) #,axis.line.x = element_line(colour = 'black',size = 0.25)
      #     panel_border(color = gray(0.5),size = 0.3)
    }
  }
  
  
  ## add background color for strip
  require(grid)
  g <- ggplot_gtable(ggplot_build(plot))
  strip_top <- which(grepl('strip-t', g$layout$name))
  fills <- colors
  k <- 1
  for (i in strip_top) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  plot <- as_ggplot(g)
  
  # output data
  if(!is.null(output_data)){
    df %>% pivot_wider() %>% write_delim(file = output_data,delim = '\t',col_names = T)
  }
  
  
  if(is.null(output_plot)){
    return(plot)
  }else{
    xleng <- 15
    yleng <- 7
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = plot,width = plot_width,height = plot_height)
  }
  
}


plot_compare_profiles_diff_strand <- function (profile1, profile2, profile_names = NULL, profile_ymax = NULL, diff_ylim = NULL, colors = NULL, condensed = FALSE) 
{
  
  ## check format: MutationType Strand Type  SubType value
  ## profile 1 and profile 2 will be the dataframe with two columns: MutationType and value
  
  ## make sure the order profile 1 and profile 2 are the same order
  profile1 <- profile1 %>% arrange(MutationType)
  profile2 <- profile2 %>% arrange(MutationType)
  
  profile1[,5] <- profile1[,5]/sum(profile1[,5])  
  profile2[,5] <- profile2[,5]/sum(profile2[,5]) 
  diff = profile1[[5]] - profile2[[5]]
  RSS = sum(diff^2)
  RSS = format(RSS, scientific = TRUE, digits = 3)
  cosine_sim = cos_sim(profile1[[5]], profile2[[5]])
  cosine_sim = round(cosine_sim, 3)
  if(colnames(profile1)[4] == colnames(profile2)[4]){
    colnames(profile1)[4] <-  "Profile1"
    colnames(profile2)[4] <-  "Profile2"
  }
  
  df <-  profile1 %>% left_join(profile2) %>% mutate(Difference=diff)
  if(is.null(profile_names)){
    profile_names <- colnames(df)[5:6]
  }
  colnames(df)[5:6] <- profile_names
  df <- df %>% pivot_longer(cols = -c(1,2,3,4)) %>% mutate(name=factor(name,levels = c(profile_names,"Difference")))
  
  if(is.null(profile_ymax)){
    profile_ymax <- max(c(profile1[[5]],profile2[[5]]))*1.1
  }
  
  if(is.null(diff_ylim)){
    diff_ylim <- range(diff)*1.1
  }
  dftmp = tibble(Type = rep(levels(as.factor(profile1$Type))[1], 4), SubType = rep(levels(as.factor(profile1$SubType))[1], 4), name = c(profile_names, "Difference", "Difference"), value = c(profile_ymax, profile_ymax, diff_ylim[1], diff_ylim[2])) %>% mutate(Strand="T")%>% mutate(name=factor(name,levels = c(profile_names,"Difference")))
  
  if (condensed) {
    plot = ggplot(data = df, aes(x = SubType, y = value,fill=Strand,group=Strand, width = 1)) + 
      geom_col(position = "dodge2")+
      #geom_bar(stat = "identity", position = "identity", colour = "black", size = 0.2) + 
      geom_point(data = dftmp, aes(x = SubType, y = value), alpha = 0,size=0) + 
      scale_fill_manual(values = rev(pal_nejm()(2)))+
      facet_grid(name ~ Type, scales = "free") +
      ylab("Relative contribution") +
      guides(fill = FALSE) + 
      labs(x="")+
      #scale_y_continuous(expand = c(0,0))+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
      theme_ipsum_rc(axis_title_just = "m",grid = "Y",axis = TRUE) + 
      ggtitle(paste("RSS = ", RSS, "; Cosine Similarity = ", cosine_sim, sep = ""))+
      theme(axis.title.y = element_text(size = 14, vjust = 1), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5), strip.text.x = element_text(size = 12,hjust = 0.5), strip.text.y = element_text(size = 14,hjust = 0.5),strip.background = element_rect(fill = "#f0f0f0"), panel.grid.major.x = element_blank(), panel.spacing.x = unit(0, "lines"),panel.spacing.y = unit(0.2, "lines"),plot.title = element_text(hjust = 0.5),axis.line.x = element_line(colour = 'black',size = 0.25),axis.line.y = element_line(colour = 'black',size = 0.25))+geom_vline(xintercept = Inf,colour = 'black',size = 0.5)
    #      panel_border(color = gray(0.5),size = 0.3)
  }
  else {
    plot = ggplot(data = df, aes(x = SubType, y = value, fill=Strand,group=Strand, width = 0.7)) + 
      #geom_bar(stat = "identity", position = "identity", colour = "black", size = 0) + 
      geom_col(position = "dodge2")+
      geom_point(data = dftmp, aes(x = SubType, y = value), alpha = 0,size=0) + 
      scale_fill_manual(values = rev(pal_nejm()(2)))+
      facet_grid(name ~ Type, scales = "free") +
      ylab("Relative contribution") +
      guides(fill = FALSE) + 
      labs(x="")+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
      #scale_y_continuous(expand = c(0,0))+
      theme_ipsum_rc(axis_title_just = "m",grid = "Y",axis = TRUE) + 
      ggtitle(paste("RSS = ", RSS, "; Cosine Similarity = ", cosine_sim, sep = ""))+
      theme(axis.title.y = element_text(size = 12, vjust = 1), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5), strip.text.x = element_text(size = 14,hjust = 0.5), strip.text.y = element_text(size = 14,hjust = 0.5),strip.background = element_rect(fill = "#f0f0f0"), panel.grid.major.x = element_blank(), panel.spacing.x = unit(0, "lines"),panel.spacing.y = unit(0.2, "lines"),plot.title = element_text(hjust = 0.5),axis.line.x = element_line(colour = 'black',size = 0.25),axis.line.y = element_line(colour = 'black',size = 0.25))+geom_vline(xintercept = Inf,colour = 'black',size = 0.5)
    # panel_border(color = gray(0.5),size = 0.3)
  }
  return(plot)
}




position_stack_and_nudge <- function(x = 0, y = 0, vjust = 1, reverse = FALSE) {
  ggproto(NULL, PositionStackAndNudge,
          x = x,
          y = y,
          vjust = vjust,
          reverse = reverse
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @noRd
PositionStackAndNudge <- ggproto("PositionStackAndNudge", PositionStack,
                                 x = 0,
                                 y = 0,
                                 
                                 setup_params = function(self, data) {
                                   c(
                                     list(x = self$x, y = self$y),
                                     ggproto_parent(PositionStack, self)$setup_params(data)
                                   )
                                 },
                                 
                                 compute_layer = function(self, data, params, panel) {
                                   # operate on the stacked positions (updated in August 2020)
                                   data = ggproto_parent(PositionStack, self)$compute_layer(data, params, panel)
                                   
                                   x_orig <- data$x
                                   y_orig <- data$y
                                   # transform only the dimensions for which non-zero nudging is requested
                                   if (any(params$x != 0)) {
                                     if (any(params$y != 0)) {
                                       data <- transform_position(data, function(x) x + params$x, function(y) y + params$y)
                                     } else {
                                       data <- transform_position(data, function(x) x + params$x, NULL)
                                     }
                                   } else if (any(params$y != 0)) {
                                     data <- transform_position(data, function(x) x, function(y) y + params$y)
                                   }
                                   data$nudge_x <- data$x
                                   data$nudge_y <- data$y
                                   data$x <- x_orig
                                   data$y <- y_orig
                                   
                                   data
                                 },
                                 
                                 compute_panel = function(self, data, params, scales) {
                                   ggproto_parent(PositionStack, self)$compute_panel(data, params, scales)
                                 }
)


# Profile mulitple profiles -----------------------------------------------

plot_compare_profiles_diff2 <- function (profile, profile_names = NULL, profile_ymax = NULL,  colors = NULL, condensed = FALSE,output_plot = NULL,plot_width=NULL, plot_height=NULL, xlab="") 
{
  ## profile 1 and profile 2 will be the dataframe with two columns: MutationType and value
  
  COLORS6 = c("#03BCEE", "#010101", "#E32926", "#CAC9C9", "#A1CE63", "#EBC6C4")
  COLORS10 = c("#03BCEE", "#0366CB", "#A1CE63", "#016601", "#FE9898", "#E32926", "#FEB166", "#FE8001", "#CB98FE", "#4C0198")
  COLORS16=c("#FCBD6F", "#FE8002", "#AFDC8A", "#36A02E", "#FCC9B4", "#FB896A", "#F04432", "#BB191A", "#CFE0F1", "#93C3DE", "#4A97C8", "#1764AA", "#E1E1EE", "#B5B5D7", "#8582BC", "#62409A")
  
  typelength = dim(profile)[1]
  if (is.null(colors)) {
    if(typelength == 96){colors = COLORS6; }
    if(typelength == 78){colors = COLORS10;}
    if(typelength == 83){colors = COLORS16;}
  }
  
  indel_short <-  FALSE
  if(typelength == 83) {indel_short = TRUE}
  profile <- profile_format_df(profile,factortype = TRUE,indel_short = indel_short) 
  
  names(colors) <- levels(profile$Type)
  
  #print(colors)
  
  profile[,4] <- profile[,4]/sum(profile[,4])  
  df <-  profile
  if(is.null(profile_names)){
    profile_names <- colnames(df)[4:5]
  }
  colnames(df)[4:5] <- profile_names
  df <- df %>% pivot_longer(cols = -c(1,2,3)) %>% mutate(name=factor(name)) %>% mutate(Type=factor(Type,levels = names(colors)))
  
  if(is.null(profile_ymax)){
    profile_ymax <- max(df$value)*1.1
  }
  
  
  #dftmp = tibble(Type = rep(levels(profile$Type)[1], 4), SubType = rep((profile$SubType)[1], 4), name = c(profile_names, "Difference", "Difference"), value = c(profile_ymax, profile_ymax, diff_ylim[1], diff_ylim[2])) %>% mutate(name=factor(name,levels = c(profile_names))) %>% mutate(Type=factor(Type,levels = names(colors)))
  
  
  if (condensed) {
    plot = ggplot(data = df, aes(x = SubType, y = value,  fill = Type, width = 1)) + 
      geom_bar(stat = "identity", position = "identity", colour = "black", size = 0.2) + 
      #geom_point(data = dftmp, aes(x = SubType, y = value), alpha = 0,size=0) + 
      scale_fill_manual(values = colors) + 
      facet_grid(name ~ Type, scales = "free",space = "free_x") +
      ylab("Relative contribution") +
      guides(fill = FALSE) + 
      labs(x=xlab)+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
      #scale_y_continuous(expand = c(0,0))+
      theme_ipsum_rc(axis_title_just = "m",grid = "Y",axis = TRUE) + 
      #ggtitle(paste("RSS = ", RSS, "; Cosine similarity = ", cosine_sim, sep = ""))+
      theme(axis.title.y = element_text(size = 14, vjust = 1), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5), strip.text.x = element_text(size = 10,hjust = 0.5,colour = "white", margin = margin(0,0,0,0,unit = 'cm')), strip.text.y = element_text(size = 10,hjust = 0.5,margin = margin(0,0,0,0,unit = 'cm')),strip.background = element_rect(fill = "#f0f0f0"), panel.grid.major.x = element_blank(), panel.spacing.x = unit(0, "lines"),panel.spacing.y = unit(0.2, "lines"),plot.title = element_text(hjust = 0.5),axis.line.y = element_line(colour = 'black',size = 0.25))+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,colour = 'black',size = 0.5)+geom_vline(xintercept = Inf,colour = 'black',size = 0.5)#axis.line.x = element_line(colour = 'black',size = 0.25),
    #panel_border(color = gray(0.5),size = 0.3)
  }
  else {
    plot = ggplot(data = df, aes(x = SubType, y = value,  fill = Type, width = 0.7)) + 
      geom_bar(stat = "identity", position = "identity", colour = "black", size = 0) + 
      #geom_point(data = dftmp, aes(x = SubType, y = value), alpha = 0,size=0) + 
      scale_fill_manual(values = colors) + 
      facet_grid(name ~ Type, scales = "free",space = "free_x") +
      ylab("Relative contribution") +
      guides(fill = FALSE) + 
      labs(x=xlab)+
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
      #scale_y_continuous(expand = c(0,0))+
      theme_ipsum_rc(axis_title_just = "m",grid = "Y",axis = TRUE) + 
      #ggtitle(paste("RSS = ", RSS, "; Cosine similarity = ", cosine_sim, sep = ""))+
      theme(axis.title.y = element_text(size = 14, vjust = 1), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5), strip.text.x = element_text(size = 10,hjust = 0.5,colour = "white",  margin = margin(0,0,0,0,unit = 'cm')), strip.text.y = element_text(size = 10,hjust = 0.5,margin = margin(0,0,0,0,unit = 'cm')),strip.background = element_rect(fill = "#f0f0f0",), panel.grid.major.x = element_blank(), panel.spacing.x = unit(0, "lines"),panel.spacing.y = unit(0.2, "lines"),plot.title = element_text(hjust = 0.5),axis.line.y = element_line(colour = 'black',size = 0.25))+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,colour = 'black',size = 0.5) +geom_vline(xintercept = Inf,colour = 'black',size = 0.5) #,axis.line.x = element_line(colour = 'black',size = 0.25)
    #     panel_border(color = gray(0.5),size = 0.3)
  }
  
  
  ## add background color for strip
  require(grid)
  g <- ggplot_gtable(ggplot_build(plot))
  strip_top <- which(grepl('strip-t', g$layout$name))
  fills <- colors
  k <- 1
  for (i in strip_top) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  plot <- as_ggplot(g)
  
  
  
  if(is.null(output_plot)){
    return(plot)
  }else{
    xleng <- 14
    yleng <- 7
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = plot,width = plot_width,height = plot_height)
  }
  
}







# Plot individual samples -------------------------------------------------
plot_individual_samples <-  function (exposure_refdata_input, signature_refsets_input, seqmatrix_refdata_input, profile_names = NULL, profile_ymax = NULL,  colors = NULL, condensed = FALSE,output_plot = NULL,plot_width=15, plot_height=10) 
{
  ## exposure_refdata_input:  Signature_name Exposure
  ## signature_refsets_input: Signature_name MutationType Contribution
  ## seqmatrix_refdata_input, MutationType xxx
  
  colnames(exposure_refdata_input) <- c("Signature_name", "Exposure")
  colnames(signature_refsets_input) <- c("Signature_name", "MutationType", "Contribution")
  colnames(seqmatrix_refdata_input) <- c("MutationType", "Contribution")
  
  present_sigs0 <- exposure_refdata_input %>% arrange(desc(Exposure)) %>% adorn_percentages(denominator = 'col') %>% filter(Exposure>0)
  present_sigs <-  present_sigs0 %>% filter(Exposure>0.01)
  label_sigs <- present_sigs %>% mutate(Exposure=percent(Exposure,accuracy = 0.1)) %>% mutate(Label=paste0(Exposure,"*",Signature_name)) %>% summarise(Label=paste0(Label,collapse = ' + ')) %>% mutate(Label=paste0('Orignal Profile = ',Label)) %>% pull(Label)
  
  orignal_profile <- seqmatrix_refdata_input %>% 
    adorn_percentages(denominator = 'col') 
  
  
  reconstructed_profile <- left_join(
    exposure_refdata_input,
    signature_refsets_input
  ) %>% mutate(Contribution=Exposure*Contribution) %>% 
    group_by(MutationType) %>% 
    summarise(Contribution=sum(Contribution)) %>% 
    adorn_percentages(denominator = 'col') 
  
  additional_profile <- signature_refsets_input %>% 
    filter(Signature_name %in% present_sigs$Signature_name) %>% 
    pivot_wider(names_from = Signature_name,values_from = Contribution) 
  
  
  
  # Barplot 
  if(str_detect(present_sigs0$Signature_name[1],'SBS')){
    sigcol <- SBScolor
  }else{
    sigcol <- c(pal_d3()(10),pal_aaas()(10),pal_npg()(10))[1:dim(present_sigs0)[1]]
  }
  
  
  p_bar <- present_sigs0 %>%
    mutate(Signature_name=factor(Signature_name),Signature_name=fct_reorder(Signature_name,Exposure)) %>% 
    ggplot(aes(x="Sample_ID",y=Exposure,fill=Signature_name,group=Signature_name)) +
    geom_bar(width = 0.2,position="stack", stat="identity")+
    geom_text_repel(aes(label=Signature_name,col=Signature_name),position = position_stack_and_nudge(vjust = 0.5,x=-0.5), direction = "y", hjust = "right")+
    theme_ipsum_rc(axis = FALSE,grid = FALSE,ticks = FALSE,axis_title_just = 'm',axis_title_size = 12)+
    theme(axis.text.x  = element_blank(),axis.text.y  = element_blank(),axis.title.y = element_blank(),legend.position = "none")+
    labs(x="")+
    scale_fill_manual(values = sigcol)+
    scale_color_manual(values = sigcol)
  
  
  p_diff <- plot_compare_profiles_diff(orignal_profile,reconstructed_profile,profile_names = c('Orignal', 'Deconstructed'),condensed = condensed)
  p_additional <- plot_compare_profiles_diff2(additional_profile,condensed = condensed,xlab=label_sigs)
  
  n_diff <- 3
  n_additional <- dim(additional_profile)[2]-1
  p_com <- plot_grid(p_diff+theme(plot.margin = unit(c(0, 0, -1, 0), "cm")),p_additional+theme(plot.margin = unit(c(-1, 0, 0, 0), "cm")),align = 'v',axis = 'lr',nrow = 2,ncol = 1,rel_heights = c(n_diff+1,n_additional))
  
  p_com_final <- plot_grid(p_bar+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),p_com+ theme(plot.margin = unit(c(0, 0, 0, -2), "cm")),align = 'h',axis = 'tb',nrow = 1,ncol = 2,rel_widths = c(1,7))
  
  if(is.null(output_plot)){
    return(p_com_final)
  }else{
    xleng <- 14
    yleng <- 7
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = p_com_final,width = plot_width,height = plot_height)
  }
  
  
}







# signature_sum_operation  --------------------------------------------------------------
## use the * and ;

signature_sum_operation <- function(sigdatabase,sigsetname,formulax,outsigname="Contribution") {
  #sigsetname <- 'COSMIC v3 Signatures (SBS)'
  #formulax <- "0.4*SBS1;0.4*SBS2;0.2*SBS13"
  formulax <- str_trim(formulax)
  inforx <- (str_split(unlist(str_split(formulax,"\\;")),"\\*"))
  #sigalltmp <- tibble(MutationType=character(),Contribution=numeric())
  for(i in 1:length(inforx)){
    tmp <- inforx[[i]]
    ratio <- as.numeric(tmp[[1]])
    signame <- tmp[[2]]
    sigtmp <- sigdatabase %>% filter(Signature_set_name==sigsetname,Signature_name==signame)
    sigtmp <- sigtmp %>% select(MutationType,tmpC=Contribution) %>% mutate(tmpC=ratio*tmpC)
    if(i==1){
      colnames(sigtmp)[2] <- 'Contribution'
      sigalltmp <- sigtmp
    }else{
      sigalltmp <- left_join(sigalltmp,sigtmp) %>% mutate(Contribution=Contribution+tmpC) %>% select(MutationType,Contribution)
    }
  }
  consum <- sum(sigalltmp$Contribution)
  sigalltmp$Contribution <- sigalltmp$Contribution/consum
  colnames(sigalltmp)[2] <- outsigname
  return(sigalltmp)
}




# Define Signature Set Colors  -------------------------------------------------------
sigsetcolor <- c(
  "Cancer Reference Signatures (RS)" = "#35978f",
  "Cancer Reference Signatures (SBS)"  = "#01665e",
  "COSMIC v2 Signatures (SBS)" = "#c994c7",
  "COSMIC v3 Signatures (DBS)" = "#df65b0",
  "COSMIC v3 Signatures (ID)" = "#e7298a",
  "COSMIC v3 Signatures (SBS)" = "#d73027",
  "Environmental Mutagen Signatures (SBS)" = "#ff7f00",
  "Organ-specific Cancer Signatures (RS)" = "#74a9cf",
  "Organ-specific Cancer Signatures (SBS)" = "#0570b0",
  "Other published signatures (ID)" = "#969696",
  "Other published signatures (SBS)" = "#525252",
  "SignatureAnalyzer PCAWG WGS 1536 Signatures (SBS)" = "#d9ef8b",
  "SignatureAnalyzer PCAWG WGS Signatures (DBS)" = "#a6d96a",
  "SignatureAnalyzer PCAWG WGS Signatures (ID)" = "#66bd63",
  "SignatureAnalyzer PCAWG WGS Signatures (SBS)" = "#1a9850",
  "SigProfiler PCAWG Strand Signatures (SBS)" = "#8073ac",
  "SigProfiler PCAWG WXS Signatures (SBS)" = "#542788"
)



# Signature pie chart -----------------------------------------------------

signature_piechart <- function(data,colset, output_plot = NULL,plot_width=NULL, plot_height=NULL){
  
  # convert nsig_data to another format for the piechart
  nsig_data_pie <- data %>%
    group_by(Profile) %>% 
    arrange(N) %>%
    mutate(
      end_angle = 2*pi*cumsum(N)/sum(N),   # ending angle for each pie slice
      start_angle = lag(end_angle, default = 0),   # starting angle for each pie slice
      mid_angle = 0.5*(start_angle + end_angle),   # middle of each pie slice, for the text label
      # horizontal and vertical justifications depend on whether we're to the left/right
      # or top/bottom of the pie
      hjust = ifelse(mid_angle > pi, 1, 0),
      vjust = ifelse(mid_angle < pi/2 | mid_angle > 3*pi/2, 0, 1)
    ) %>% 
    ungroup()
  
  # radius of the pie and radius for outside and inside labels
  rpie <- 1
  rlabel_out <- 1.05 * rpie
  rlabel_in <- 0.6 * rpie
  
  
  
  p <- ggplot(nsig_data_pie) +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = rpie, start = start_angle, end = end_angle, fill = Signature_set_name)) +
    geom_text(aes(x = rlabel_in * sin(mid_angle), y = rlabel_in * cos(mid_angle), label = N2 ), size = 14/.pt)+
    facet_wrap(~Profile,nrow = 2)+
    coord_fixed()+
    labs(x="",y="")+
    scale_fill_manual("Signature Set Name",values = colset)+ #breaks = names(colset)
    theme_ipsum_rc(axis = FALSE, grid = FALSE)+
    theme(axis.text.y = element_blank(),axis.text.x = element_blank(),strip.text = element_text(size = 14,hjust = 0.5,face = "bold"))
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 16
    yleng <- 7
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
  
} 



## TMB plot ###
TMBplot <- function(data,output_plot = NULL,plot_width=NULL, plot_height=NULL,addnote= NULL,ylab="Number of Mutations per Megabase (log10)\n"){
  # Cancer_Type,     Sample,   Burden
  # group Cancer Type
  if(dim(data)[1]<5){
    stop('Error: Not enough data to show the TMB plot')
  }
  
  data <- data %>%mutate(Burden=as.numeric(Burden)) %>%  mutate(Burden=if_else(is.infinite(Burden),NA_real_,Burden))
  mdata <- suppressMessages(data %>% group_by(Cancer_Type) %>% summarise(median=median(Burden,na.rm = TRUE),total=n(),nsample=sum(!is.na(Burden),na.rm=TRUE)) %>% arrange(median) %>% mutate(Seq=seq_along(Cancer_Type)*10-10))
  ## remove na cancer type
  natype <- mdata %>% filter(is.na(median)) %>% pull(Cancer_Type)
  mdata <- mdata %>% filter(!is.na(median))
  data <- data %>% filter(!(Cancer_Type %in% natype))
  
  data <- suppressMessages(data %>% mutate(Cancer_Type=factor(Cancer_Type,levels = mdata$Cancer_Type)) %>% arrange(Cancer_Type,Burden) %>% left_join(mdata))
  
  data <- data %>% filter(!is.na(Burden))%>% group_by(Cancer_Type) %>% mutate(Subseq=seq_along(Sample)) %>% ungroup() %>% mutate(score=Seq+1+8*Subseq/nsample)
  mdata <- mdata %>% mutate(Seq2=Seq+10) %>% mutate(score=(Seq+Seq2)/2) %>% 
    mutate(Burden=median,Type=if_else(seq_along(Cancer_Type) %% 2 ==0,"A","B"),Label=paste0(nsample,'<br /> - <br />',total)) %>% 
    mutate(Label_col = glue::glue("<span style='color:#377eb8'>{nsample}</span><br />\u2015<br /><span style='color:#4daf4a'>{total}</span>"))
  # "<span style='color:#377eb8'>{nsample}</span><br />\u2501<br /><span style='color:#4daf4a'>{total}</span>"
  f <- function(y) seq(floor(min(y,na.rm = TRUE)), ceiling(max(y,na.rm = TRUE)))
  
  ymint <- floor(min(data$Burden,na.rm = TRUE))
  ymint <- if_else(is.infinite(ymint),-4,ymint)
  ymaxt <- ceiling(max(data$Burden,na.rm = TRUE))
  
  print(data)
  p <-  data %>% filter(!is.na(Burden)) %>% 
    ggplot(aes(score,(Burden),group=Cancer_Type))+
    geom_rect(data=mdata,aes(xmin=Seq,xmax=Seq2,ymin=-Inf,ymax=Inf,fill=Type),alpha=0.5)+
    geom_point(pch=21,size=1,stroke=0.5)+
    geom_segment(data=mdata,aes(x=Seq+1,y=(median),xend=Seq2-1,yend=(median)),col="#e41a1c")+
    scale_fill_manual(values=c('#F2F2F2','#D4D4D4'))+
    scale_x_continuous(limits = c(min(mdata$Seq),max(mdata$Seq2)),expand = c(0,0),breaks = mdata$score,labels = mdata$Label_col,sec.axis = dup_axis(labels = mdata$Cancer_Type))+
    scale_y_continuous(limits = c(ymint,ymaxt),expand = c(0,0),breaks = f)+
    labs(x="",y=ylab)+
    theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 10,axis_text_size = 8,grid = "Y",ticks = TRUE)+
    theme(legend.position = "None",axis.text.x.top = element_text(size = 8,angle = -45,hjust = 1,vjust = 0),axis.text.x.bottom = element_markdown(size = 8),panel.border = element_rect(linewidth = 0.35,colour = "black",fill = NA))+
    coord_cartesian(clip = 'off')
  
  if(!is.null(addnote)){
    p <- p+annotate("text",x=-Inf, y = Inf, label = addnote, vjust=2, hjust=-0.2,size=4)
  }
  
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 2.5+length(unique(mdata$Cancer_Type))*0.5
    xleng <- if_else(xleng>18,18,xleng)
    yleng <- 6
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
  
}

TMBplot2 <- function(data,output_plot = NULL,plot_width=NULL,sgroupcol = pal_d3(palette = 'category20'), plot_height=NULL,addnote= NULL,ylab="Number of Mutations per Megabase (log10)\n",cancer_type_levels=NULL){
  # Cancer_Type,     Sample,   Burden, Sample_Group
  # group Cancer Type
  if(dim(data)[1]<5){
    stop('Error: Not enough data to show the TMB plot')
  }
  
  data <- data %>%mutate(Burden=as.numeric(Burden)) %>%  mutate(Burden=if_else(is.infinite(Burden),NA_real_,Burden))
  
  mdata <- suppressMessages(data %>% group_by(Cancer_Type) %>% summarise(median=median(Burden,na.rm = TRUE),total=n(),nsample=sum(!is.na(Burden),na.rm=TRUE)) %>% arrange(median) %>% mutate(Seq=seq_along(Cancer_Type)*10-10))
  
  if(!is.null(cancer_type_levels)){
    mdata <- suppressMessages(data %>% group_by(Cancer_Type) %>% summarise(median=median(Burden,na.rm = TRUE),total=n(),nsample=sum(!is.na(Burden),na.rm=TRUE)) %>% mutate(Cancer_Type=factor(Cancer_Type,levels=cancer_type_levels))%>% arrange(Cancer_Type) %>% mutate(Seq=seq_along(Cancer_Type)*10-10))
  }
  
  
  ## remove na cancer type
  natype <- mdata %>% filter(is.na(median)) %>% pull(Cancer_Type)
  mdata <- mdata %>% filter(!is.na(median))
  data <- data %>% filter(!(Cancer_Type %in% natype))
  
  data <- suppressMessages(data %>% mutate(Cancer_Type=factor(Cancer_Type,levels = mdata$Cancer_Type)) %>% arrange(Cancer_Type,Burden) %>% left_join(mdata))
  
  data <- data %>% filter(!is.na(Burden))%>% group_by(Cancer_Type) %>% mutate(Subseq=seq_along(Sample)) %>% ungroup() %>% mutate(score=Seq+1+8*Subseq/nsample)
  mdata <- mdata %>% mutate(Seq2=Seq+10) %>% mutate(score=(Seq+Seq2)/2) %>% 
    mutate(Burden=median,Type=if_else(seq_along(Cancer_Type) %% 2 ==0,"A","B"),Label=paste0(nsample,'<br /> - <br />',total)) %>% 
    mutate(Label_col = glue::glue("<span style='color:#377eb8'>{nsample}</span><br />\u2015<br /><span style='color:#4daf4a'>{total}</span>"))
  # "<span style='color:#377eb8'>{nsample}</span><br />\u2501<br /><span style='color:#4daf4a'>{total}</span>"
  f <- function(y) seq(floor(min(y,na.rm = TRUE)), ceiling(max(y,na.rm = TRUE)))
  
  ymint <- floor(min(data$Burden,na.rm = TRUE))
  ymint <- if_else(is.infinite(ymint),-4,ymint)
  ymaxt <- ceiling(max(data$Burden,na.rm = TRUE))
  
  print(data)
  p <-  data %>% filter(!is.na(Burden)) %>% 
    ggplot(aes(score,(Burden),group=Cancer_Type))+
    geom_rect(data=mdata,aes(xmin=Seq,xmax=Seq2,ymin=-Inf,ymax=Inf,fill=Type),alpha=0.5)+
    geom_point(aes(col=Sample_Group), pch=21,size=1,stroke=0.5)+
    geom_segment(data=mdata,aes(x=Seq+1,y=(median),xend=Seq2-1,yend=(median)),col="#542788")+#"#e41a1c"
    scale_fill_manual(values=c('#F2F2F2','#D4D4D4'))+
    scale_color_manual(values = sgroupcol)+
    scale_x_continuous(limits = c(min(mdata$Seq),max(mdata$Seq2)),expand = c(0,0),breaks = mdata$score,labels = mdata$Cancer_Type,sec.axis = dup_axis(labels = mdata$Label_col))+
    scale_y_continuous(limits = c(ymint,ymaxt),expand = c(0,0),breaks = f)+
    labs(x="",y=ylab)+
    theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 10,axis_text_size = 8,grid = "Y",ticks = TRUE,plot_margin=margin(5.5,5.5,5.5,5.5))+
    theme(legend.position = "None",axis.text.x.bottom = element_text(size = 8,angle = 45,hjust = 0.5,vjust =0.5 ),axis.text.x.top = element_markdown(size = 8),panel.border = element_rect(linewidth = 0.2,colour = "black",fill = NA))+
    coord_cartesian(clip = 'off')
  
  if(!is.null(addnote)){
    p <- p+annotate("text",x=-Inf, y = Inf, label = addnote, vjust=2, hjust=-0.2,size=4)
  }
  
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 2.5+length(unique(mdata$Cancer_Type))*0.5
    xleng <- if_else(xleng>18,18,xleng)
    yleng <- 6
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
  
}





decompsite_distribution <- function(decompsite,output_plot = NULL,plot_width=NULL, plot_height=NULL){
  #Cancer_Type   Sample Total_Mutations Cosine_similarity .... 
  fealist <- c('Cancer_Type','Sample','Cosine_similarity','100-L1_Norm_%','100-L2_Norm_%','KL_Divergence','Correlation')
  decompsite2 <- decompsite %>% select(one_of(fealist)) %>% pivot_longer(cols = -c(Cancer_Type,Sample))
  mtmp <- decompsite %>% group_by(Cancer_Type) %>% summarise(m=median(Cosine_similarity,na.rm = TRUE)) 
  mtmp2 <- decompsite2 %>% group_by(Cancer_Type,name) %>% summarise(m=median(value,na.rm = TRUE)) %>% ungroup() %>% mutate(m=if_else(name=="KL_Divergence",1-m,m))%>% group_by(name) %>% arrange(desc(m)) %>% mutate(Seq=seq_along(Cancer_Type)) %>% ungroup()
  ntmp <- decompsite %>% count(Cancer_Type) %>% mutate(Caner_type2=paste0(Cancer_Type," (",n,")")) %>% left_join(mtmp) %>% arrange(m)
  
  library(ggridges)
  library(scales)
  p <- decompsite2 %>% 
    left_join(mtmp2) %>% 
    left_join(ntmp) %>% 
    mutate(Cancer_Type=factor(Cancer_Type,levels = ntmp$Cancer_Type,labels =ntmp$Caner_type2 )) %>% 
    mutate(name=factor(name,levels = c('Cosine_similarity','100-L1_Norm_%','100-L2_Norm_%','KL_Divergence','Correlation'))) %>% 
    mutate(value=if_else(value<0,0,value)) %>% 
    ggplot(aes(value,y=Cancer_Type,fill=Cancer_Type,height = ..ndensity..))+
    geom_density_ridges(color="black")+
    facet_wrap(~name,scales = 'free_x',nrow = 1)+
    theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = "X",ticks = TRUE,base_size = 12,grid = FALSE,strip_text_size = 12)+
    scale_x_continuous(name="",expand = c(0,0),breaks = pretty_breaks())+
    scale_y_discrete(name = "",expand = c(0,0))+
    #scale_fill_viridis_c(direction = -1)+
    scale_fill_aaas()+
    theme(strip.text.x = element_text(hjust = 0.5),legend.position = "none",panel.grid.major.x = element_line(colour = gray(0.85),size = 0.5,linetype = 2),panel.spacing = unit(0.8, "lines"),strip.background = element_rect(colour = '#cccccc',fill = c('#a1d99b')))+
    #panel_border(color = "black",size = 0.8)+
    guides(fill=guide_legend(title = "",nrow = 4))
  
  if(is.null(output_plot)){
    return(p)
  }else{
    yleng <- 2.5+length(unique(ntmp$Cancer_Type))*0.5
    yleng <- if_else(yleng>12,12,yleng)
    xleng <- 16
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
  
  
}


### Landscape of mutational signature cluster function ####
Exposure_Clustering <- function(sigdata,sigcolor=NULL,studydata=NULL,studydata_cat = TRUE, studycolor=NULL,study_cutoff = -Inf,studyname = NULL,puritydata=NULL,puritydata_cat=FALSE,puritycol=NULL,purity_cutoff=-Inf,purityname = NULL,cosinedata=NULL,cosine_cutoff=0,highlight=NULL,legendnrow=NULL,sampletext_size=6,output_plot=NULL,plot_height=NULL,plot_width=NULL,clustern=NULL,kcolors=NULL,hc_func='hclust',hc_metric = 'euclidean',hc_method = 'ward.D2',stand=TRUE){
  require(tidyverse)
  require(scales)
  require(janitor)
  require(ggsci)
  require(factoextra)
  require(cowplot)
  require(ggpubr)
  require(hrbrthemes)
  
  # remove signature name with 0 contribution
  sigdata <- sigdata %>% select(where(~ is.character(.x) || sum(.x) !=0 ))
  
  
  # define color for categoly variable
  colset <- c(pal_npg()(10),pal_jama()(7),pal_igv()(51))
  #colset <- pal_igv()(51)
  colnames(sigdata)[1] <- 'Samples'
  
  ## define color for signature
  if(is.null(sigcolor)){
    uvalues <- colnames(sigdata)[-1] #sort
    if(unique(uvalues %in% names(Subscolor)) == TRUE && length(unique(uvalues %in% names(Subscolor))) == 1){
      sigcolorindex <- as.character(Subscolor[uvalues])
      names(sigcolorindex) <- uvalues
    } else if(unique(uvalues %in% names(SBScolor)) == TRUE && length(unique(uvalues %in% names(SBScolor))) == 1){
      sigcolorindex <- as.character(SBScolor[uvalues])
      names(sigcolorindex) <- uvalues
    } else {
      sigcolorindex <- as.character(SBScolor[1:length(uvalues)])
      names(sigcolorindex) <- uvalues
    }
  }else{
    sigcolorindex <- sigcolor
  }
  
  # cluster data
  tmp=sigdata %>% adorn_percentages('row') 
  #%>% mutate_if(is.numeric, funs(replace_na(., 0)))
  mdata=as.matrix(tmp[,-1])
  rownames(mdata) <- tmp$Samples
  
  #fviz_nbclust(mdata, kmeans, method = "gap_stat")
  if(is.null(clustern)){
    clustern=2
    kcolors="black"
  }else{
    clustern <- if_else(dim(mdata)[1]<clustern+5,2L,as.integer(clustern))
    if(is.null(kcolors)) {kcolors="black"}
  }
  
  res <- hcut(mdata,k = clustern,hc_func = hc_func,hc_metric = hc_metric,hc_method = hc_method,stand=stand)
  
  # p_cluster
  p_cluster <- fviz_dend(res, rect = TRUE, cex = 0.5,k_colors = kcolors,lwd = 0.3,show_labels = F)+scale_x_discrete(expand = c(0,0))+ theme(plot.margin=margin(b=-0.7,unit="cm"),title = element_blank())
  
  ## p_mutation
  p_mutation <- sigdata %>% gather(Signature,Weight,-Samples) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,(Weight),fill=factor(Signature,levels = names(sigcolorindex))))+geom_bar(stat="identity",col="gray95",width=1,size=0.1)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = "xy")+theme(legend.title = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_blank(),panel.grid.major.x=element_blank(),legend.position = "none",plot.margin=margin(b=-1,t = 1,unit="cm"))+scale_y_continuous(expand = c(0, 0),breaks = pretty_breaks(),labels = comma)+scale_fill_manual(values = sigcolorindex,drop=FALSE)+guides(fill=guide_legend(nrow=2,byrow=TRUE))+xlab("")+ylab("Number of mutations \n")
  
  # p_proportion plot
  if(is.null(legendnrow)){
    nsize <- dim(sigdata)[2]-1
    if(nsize > 20){
      legendnrow <- 2
    }else {
      legendnrow <- 1
    }
  }
  
  
  # define the legend size
  legend_size_cat <- theme(legend.box.background = element_blank(),legend.key.size = unit(0.5, "cm"), legend.key.height = unit(0.3, "cm"), legend.key.width =unit(0.4, "cm"), legend.position=c(0,1),legend.justification = c(0,1))
  legend_size_num <-  theme(legend.box.background = element_blank(),legend.key.width =unit(0.3,"cm"),legend.key.height = unit(1,"cm"),legend.position=c(0,1),legend.justification = c(0,1))
  
  p_proportion <- sigdata %>% gather(Signature,Weight,-Samples) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,Weight,fill=factor(Signature,levels = names(sigcolorindex))))+geom_bar(stat="identity",position="fill",col="gray95",width = 1,size=0.1)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,grid = FALSE)+theme(panel.background = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_text(size = sampletext_size, angle = 90, vjust = 0.5, hjust = 1),panel.grid.major=element_line(),axis.ticks.y = element_line(colour = "black"))+scale_y_continuous(expand = c(0, 0),breaks = pretty_breaks())+xlab("")+scale_fill_manual("Sigantures",values = sigcolorindex,drop=FALSE)+ylab("Signature contribution\n")+guides(fill=guide_legend(ncol=legendnrow,byrow=FALSE))+legend_size_cat
  
  p_proportion_legend <- get_legend(p_proportion)
  p_proportion <- p_proportion + theme(legend.position = "none")
  
  if(!is.null(highlight)){
    samhigh <- sigdata %>% mutate(Samples_high=if_else(Samples %in% highlight,paste0("*",Samples),Samples))
    p_proportion <-  p_proportion+scale_x_discrete(breaks=samhigh$Samples,labels=samhigh$Samples_high)
  }
  
  if(!is.null(cosinedata)){
    p_cosine <- cosinedata %>%  mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% mutate(Similarity=if_else(Similarity<cosine_cutoff,NA_real_,Similarity)) %>% ggplot(aes(Samples,1,fill=Similarity))+geom_tile(col="black")+scale_fill_viridis_c("Cosine\nSimilarity",na.value = "#cccccc",option = 'C',limits = c(0.6, 1), oob = scales::squish)+theme_minimal()+theme(panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+ylim(c(0,2))+legend_size_num
    p_cosine_legend <- get_legend(p_cosine)
    p_cosine <- p_cosine + theme(legend.position = "none",plot.margin=margin(b=-1.2,t = 0.5,unit="cm"))
    
  }else{
    p_cosine <- NULL
    p_cosine_legend <- NULL
  }
  
  
  # study color bar
  if(!is.null(studydata)){
    colnames(studydata) <- c('Samples','Study')
    studydata <- studydata %>% filter(Samples %in% res$labels) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) 
    studyname <- if_else(is.null(studyname),"",studyname)
    if(!studydata_cat) {
      #studyname <- paste0(studyname,"\n")
      studydata <- studydata %>% mutate(Study=if_else(Study<study_cutoff,NA_real_,Study))
      p_study <- studydata %>% ggplot(aes(Samples,1,fill=Study))+geom_tile(col="black")+scale_fill_viridis_c(studyname,na.value = "#cccccc",option = 'D',breaks=pretty_breaks(5))+theme_minimal()+theme(panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+ylim(c(0,2))+legend_size_num
      p_study_legend <- get_legend(p_study)
      p_study <- p_study + theme(legend.position = "none",plot.margin=margin(b=-1,t = 0.5,unit="cm"))
      
    }else { 
      if(is.null(studycolor)){
        studycolor <-  colset[1:length(unique(studydata$Study))]
        names(studycolor) <- unique(studydata$Study)
      }
      p_study <- studydata %>% ggplot(aes(Samples,1,fill=factor(Study,levels = names(studycolor))))+geom_tile(col="black")+scale_fill_manual(studyname,values =studycolor)+theme_minimal()+theme(panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+ylim(c(0,2))+legend_size_cat
      p_study_legend <- get_legend(p_study)
      p_study <- p_study + theme(legend.position = "none",plot.margin=margin(b=-1,t=0.5,unit="cm"))
    }
  }else {
    p_study <- NULL
    p_study_legend <- NULL
  }
  
  if(!is.null(puritydata)){
    colnames(puritydata) <- c('Samples','Purity')
    puritydata <- puritydata %>% filter(Samples %in% res$labels) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) 
    purityname <- if_else(is.null(purityname),"",purityname)
    if(!puritydata_cat) {
      #purityname <-  paste0(purityname,"\n")
      puritydata <- puritydata %>% mutate(Purity=if_else(Purity<purity_cutoff,NA_real_,Purity))
      p_purity <- puritydata %>% ggplot(aes(Samples,1,fill=Purity))+geom_tile(col="black")+scale_fill_viridis_c(purityname,na.value = "#cccccc",option = 'D',breaks=pretty_breaks(5))+theme_minimal()+theme(panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+ylim(c(0,2))+legend_size_num
      p_purity_legend <- get_legend(p_purity)
      p_purity <- p_purity + theme(legend.position = "none",plot.margin=margin(b=-1,t = 0.5,unit="cm"))
    }else {
      #purityname <-  paste0(purityname,"\n")
      if(is.null(puritycol)){
        puritycol <- colset[1:length(unique(puritydata$Purity))]
        names(puritycol) <- unique(puritydata$Purity)
      }
      p_purity <- puritydata %>% ggplot(aes(Samples,1,fill=factor(Purity,levels = names(puritycol))))+geom_tile(col="black")+scale_fill_manual(purityname,values =puritycol)+theme_minimal()+theme(panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+ylim(c(0,2))+legend_size_cat
      p_purity_legend <- get_legend(p_purity)
      p_purity <- p_purity + theme(legend.position = "none",plot.margin=margin(b=-1,t=0.5,unit="cm"))
      
    }
  }else {
    p_purity <- NULL
    p_purity_legend <- NULL
  }
  
  h_study <- if_else(is.null(p_study),0,0.14)
  h_purity <- if_else(is.null(p_study),0,0.11)
  
  xleng <- 2.5+dim(sigdata)[1]*0.1
  xleng <- if_else(xleng>20,20,xleng)
  yleng <- 12
  if(is.null(plot_width)){ plot_width <-  xleng}
  if(is.null(plot_height)){ plot_height <-  yleng}
  
  
  p_legends <- plot_grid(p_study_legend,p_purity_legend, p_cosine_legend, p_proportion_legend, align = 'v',axis = 'l',ncol = 1,rel_heights = c(1,1,1,1.5))
  p_main <- plot_grid(p_cluster,p_study,p_purity,p_mutation,p_cosine,p_proportion,align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2,h_study,h_purity,4,0.05,7))
  p_all <- plot_grid(p_main,p_legends+theme(plot.margin = margin(l = -0.5,unit = 'cm')),ncol = 2,rel_widths = c(10,2))
  
  if(is.null(output_plot)){
    return(pall)
  }else{
    ggsave(filename = output_plot,plot = p_all,width = plot_width,height = plot_height)
  }
  
}


### Landscape of mutational signature cluster function ####
Exposure_Clustering_old <- function(sigdata,sigcolor=NULL,studydata=NULL,studydata_cat = TRUE, studycolor=NULL,study_cutoff = -Inf,studyname = NULL,puritydata=NULL,puritydata_cat=FALSE,puritycol=NULL,purity_cutoff=-Inf,purityname = NULL,cosinedata=NULL,cosine_cutoff=0,highlight=NULL,legendnrow=NULL,sampletext_size=6,output_plot=NULL,plot_height=NULL,plot_width=NULL,clustern=NULL,kcolors=NULL,hc_func='hclust',hc_metric = 'euclidean',hc_method = 'ward.D2',stand=TRUE){
  require(tidyverse)
  require(scales)
  require(janitor)
  require(ggsci)
  require(ggpubr)
  require(factoextra)
  require(cowplot)
  
  # define color for categoly variable
  colset <- c(pal_npg()(10),pal_jama()(7),pal_igv()(51))
  #colset <- pal_igv()(51)
  
  colnames(sigdata)[1] <- 'Samples'
  
  ## define color for signature
  if(is.null(sigcolor)){
    uvalues <- colnames(sigdata)[-1] #sort
    if(unique(uvalues %in% names(Subscolor)) == TRUE && length(unique(uvalues %in% names(Subscolor))) == 1){
      sigcolorindex <- as.character(Subscolor[uvalues])
      names(sigcolorindex) <- uvalues
    } else if(unique(uvalues %in% names(SBScolor)) == TRUE && length(unique(uvalues %in% names(SBScolor))) == 1){
      sigcolorindex <- as.character(SBScolor[uvalues])
      names(sigcolorindex) <- uvalues
    } else {
      sigcolorindex <- as.character(SBScolor[1:length(uvalues)])
      names(sigcolorindex) <- uvalues
    }
  }else{
    sigcolorindex <- sigcolor
  }
  
  # cluster data
  tmp=sigdata %>% adorn_percentages('row') 
  #%>% mutate_if(is.numeric, funs(replace_na(., 0)))
  mdata=as.matrix(tmp[,-1])
  rownames(mdata) <- tmp$Samples
  
  #fviz_nbclust(mdata, kmeans, method = "gap_stat")
  if(is.null(clustern)){
    clustern=2
    kcolors="black"
  }else{
    clustern <- if_else(dim(mdata)[1]<clustern+5,2L,as.integer(clustern))
    if(is.null(kcolors)) {kcolors="black"}
  }
  
  res <- hcut(mdata,k = clustern,hc_func = hc_func,hc_metric = hc_metric,hc_method = hc_method,stand=stand)
  
  # p_cluster
  p_cluster <- fviz_dend(res, rect = TRUE, cex = 0.5,k_colors = kcolors,lwd = 0.3,show_labels = F)+scale_x_discrete(expand = c(0,0))+ theme(plot.margin=margin(b=-0.7,unit="cm"),title = element_blank())
  #plot.margin=margin(b=-1,unit="cm")
  
  ## p_mutation
  p_mutation <- sigdata %>% gather(Signature,Weight,-Samples) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,(Weight),fill=factor(Signature,levels = names(sigcolorindex))))+geom_bar(stat="identity",col="gray95",width=1,size=0.1)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,axis = "xy")+theme(legend.title = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_blank(),panel.grid.major.x=element_blank(),legend.position = "none",legend.box.spacing = unit(0,"cm"),plot.margin=margin(b=-1,t = 1,unit="cm"))+scale_y_continuous(expand = c(0, 0),breaks = pretty_breaks(),labels = comma)+scale_fill_manual(values = sigcolorindex,drop=FALSE)+guides(fill=guide_legend(nrow=2,byrow=TRUE))+xlab("")+ylab("Number of mutations \n")
  #+theme(plot.margin=margin(b=4,unit="pt"))
  
  # p_proportion plot
  if(is.null(legendnrow)){
    nsize <- dim(sigdata)[2]-1
    if(nsize > 15){
      legendnrow <- 2
    }else {
      legendnrow <- 1
    }
  }
  
  p_proportion <- sigdata %>% gather(Signature,Weight,-Samples) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% ggplot(aes(Samples,Weight,fill=factor(Signature,levels = names(sigcolorindex))))+geom_bar(stat="identity",position="fill",col="gray95",width = 1,size=0.1)+theme_ipsum_rc(axis_title_just = "m",axis_title_size = 12,grid = FALSE)+theme(panel.background = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_text(size = sampletext_size, angle = 90, vjust = 0.5, hjust = 1),panel.grid.major=element_line(),legend.position = "bottom",legend.box.background = element_blank(),legend.box.spacing = unit(-0.5,"cm"),legend.key = element_rect(size = 0),axis.ticks.y = element_line(colour = "black"),legend.key.size = unit(0.25, "cm"),legend.key.width =unit(1.5, "cm"))+scale_y_continuous(expand = c(0, 0),breaks = pretty_breaks())+xlab("")+scale_fill_manual("Mutational sigantures",values = sigcolorindex,drop=FALSE)+guides(fill=guide_legend(nrow=legendnrow,byrow=TRUE,label.position = "bottom"))+ylab("Signature contribution\n")
  #+theme(plot.margin=margin(t=4,unit="pt"))
  #legend.title = element_blank(),
  
  p_proportion_legend <- as_ggplot(get_legend(p_proportion))
  #p_proportion_legend <- p_proportion_legend+theme(plot.margin = margin(b = 0))
  p_proportion <- p_proportion + theme(legend.position = "none")
  
  if(!is.null(highlight)){
    samhigh <- sigdata %>% mutate(Samples_high=if_else(Samples %in% highlight,paste0("*",Samples),Samples))
    p_proportion <-  p_proportion+scale_x_discrete(breaks=samhigh$Samples,labels=samhigh$Samples_high)
  }
  
  
  #p3 <- flush_ticks(p3,flush = "Y",plot = FALSE)
  #p4 <- flush_ticks(p4,flush = "Y",plot = FALSE)
  
  if(!is.null(cosinedata)){
    # cosine similarity bar
    p_cosine <- cosinedata %>%  mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% mutate(Similarity=if_else(Similarity<cosine_cutoff,NA_real_,Similarity)) %>% ggplot(aes(Samples,1,fill=Similarity))+geom_tile(col="black")+scale_fill_viridis_c("Cosine similarity\n",na.value = "#cccccc",option = 'C',limits = c(0.6, 1), oob = scales::squish)+theme_minimal()+theme(legend.position = "bottom",legend.key.width =unit(2,"cm"),legend.key.height = unit(0.3,"cm"), panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+ylim(c(0,2))
    #+theme(title = element_blank())
    ## p_cosine_legend = cosine similarity bar (legend)
    p_cosine_legend <- as_ggplot(get_legend(p_cosine))
    #p_cosine_legend <- p_cosine_legend+theme(plot.margin = margin(b = 0))
    p_cosine <- cosinedata %>%  mutate(Samples=factor(Samples,levels=res$labels[res$order])) %>% mutate(Similarity=if_else(Similarity<cosine_cutoff,NA_real_,Similarity)) %>% ggplot(aes(Samples,1,fill=Similarity))+geom_tile(col="black")+scale_fill_viridis_c(na.value = "#cccccc",option = 'C',limits = c(0.6, 1), oob = scales::squish)+theme_minimal()+theme(legend.position = "none", panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1.2,t = 0.5,unit="cm"))+ylim(c(0,2))
    #,title = element_blank()
  }else{
    p_cosine <- NULL
    p_cosine_legend <- NULL
  }
  
  legend_size <- theme(legend.position = "bottom",legend.direction = "horizontal", legend.box.background = element_blank(),legend.key = element_rect(size = 0),legend.key.size = unit(0.25, "cm"),legend.key.width =unit(0.6, "cm"))
  #legend.box.spacing = unit(-0.5,"cm"),
  
  # study color bar
  if(!is.null(studydata)){
    colnames(studydata) <- c('Samples','Study')
    studydata <- studydata %>% filter(Samples %in% res$labels) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) 
    studyname <- if_else(is.null(studyname),"",studyname)
    if(!studydata_cat) {
      studyname <- paste0(studyname,"\n")
      studydata <- studydata %>% mutate(Study=if_else(Study<study_cutoff,NA_real_,Study))
      p_study <- studydata  %>% ggplot(aes(Samples,1,fill=Study))+geom_tile(col="black")+scale_fill_viridis_c(studyname,na.value = "#cccccc",option = 'D',breaks=pretty_breaks(5))+theme_minimal()+theme(legend.position = "bottom",legend.key.width =unit(2, "cm"),legend.key.height = unit(0.3,"cm"), panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+ylim(c(0,2))
      p_study_legend <- as_ggplot(get_legend(p_study+legend_size))
      #p_study_legend <- p_study_legend+theme(plot.margin = margin(b = 0))
      p_study <- studydata %>% ggplot(aes(Samples,1,fill=Study))+geom_tile(col="black")+scale_fill_viridis_c(studyname,na.value = "#cccccc",option = 'D',breaks=pretty_breaks(5))+theme_minimal()+theme(legend.position = "none", panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t = 0.5,unit="cm"))+ylim(c(0,2))
      
    }else { 
      if(is.null(studycolor)){
        studycolor <-  colset[1:length(unique(studydata$Study))]
        names(studycolor) <- unique(studydata$Study)
      }
      p_study <- studydata %>% ggplot(aes(Samples,1,fill=factor(Study,levels = names(studycolor))))+geom_tile(col="black")+scale_fill_manual(studyname,values =studycolor)+theme_minimal()+theme(legend.position = "bottom",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t=0.5,unit="cm"))+ylim(c(0,2))+guides(fill = guide_legend(nrow = 1))
      p_study_legend <- as_ggplot(get_legend(p_study+legend_size))
      #p_study_legend <- p_study_legend+theme(plot.margin = margin(b = 0))
      p_study <- studydata %>% ggplot(aes(Samples,1,fill=factor(Study,levels = names(studycolor))))+geom_tile(col="black")+scale_fill_manual(studyname,values =studycolor)+theme_minimal()+theme(legend.position = "none",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t=0.5,unit="cm"))+ylim(c(0,2))
    }
  }else {
    p_study <- NULL
    p_study_legend <- NULL
  }
  
  # purity color bar
  # if(is.null(purityname)){
  #   purityname <-  "              "
  # }
  
  if(!is.null(puritydata)){
    colnames(puritydata) <- c('Samples','Purity')
    puritydata <- puritydata %>% filter(Samples %in% res$labels) %>% mutate(Samples=factor(Samples,levels=res$labels[res$order])) 
    purityname <- if_else(is.null(purityname),"",purityname)
    if(!puritydata_cat) {
      purityname <-  paste0(purityname,"\n")
      puritydata <- puritydata %>% mutate(Purity=if_else(Purity<purity_cutoff,NA_real_,Purity))
      p_purity <- puritydata  %>% ggplot(aes(Samples,1,fill=Purity))+geom_tile(col="black")+scale_fill_viridis_c(purityname,na.value = "#cccccc",option = 'D',breaks=pretty_breaks(5))+theme_minimal()+theme(legend.position = "bottom",legend.key.width =unit(2, "cm"),legend.key.height = unit(0.3,"cm"), panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+ylim(c(0,2))
      p_purity_legend <- as_ggplot(get_legend(p_purity+legend_size))
      #p_purity_legend <- p_purity_legend+theme(plot.margin = margin(b = 0))
      p_purity <- puritydata %>% ggplot(aes(Samples,1,fill=Purity))+geom_tile(col="black")+scale_fill_viridis_c(purityname,na.value = "#cccccc",option = 'D',breaks=pretty_breaks(5))+theme_minimal()+theme(legend.position = "none", panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t = 0.5,unit="cm"))+ylim(c(0,2))
    }else {
      purityname <-  paste0(purityname,"\n")
      if(is.null(puritycol)){
        puritycol <- colset[1:length(unique(puritydata$Purity))]
        names(puritycol) <- unique(puritydata$Purity)
      }
      p_purity <- puritydata %>% ggplot(aes(Samples,1,fill=factor(Purity,levels = names(puritycol))))+geom_tile(col="black")+scale_fill_manual(purityname,values =puritycol)+theme_minimal()+theme(legend.position = "bottom",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t=0.5,unit="cm"))+ylim(c(0,2))+guides(fill = guide_legend(nrow = 1))
      p_purity_legend <- as_ggplot(get_legend(p_purity+legend_size))
      #p_purity_legend <- p_purity_legend+theme(plot.margin = margin(b = 0))
      p_purity <- puritydata %>% ggplot(aes(Samples,1,fill=factor(Purity,levels = names(puritycol))))+geom_tile(col="black")+scale_fill_manual(purityname,values =puritycol)+theme_minimal()+theme(legend.position = "none",panel.background = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),panel.grid = element_blank())+theme(plot.margin=margin(b=-1,t=0.5,unit="cm"))+ylim(c(0,2))
    }
  }else {
    p_purity <- NULL
    p_purity_legend <- NULL
  }
  
  h_study <- if_else(is.null(p_study),0,0.14)
  h_purity <- if_else(is.null(p_study),0,0.11)
  
  p_legends <- plot_grid(p_study_legend+theme(plot.margin = margin(b = -2,t=0,unit = 'cm')), p_purity_legend+theme(plot.margin = margin(b = -2,t=0,unit = 'cm')), p_proportion_legend+theme(plot.margin = margin(t = -1,b=-0.5,unit = 'cm')), p_cosine_legend+theme(plot.margin = margin(t = -1,b=-0.5,unit = 'cm')),nrow = 2,ncol = 2)
  
  if(!is.null(p_purity) & !is.null(p_study) ){
    pall <- plot_grid(p_cluster,p_study,p_purity,p_mutation,p_cosine,p_proportion,p_legends+theme(plot.margin = margin(t=-3,unit = 'cm')),align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2,h_study,h_purity,4,0.05,7,1))
  }else if(is.null(p_purity) & is.null(p_study)){
    p_legends <- plot_grid(p_proportion_legend+theme(plot.margin = margin(t = 2,b=-0.5,unit = 'cm')), p_cosine_legend+theme(plot.margin = margin(t = -1,b=-0.5,unit = 'cm')),nrow=2,ncol = 1)
    pall <- plot_grid(p_cluster,p_mutation,p_cosine,p_proportion,p_legends+theme(plot.margin = margin(t=-3,unit = 'cm')),align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2,4,0.05,7,1))
  } else if(is.null(p_purity)){
    pall <- plot_grid(p_cluster,p_study,p_mutation,p_cosine,p_proportion,p_legends+theme(plot.margin = margin(t=-3,unit = 'cm')),align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2,h_study,4,0.05,7,1))
  }else{
    pall <- plot_grid(p_cluster,p_purity,p_mutation,p_cosine,p_proportion,p_legends+theme(plot.margin = margin(t=-3,unit = 'cm')),align = 'v',axis = 'lr',ncol = 1,rel_heights = c(2,h_purity,4,0.05,7,1))
  }
  
  if(is.null(output_plot)){
    return(pall)
  }else{
    xleng <- 2.5+dim(sigdata)[1]*0.1
    xleng <- if_else(xleng>20,20,xleng)
    yleng <- 12
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = pall,width = plot_width,height = plot_height)
  }
  
}


## piechar common
piechart_plot <- function(data,colset=NULL,keep_legend = TRUE,legend_name=NULL, output_plot = NULL,plot_width=NULL, plot_height=NULL){
  
  ## data format: Type,Catelogy,Freq,Label
  colnames(data) <- c("Type","Catelogy",'Frequency','Label')
  
  if(is.null(colset)){
    uvalues <- unique(data$Catelogy)
    if((unique(uvalues %in% names(Subscolor))) == TRUE && length(unique(uvalues %in% names(Subscolor))) == 1){
      sigcolorindex <- as.character(Subscolor[uvalues])
      names(sigcolorindex) <- uvalues
    } else if((unique(uvalues %in% names(SBScolor))) == TRUE && length(unique(uvalues %in% names(SBScolor))) == 1){
      sigcolorindex <- as.character(SBScolor[uvalues])
      names(sigcolorindex) <- uvalues
    } else {
      sigcolorindex <- as.character(SBScolor[1:length(uvalues)])
      names(sigcolorindex) <- uvalues
    }
  }else{
    sigcolorindex <- colset
  }
  
  
  # convert nsig_data to another format for the piechart
  data_pie <- data %>%
    group_by(Type) %>% 
    arrange(Frequency) %>%
    mutate(
      end_angle = 2*pi*cumsum(Frequency)/sum(Frequency),   # ending angle for each pie slice
      start_angle = lag(end_angle, default = 0),   # starting angle for each pie slice
      mid_angle = 0.5*(start_angle + end_angle),   # middle of each pie slice, for the text label
      # horizontal and vertical justifications depend on whether we're to the left/right
      # or top/bottom of the pie
      hjust = ifelse(mid_angle > pi, 1, 0),
      vjust = ifelse(mid_angle < pi/2 | mid_angle > 3*pi/2, 0, 1)
    ) %>% 
    ungroup()
  
  # radius of the pie and radius for outside and inside labels
  rpie <- 1
  rlabel_out <- 1.05 * rpie
  rlabel_in <- 0.6 * rpie
  if(is.null(legend_name)){
    legend_name <- ''
  }
  #legend_name <- if_else(is.null(legend_name),'',legend_name)
  
  p <- ggplot(data_pie) +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = rpie, start = start_angle, end = end_angle, fill = Catelogy)) +
    geom_text(aes(x = rlabel_in * sin(mid_angle), y = rlabel_in * cos(mid_angle), label = Label ), size = 4)+
    facet_wrap(~Type,nrow = 2)+
    coord_fixed()+
    labs(x="",y="")+
    scale_fill_manual(legend_name,values = sigcolorindex)+
    theme_ipsum_rc(axis = FALSE, grid = FALSE)+
    theme(axis.text.y = element_blank(),axis.text.x = element_blank(),strip.text = element_text(size = 14,hjust = 0.5,face = "bold"))
  
  if(!(keep_legend)){
    p <- p+theme(legend.position = 'none')
  }
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 7
    yleng <- 7
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
  
} 



## bar plot 
barchart_plot <- function(data,colset=NULL,keep_legend = TRUE,legend_name=NULL, output_plot = NULL,plot_width=NULL, plot_height=NULL){
  
  ## data format: Type,Catelogy,Freq,Label
  colnames(data) <- c("Type","Catelogy",'Frequency')
  
  if(is.null(colset)){
    uvalues <- unique(data$Catelogy)
    if((unique(uvalues %in% names(Subscolor))) == TRUE && length(unique(uvalues %in% names(Subscolor))) == 1){
      sigcolorindex <- as.character(Subscolor[uvalues])
      names(sigcolorindex) <- uvalues
    } else if((unique(uvalues %in% names(SBScolor))) == TRUE && length(unique(uvalues %in% names(SBScolor))) == 1){
      sigcolorindex <- as.character(SBScolor[uvalues])
      names(sigcolorindex) <- uvalues
    } else {
      sigcolorindex <- as.character(SBScolor[1:length(uvalues)])
      names(sigcolorindex) <- uvalues
    }
  }else{
    sigcolorindex <- colset
  }
  legend_name <- if_else(is.null(legend_name),'',legend_name)
  
  p <- data %>% 
    mutate(Catelogy=fct_reorder(Catelogy,Frequency,.desc = TRUE)) %>% 
    ggplot(aes(Catelogy,Frequency,fill=Catelogy))+
    geom_col(width = 0.8)+
    geom_text(aes(label = percent(Frequency,accuracy = 0.1)), vjust = -.5,size=3.5)+
    facet_wrap(~Type,nrow = 2)+
    labs(x="",y="Frequency (%)")+
    scale_y_percent(breaks = pretty_breaks(),limits = c(0,1),expand = expand_scale(mult = 0.1))+
    scale_fill_manual(legend_name,values = sigcolorindex,breaks = names(sigcolorindex))+
    theme_ipsum_rc(axis_title_size = 12,axis_title_just = "m",axis = TRUE, grid = "Yy")+
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),strip.text = element_text(size = 14,hjust = 0.5,face = "bold"),axis.line.x = element_line(colour = 'black',size = 0.25),axis.line.y = element_line(colour = 'black',size = 0.25),legend.position = 'bottom')+
    guides(fill=guide_legend(nrow=1,byrow=TRUE))
  
  
  if(!(keep_legend)){
    p <- p+theme(legend.position = 'none')
  }
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 2.5+length(uvalues)*0.5
    xleng <- if_else(xleng>15,15,xleng)
    yleng <- 7
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
  
} 

barchart_plot2 <- function(data, output_plot = NULL,plot_width=NULL, plot_height=NULL){
  
  ## data format: Type,Catelogy,Freq,Label
  colnames(data) <- c("Type","Catelogy",'Frequency')
  uvalues <- sort(unique(data$Catelogy))
  p <- data %>% 
    mutate(Catelogy=fct_reorder(Catelogy,Frequency,.desc = TRUE)) %>% 
    ggplot(aes(Catelogy,Frequency,fill=Frequency))+
    geom_col(width = 0.8)+
    #geom_text(aes(label = percent(Frequency,accuracy = 0.1)), vjust = -.5,size=3.5)+
    facet_wrap(~Type,nrow = 2)+
    labs(x="",y="Frequency")+
    scale_y_continuous(breaks = pretty_breaks(),expand = expansion(mult = c(0,0.1)))+
    #scale_fill_manual(legend_name,values = sigcolorindex,breaks = names(sigcolorindex))+
    theme_ipsum_rc(axis_title_size = 12,axis_title_just = "m",axis = TRUE, grid = "Yy",)+
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size = 8),strip.text = element_text(size = 14,hjust = 0.5,face = "bold"),axis.line.x = element_line(colour = 'black',size = 0.25),axis.line.y = element_line(colour = 'black',size = 0.25),legend.position = 'none')
  #    guides(fill=guide_legend(nrow=1,byrow=TRUE))
  
  if(is.null(output_plot)){
    return(p)
  }else{
    xleng <- 2.5+length(uvalues)*0.5
    xleng <- if_else(xleng>15,15,xleng)
    yleng <- 4
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
  
} 



prevalence_plot <- function(sigdata,nmutation = 0, legend_name="Sigantures", colset = NULL, output_plot = NULL,plot_width=NULL, plot_height=NULL,rel_widths=NULL){
  require(janitor)
  pie_input <- sigdata %>% 
    summarise(across(where(is.numeric),~sum(.x,na.rm=TRUE))) %>% 
    mutate(Type="Prevalence by mutations") %>%
    select(Type,everything()) %>%
    adorn_percentages(denominator = 'row') %>% 
    pivot_longer(cols = -Type) %>% 
    mutate(lab=if_else(value>0.05,percent(value,accuracy = 0.1),''))
  
  nfilter <- sigdata %>% pivot_longer(cols=-Samples) %>% filter(value>nmutation) %>% dim() %>% .[[1]]
  
  if(nfilter==0){
    bar_input <- NULL
  }else{
    bar_input <- sigdata %>% pivot_longer(cols=-Samples) %>% #,names_transform = list(key = forcats::fct_inorder)
      filter(value>nmutation) %>%
      count(name) %>%
      mutate(value=n/dim(sigdata)[1]) %>% 
      mutate(Type="Prevalence by samples") %>% 
      select(Type,Catelogy=name,Frequency=value) %>% 
      mutate(Catelogy=factor(Catelogy,levels = c(colnames(sigdata)[-1]))) %>% 
      arrange(Catelogy) %>% 
      mutate(Catelogy=as.character(Catelogy))
  }
  
  p_barchart <- NULL
  
  if(!is.null(bar_input)) {
    p_piechart <- piechart_plot(data = pie_input,keep_legend = FALSE,colset = colset)
    p_barchart <- barchart_plot(data = bar_input,legend_name = legend_name,keep_legend = FALSE,colset = colset)
    uvalues <- sort(unique(bar_input$Catelogy))
    if(is.null(rel_widths)){rel_widths <- c(3,length(uvalues))}
    pall <- plot_grid(p_piechart+theme(plot.margin = margin(r = -2)),p_barchart,align = 'h',nrow = 1,rel_widths = rel_widths)
  }else{
    uvalues <- 3
    p_piechart <- piechart_plot(data = pie_input,keep_legend = TRUE)
    if(is.null(rel_widths)){rel_widths <- c(3,length(uvalues))}
    pall <- plot_grid(p_piechart+theme(plot.margin = margin(r = -2)),p_barchart,align = 'h',nrow = 1,rel_widths = c(3,length(uvalues)))
    
  }
  
  if(!is.null(output_plot)){
    uvalues <- sort(unique(bar_input$Catelogy))
    xleng <- 6+length(uvalues)*0.5
    xleng <- if_else(xleng>15,15,xleng)
    yleng <- 5
    if(is.null(plot_width)){ plot_width <-  xleng}
    if(is.null(plot_height)){ plot_height <-  yleng}
    ggsave(filename = output_plot,plot = pall,width = plot_width,height = plot_height)
  }
  
  pie_input <- pie_input %>% mutate(Percentage=percent(value,accuracy = 0.1)) %>% select(Type,Signature=name,Value=value,Percentage)
  if(!is.null(bar_input)){
    bar_input <- bar_input %>% select(Type,Signature=Catelogy,Value=Frequency)%>% mutate(Percentage=percent(Value,accuracy = 0.1))
    freq_data <- bind_rows(pie_input,bar_input)
  }else{
    freq_data <- pie_input
  }
  
  list_return <- list(freq_data=freq_data,pie_plot=p_piechart,bar_plot=p_barchart)
  
  return(list_return)
  
}



iupac <- function(base){
  base2 <- NA_character_
  if(base == "A"){base2 = "A"}
  if(base == "T"){base2 = "T"}
  if(base == "C"){base2 = "C"}
  if(base == "G"){base2 = "G"}
  if(base == "R"){base2 = c("A","G")}
  if(base == "Y"){base2 = c("C", "T")}
  if(base == "S"){base2 = c("G","C")}
  if(base == "W"){base2 = c("A","T")}
  if(base == "K"){base2 = c("G","T")}
  if(base == "M"){base2 = c("A","C")}
  if(base == "B"){base2 = c("C","G","T")}
  if(base == "D"){base2 = c("A","G","T")}
  if(base == "H"){base2 = c("A","C","T")}
  if(base == "V"){base2 = c("A","C","G")}
  if(base == "N"){base2 = c("A","T","C","G")}
  return(base2)
}

iupac_list <- c('A','T','C','G','R','Y','S','W','K','M','B','D','H','V','N')

context_plot <- function(data,pattern,data_return = FALSE,output_plot = NULL,plot_width=14, plot_height=9){
  # type <- 'C>T'
  # subtype1 <- 'N'
  # subtype2 <- 'G'
  #pattern <- paste0(subtype1,str_sub(type,1,1),subtype2,">",subtype1,str_sub(type,3,3),subtype2)
  type <- paste0(str_sub(pattern,2,2),str_sub(pattern,4,4),str_sub(pattern,6,6))
  subtype1 <- str_sub(pattern,1,1)
  subtype2 <- str_sub(pattern,3,3)
  pattern1 <- paste0("\n",subtype1,str_sub(type,1,1),subtype2,">",subtype1,str_sub(type,3,3),subtype2," context")
  pattern2 <- paste0(type," other contexts\n")
  
  tmpdata0 <- suppressMessages(data %>% group_by(Study,Sample) %>% summarise(Total=sum(Mutations,na.rm = TRUE)) %>% ungroup())
  tmpdata1 <- suppressMessages(data %>% filter(Type == type) %>% group_by(Study,Sample,Type) %>% summarise(N0=sum(Mutations,na.rm = TRUE)) %>% ungroup())
  tmpdata2 <- suppressMessages(data %>% filter(Type == type,SubType1 %in% iupac(subtype1),SubType2 %in% iupac(subtype2)) %>% group_by(Study,Sample,Type) %>% summarise(N1=sum(Mutations,na.rm = TRUE)) %>% ungroup())
  
  data <- suppressMessages(left_join(tmpdata0,tmpdata1) %>% 
                             left_join(tmpdata2) %>% 
                             mutate(N2=N0-N1) %>% 
                             mutate(N1=N1/Total,N2=N2/Total,Type=pattern) %>% 
                             rename(Pattern=Type) %>% 
                             filter(Total>200))
  
  if(data_return){
    return(data)
  }
  p <- data %>% 
    ggplot(aes(N1,N2,size=(Total),fill=Study))+
    geom_point(pch=21,stroke=0.25,col="black")+
    labs(x=pattern1,y=pattern2,size="Number of mutations")+
    #guides(fill=FALSE)+
    guides(fill = guide_legend(override.aes = list(size=4)))+
    xlim(c(0,1))+
    ylim(c(0,1))+
    scale_fill_manual(values = as.character(SBScolor))+
    theme_ipsum_rc(axis_title_size = 14,axis_title_just = "m",axis = TRUE, grid = "XxYy")+
    theme(axis.line.x = element_line(colour = 'black',size = 0),axis.line.y = element_line(colour = 'black',size = 0),legend.position = 'right')+
    scale_size_continuous(trans = "log10",labels = comma_format())
  
  if(is.null(output_plot)){
    return(p)
  }else{
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
  }
}


content_extraction <- function(data){
  content_data_all <- tibble(Study=character(),Sample=character(),Total=numeric(),Pattern=character(),N0=numeric(),N1=numeric(),N2=numeric())
  for(i in c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")){
    for(j1 in iupac_list){
      for(j2 in iupac_list){
        if(j1 == j2 & j1== "N"){ next }
        pattern <- paste0(j1,str_sub(i,1,1),j2,">",j1,str_sub(i,3,3),j2)
        tmp <- context_plot(data = data,pattern = pattern,data_return = TRUE)
        content_data_all <- bind_rows(content_data_all,tmp)
      }
    }
  }
  return(content_data_all)
}



signature_association <- function(data,signature_name_input1="Siganture 1",signature_name_input2="Siganture 2", signature_both=FALSE,output_plot = NULL,plot_width=10, plot_height=10){
  # cancer_type_only
  # if(!is.null(cancer_type_input)){
  #   data <- data %>% filter(Cancer_Type==cancer_type_input)
  # }
  
  if(signature_both){
    data <- data %>% filter(Exposure1>0,Exposure2>0)
  }
  
  ntotal <- data %>% filter(!is.na(Exposure1),!is.na(Exposure2)) %>% dim() %>% .[[1]]
  
  if(ntotal==0){
    errinfo <- "ERROR: No sample assigned to both signatures in selected study"
    return(errinfo)
  }else {
    
    # data %>% 
    #   ggplot(aes(log10(Exposure1+1),log10(Exposure2+1)))+
    #   geom_point()+
    #   geom_smooth(method = "lm",se = TRUE)+
    #   labs(x=paste0('Nubmer of mutations in ',signature_name_input1, ' (log10)'),y=paste0('Nubmer of mutations in ',signature_name_input2,' (log10)'))+
    #   theme_ipsum_rc(axis_title_size = 12,axis_title_just = 'm',axis_col = 'black',ticks = T)
    #   
    
    ## generated another Rplot blank file???
    p <- ggstatsplot::ggscatterstats(
      data=data %>% mutate(Exposure1=log10(Exposure1+1),Exposure2=log10(Exposure2+1)),
      x=Exposure1,
      y=Exposure2,
      xlab=paste0('Number of mutations in ',signature_name_input1,' (log10)'),
      ylab=paste0('Number of mutations in ',signature_name_input2,' (log10)'),
      marginal.type = "density",
      messages=FALSE,
    )
    
    if(is.null(output_plot)){
      return(p)
    }else{
      ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
    }
    
  }
  
  
}


# Genome2Size  ------------------------------------------------------------
genome2size <- function(genome){
  genomesize <- case_when(
    genome == "hg18" ~ 3080436051/10^6, 
    genome == "GRCh38" ~ 3217346917/10^6, 
    genome == "GRCh37" ~ 3101976562/10^6,
    genome == "hg19" ~ 3101976562/10^6,
    genome == "mmc10" ~ 2725537669/10^6,
    genome == "mmc9" ~ 2654911517/10^6,
    TRUE ~ NA_real_
  )
  return(genomesize)
}




# Rainfall_Plot -----------------------------------------------------------

# mutSNP ----------------------------------------------------------------
mutSNP.input <- function(mut.data, chr = "chr", pos = "pos", ref = "ref", 
                         alt = "alt", build = NULL, k = 10) {
  if (exists("mut.data", mode = "list")) {
    mut.data <- mut.data
  } else {
    if (file.exists(mut.data)) {
      mut.data <- utils::read.table(mut.data, sep = "\t", header = TRUE, 
                                    as.is = FALSE, check.names = FALSE)
    } else {
      stop("mut.data is neither a file nor a loaded data frame")
    }
  }
  mut.data <- mut.data[, c(chr, pos, ref, alt)]
  genome.opts = c("hg19", "hg18", "hg38")
  if (!build %in% genome.opts || is.null(build)) {
    stop("Available reference builds: hg18, hg19, hg38")
  }
  if (build == "hg19") {
    chr.lens = c(249250621, 243199373, 198022430, 191154276, 180915260, 
                 171115067, 159138663, 146364022, 141213431, 135534747, 
                 135006516, 133851895, 115169878, 107349540, 102531392, 
                 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 
                 51304566, 155270560, 59373566)
    bsg = BSgenome.Hsapiens.UCSC.hg19
  } else if (build == "hg18") {
    chr.lens = c(247249719, 242951149, 199501827, 191273063, 180857866, 
                 170899992, 158821424, 146274826, 140273252, 135374737, 
                 134452384, 132349534, 114142980, 106368585, 100338915, 
                 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 
                 49691432, 154913754, 57772954)
    bsg = BSgenome.Hsapiens.UCSC.hg18
  } else if (build == "hg38") {
    chr.lens = c(248956422, 242193529, 198295559, 190214555, 181538259, 
                 170805979, 159345973, 145138636, 138394717, 133797422, 
                 135086622, 133275309, 114364328, 107043718, 101991189, 
                 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 
                 50818468, 156040895, 57227415)
    bsg = BSgenome.Hsapiens.UCSC.hg38
  } else {
    stop("Available reference builds: hg18, hg19, hg38")
  }
  mut.data$build = build
  if (!all(mut.data$ref %in% DNA_BASES & mut.data$alt %in% DNA_BASES)) {
    stop("Only SNV substitutions are currently supported.")
  }
  ref_base = DNAStringSet(mut.data$ref)
  alt_base = DNAStringSet(mut.data$alt)
  conv.start = mut.data$pos - k
  conv.end = mut.data$pos + k
  context = getSeq(bsg, mut.data$chr, start = conv.start, end = conv.end)
  if (TRUE) {
    idx = mut.data$ref %in% c("A", "G")
    context[idx] = reverseComplement(context[idx])
    ref_base[idx] = reverseComplement(ref_base[idx])
    alt_base[idx] = reverseComplement(alt_base[idx])
  }
  mut.data$alteration = paste(ref_base, alt_base, sep = ">")
  mut.data$context = context
  # Replace chr x and y with numeric value (23 and 24) for better
  # ordering
  seq = gsub(pattern = "chr", replacement = "", x = mut.data$chr, 
             fixed = TRUE)
  seq = gsub(pattern = "X", replacement = "23", x = seq, fixed = TRUE)
  seq = gsub(pattern = "Y", replacement = "24", x = seq, fixed = TRUE)
  mut.data$seq = as.numeric(seq)
  mut.data = mut.data[order(mut.data$seq, mut.data$pos), ]
  chr.lens.sum = cumsum(chr.lens)
  chr.lens.sum = c(0, chr.lens.sum)
  mut.data$dis = c(mut.data$pos[1], diff(mut.data$pos + chr.lens.sum[mut.data$seq]))
  return(mut.data)
}




# Kataegis point ----------------------------------------------------------
katPoint <- function(data, sample = "sample", min.mut = 5, max.dis = 1000, 
                     txdb = NULL) {
  build = data$build[1]
  genome.opts = c("hg19", "hg18", "hg38")
  if (!build %in% genome.opts) {
    stop("Available reference builds: hg18, hg19, hg38")
  }
  if (build == "hg19") {
    chr.arm = c(1.25e+08, 93300000, 9.1e+07, 50400000, 48400000, 
                6.1e+07, 59900000, 45600000, 4.9e+07, 40200000, 53700000, 
                35800000, 17900000, 17600000, 1.9e+07, 36600000, 2.4e+07, 
                17200000, 26500000, 27500000, 13200000, 14700000, 60600000, 
                12500000)
  } else if (build == "hg18") {
    chr.arm = c(124300000, 93300000, 91700000, 50700000, 47700000, 
                60500000, 59100000, 45200000, 51800000, 40300000, 52900000, 
                35400000, 1.6e+07, 15600000, 1.7e+07, 38200000, 22200000, 
                16100000, 28500000, 27100000, 12300000, 11800000, 59500000, 
                11300000)
  } else if (build == "hg38") {
    chr.arm = c(123400000, 93900000, 90900000, 5e+07, 48800000, 
                59800000, 60100000, 45200000, 4.3e+07, 39800000, 53400000, 
                35500000, 17700000, 17200000, 1.9e+07, 36800000, 25100000, 
                18500000, 26200000, 28100000, 1.2e+07, 1.5e+07, 6.1e+07, 
                10400000)
  } else {
    stop("Available reference builds: hg18, hg19, hg38")
  }
  num = dim(data)[1] - 5
  katPoint <- matrix(nrow = num, ncol = 8)
  i = 1
  mutnum = 1
  Cmutnum = 0
  for (i in 1:num) {
    if (data$ref[i] %in% c("C", "G")){
      Cmutnum = Cmutnum + 1
    }
    if (data$dis[i + 1] <= max.dis) {
      mutnum = mutnum + 1
    } else {
      if (mutnum >= min.mut) {
        len = data$pos[i] - data$pos[i - mutnum + 1] + 1
        chr.n = gsub(pattern = "chr", replacement = "", x = data$chr[i], 
                     fixed = TRUE)
        chr.n = gsub(pattern = "X", replacement = "23", x = chr.n, 
                     fixed = TRUE)
        chr.n = gsub(pattern = "Y", replacement = "24", x = chr.n, 
                     fixed = TRUE)
        chr.n = as.numeric(chr.n)
        if (data$pos[i] <= chr.arm[chr.n]) {
          arm = paste(chr.n, "p", sep = "")
        } else if (data$pos[i - mutnum + 1] >= chr.arm[chr.n]) {
          arm = paste(chr.n, "q", sep = "")
        } else {
          arm = paste(chr.n, "p, ", chr.n, "q", sep = "")
        }
        katPoint[i, 1:8] = c(sample, data$chr[i], data$pos[i - mutnum + 1], data$pos[i], arm, len, mutnum, round(Cmutnum/mutnum,3))
      }
      mutnum = 1
      Cmutnum = 0
    }
  }
  katPoint.out = data.frame(na.omit(katPoint))
  names(katPoint.out) = c("sample", "chrom", "start", "end", "chrom.arm", "length", "number.mut", 
                          "weight.C>X")
  for (i in 1:dim(katPoint.out)[1]) {
    if (as.numeric(as.character(katPoint.out$"weight.C>X"[i])) < 0.8) {
      katPoint.out$confidence[i] = 0
    } else {
      katPoint.out$confidence[i] <- length(which(subset(katPoint.out,as.numeric(as.character(katPoint.out$"weight.C>X")) >= 0.8)$chrom == katPoint.out$chrom[i]))
      if (katPoint.out$confidence[i] > 3) {
        katPoint.out$confidence[i] = 3
      }
    }
  }
  if (!is.null(txdb)) {
    gr <- GRanges(seqnames = Rle(katPoint.out$chrom), ranges=IRanges(start = as.numeric(as.character(katPoint.out$start)), end =as.numeric(as.character(katPoint.out$end))))
    peakAnno <- ChIPseeker::annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
    katPoint.out$annotation <- peakAnno@anno$annotation
    katPoint.out$distanceToTSS <- peakAnno@anno$distanceToTSS
    katPoint.out$geneName <- peakAnno@anno$SYMBOL
    katPoint.out$geneID <- peakAnno@anno$geneId
  } 
  message(paste(dim(katPoint.out)[1], "potential kataegis events identified", 
                sep = " "))
  return(katPoint.out)
}


# Kataegis plot -----------------------------------------------------------

kataegis_rainfall_plot <- function(mutdata,sample_name="sample",genome_build = "hg19",reference_data_folder=NULL, chromsome=NULL,kataegis_highligh=FALSE,min.mut = 5,max.dis = 1000,min.dis = 10, filename=NULL){
  require(tidyverse)
  require(ggsci)
  require(hrbrthemes)
  require(cowplot)
  #require(ChIPseeker)
  
  genome_build <- str_to_lower(genome_build)
  
  if(genome_build %in% c('hg19','grch37')){
    require(BSgenome.Hsapiens.UCSC.hg19)
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    genome_build <- 'hg19'
    if(!is.null(reference_data_folder)){
      ref_file <- paste0(reference_data_folder,"/hg19_ref.RData")
      load(ref_file)
    }
    hgref <- hg19
    TxDb.Hsapiens <- TxDb.Hsapiens.UCSC.hg19.knownGene
  }
  
  if(genome_build %in% c('hg38','grch38')){
    require(BSgenome.Hsapiens.UCSC.hg38)
    require(TxDb.Hsapiens.UCSC.hg38.knownGene)
    genome_build <- 'hg38'
    if(!is.null(reference_data_folder)){
      ref_file <- paste0(reference_data_folder,"/hg38_ref.RData")
      load(ref_file)
    }
    hgref <- hg38
    TxDb.Hsapiens <- TxDb.Hsapiens.UCSC.hg38.knownGene
  }
  
  
  #mutdata format:chr,pos,ref,alt
  mutdata <- mutdata %>% 
    mutate(chr=str_remove(chr,"^chr")) %>% 
    dplyr::filter(str_length(ref)==1,str_length(alt)==1,!str_detect(alt,"-"),!str_detect(ref,"-"),chr %in% c(1:22,"X","Y"))%>% 
    mutate(chr=paste0("chr",chr)) 
  
  if(!is.null(chromsome)){
    hgref <- hgref %>% dplyr::filter(chr==chromsome) %>% mutate(start2=1,end2=len)
    mutdata <- mutdata %>% dplyr::filter(chr==chromsome)
  }
  
  if(dim(mutdata)[1] < min.mut){
    stop("ERROR: Not enough mutations for kataegis identification!")
  }
  
  mutSNP = mutSNP.input(mut.data = as.data.frame(mutdata), chr = "chr", pos = "pos", ref = "ref", alt = "alt", build = genome_build)
  mutSNP$context <- ""
  katdata <- NA
  try({
    katdata <- katPoint(mutSNP %>% filter(dis > min.dis),txdb = TxDb.Hsapiens,sample = sample_name,min.mut = min.mut,max.dis = max.dis)
  },silent = FALSE)
  #print(katdata)
  
  # extend the min.dist for capture the DBS or MBS; number.mut will still be the mutations > min.dist
  if(is.data.frame(katdata)){
    if(dim(katdata)[1]>0){
      katdata <- katdata %>% mutate(start=as.integer(start),end=as.integer(end)) %>% mutate(start=start-min.dis, end=end+min.dis) %>% left_join(hgref,by=c('chrom'='chr'))
      for(i in 1:dim(katdata)[1]){
        kchr <- katdata$chrom[i]
        kstart <- katdata$start[i]
        kend <- katdata$end[i]
        mutSNP$context[mutSNP$chr==kchr & mutSNP$pos>=kstart & mutSNP$pos<=kend] <- 'Kataegsis'
      }
    }
  }
  
  #color themes for the mutation subtype
  SNVcolor <- pal_d3()(6)
  names(SNVcolor) = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  addcolor <- "gray80"
  names(addcolor) <- "Non-kataegis mutations"
  
  mutSNP <- mutSNP %>% 
    left_join(hgref %>% dplyr::select(chr,start2),by=c('chr'='chr')) %>% 
    mutate(pos=pos+start2-1) %>% 
    dplyr::rename(dist=dis) 
  
  if(kataegis_highligh){
    mutSNP <- mutSNP %>% mutate(SNV=if_else(context=="","Non-kataegis mutations",alteration))
    SNVcolor <- c(SNVcolor,addcolor)
  }else{
    mutSNP <- mutSNP %>% mutate(SNV=alteration) 
  }
  
  mutSNP <- mutSNP %>% mutate(SNV=factor(SNV,level=names(SNVcolor))) 
  
  p <- mutSNP %>% 
    ggplot(aes(pos,log10(dist),fill=SNV))+
    geom_vline(xintercept = c(hgref$start2,hgref$end2),col="#cccccc",lwd=0.4)+
    geom_hline(yintercept = c(2,3),col=c("red","blue"),lwd=0.4,linetype=2)+
    geom_point(shape=21,stroke=0.1,size=1.5)+
    #labs(y="Intermutation distance (bp,log10)\n",x = "Chromosome")+
    #scale_x_continuous(breaks = (hgref$start2+hgref$end2)/2,labels = str_remove(hgref$chr,'chr'),expand = c(0.005,0.005),limits = c(0,tail(hgref$end2,1)),guide = guide_axis(n.dodge=1))+
    scale_y_continuous(breaks = pretty_breaks())+
    scale_fill_manual("Subtype",values = SNVcolor,drop=FALSE,na.value="gray90")+
    theme_ipsum_rc(axis_title_just = "m",axis_title_size = 14,grid="Y",ticks = T,axis = FALSE)+
    theme(plot.title = element_text(hjust = 0.5,face = "plain"),legend.text = element_text(size = 12),legend.position = "bottom")+
    guides(fill = guide_legend(override.aes = list(size = 3),nrow = 1))+
    panel_border(color = 'black',size = 0.5)+
    coord_cartesian(clip = 'off')
  
  if(!is.null(chromsome)){
    p <- p+scale_x_continuous(breaks = pretty_breaks(),labels = label_comma(),expand = c(0.005,0.005))+labs(y="Intermutation distance (bp,log10)\n",x = paste0("Chromosome ",extract_numeric(chromsome)))
  }else{
    p <- p+scale_x_continuous(breaks = (hgref$start2+hgref$end2)/2,labels = str_remove(hgref$chr,'chr'),expand = c(0.005,0.005),limits = c(0,tail(hgref$end2,1)),guide = guide_axis(n.dodge=1))+labs(y="Intermutation distance (bp,log10)\n",x = "Chromosome")
  }
  
  
  if(is.data.frame(katdata)){
    if(dim(katdata)[1]>0){
      p <- p+annotate("segment", x = katdata$start2+(katdata$start+katdata$end)/2, xend = katdata$start2+(katdata$start+katdata$end)/2, y = 0.5, yend = 1,size=0.3, colour = "black",arrow=arrow(length = unit(0.25,"cm")))
      katdata <- katdata %>% dplyr::select(sample:geneID)
    }}
  ggsave(filename = filename,plot = p,width = 15,height = 7)
  
  return(katdata)
}





# seqmatrix_public_download  ---------------------------------------------------------------

seqmatrix_public_download <- function(seqmatrix_refdata_public,tmpfolder= "./seqmatrix_public_testdownload"){
  
  #tmpfolder <- '~/Downloads/vistmp'
  tmpfile <- paste0(tmpfolder,'/',study,"_",cancer_type)
  dir.create(path = tmpfolder,recursive = TRUE)
  profile_list <- seqmatrix_refdata_public %>% count(Profile) %>% pull(Profile)
  
  for(profiletmp in profile_list){
    dupsample <- seqmatrix_refdata_public %>% select(Study,Cancer_Type,Sample) %>% unique() %>% count(Sample) %>% filter(n>1) %>% dim() %>% .[[1]]
    if(dupsample == 0){
      seqmatrix_refdata_public %>%
        filter(Profile==profiletmp) %>% 
        select(Sample,MutationType,Mutations) %>% 
        pivot_wider(id_cols = MutationType,names_from=Sample,values_from=Mutations) %>% 
        write_delim(paste0(tmpfile,'.',profiletmp,'.txt'),delim = '\t',col_names = T)
    }else{
      seqmatrix_refdata_public %>%
        filter(Profile==profiletmp) %>% 
        mutate(Sample=paste0(Study,"_",Cancer_Type,"_",Sample)) %>% 
        select(Sample,MutationType,Mutations) %>% 
        pivot_wider(id_cols = MutationType,names_from=Sample,values_from=Mutations) %>% 
        write_delim(paste0(tmpfile,'.',profiletmp,'.txt'),delim = '\t',col_names = T)
    }
  }
  
  
  cmdfile1 <- paste0('tar -C ',str_remove(tmpfolder,"/[^/]*$"),' -zcvf ', tmpfolder,".tar.gz ",str_remove(tmpfolder,".*/"))
  system(cmdfile1)
  cmdfile2 <- paste0('rm -rf ',tmpfolder)
  system(cmdfile2)
}



### Association function ### 

#require(hablar)
change_data_type <- function(data){
  if(is.numeric(data)){
    data <- as.character(data)
  }else if(is.character(data)){
    data <- as.numeric(data)
  }
  return(data)
}

validate_vardf <- function(data, forces=NULL, nachars=c("","na","Na","NA","nan","NAN","Nan"), nacode='NA',Nmin=5, Nmin_drop=FALSE,excludes=NULL, lump=FALSE){
  #names_oringal <- colnames(data)
  # force data type
  
  # for eqaul number of unique value
  equal_chars <- NULL
  chrlist <- colnames(data %>% select(where(~is.character(.) & (n_distinct(.)>1))))
  
  if(!is.null(excludes)){ chrlist <- chrlist[!chrlist %in% excludes]}
  
  if(length(chrlist)>0){
    
    for(chr in chrlist){
      tmpn <- data %>% select(varname=one_of(chr)) %>% count(varname) %>% count(n) %>% dim() %>% .[[1]]
      if(tmpn==1){
        equal_chars <- c(equal_chars,chr)
      }
    }
    
    data <- data %>% mutate(across(one_of(equal_chars),as.factor))
    excludes <- c(excludes,equal_chars)
  }
  
  
  if(!is.null(excludes)){
    names_all <- colnames(data)
    data_excludes <- data  %>% dplyr::select(one_of(excludes))
    data <- data %>% dplyr::select(-one_of(excludes))
  }
  
  if(!is.null(forces)){ 
    data <- data %>% mutate(across(one_of(forces),change_data_type))  
  }
  
  # process na
  data <- data %>% mutate(across(where(is.character), ~  if_else(.x %in% nachars, NA_character_, .x)))
  data <- data %>% mutate(across(where(is.character), ~  if_else(is.na(.x), nacode, .x)))
  
  # process the integer
  data <- data %>% mutate(across(where(is.integer), as.numeric))
  
  # convert numbers to characters if less than xx unique value
  #data %>% mutate(across(where(is.numeric), ~ if_else(n_distinct(.) < Nmin, .x, .x)))
  Nmin_names <- data %>% summarise(across(where(is.numeric),n_distinct)) %>% dplyr::select_if(function(x) x<Nmin) %>% colnames()
  if(length(Nmin_names)>0){
    if(Nmin_drop){
      data <- data %>% dplyr::select(-one_of(Nmin_names))
      #names_oringal <- names_oringal[!(names_oringal %in% Nmin_names)]
    }
    else{
      data <- data %>% mutate(across(one_of(Nmin_names),change_data_type)) 
    }
  }
  
  ## make the factors and reorder the levels
  ## exclude for the all unique value columns
  #data <- data %>% mutate(across(where(is.character),~ fct_lump(fct_infreq(as.factor(.x)),prop = 0.2)))
  if(lump){
    nrows <- dim(data)[1]
    data <- data %>% mutate(across(where(~ is.character(.) & (n_distinct(.)>1) & (n_distinct(.) < 0.8*nrows)),~ fct_lump(fct_infreq(as.factor(.x)),prop = 0.2)))
  }
  # if change the order
  #data <- data %>% select(names_oringal)
  
  if(!is.null(excludes)){
    names_keep <- c(colnames(data),colnames(data_excludes))
    names_all <- names_all[names_all %in% names_keep]
    data <- bind_cols(data_excludes,data) %>% select(one_of(names_all))
  }
  
  return(as_tibble(data))
}


### other funcitons ###
rowAny <- function(x) rowSums(x) > 0


# mSigPortal_associaiton --------------------------------------------------
mSigPortal_associaiton <- function(data, Var1, Var2, regression=FALSE, formula=NULL, xlab="Variable1", ylab="Variable2",filter1=NULL, filter2=NULL,log1=FALSE,log2=FALSE, type="parametric", collapse_var1=NULL, collapse_var2=NULL, output_plot=NULL,plot_width=12,plot_height=8) {
  
  data <- validate_vardf(data,lump = T)
  
  if(regression){
    ## for regression module
    supported_types <- c("lm", "glm")
    
    if(!str_detect(formula,"~")){
      stop("Please check your formula for regression, for example, lm( mpg ~ vs + gear")
    }
    
    input_formula <- paste0("mod <- data %>% ",type, "(", formula,", data=.)")
    eval(parse(text=input_formula))
    
    p <- ggstatsplot::ggcoefstats(
      x= mod,
      point.args = list(color = "red", size = 3, shape = 15),
      exclude.intercept = TRUE,
      #title = formula,
      ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14)
    ) + # note the order in which the labels are entered
      ggplot2::labs(x = "Regression Coefficient", y = NULL)
  }else{
    
    ## subset data
    data <- data %>% select(one_of(c(Var1,Var2)))
    colnames(data) <- c("Var1","Var2")
    var1_type <- if_else(is.factor(data[[1]]),"categorical", if_else(is.numeric(data[[1]]),"continuous",NA_character_))
    var2_type <- if_else(is.factor(data[[2]]),"categorical", if_else(is.numeric(data[[2]]),"continuous",NA_character_))
    
    if(is.na(var1_type)|is.na(var2_type)){
      stop("Please check your data type of these two selected variables")
    }
    
    # process data or filtering data
    if(!is.null(filter1) & var1_type == 'continuous') {
      filter1 <-  as.numeric(filter1)
      if(!is.na(filter1)){
        data <- data %>% filter(Var1 > filter1)
      }
      
    }
    
    if(!is.null(filter2) & var2_type == 'continuous') {
      filter2 <-  as.numeric(filter2)
      if(!is.na(filter2)){
        data <- data %>% filter(Var2 > filter2)
      }
    }
    
    
    if(log1 & var1_type == 'continuous') {
      data <- data %>% filter(Var1>0) %>% mutate(Var1 = log2(Var1))
    }
    
    if(log2 & var2_type == 'continuous') {
      data <- data %>% filter(Var2>0) %>% mutate(Var2 = log2(Var2))
    }
    
    if(var1_type =="categorical" && !is.null(collapse_var1)){
      if(! (collapse_var1 %in% data$Var1)){ 
        print("Warning: categorical value does not exist in data, please input the correct level of the categorical variables for variable1.")
      }else{
        data$Var1 <- fct_other(data$Var1,keep = collapse_var1)
      }
    }
    
    if(var2_type =="categorical" && !is.null(collapse_var2)){
      if(! (collapse_var2 %in% data$Var2)){ 
        print("Warning: categorical value does not exist in data, please input the correct level of the categorical variables for variable1.")
      }else{
        data$Var2 <- fct_other(data$Var2,keep = collapse_var2)
      }
    }
    
    ## association test based on the types
    
    # for continues vs continues
    if(var1_type == 'continuous' & var2_type == 'continuous'){
      # supported types: "parametric", "nonparametric", "robust", "bayes", "skit
      supported_types <- c("parametric", "nonparametric", "robust", "bayes", "skit")
      if(!(type %in% supported_types)){
        print("Warning: selected statistic parameter type does not supported for selected data types, use the default parameteric type")
        type = "parametric"
      }
      
      if(type == "skit"){
        skit_res <- SKIT::skit(data$Var1,data$Var2,nboot = 1000)
        skit_res <- skit_res$pvalues[1]
        skit_lab <- paste0("P-value by SKIT test: ",if_else(skit_res <= 1/1000, "<1e-03",as.character(scientific(skit_res,digits = 3))))
        xlab=paste0(xlab,"\n",skit_lab)
        type = "parametric"
      }
      
      p <-  ggstatsplot::ggscatterstats(
        data = data,
        x = Var1,
        y = Var2, 
        xlab= xlab,
        ylab = ylab,
        marginal.type = "boxplot",
        xfill = "#009E73",
        yfill = "#D55E00",
        ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14),
        type=type
      )
      
    }
    
    # for categorical vs categorical
    if(var1_type == 'categorical' & var2_type == 'categorical'){
      # supported types: "parametric", "nonparametric", "robust", "bayes"
      supported_types <- c("parametric", "nonparametric", "robust", "bayes","fisher")
      if(!(type %in% supported_types)){
        print("Warning: selected statistic parameter type does not supported for selected data types, use the default parameteric type")
        type = "parametric"
      }
      
      if(type == "fisher"){
        fisher_res <- tidy(fisher.test(data$Var1,data$Var2))
        if("estimate" %in% colnames(fisher_res)){
          fisher_lab <- paste0("Fisher Exact Test: P-value =", scientific(fisher_res$p.value,digits = 3),", OR = ",number_format(accuracy = 0.01)(fisher_res$estimate))
        }else{
          fisher_lab <- paste0("Fisher Exact Test: P-value =", scientific(fisher_res$p.value,digits = 3))
        }
        xlab=paste0(xlab,"\n",fisher_lab)
        type = "parametric"
      }
      
      p <-  ggstatsplot::ggbarstats(
        data = data,
        x = Var1,
        y = Var2, 
        xlab= xlab,
        legend.title = ylab,
        ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14),
        type=type
      )
    }
    
    # for categorical vs continues 
    if(var1_type != var2_type ){
      # supported types: "parametric", "nonparametric", "robust", "bayes"
      # switch the name if Var2 is categorical
      if(var2_type == 'categorical'){
        p <-  ggstatsplot::ggbetweenstats(
          data = data,
          x = Var2,
          y = Var1, 
          xlab= ylab,
          ylab = xlab,
          ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14),
          type=type
        )
        
      }else {
        p <-  ggstatsplot::ggbetweenstats(
          data = data,
          x = Var1,
          y = Var2, 
          xlab= xlab,
          ylab = ylab,
          ggtheme = hrbrthemes::theme_ipsum_rc(axis_title_just = 'm',axis_title_size = 14),
          type=type
        )
        
      }
    }
    
  } 
  
  if(is.null(output_plot)){
    return(p)
  }else{
    ggsave(filename = output_plot,plot = p,width = plot_width,height = plot_height)
    #return(p)
  }
}



mSigPortal_associaiton_group <- function(data, Var1, Var2, Group_Var, regression=FALSE, formula=NULL, filter1=NULL, filter2=NULL,log1=FALSE,log2=FALSE,type="parametric", collapse_var1=NULL, collapse_var2=NULL) {
  
  data <- validate_vardf(data,excludes = Group_Var)
  
  if(regression){
    ## for regression module
    supported_types <- c("lm", "glm")
    
    if(is.null(formula)|!str_detect(formula,"~")){
      stop("Please check your formula for regression, for example, lm( mpg ~ vs + gear")
    }
    
    colnames(data)[colnames(data) == Group_Var] <- 'Group'
    
    if(type == "lm"){
      result <- data %>% group_by(Group) %>% do(tidy(lm(formula,data=.))) %>% ungroup() %>% filter(term!="(Intercept)") %>% filter(!is.na(p.value)) %>% ungroup()
    }
    
    if(type == "glm"){
      result <- data %>% group_by(Group) %>% do(tidy(glm(formula,data=.))) %>% ungroup() %>% filter(term!="(Intercept)") %>% filter(!is.na(p.value)) %>% ungroup()
    } 
    
    colnames(result)[1] <- tolower(Group_Var)
    result <- result %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) #%>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH") #%>% mutate(formula=formula)
    
  }else{
    
    ## subset data
    data <- data %>% select(one_of(c(Group_Var,Var1,Var2)))
    colnames(data) <- c("Group","Var1","Var2")
    var1_type <- if_else(is.factor(data[["Var1"]]),"categorical", if_else(is.numeric(data[["Var1"]]),"continuous",NA_character_))
    var2_type <- if_else(is.factor(data[["Var2"]]),"categorical", if_else(is.numeric(data[["Var2"]]),"continuous",NA_character_))
    
    if(is.na(var1_type)|is.na(var2_type)){
      stop("Please check your data type of these two selected variables")
    }
    
    # process data or filtering data
    if(!is.null(filter1) & var1_type == 'continuous') {
      filter1 <-  as.numeric(filter1)
      if(!is.na(filter1)){
        data <- data %>% filter(Var1 > filter1)
      }
      
    }
    
    if(!is.null(filter2) & var2_type == 'continuous') {
      filter2 <-  as.numeric(filter2)
      if(!is.na(filter2)){
        data <- data %>% filter(Var2 > filter2)
      }
    }
    
    if(log1 & var1_type == 'continuous') {
      data <- data %>% filter(Var1>0) %>% mutate(Var1 = log2(Var1))
    }
    
    if(log2 & var2_type == 'continuous') {
      data <- data %>% filter(Var2>0) %>% mutate(Var2 = log2(Var2))
    }
    
    if(var1_type =="categorical" && !is.null(collapse_var1)){
      if(! (collapse_var1 %in% data$Var1)){ 
        print("Warning: categorical value does not exist in data, please input the correct level of the categorical variables for variable1.")
      }else{
        data$Var1 <- fct_other(data$Var1,keep = collapse_var1)
      }
    }
    
    if(var2_type =="categorical" && !is.null(collapse_var2)){
      if(! (collapse_var2 %in% data$Var2)){ 
        print("Warning: categorical value does not exist in data, please input the correct level of the categorical variables for variable1.")
      }else{
        data$Var2 <- fct_other(data$Var2,keep = collapse_var2)
      }
    }
    
    ## association test based on the types
    
    # for continues vs continues
    if(var1_type == 'continuous' & var2_type == 'continuous'){
      # supported types: "parametric", "nonparametric", "robust", "bayes", "skit
      supported_types <- c("parametric", "nonparametric", "robust", "bayes", "skit")
      if(!(type %in% supported_types)){
        print("Warning: selected statistic parameter type does not supported for selected data types, use the default parameteric type")
        type = "parametric"
      }
      
      if(type == "skit"){
        
        result <- tibble(Group=character(),parameter1=character(),parameter2=character(),p.value=numeric(), method=character(), n.obs=integer())
        for(sig in unique(data$Group)){
          tmp <- data %>% filter(Group==sig)
          skit_res <- SKIT::skit(tmp$Var1,tmp$Var2,nboot = 1000)
          skit_res <- as.numeric(skit_res$pvalues[1])
          nobs <- tmp %>% filter(!is.na(Var1),!is.na(Var2)) %>% dim() %>% .[[1]]
          result <- tibble(Group=sig,parameter1="Var1",parameter2="Var2",p.value=skit_res, method="SKIT test", n.obs=nobs) %>% 
            bind_rows(result)
        }
        
      }else{
        result <- data %>% group_by(Group) %>% group_modify(~statsExpressions::corr_test(data = .,x=Var1,y=Var2,type=type) %>% select(-expression)) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
      }
      
      result$parameter1 <- Var1
      result$parameter2 <- Var2
      colnames(result)[1:3] <- c(tolower(Group_Var),"variable_name1","varible_name2")
      
    }
    
    # for categorical vs categorical
    if(var1_type == 'categorical' & var2_type == 'categorical'){
      # supported types: "parametric", "nonparametric", "robust", "bayes"
      supported_types <- c("parametric", "nonparametric", "robust", "bayes","fisher")
      if(!(type %in% supported_types)){
        print("Warning: selected statistic parameter type does not supported for selected data types, use the default parameteric type")
        type = "parametric"
      }
      
      tmp <- data %>% count(Group,Var1,Var2) %>% count(Group) %>% filter(n>2) %>% pull(Group)
      
      if(type == "fisher"){
        result <-  data %>% filter(Group %in% tmp) %>%  nest_by(Group) %>% mutate(test=list(fisher.test(data$Var1,data$Var2))) %>% summarise(tidy(test)) %>% arrange(p.value) %>% ungroup() %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% ungroup() 
      }else{
        result <- data %>% filter(Group %in% tmp) %>%  group_by(Group) %>% group_modify(~tryCatch(expr = statsExpressions::contingency_table(data = .,x=Var1,y=Var2,type=type), error = function(e) NULL)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
      }
      result <- result %>% mutate(variable_name1=Var1, variable_name2 = Var2) %>% select(Group,variable_name1,variable_name2,everything())
      colnames(result)[1] <- c(tolower(Group_Var))
    }
    
    # for categorical vs continues 
    if(var1_type != var2_type ){
      # supported types: "parametric", "nonparametric", "robust", "bayes"
      # switch the name if Var2 is categorical
      ## remove unique value
      
      ## decide two sample test or oneway_annovar 
      
      if(var1_type=="categorical"){
        if(length(levels(data$Var1))==2){
          
          tmp <- data %>% group_by(Group) %>% filter(!is.na(Var1),!is.na(Var2))%>%  summarise(n1=n_distinct(Var1),n2=n_distinct(Var2)) %>% filter(n1!=2|n2==1) %>% pull(Group)
          result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::two_sample_test(data = .,x=Var1,y=Var2,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
        }
        
        if(length(levels(data$Var1))>2){
          tmp <- data %>% filter(!is.na(Var1),!is.na(Var2)) %>% group_by(Group,Var1) %>% summarise(SD=sd(Var2)) %>% filter(SD==0) %>% select(Group) %>% unique() %>% pull(Group) %>% as.character()
          result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::oneway_anova(data = .,x=Var1,y=Var2,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
        }
        
        result$parameter1 <- Var1
        result$parameter2 <- Var2
        colnames(result)[1:3] <- c(tolower(Group_Var),"variable_name1","varible_name2")
      }else{
        
        # 
        # result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::two_sample_test(data = .,x=Var2,y=Var1,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")  
        # 
        if(length(levels(data$Var2))==2){
          
          tmp <- data %>% group_by(Group) %>% filter(!is.na(Var1),!is.na(Var2))%>%  summarise(n1=n_distinct(Var2),n2=n_distinct(Var1)) %>% filter(n1!=2|n2==1) %>% pull(Group)
          result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::two_sample_test(data = .,x=Var2,y=Var1,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
        }
        
        if(length(levels(data$Var2))>2){
          tmp <- data %>% filter(!is.na(Var1),!is.na(Var2)) %>% group_by(Group,Var2) %>% summarise(SD=sd(Var1)) %>% filter(SD==0) %>% select(Group) %>% unique() %>% pull(Group) %>% as.character()
          result <- data  %>% filter(!Group %in% tmp) %>%   group_by(Group) %>% group_modify(~statsExpressions::oneway_anova(data = .,x=Var2,y=Var1,type=type)) %>% select(-expression) %>% ungroup() %>% filter(!is.na(p.value)) %>% arrange(p.value) %>% mutate(fdr=p.adjust(p.value,method = 'BH')) %>% mutate(fdr.method="BH")
        }
        
        result$parameter1 <- Var2
        result$parameter2 <- Var1
        colnames(result)[1:3] <- c(tolower(Group_Var),"variable_name1","varible_name2")
      }
    }
    
  } 
  
  return(result)
}





multivariable_inputs <- function(data,listpars) {
  
  result <- NULL
  
  for(i in 1:length(listpars)){
    Var <- listpars[[i]]
    # data format: data_source     data_type          variable_name variable_value variable_value_type
    # parameters: c('quality control', 'sequencing metrics', 'FWHM_Normal',NULL, FALSE, NULL)
    datav <- data %>%
      filter(data_source == Var$source, data_type == Var$type, variable_name == Var$name)
    
    var_type <- if_else(unique(datav$variable_value_type) == "character","categorical", if_else(unique(datav$variable_value_type) == "numeric","continuous",NA_character_))
    
    # change data type for regression
    if(var_type == 'continuous') { datav$variable_value = as.numeric(datav$variable_value)}
    if(var_type == 'categorical') { datav$variable_value = as.factor(datav$variable_value)}
    
    # process data or filtering data
    if(!is.null(Var$filter) & var_type == 'continuous') {
      Var$filter <-  as.numeric(Var$filter)
      if(!is.na(Var$filter)){
        datav <- datav %>% filter(Var1 > filter1)
      }
      
    }
    
    ## Var[5] shown as character
    
    if(Var$log2== "TRUE" & var_type == 'continuous') {
      datav <- datav %>% filter(variable_value>0) %>% mutate(variable_value = log2(variable_value))
    }
    
    if(var_type =="categorical" && !is.null(Var$collapse)){
      if(! (Var$collapse %in% datav$variable_value)){ 
        print("Warning: categorical value does not exist in data, please input the correct level of the categorical variables for variable1.")
      }else{
        datav$variable_value <- fct_other(datav$variable_value,keep = Var$collapse)
      }
    }
    
    datav <- datav %>% select(Sample,variable_value) 
    #colnames(datav)[2] <- paste0('Variable',i)
    colnames(datav)[2] <- Var$name
    
    if(is.null(result)){
      result <-  datav
    }else{
      result <- left_join(result,datav)
    }
    
  }
  
  return(result)
}


