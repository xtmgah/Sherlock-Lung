
# uniform funciton for data visualization


barplot_fisher <- function(mdata0,var1name,var2name,samplelist=NULL,filename=NULL, var1lab=NULL, var2lab=NULL,width = 5,height = 8,fisher=TRUE,textcol='white'){
  #pdfhr2()
  
  if(is.null(var1lab)){ var1lab = var1name}
  if(is.null(var2lab)){ var1lab = var2name}
  print(var1lab)
  
  mdata <- mdata0
  
  if(!is.null(samplelist)){
    mdata <- mdata0 %>% filter(Tumor_Barcode %in% samplelist)
  }
  
  mdata <- mdata %>% select(Tumor_Barcode,one_of(c(var1name,var2name))) %>% drop_na() 
  
  if(fisher){
    tresult <- mdata %>% select(-Tumor_Barcode) %>% table() %>% fisher.test() %>% tidy()
    print(tresult)
    titletext <- paste0("OR = ",round(tresult$estimate,2),"; P = ",scientific(tresult$p.value,3))
  }else{
    titletext <- ""
  }
  

  colnames(mdata)[2] <- 'Var1'
  colnames(mdata)[3] <- 'Var2'
  
  vartmp1 <- mdata %>% group_by(Var1) %>% tally() %>% dplyr::mutate(percent=percent_format(accuracy = 0.1)(n/sum(n))) %>% mutate(Var1_Lab=paste0(Var1,' (',percent,')')) %>% select(Var1,Var1_Lab)
  vartmp2 <- mdata %>% group_by(Var2) %>% tally() %>% dplyr::mutate(percent=percent_format(accuracy = 0.1)(n/sum(n))) %>% mutate(Var2_Lab=paste0(Var2,' (',percent,')')) %>% select(Var2,Var2_Lab)
  
  p <-  mdata %>% 
    group_by(Var1,Var2) %>% 
    dplyr::tally()%>%
    dplyr::mutate(percent=n/sum(n)) %>% 
    ungroup() %>% 
    arrange(Var1,Var2) %>% 
    left_join(vartmp1) %>% 
    left_join(vartmp2) %>%
    mutate(Var1_Lab = fct_inorder(Var1_Lab)) %>% 
    mutate(Var2_Lab = fct_inorder(Var2_Lab)) %>% 
    ggplot(aes(x=Var1_Lab, y=n, fill=Var2_Lab))+
    geom_bar(stat="identity", position ="fill",width = 0.8)+
    geom_text(aes(label=paste0(n,'\n',sprintf("%1.1f", percent*100),"%")), position=position_fill(vjust=0.5), colour=textcol,family='Roboto Condensed')+
    scale_y_continuous(breaks = pretty_breaks(),labels = percent_format(),expand = c(0,0))+
    scale_fill_jama()+
    theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,axis_text_size = 12,grid = 'Yy',ticks = FALSE)+
    labs(title = titletext,fill=var2lab,x=var1lab, y = 'Percentage')+
    theme(text = element_text(family = 'Roboto Condensed'),
          plot.title = element_text(size = 15,hjust = 0.5,face = 'plain',family = 'Roboto Condensed'),
          axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1))+
    coord_cartesian(clip = 'off')
  
  if(!is.null(filename)){
    ggsave(filename = filename,width = width,height = height,device = cairo_pdf )
  }else{
    return(p)
  }
  
}





