library("survminer")
library("survival")

# Function ----------------------------------------------------------------


SurvZTW <- function(suvdata_input,plot=TRUE,legend.labs=NULL,keyname="Strata",filename=NULL,width = 7,height = 5.5,keyterm="KeyY"){
  fit <- coxph(Surv(Survival_Month, Death) ~ Key+Age+Gender+Stage+Smoking+Histology, data = suvdata_input) ### overall 
  suvpvalue <- broom::tidy(fit,conf.int=T) %>% filter(term==keyterm) %>% mutate(lab=paste0("P = ",format(p.value,digits = 2,scientific = FALSE),"\nHR = ",format(exp(estimate),digits = 2,scientific = FALSE)," (95% CI: ",format(exp(conf.low),digits = 2,scientific = FALSE),"-",format(exp(conf.high),digits = 2,scientific = FALSE),")")) %>% select(term,p.value,lab)
  if(plot){
    fit2 <- coxph(Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Smoking+Histology, data = suvdata_input) ## strata curve ##
    ##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata_input)  ## all strata, wrong way
    suvfit <- survfit(fit2)
    #plot(suvfit)
    suvfit$call$formula <- as.formula("Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Smoking+Histology")
    
    if(is.null(legend.labs)) 
    {
      legend.labs <- c(paste0("N (",suvfit$n[1],")"),paste0("Y (",suvfit$n[2],")"))
    }
    keyname=keyname
    ggsurv <- ggsurvplot(
      suvfit,                     # survfit object with calculated statistics.
      data = suvdata_input,             # data used to fit survival curves.
      risk.table = FALSE,       # show risk table.
      pval = suvpvalue$lab,             # show p-value of log-rank test.
      pval.size=4,
      #pval = TRUE,
      #pval.method
      #conf.int = TRUE,         # show confidence intervals for 
      # point estimates of survival curves.
      palette = pal_jama()(2),
      #xlim = c(0,210),         # present narrower X axis, but not affect
      # survival estimates.
      xlab = "Time in months",   # customize X axis label.
      break.time.by = 50,     # break X axis in time intervals by 500.
      ggtheme = theme_ipsum_rc(axis_title_size = 14,base_size = 12,axis_title_just = "m",axis = "XY",axis_col = "black"), # customize plot and risk table with a theme.
      font.family="Roboto Condensed" ,
      risk.table.y.text.col = T,# colour risk table text annotations.
      risk.table.height = 0.25, # the height of the risk table
      risk.table.y.text = FALSE,# show bars instead of names in text annotations
      # in legend of risk table.
      censor=TRUE,
      ncensor.plot = FALSE,      # plot the number of censored subjects at time t
      ncensor.plot.height = 0.25,
      conf.int.style = "step",  # customize style of confidence intervals
      surv.median.line = "none",
      legend.labs = legend.labs
    )
    # ggsurv$plot <- ggsurv$plot + labs(
    #   title    = "Survival curve",                     
    #   subtitle = "Based on Kaplan-Meier estimates"
    # )
    ggsurv <- ggpar(
      ggsurv,
      # font.title    = c(15, "bold", "darkblue"),         
      # font.subtitle = c(14, "bold.italic", "purple"), 
      # font.caption  = c(14, "plain", "orange"),        
      # font.x        = c(14, "bold.italic", "red"),          
      # font.y        = c(14, "bold.italic", "darkred"),      
      # font.xtickslab = c(12, "plain", "darkgreen"),
      #legend = c(0.8,1.1)
      legend = "top"
    )
    ggsurv$plot <- ggsurv$plot+
      scale_x_continuous(breaks = pretty_breaks(6),expand = c(0,0))+
      scale_y_continuous(breaks = pretty_breaks(6),expand = c(0,0))+
      theme(legend.direction = "horizontal")+
      labs(color=keyname)+
      coord_cartesian(clip = 'off')
    ggsurv$plot <- flush_ticks(ggsurv$plot)
    ggsurv
    
    if(is.null(filename)){
      filename=paste0(keyname,"-Survival.pdf")
    }
    
    ggsave(filename,width = width,height = height,device = cairo_pdf)
    
    
  }
  return(c(suvpvalue$p.value,suvpvalue$lab))
}

SurvZTWm <- function(suvdata_input,plot=TRUE,legend.labs=NULL,keyname="Strata",filename=NULL,width = 7,height = 5.5,pvalsize=4,pvalab=NULL){
  
  fit <- coxph(Surv(Survival_Month, Death) ~ Key+Age+Gender+Stage+Smoking+Histology, data = suvdata_input) ### overall 
  suvpvalue <- tidy(fit,conf.int = T) %>% filter(str_detect(term,"Key")) %>% mutate(lab=paste0("P = ",format(p.value,digits = 2,scientific = FALSE),"\nHR = ",format(exp(estimate),digits = 2,scientific = FALSE)," (95% CI: ",format(exp(conf.low),digits = 2,scientific = FALSE),"-",format(exp(conf.high),digits = 2,scientific = FALSE),")")) %>% select(term,p.value,lab) %>% 
    mutate(term=str_remove(term,"Key")) %>% mutate(lab=paste0(term,":\n",lab))
  if(is.null(pvalab)){
    pvalab <- paste0(suvpvalue$lab,collapse = "\n")
  }
  
  fit2 <- coxph(Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Smoking+Histology, data = suvdata_input) ## strata curve ##
  ##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata_input)  ## all strata, wrong way
  suvfit <- survfit(fit2)
  #plot(suvfit)
  suvfit$call$formula <- as.formula("Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Smoking+Histology")
  
  if(is.null(legend.labs)) 
  {
    legend.labs <- paste0(levels(suvdata_input$Key)," (", suvfit$n,")")
  }
  keyname=keyname
  ggsurv <- ggsurvplot(
    suvfit,                     # survfit object with calculated statistics.
    data = suvdata_input,             # data used to fit survival curves.
    risk.table = FALSE,       # show risk table.
    pval = pvalab,             # show p-value of log-rank test.
    pval.size=pvalsize,
    #pval = TRUE,
    #pval.method
    #conf.int = TRUE,         # show confidence intervals for 
    # point estimates of survival curves.
    palette = c("black",pal_d3()(length(levels(suvdata_input$Key))-1)),
    #xlim = c(0,210),         # present narrower X axis, but not affect
    # survival estimates.
    xlab = "Time in months",   # customize X axis label.
    break.time.by = 50,     # break X axis in time intervals by 500.
    ggtheme = theme_ipsum_rc(axis_title_size = 14,base_size = 12,axis_title_just = "m",axis = "XY",axis_col = "black"), # customize plot and risk table with a theme.
    font.family="Roboto Condensed" ,
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.height = 0.25, # the height of the risk table
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    # in legend of risk table.
    censor=TRUE,
    ncensor.plot = FALSE,      # plot the number of censored subjects at time t
    ncensor.plot.height = 0.25,
    conf.int.style = "step",  # customize style of confidence intervals
    surv.median.line = "none",
    legend.labs = legend.labs
  )
  # ggsurv$plot <- ggsurv$plot + labs(
  #   title    = "Survival curve",                     
  #   subtitle = "Based on Kaplan-Meier estimates"
  # )
  ggsurv <- ggpar(
    ggsurv,
    # font.title    = c(15, "bold", "darkblue"),         
    # font.subtitle = c(14, "bold.italic", "purple"), 
    # font.caption  = c(14, "plain", "orange"),        
    # font.x        = c(14, "bold.italic", "red"),          
    # font.y        = c(14, "bold.italic", "darkred"),      
    # font.xtickslab = c(12, "plain", "darkgreen"),
    #legend = c(0.8,1.1)
    legend = "top"
  )
  ggsurv$plot <- ggsurv$plot+
    scale_x_continuous(breaks = pretty_breaks(6),expand = c(0,0))+
    scale_y_continuous(breaks = pretty_breaks(6),expand = c(0,0))+
    theme(legend.direction = "horizontal")+
    labs(color=keyname)+
    coord_cartesian(clip = 'off')
  ggsurv$plot <- flush_ticks(ggsurv$plot)
  ggsurv
  
  if(is.null(filename)){
    filename=paste0(keyname,"-Survival.pdf")
  }
  
  ggsave(filename,width = width,height = height,device = cairo_pdf)
  
  return(suvpvalue)
}

SurvZTWms <- function(suvdata_input,plot=TRUE,legend.labs=NULL,keyname="Strata",filename=NULL,gcolors=NULL,width = 7,height = 5.5,pvalsize=4,pvalab=NULL){
  
  fit <- coxph(Surv(Survival_Month, Death) ~ Key+Age+Gender+Stage+Histology, data = suvdata_input) ### overall 
  suvpvalue <- tidy(fit,conf.int = T) %>% filter(str_detect(term,"Key")) %>% mutate(lab=paste0("P = ",format(p.value,digits = 2,scientific = FALSE),"\nHR = ",format(exp(estimate),digits = 2,scientific = FALSE)," (95% CI: ",format(exp(conf.low),digits = 2,scientific = FALSE),"-",format(exp(conf.high),digits = 2,scientific = FALSE),")")) %>% select(term,p.value,lab) %>% 
    mutate(term=str_remove(term,"Key")) %>% mutate(lab=paste0(term,":\n",lab))
  
  print(suvpvalue)
  
  suvdata_input <- suvdata_input %>% mutate(Survival_Month = if_else(Survival_Month > 60, NA, Survival_Month))
  
  if(is.null(pvalab)){
    pvalab <- paste0(suvpvalue$lab,collapse = "\n")
  }
  
  fit2 <- coxph(Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Histology, data = suvdata_input) ## strata curve ##
  ##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata_input)  ## all strata, wrong way
  suvfit <- survfit(fit2)
  #plot(suvfit)
  suvfit$call$formula <- as.formula("Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Histology")
  
  if(is.null(legend.labs)) 
  {
    legend.labs <- paste0(levels(suvdata_input$Key)," (", suvfit$n,")")
  }
  if(is.null(gcolors)){
    gcolors <- pal_jama()(length(levels(suvdata_input$Key)))
  }
  
  keyname=keyname
  ggsurv <- ggsurvplot(
    suvfit,                     # survfit object with calculated statistics.
    data = suvdata_input,             # data used to fit survival curves.
    risk.table = FALSE,       # show risk table.
    pval = pvalab,             # show p-value of log-rank test.
    pval.size=pvalsize,
    #pval = TRUE,
    #pval.method
    #conf.int = TRUE,         # show confidence intervals for 
    # point estimates of survival curves.
    palette = gcolors,
    #xlim = c(0,210),         # present narrower X axis, but not affect
    # survival estimates.
    xlab = "Time in months",   # customize X axis label.
    break.time.by = 50,     # break X axis in time intervals by 500.
    ggtheme = theme_ipsum_rc(axis_title_size = 14,base_size = 12,axis_title_just = "m",axis = "XY",axis_col = "black"), # customize plot and risk table with a theme.
    font.family="Roboto Condensed" ,
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.height = 0.25, # the height of the risk table
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    # in legend of risk table.
    censor=TRUE,
    ncensor.plot = FALSE,      # plot the number of censored subjects at time t
    ncensor.plot.height = 0.25,
    conf.int.style = "step",  # customize style of confidence intervals
    surv.median.line = "none",
    legend.labs = legend.labs
  )
  # ggsurv$plot <- ggsurv$plot + labs(
  #   title    = "Survival curve",                     
  #   subtitle = "Based on Kaplan-Meier estimates"
  # )
  ggsurv <- ggpar(
    ggsurv,
    # font.title    = c(15, "bold", "darkblue"),         
    # font.subtitle = c(14, "bold.italic", "purple"), 
    # font.caption  = c(14, "plain", "orange"),        
    # font.x        = c(14, "bold.italic", "red"),          
    # font.y        = c(14, "bold.italic", "darkred"),      
    # font.xtickslab = c(12, "plain", "darkgreen"),
    #legend = c(0.8,1.1)
    legend = "top"
  )
  ggsurv$plot <- ggsurv$plot+
    scale_x_continuous(breaks = pretty_breaks(6),expand = c(0,0))+
    scale_y_continuous(breaks = pretty_breaks(6),expand = c(0,0))+
    theme(legend.direction = "horizontal")+
    labs(color=keyname)+
    coord_cartesian(clip = 'off')
  ggsurv$plot <- flush_ticks(ggsurv$plot)
  ggsurv
  
  if(is.null(filename)){
    filename=paste0(keyname,"-Survival.pdf")
  }
  
  ggsave(filename,width = width,height = height,device = cairo_pdf)
  
  return(suvpvalue)
}


SurvZTWmsH <- function(suvdata_input,plot=TRUE,legend.labs=NULL,keyname="Strata",filename=NULL,gcolors=NULL,width = 7,height = 5.5,pvalsize=4,pvalab=NULL){
  
  fit <- coxph(Surv(Survival_Month, Death) ~ Key+Age+Gender+Stage, data = suvdata_input) ### overall 
  suvpvalue <- tidy(fit,conf.int = T) %>% filter(str_detect(term,"Key")) %>% mutate(lab=paste0("P = ",format(p.value,digits = 2,scientific = FALSE),"\nHR = ",format(exp(estimate),digits = 2,scientific = FALSE)," (95% CI: ",format(exp(conf.low),digits = 2,scientific = FALSE),"-",format(exp(conf.high),digits = 2,scientific = FALSE),")")) %>% select(term,p.value,lab) %>% 
    mutate(term=str_remove(term,"Key")) %>% mutate(lab=paste0(term,":\n",lab))
  
  print(suvpvalue)
  
  if(is.null(pvalab)){
    pvalab <- paste0(suvpvalue$lab,collapse = "\n")
  }
  
  fit2 <- coxph(Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage, data = suvdata_input) ## strata curve ##
  ##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata_input)  ## all strata, wrong way
  suvfit <- survfit(fit2)
  #plot(suvfit)
  suvfit$call$formula <- as.formula("Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage")
  
  if(is.null(legend.labs)) 
  {
    legend.labs <- paste0(levels(suvdata_input$Key)," (", suvfit$n,")")
  }
  if(is.null(gcolors)){
    gcolors <- pal_jama()(length(levels(suvdata_input$Key)))
  }
  
  keyname=keyname
  ggsurv <- ggsurvplot(
    suvfit,                     # survfit object with calculated statistics.
    data = suvdata_input,             # data used to fit survival curves.
    risk.table = FALSE,       # show risk table.
    pval = pvalab,             # show p-value of log-rank test.
    pval.size=pvalsize,
    #pval = TRUE,
    #pval.method
    #conf.int = TRUE,         # show confidence intervals for 
    # point estimates of survival curves.
    palette = gcolors,
    #xlim = c(0,210),         # present narrower X axis, but not affect
    # survival estimates.
    xlab = "Time in months",   # customize X axis label.
    break.time.by = 50,     # break X axis in time intervals by 500.
    ggtheme = theme_ipsum_rc(axis_title_size = 14,base_size = 12,axis_title_just = "m",axis = "XY",axis_col = "black"), # customize plot and risk table with a theme.
    font.family="Roboto Condensed" ,
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.height = 0.25, # the height of the risk table
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    # in legend of risk table.
    censor=TRUE,
    ncensor.plot = FALSE,      # plot the number of censored subjects at time t
    ncensor.plot.height = 0.25,
    conf.int.style = "step",  # customize style of confidence intervals
    surv.median.line = "none",
    legend.labs = legend.labs
  )
  # ggsurv$plot <- ggsurv$plot + labs(
  #   title    = "Survival curve",                     
  #   subtitle = "Based on Kaplan-Meier estimates"
  # )
  ggsurv <- ggpar(
    ggsurv,
    # font.title    = c(15, "bold", "darkblue"),         
    # font.subtitle = c(14, "bold.italic", "purple"), 
    # font.caption  = c(14, "plain", "orange"),        
    # font.x        = c(14, "bold.italic", "red"),          
    # font.y        = c(14, "bold.italic", "darkred"),      
    # font.xtickslab = c(12, "plain", "darkgreen"),
    #legend = c(0.8,1.1)
    legend = "top"
  )
  ggsurv$plot <- ggsurv$plot+
    scale_x_continuous(breaks = pretty_breaks(6),expand = c(0,0))+
    scale_y_continuous(breaks = pretty_breaks(6),expand = c(0,0))+
    theme(legend.direction = "horizontal")+
    labs(color=keyname)+
    coord_cartesian(clip = 'off')
  ggsurv$plot <- flush_ticks(ggsurv$plot)
  ggsurv
  
  if(is.null(filename)){
    filename=paste0(keyname,"-Survival.pdf")
  }
  
  ggsave(filename,width = width,height = height,device = cairo_pdf)
  
  return(suvpvalue)
}




SurvZTWms_tp53 <- function(suvdata_input,plot=TRUE,legend.labs=NULL,keyname="Strata",gcolors=NULL,filename=NULL,width = 7,height = 5.5,pvalsize=4,pvalab=NULL){
  
  fit <- coxph(Surv(Survival_Month, Death) ~ Key+Age+Gender+Stage+Histology+TP53_Status, data = suvdata_input) ### overall 
  suvpvalue <- tidy(fit,conf.int = T) %>% filter(str_detect(term,"Key")) %>% mutate(lab=paste0("P = ",format(p.value,digits = 2,scientific = FALSE),"\nHR = ",format(exp(estimate),digits = 2,scientific = FALSE)," (95% CI: ",format(exp(conf.low),digits = 2,scientific = FALSE),"-",format(exp(conf.high),digits = 2,scientific = FALSE),")")) %>% select(term,p.value,lab) %>% 
    mutate(term=str_remove(term,"Key")) %>% mutate(lab=paste0(term,":\n",lab))
  if(is.null(pvalab)){
    pvalab <- paste0(suvpvalue$lab,collapse = "\n")
  }
  
  suvdata_input <- suvdata_input %>% mutate(Survival_Month = if_else(Survival_Month > 60, NA, Survival_Month))
  
  fit2 <- coxph(Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Histology+TP53_Status, data = suvdata_input) ## strata curve ##
  ##suvfit <- survfit(Surv(survival_months, death) ~ group+age+sex+stage, data = suvdata_input)  ## all strata, wrong way
  suvfit <- survfit(fit2)
  #plot(suvfit)
  suvfit$call$formula <- as.formula("Surv(Survival_Month, Death) ~ strata(Key)+Age+Gender+Stage+Histology+TP53_Status")
  
  if(is.null(legend.labs)) 
  {
    legend.labs <- paste0(levels(suvdata_input$Key)," (", suvfit$n,")")
  }
  if(is.null(gcolors)){
    gcolors <- pal_jama()(length(levels(suvdata_input$Key)))
  }
  keyname=keyname
  ggsurv <- ggsurvplot(
    suvfit,                     # survfit object with calculated statistics.
    data = suvdata_input,             # data used to fit survival curves.
    risk.table = FALSE,       # show risk table.
    pval = pvalab,             # show p-value of log-rank test.
    pval.size=pvalsize,
    #pval = TRUE,
    #pval.method
    #conf.int = TRUE,         # show confidence intervals for 
    # point estimates of survival curves.
    palette = gcolors,
    #xlim = c(0,210),         # present narrower X axis, but not affect
    # survival estimates.
    xlab = "Time in months",   # customize X axis label.
    break.time.by = 50,     # break X axis in time intervals by 500.
    ggtheme = theme_ipsum_rc(axis_title_size = 14,base_size = 12,axis_title_just = "m",axis = "XY",axis_col = "black"), # customize plot and risk table with a theme.
    font.family="Roboto Condensed" ,
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.height = 0.25, # the height of the risk table
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    # in legend of risk table.
    censor=TRUE,
    ncensor.plot = FALSE,      # plot the number of censored subjects at time t
    ncensor.plot.height = 0.25,
    conf.int.style = "step",  # customize style of confidence intervals
    surv.median.line = "none",
    legend.labs = legend.labs
  )
  # ggsurv$plot <- ggsurv$plot + labs(
  #   title    = "Survival curve",                     
  #   subtitle = "Based on Kaplan-Meier estimates"
  # )
  ggsurv <- ggpar(
    ggsurv,
    # font.title    = c(15, "bold", "darkblue"),         
    # font.subtitle = c(14, "bold.italic", "purple"), 
    # font.caption  = c(14, "plain", "orange"),        
    # font.x        = c(14, "bold.italic", "red"),          
    # font.y        = c(14, "bold.italic", "darkred"),      
    # font.xtickslab = c(12, "plain", "darkgreen"),
    #legend = c(0.8,1.1)
    legend = "top"
  )
  ggsurv$plot <- ggsurv$plot+
    scale_x_continuous(breaks = pretty_breaks(6),expand = c(0,0))+
    scale_y_continuous(breaks = pretty_breaks(6),expand = c(0,0))+
    theme(legend.direction = "horizontal")+
    labs(color=keyname)+
    coord_cartesian(clip = 'off')
  ggsurv$plot <- flush_ticks(ggsurv$plot)
  ggsurv
  
  if(is.null(filename)){
    filename=paste0(keyname,"-Survival.pdf")
  }
  
  ggsave(filename,width = width,height = height,device = cairo_pdf)
  
  return(suvpvalue)
}










