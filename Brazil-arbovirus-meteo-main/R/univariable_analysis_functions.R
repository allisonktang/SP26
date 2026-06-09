#### univariable analysis

run_univar_analysis <- function(data, variables, save_progress=TRUE, pathogen){
  
  options(mc.cores = parallel::detectCores())
  
  list <- list()
  models <- list()
  var_names <- variables
  
  graph <- inla.read.graph(GRAPH)
  
  for(i in 1:length(var_names)){
    
    if(var_names[i] %in% c("flood_0","flood_1","flood_2" ,"flood_3","flood_01","flood_02",           
                           "flood_03")){ # categorical variable
      formula <- CASES ~ 1 + lag_log_cases + 
        f(WEEK_NO, replicate=STATE_NO, model="rw1", 
          cyclic=TRUE) + 
        f(GEOM, model="bym2", 
          graph=graph, 
          replicate=YEAR_NO)+
        f(data[,var_names[i]],model="iid")
    
    }else{
      formula <- CASES ~ 1 + lag_log_cases + 
        f(WEEK_NO, replicate=STATE_NO, model="rw1", 
          cyclic=TRUE) + 
        f(GEOM, model="bym2", 
          graph=graph, 
          replicate=YEAR_NO)+
        data[,var_names[i]]
    }
    
    models[[i]] <- inla(formula = formula,
                        offset=log(pop),
                        data = data, 
                        family="nbinomial",
                        control.family=list(link='log'),
                        control.compute = list(waic = TRUE),
                        control.inla=list(int.strategy = "eb"), 
                        verbose = FALSE,
                        safe = TRUE)
    
    waic<- models[[i]]$waic$waic
    beta<-models[[i]]$summary.fixed[3,]
    list[[i]]<-data.frame(var = variables[i], 
                          waic = waic,
                          beta_mean = beta$mean, 
                          beta_low = beta$`0.025quant`, 
                          beta_upp = beta$`0.975quant`)
    models[[i]]<-NA
    
    print(paste0("finished model ", i, " out of ", length(variables)))
  }
  result<- bind_rows(list)
  
  if(save_progress){
   saveRDS(result, paste0("outputs/univar/temp/",pathogen,"_", i,".RDS")) 
  }
  return(result)
  
  }
  

## bind univariable model results together
binding_univar_results <-function(pathogen, file_path="outputs/univar/temp/"){
  filenames <- list.files(file_path, pattern=pathogen, full.names=TRUE)
  ldf <- lapply(filenames, readRDS)
  my.fcn <- function(x){
    ifelse("cat" %in% colnames(x),x$cat<-as.character(x$cat), x$cat<-NA)
    return(x)
  }
  test<-lapply(ldf,my.fcn)
  res<- bind_rows(test, id=NULL)
  res<- res[!duplicated(res$var),]

  res <- process_univar_names(res)
  
  return(res)
}


## choose the best univariable model
best_univar <- function(univar_results_table){
  univar_results_table %>% filter(beta_low < 0 & beta_upp <0 | beta_low > 0 & beta_upp >0) -> sig
  sig$var[sig$waic==min(sig$waic)] -> var
  waic<-min(sig$waic)
  return(c(var,waic))
}


## plot univariable beta results
plot_univar_betas <- function(choice = "temperature"){
  
  if(choice!="temperature" & choice !="humidity" & choice !="rain" & choice !="socioeconomic"){
    warning("choice must be temperature or humidity or rain or socioeconomic")
  }
  
  data <- process_univar_results_for_plot()
  plots <- list()
  
  if(choice=="temperature"){
  sub_tot <- data |> filter(str_detect(group, "Tmin") | str_detect(group, "Tmax") | str_detect(group, "Tmean") | str_detect(group, "TR"))
  }
  if(choice=="humidity"){
  sub_tot <- data |> filter(str_detect(group, "RH"))
  }
  if(choice=="rain"){
  sub_tot <- data |> filter(str_detect(group, "rf") | str_detect(group, "SPEI") | str_detect(group, "ENSO") | str_detect(group, "days"))
  }
  if(choice=="socioeconomic"){
    sub_tot <- data |> filter(str_detect(group, "cent") | str_detect(group, "sanitation") | str_detect(group, "income") | str_detect(group, "illiteracy"))
  }
  for(i in 1:length(unique(sub_tot$group))){
    sub<-sub_tot[sub_tot$group == unique(sub_tot$group)[i],]
    lim1<- ifelse(min(sub$beta_low)< -max(sub$beta_upp),
                  min(sub$beta_low), - max(sub$beta_upp))
    lim2<-ifelse(min(sub$beta_low)< -max(sub$beta_upp),
                 - min(sub$beta_low), max(sub$beta_upp))
    plots[[i]]<-ggplot(sub)+geom_hline(aes(yintercept = 0), linetype="dashed") +
      geom_point(aes(x=var, y=beta_mean, col=path), 
                 size=2,position=position_dodge(width=0.2)) +
      geom_pointrange(aes(x=var, y=beta_mean, 
                          ymin=beta_low, 
                          ymax=beta_upp, col=path), 
                      size=0.3,position=position_dodge(width=0.2)) + 
      labs(x = NULL, y = "Beta value", title = unique(sub_tot$group)[i]) +
      coord_flip(ylim=c(lim1,lim2)) + #makes horizontal
      theme_classic() +
      theme(legend.title=element_blank(),legend.position="none",
            plot.title = element_text(hjust=0.5),
            plot.margin = margin(10, 12, 10, 10))#, 
    
  }
  plot1 <- sjPlot::plot_grid(plots)
  
  return(plot1)
  }
