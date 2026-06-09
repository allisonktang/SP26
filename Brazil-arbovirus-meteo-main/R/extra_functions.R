## extra functions

data_processing <- function(df){
  # filter out municipalities with no local transmission
  df |> 
    select(MUNICIP, CASES, YEAR)|>
    group_by(MUNICIP, YEAR) |>
    summarise(CASES=sum(CASES)) |> 
    mutate(bin=ifelse(CASES>=5,1,0), # continuous low level transmission
           bin2=ifelse(CASES>20,1,0)) |> # at least one year with many cases
    group_by(MUNICIP) |>
    summarise(years_n=sum(bin), years_n2=sum(bin2))|>
    filter(years_n==5|years_n2>=1)-> df2
  df |> 
    filter(MUNICIP %in% unique(df2$MUNICIP))-> df
  
  return(df)
}

run_baseline <- function(data, GRAPH){
  
  graph<-inla.read.graph(GRAPH)
  
  formula <- CASES ~ 1 + lag_log_cases + 
    f(WEEK_NO, replicate=STATE_NO, model="rw1", 
      cyclic=TRUE) + 
    f(GEOM, model="bym2", 
      graph=graph, 
      replicate=YEAR_NO)
  
  model <- inla(formula = formula,
                offset=log(pop),
                data = data, 
                family="nbinomial",
                control.family=list(link='log'),
                control.compute = list(waic = TRUE, cpo = TRUE,
                                       config = TRUE, return.marginals.predictor = TRUE),
                verbose = F)
  
  return(model)
  
}

select_variables <- function(data, type="all"){
  
variables <- readRDS("data/variables.RDS")
if(type=="temperature") variables <- variables[(str_detect(variables,"Tmin") | str_detect(variables,"Tmax") | str_detect(variables,"TR"))]
if(type=="subset") variables <- variables[c(8,11,27)] 
return(variables)

}


process_univar_names <- function(res){
  
  res |> 
    mutate(group = str_replace_all(var, "(_|\\d)", " ")) |>
    mutate(group = case_when(
      var == "ENSO_anom" ~ "ENSO",
      substring(var,1,4) =="cent" ~ "centrality",
      substring(var,1,4) =="spei" ~ "SPEI", .default=group)) -> res
  
  return(res)
}

process_univar_results_for_plot <- function(){

  a<-b<-c<-NULL
  
  if(file.exists("outputs/univar/chikv_univar_all.RDS")){
  a<-readRDS("outputs/univar/chikv_univar_all.RDS")
  a$path<-"chikv"
  a<-a[order(a$waic),]
  a$ORDER<-1:nrow(a)
  }
  
  if(file.exists("outputs/univar/denv_univar_all.RDS")){
  b<-readRDS("outputs/univar/denv_univar_all.RDS")
  b$path<-"denv"
  b<-b[order(b$waic),]
  b$ORDER<-1:nrow(b)
  }
  
  if(file.exists("outputs/univar/zikv_univar_all.RDS")){
  c<-readRDS("outputs/univar/zikv_univar_all.RDS")
  c$path<-"zikv"
  c<-c[order(c$waic),]
  c$ORDER<-1:nrow(c)
  }
  
  data <- rbind(a,b,c)
  
  data |> 
    mutate(var = case_when(str_detect(var,"_03") ~ "Previous 0-3 weeks", 
                           str_detect(var,"_02") ~ "Previous 0-2 weeks", 
                           str_detect(var,"_01") ~ "Previous 0-1 weeks", 
                           str_detect(var,"_3") ~ "Week t -3 ", 
                           str_detect(var,"_2") ~ "Week t -2", 
                           str_detect(var,"_1") ~ "Week t -1", 
                           str_detect(var,"_0") ~ "Week t", 
                           .default = var)) |>
    mutate(var = case_when(var %in% c("raindays_cons", "drydays_cons", "raindays", "raindays_heavy",
                                      "raindays_heavy+", "TR", "RH_min_min", "Tmin_min", "RH_max_max",
                                      "Tmax_max", "Tmin", "Tmax", "Tmean", "RH_min", "RH_max", "RH_mean",
                                      "Rf_mean", "Rf_cum") ~ "Week t", 
                           .default = var)) -> data
  data$var <-sub("ENSO_anom","ENSO", data$var)
  data$var <-sub("MIR","Mid infra-red index", data$var)
  data$var <-sub("income","Median household\nincome", data$var)
  data$var <-sub("illiteracy","Rate of adult\nilliteracy", data$var)
  data$var <-sub("sanitation","Proportion of\nhouseholds with\naccess to adequate\nsanitation", data$var)
  data$var <-sub("cent_between","Between centrality", data$var)
  data$var <-sub("cent_alpha","Alpha centrality", data$var)
  data$var <-sub("cent_degree","Degree centrality", data$var)
  data$var <-sub("spei12","SPEI 12 months", data$var)
  data$var <-sub("spei1","SPEI 1 month", data$var)
  data$var <-sub("spei3","SPEI 3 months", data$var)
  data$var <-sub("spei6","SPEI 6 months", data$var)
  
  data$path <-toupper(data$path)

return(data)
}
