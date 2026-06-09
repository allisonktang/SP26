### plotting functions


## data visualisation

## temporal heterogeneity - average cases per week number per state
temporal_cases_state <- function(data_type){
list<-list()
for(i in 1:3){
  path<-unique(c("zikv","chikv", "denv"))[i]
  df<-readRDS(paste0("data/", path, "_", data_type,".RDS"))
  df |> mutate(Inc = ((CASES+1e-05)/pop)*100000) -> df
  df |> select(STATE,WEEK_NO,CASES,YEAR, Inc) |> 
    group_by(STATE,WEEK_NO) |>
    reframe(Inc = mean((Inc)),
            cases_50 = mean((CASES))) -> df
  df$path <- path
  list[[i]]<-df
  remove("df")
}

combo <- bind_rows(list)
GRID <- br_states_grid1
combo$path <-toupper(combo$path)
combo <- combo[combo$WEEK_NO<53,]

week_RE <- ggplot(combo) + facet_geo(facets = "STATE", 
                                   grid = br_states_grid1,scales="free")+
  geom_line(aes(x = WEEK_NO, y = log(Inc), col = path), linewidth=1) +
  xlab("\nWeek number\n") +
  ylab("Log(Incidence)") +
  scale_y_continuous(breaks=scales::breaks_pretty(n=4))+
  scale_x_continuous(breaks=seq(from=1,to=53, by=10), labels=seq(from=1,to=53, by=10))+
  theme(legend.position = c(0.8,0.1),
        panel.background = element_rect(fill = NA), 
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.key  = element_blank(),
        panel.border = element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title =element_text(size=20),
        strip.text = element_text(size=12),
        strip.background = element_blank())

return(week_RE)

}


## spatial heterogeneity - average cases municipality
spatial_cases_municip <- function(data_type){
list<-list()
for(i in 1:3){
  path<-unique(c("zikv","chikv", "denv"))[i]
  
  df<-readRDS(paste0("data/", path, "_", data_type,".RDS"))
  df |> mutate(Inc = ((CASES+1e-05)/pop)*100000) -> df
  
  df |> select(STATE,MUNICIP,WEEK_NO,CASES,YEAR, Inc) |>
    group_by(MUNICIP) |>
    reframe(Inc = mean((Inc)),
            cases_50 = mean((CASES))) -> df
  
  df$path <- toupper(path)
  list[[i]]<-df
  remove("df")
}
combo <- bind_rows(list)

combo$loginc<-log(combo$Inc)
space <- left_join(shapefile, combo, by = c("code_muni" = "MUNICIP"))
space<-space[is.na(space$path)==F,]

space <- ggplot() + facet_wrap("path",ncol=2)+
  geom_sf(data = space, aes(fill = loginc),
          lwd = 0, color = NA) +
  scale_fill_distiller(palette = "RdYlBu",
                       direction = -1
  ) +
  labs(fill="Log(Incidence)", col="") +
  theme_void() +
  theme(strip.text = element_text(size=15),
        legend.title =  element_text(size=12),
        legend.text =  element_text(size=10),
        legend.position = "bottom")+
  guides(fill=guide_colourbar(title.position="top",
                              title.hjust=0.5,
                              barwidth=8)) 

return(space)
}


## spatiotemporal heterogeneity - heatmaps
heatmaps <- function(shapefile = shapefile, pathogen = "denv", data_type = "exp_resid"){
  
  if(pathogen!="denv" & pathogen !="zikv"& pathogen !="chikv"){
    warning("pathogen must be chikv, denv, or zikv")
  }
  
  df <- readRDS(paste0("data/", pathogen, "_", data_type,".RDS"))
  df$WEEK <- (paste0(df$YEAR,"-W",sprintf("%02d",(df$WEEK_NO)),"-1"))
  df |> select(MUNICIP, WEEK, CASES, pop, GEOM, YEAR) |> mutate(Inc = ((CASES+1e-05)/pop)*100000) -> df
  
  cent <- st_centroid(shapefile)
  cent.coords <- st_coordinates(cent)
  y <- data.frame(MUNICIP = shapefile$code_muni, 
                  NAME = shapefile$name_muni, 
                  x = cent.coords[,1],
                  y = cent.coords[,2])
  y$order <- 1:nrow(y)
  
  df$logInc <- log(df$Inc)
  df$date<- ISOweek2date(df$WEEK)
  
  
  # geom
  zD <- merge(y, df, by = "MUNICIP", all.x = F, all.y = F)
    zD <- zD[order(zD$x), ] # ordering geographically
    zD$MUNICIP<-as.character(zD$MUNICIP)
    zD <- zD |>
      mutate(MUNICIP = fct_reorder(MUNICIP, desc(x))) 
  
    min_year<-min(zD$YEAR)
    max_year<-max(zD$YEAR)

plot <- 
  ggplot(zD,aes(x=date, y=MUNICIP, fill=logInc))  +
  geom_tile()+
  labs(title = toupper(pathogen),
       fill = "Log(Incidence)",
       x = NULL,
       y = "Municipality") +
  theme(axis.text.x = element_text(angle=90, hjust=1, size=10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1.2, 'cm'),
        panel.background = element_rect(fill = NA), 
        panel.border = element_rect(colour = "black", fill=NA),
        plot.title = element_text(hjust = 0.5, size=15),
        axis.title  = element_text(size=15))+
  scale_fill_distiller(palette = "RdYlBu", 
                       direction = -1,
                       na.value="white" ,
                       guide = "colourbar",
                       limits =c(-18.4,9))+
  scale_x_date(date_breaks="1 year",date_labels = "%Y") +
  guides(fill=guide_colourbar(title.position="top",
                              title.hjust=0.5,
                              barwidth=8)) 

return(plot)
}


## plotting multivariable model random effects - temporal 
plotting_RE_week_number <- function(model, data, pathogen){
  
    random <- model$summary.random$WEEK_NO
    colnames(random) <- c("ID","mean","sd", "quant2.5", "0.5quant" ,  "quant97.5", "mode", "kld"   )
    
    week_effects <- data.table(cbind(rep(unique(data$STATE), 
                                          each = length(unique(data$WEEK_NO))),
                                      model$summary.random$WEEK_NO
    ))
    names(week_effects)[1:2] <- c("state_code", "Week")
    week_effects$path <- toupper(pathogen)
    
    GRID <- br_states_grid1
  
    plot<-ggplot(week_effects) + facet_geo(facets = "state_code", 
                                  grid = br_states_grid1)+
    geom_ribbon(aes(x = Week, ymin = `0.025quant`, ymax = `0.975quant`, 
                    fill = path), alpha = 0.5) + 
    geom_line(aes(x = Week, y = `mean`, col = path)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    xlab("Week number") +
    ylab("Contribution to the linear predictor") +
    coord_cartesian(ylim=c(-1,1.2))+
    scale_y_continuous()+
    scale_x_continuous(breaks=seq(from=1,to=53, by=10))+
    theme(#legend.position = "none",
      panel.background = element_rect(fill = NA), 
      legend.title = element_blank(),
      legend.text = element_text(size=10),
      panel.border = element_rect(colour = "black", fill=NA),
      #plot.title = element_text(hjust = 0.5),
      axis.text.x=element_text(size=8),
      axis.text.y=element_text(size=8),
      strip.text = element_text(size=10),
      strip.background = element_blank())#+
  
  return(plot)
}

## plotting multivariable model random effects - spatial 
plotting_RE_municip <- function(model, data, shapefile){
  
  FILTER<-unique(data$MUNICIP)
  
  space <- data.table(model$summary.random$GEOM)
  space$year <-rep(min(data$YEAR):max(data$YEAR), 
                    each = 2*5572)
  space$re <- rep(c(rep("iid_adj",5572),rep("adj",5572)),length(unique(data$YEAR_NO)))
  
  space <- space[space$re == "iid_adj",]
  space$micro_code <- rep(unique(shapefile$code_muni), length(seq(min(data$YEAR_NO):max(data$YEAR_NO))))
  
  # Add the map geometry to the space dataframe
  space <- left_join(shapefile, space, by = c("code_muni" = "micro_code"))
  
  space |> mutate(mean=ifelse(code_muni%in%FILTER, mean, NA))->space
  
  space_effects <- ggplot() + 
    geom_sf(data = space, aes(fill = mean), lwd = 0, color = NA) +
    scale_fill_distiller(palette = "RdYlBu",
                         na.value="lightgrey",
                         direction = -1) +
    labs(fill="", col="",title =paste0(toupper(pathogen), "\n \n")) +
    theme_void() +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust=0.5))+
    facet_wrap(~space$year)
  
  return(space_effects)
}


## plotting model fit
plot_model_fit <- function(data, pathogen, nsamp, samples, level="national"){
  
  data$date <- ISOweek2date(paste0(data$YEAR,"-W",sprintf("%02d", data$WEEK_NO),"-1"))
  
  if(level=="national"){
    
  store<-matrix(NA, length(unique(data$date)), nsamp)
  for(n in 1:nsamp){
    new<-data
    new$est_new<-NA
    new$est_new<-exp(samples[,n]) 
    
    new |> select(date, MUNICIP, STATE, CASES, YEAR, WEEK_NO, pop, est_new)|>
      group_by(date) |>
      summarise(est_new = (sum(est_new)/sum(pop))*100000, .groups="keep")->new
    store[,n]<-new$est_new
  }
  
  data |> select(date, MUNICIP, STATE, CASES, YEAR, WEEK_NO, pop)|>
    group_by(date) |>
    summarise(CASES=(sum(CASES)/sum(pop))*100000, .groups="keep") -> df_nat
  
  df_nat$median <- apply(store, 1, median)
  df_nat$lci <- apply(store, 1, quantile, probs = c(0.025))
  df_nat$uci <- apply(store, 1, quantile, probs = c(0.975))
  
  
  plot<- ggplot(df_nat)+ 
    geom_ribbon(aes(x = date, ymin = (lci), ymax = (uci)), 
                fill = "#D37295", alpha = 0.5) + 
    geom_line(aes(x = date, y = (median), col = "#D37295")) +
    geom_line(aes(x = date, y = (CASES), col = "#499894")) +
    xlab("Time") +
    ylab("Incidence per 100000") +
    scale_colour_identity(name = "",
                          breaks = c("#499894", "#D37295"),
                          labels = c("Observed", "Posterior predicted"),
                          guide = "legend") +
    theme_bw() + 
    theme(axis.text = element_text(size = 12), 
          legend.position = "bottom",
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 15),
          axis.title = element_text(size =15))
  
  }
  
  if(level=="state"){
    data |> select(date,MUNICIP, STATE, CASES, YEAR, WEEK_NO, pop)|>
      group_by(date, STATE) |>
      summarise(CASES=(sum(CASES)/sum(pop))*100000, .groups="keep") -> df_st
    
    store2<-matrix(NA, nrow(df_st), nsamp)
    for(n in 1:nsamp){
      new<-data
      new$est_new<-NA
      new$est_new<-exp(samples[,n]) 
      
      new |> select(date,MUNICIP, STATE, CASES, YEAR, WEEK_NO, pop, est_new)|>
        group_by(date, STATE) |>
        summarise(est_new = (sum(est_new)/sum(pop))*100000, .groups="keep") -> new
      store2[,n]<-new$est_new
      
    }
    
    df_st$median <- apply(store2, 1, median)
    df_st$lci <- apply(store2, 1, quantile, probs = c(0.025))
    df_st$uci <- apply(store2, 1, quantile, probs = c(0.975))
    
    plot <- ggplot(df_st)+ facet_geo("STATE",grid=br_states_grid1, scales="free_y")+
      geom_ribbon(aes(x = date, ymin = (lci), ymax = (uci)), 
                  fill = "#D37295", alpha = 0.5) + 
      geom_line(aes(x = date, y = (median), col = "#D37295")) +
      geom_line(aes(x = date, y = (CASES), col = "#499894")) +
      xlab("Time") +
      ylab("Incidence per 100000") +
      scale_colour_identity(name = "",
                            breaks = c("#499894", "#D37295"),
                            labels = c("Observed", "Posterior predicted"),
                            guide = "legend") +
      theme_bw() + 
      theme(panel.grid = element_blank(),
            axis.text = element_text(size = 12), 
            legend.position = "bottom",
            legend.text = element_text(size = 15),
            strip.text = element_text(size = 15),
            
            axis.text.x = element_text(angle=25),
            axis.title = element_text(size=15),
            strip.background = element_blank())+
      guides(col = guide_legend(override.aes = list(linewidth = 10))) 
    
  
    }
  return(plot)
  
}


## plotting model error
plot_map_MAE <- function(data, model, pathogen, shapefile){
  
  
  data$estimate<-exp(model$summary.linear.predictor$'0.5quant')
  data$estimate_low<-exp(model$summary.linear.predictor$'0.025quant')
  data$estimate_upp<-exp(model$summary.linear.predictor$'0.975quant')
    
  data |> select(MUNICIP, STATE, CASES, YEAR, WEEK_NO, pop, estimate)|>
      mutate(CASES=(CASES/pop)*100000,
             estimate= (estimate/pop)*100000)|>
      mutate(error = abs((CASES - estimate))) |> 
      group_by(MUNICIP) |> 
      summarise(error_t = sum(error)/n()) -> d
    
    ## error map
    colnames(shapefile)[1]<-"MUNICIP"
    
   MAP <- dplyr::left_join(shapefile, d, by = "MUNICIP")
   MAP <-MAP[is.na(MAP$error_t)==F,]
  
  pal<-c("lightblue", "blue", "red")
  sf_gadm<-read_sf("data/gadm41_BRA_shp/gadm41_BRA_0.shp")
  
  
  A<-ggplot(MAP) + 
    geom_sf(data=sf_gadm$geometry,col="black",fill=NA)+
    geom_sf(aes(fill = error_t),color=NA) +
    scale_fill_gradientn(colours=pal, limits=c(0,200))+
    labs(fill = "Mean absolute\nerror", title=toupper(pathogen)) +
    theme_void()+
    theme(plot.title = element_text(hjust=0.5, size=15)
    )
 
  return(A)
}
