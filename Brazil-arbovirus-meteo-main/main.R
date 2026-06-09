##############################################
###### Spatiotemporal relationships ##########
###### between extreme meteo and #############
###### arbovirus infections in Brazil ########
##############################################

# The code here recreates the analysis in the paper. 
# Data must first be downloaded from zenodo and saved in the data folder in this repository.

install.packages("INLA")
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
## load packages
library(INLA)
library(data.table)
library(ggplot2)
library(geofacet)
library(dplyr)
library(forcats)
library(cowplot)
library(sf)
library(ISOweek)
library(ggpubr)
library(ggrepel)
library(beepr)
library(sjPlot)
library(stringr)

## source functions
source("R/plotting_functions.R")
source("R/extra_functions.R")
source("R/univariable_analysis_functions.R")
source("R/multivariable_analysis_functions.R")

## read in shapefile
shapefile <- readRDS("data/Brazil_2_geobr.RDS")

## path to read in graph for INLA models
GRAPH <- "data/Brazil2.graph"


####################################################################################
###################### Disease data visualisation ##################################
####################################################################################

## temporal heterogeneity - average cases per week number per state
temporal_cases_state(data_type = "exp_resid")

## spatial heterogeneity - average cases per municipality
spatial_cases_municip(data_type = "exp_resid")

## spatiotemporal heterogeneity - heatmaps per pathogen
heatmaps(shapefile = shapefile, pathogen = "zikv", data_type = "exp_resid")


####################################################################################
########################## Run baseline models  ####################################
####################################################################################

## choose the pathogen (zikv, chikv, denv)
pathogen <- "zikv"
## choose data type (exp_resid, exp_report, resid, report which specify how the municipality is assigned - see README)
data_type <- "exp_resid"

## read in data
data_tot <- readRDS(paste0("data/", pathogen, "_", data_type, ".RDS"))
data_tot <- data_tot |> 
            data_processing()  

## run model 
    # note that the models take a long time to run
    # beep will sound when this step is finished
baseline <- run_baseline(data = data_tot, GRAPH = GRAPH)
beep()

## save model
summary(baseline)
dir.create("outputs")
saveRDS(baseline, paste0("outputs/baseline_model_", pathogen, ".RDS"))


####################################################################################
########################## Univariable analysis ####################################
####################################################################################

## select variables
all_variables <- select_variables(data_tot, "all")
temperature_variables <- select_variables(data_tot, "temperature")
testing_variables <- select_variables(data_tot, "subset")

## run univariable models
  # the models take a long time to run without HPC 
  # for testing - run for Zika because it is a smaller dataset, and use only a subset of variables 
dir.create("outputs/univar")
dir.create("outputs/univar/temp")
univar_results_table <- run_univar_analysis(data=data_tot, 
                                            variables=testing_variables, # suggest running in parallel for sets of variables
                                            save_progress=TRUE, # save in outputs/temp/
                                            pathogen=pathogen)
beep()
univar_results_table

## bind all saved univariable results files together 
univar_results_table_full <- binding_univar_results(pathogen, file_path="outputs/univar/temp/")
saveRDS(univar_results_table_full, paste0("outputs/univar/", pathogen, "_univar_all.RDS"))


## best univariable model
best_univariable_result <- best_univar(univar_results_table_full) # returns variable and WAIC score
best_univariable_result

## plot univariable model betas 
plot_univar_betas(choice = "temperature")
plot_univar_betas(choice = "humidity")
plot_univar_betas(choice = "rain")
plot_univar_betas(choice = "socioeconomic")


####################################################################################
########################## Multivariable analysis ##################################
####################################################################################

##### run multivariable analysis - forward selection to choose variables ######
dir.create("outputs/multivar")
dir.create("outputs/multivar/temp")

      #### stage 1
       
      # remove variables which are collinear with best variable from univariable analysis
      to_remove <- collinarity_check(data_tot, all_variables, variable_chosen = best_univariable_result[[1]])
      step_variables <- all_variables[!all_variables %in% to_remove]
      
      # run step 1 of forward selection process
        # note - these steps will be slow 
        # suggest splitting step_variables into parallel runs
      step1 <- run_multivar_analysis(data = data_tot, GRAPH = GRAPH, variables_to_try = step_variables, 
                                   variables_chosen_already = best_univariable_result[[1]], stage = 1)
      beep()
      
      # check waic decrease threshold reached
      (round(as.numeric(best_univariable_result[[2]])) - round(as.numeric(step1[[2]]))) > 4
      
      #### stage 2
      
      # remove variables collinear with best variable from step 1
      to_remove <- collinarity_check(data_tot, all_variables, variable_chosen = step1[[1]])
      step_variables <- step_variables[!step_variables %in% to_remove]
      
      # run step 2 of forward selection process
      step2 <- run_multivar_analysis(data=data_tot,GRAPH=GRAPH, variables_to_try = step_variables, 
                                     variables_chosen_already = c(best_univariable_result[[1]], step1[[1]]), 
                                     stage = 2)
      beep()
      
      # check waic decrease threshold reached
      (round(as.numeric(step1[[2]])) - round(as.numeric(step2[[2]]))) > 4
      
      #### stage 3
      
      # remove variables collinear with best variable from step 2
      to_remove <- collinarity_check(data_tot, all_variables, variable_chosen = step2[[1]])
      step_variables <- step_variables[!step_variables %in% to_remove]
      
      # run step 3 of forward selection process
      step3 <- run_multivar_analysis(data=data_tot,GRAPH=GRAPH, variables_to_try = step_variables, 
                                     variables_chosen_already = c(best_univariable_result[[1]], step1[[1]], step2[[1]]), 
                                     stage = 3)
      beep()
      
      # check waic decrease threshold reached
      (round(as.numeric(step2[[2]])) - round(as.numeric(step3[[2]]))) > 4
      
      #### stage 4
      
      # remove variables collinear with best variable from step 3
      to_remove <- collinarity_check(data_tot, all_variables, variable_chosen = step3[[1]])
      step_variables <- step_variables[!step_variables %in% to_remove]
      
      # run step 4 of forward selection process
      step4 <- run_multivar_analysis(data=data_tot,GRAPH=GRAPH, variables_to_try = step_variables, 
                                     variables_chosen_already = c(best_univariable_result[[1]], step1[[1]], step2[[1]], step3[[1]]), 
                                     stage = 4)
      beep()
      
      # check waic decrease threshold reached
      (round(as.numeric(step3[[2]])) - round(as.numeric(step4[[2]]))) > 4
      
      #### stage 5
      
      # remove variables collinear with best variable from step 4
      to_remove <- collinarity_check(data_tot, all_variables, variable_chosen = step4[[1]])
      step_variables <- step_variables[!step_variables %in% to_remove]
      
      # run step 5 of forward selection process
      step5 <- run_multivar_analysis(data=data_tot,GRAPH=GRAPH, variables_to_try = step_variables, 
                                     variables_chosen_already = c(best_univariable_result[[1]], step1[[1]], step2[[1]], step3[[1]], step4[[1]]), 
                                     stage = 5)
      beep()
      
      # check waic decrease threshold reached
      (round(as.numeric(step4[[2]])) - round(as.numeric(step5[[2]]))) > 4
      
      #### stage 6
      
      # remove variables collinear with best variable from step 5
      to_remove <- collinarity_check(data_tot, all_variables, variable_chosen = step5[[1]])
      step_variables <- step_variables[!step_variables %in% to_remove]
      
      # run step 6 of forward selection process
      step6 <- run_multivar_analysis(data=data_tot,GRAPH=GRAPH, variables_to_try = step_variables, 
                                     variables_chosen_already = c(best_univariable_result[[1]], step1[[1]], 
                                                                  step2[[1]], step3[[1]], step4[[1]], step5[[1]]), 
                                     stage = 6)
      beep()
      
      # check waic decrease threshold reached
      (round(as.numeric(step5[[2]])) - round(as.numeric(step6[[2]]))) > 4
      
      #### stage 7
      
      # remove variables collinear with best variable from step 6
      to_remove <- collinarity_check(data_tot, all_variables, variable_chosen = step6[[1]])
      step_variables <- step_variables[!step_variables %in% to_remove]
      
      # run step 7 of forward selection process
      step7 <- run_multivar_analysis(data=data_tot,GRAPH=GRAPH, variables_to_try = step_variables, 
                                     variables_chosen_already = c(best_univariable_result[[1]], step1[[1]], 
                                                                  step2[[1]], step3[[1]], step4[[1]], step5[[1]], step6[[1]]), 
                                     stage = 7)
      beep()
      
      # check waic decrease threshold reached
      (round(as.numeric(step6[[2]])) - round(as.numeric(step7[[2]]))) > 4
      
      #### stage 8
      
      # remove variables collinear with best variable from step 7
      to_remove <- collinarity_check(data_tot, all_variables, variable_chosen = step7[[1]])
      step_variables <- step_variables[!step_variables %in% to_remove]
      
      # run step 8 of forward selection process
      step8 <- run_multivar_analysis(data=data_tot,GRAPH=GRAPH, variables_to_try = step_variables, 
                                     variables_chosen_already = c(best_univariable_result[[1]], step1[[1]], 
                                                                  step2[[1]], step3[[1]], step4[[1]], step5[[1]], 
                                                                  step6[[1]], step7[[1]]), 
                                     stage = 8)
      beep()
      
      # check waic decrease threshold reached
      (round(as.numeric(step7[[2]])) - round(as.numeric(step8[[2]]))) > 4
      
      #### stage 9
      
      # remove variables collinear with best variable from step 8
      to_remove <- collinarity_check(data_tot, all_variables, variable_chosen = step8[[1]])
      step_variables <- step_variables[!step_variables %in% to_remove]
      
      # run step 9 of forward selection process
      step9 <- run_multivar_analysis(data=data_tot,GRAPH=GRAPH, variables_to_try = step_variables, 
                                     variables_chosen_already = c(best_univariable_result[[1]], step1[[1]], 
                                                                  step2[[1]], step3[[1]], step4[[1]], step5[[1]], 
                                                                  step6[[1]], step7[[1]], step8[[1]]), 
                                     stage = 9)
      beep()
      
      # check waic decrease threshold reached
      (round(as.numeric(step8[[2]])) - round(as.numeric(step9[[2]]))) > 4
      


############# run multivariable analysis - final models ##############

# run final multivariable model
final_model<- run_final_model(data=data_tot, graph=GRAPH, pathogen=pathogen) # this function includes the variables chosen in the models from the paper
beep()
summary(final_model)

saveRDS(final_model, paste0("outputs/", pathogen, "_final_model.RDS"))


### random effect plots
plotting_RE_week_number(model = final_model, data = data_tot, pathogen = pathogen)

plotting_RE_municip(model = final_model, data = data_tot, shapefile = shapefile) 

    
### sampling
n <- 1000
index <- 1:nrow(data_tot)
samples <- sampling_models(model = final_model, index = index, nsamp = n)
beep()
saveRDS(samples, paste0("outputs/samples_", pathogen, ".RDS"))


### fit plots
plot_model_fit(data=data_tot, pathogen=pathogen, nsamp=n, samples=samples, level="national")
plot_model_fit(data=data_tot, pathogen=pathogen, nsamp=n, samples=samples, level="state")


### error maps
plot_map_MAE(data=data_tot, model=final_model, pathogen=pathogen, shapefile=shapefile)



####################################################################################
########################## Sensitivity analysis ####################################
####################################################################################

data_type <-"exp_resid" 
      # response must be report, exp_report, resid, or exp_resid - main analysis used exp_resid
      # report = municipality where a case was reported
      # resid = municipality where the person reporting a case lives
      # exp_report = municipality where a case was suspected infected (where this information is available), and otherwise report
      # exp_resid = municipality where a case was suspected infected (where this information is available), and otherwise resid 

data_sens <- readRDS(paste0("data/", pathogen,"_", data_type,".RDS"))
data_sens <- data_processing(data_sens)

sensitivity <- run_final_model(data=data_sens, graph=GRAPH, pathogen=pathogen)
beep()
summary(sensitivity)
