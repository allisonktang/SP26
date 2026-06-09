# calculate zika incidence by municipality, year

# output columns:
# MUNI  year cases_tot pop_tot incidence_tot 
# Incidence calculated per 100k

rm(list = ls())

library(dplyr)
library(tidyr)
library(readr)

# 1. load datasets
cases_raw <- read_csv("combined_zika_2016_2025.csv")
cases_raw <- cases_raw %>%
  filter(NU_ANO >= 2016 & NU_ANO <= 2025)
pop <- read_csv("data/Pop_size_estimative_2010_2025_IBGE_collected_in_April_14_2026_clean.csv")

# 2. aggregate total cases per muni-date
total_cases <- cases_raw %>%
  filter(CLASSI_FIN == 1) %>%                           # keep confirmed case
  transmute(
    muni = sprintf("%07d", as.numeric(MUNI)),
    year = NU_ANO
  ) %>%
  group_by(muni, year) %>%    
  summarise(total_cases = n(), .groups = "drop")
  
# 3. prep population data
pop_total <- pop %>%
  select(code_muni, Year, Pop_estimated, Sigla) %>%
  transmute(
    muni = sprintf("%07d", as.numeric(code_muni)),
    year = Year,
    pop_tot = as.numeric(Pop_estimated)
  )

# 4. join cases with population
incidence_data <- total_cases %>%
  left_join(pop_total, by = c("muni", "year")
  )

# 6. calculate incidence per 100k
incidence_data <- incidence_data %>%
  rename(
    cases_tot = total_cases
  ) %>%
  mutate(
    incidence = round(cases_tot / pop_tot * 100000, 3)
  )

# 7. save output
write_csv(incidence_data, "/Users/allison/Desktop/Research_BentoLab/SP26/zika/zika_incidence_2016_2025.csv")

