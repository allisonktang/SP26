# calculate dengue incidence by municipality, year

# output columns:
# MUNI  year  cases_tot cases_m cases_f pop_tot pop_m pop_f incidence_tot incidence_m incidence_f
# Incidence calculated per 100k

rm(list = ls())

library(dplyr)
library(tidyr)
library(readr)

# 1. load datasets
cases_raw <- read_csv("data/dengue/deng_data_2016_2025.csv")
cases_raw <- cases_raw %>%
  mutate(NU_ANO = as.numeric(NU_ANO)) %>%   # convert to numeric
  filter(NU_ANO >= 2016 & NU_ANO <= 2025)
pop <- read_csv("data/Pop_size_estimative_2010_2025_IBGE_collected_in_April_14_2026_clean.csv") 
pop <- pop %>%
  filter(!Municipio_age %in% "110000") # no population data for state department

# removing these
cases_raw <- cases_raw %>%
  filter(!MUNI %in% c("110000", "35", "47", "35050", "35503", NA)) # possible typos, aggregate data? but no population data for this

# 2. aggregate total cases per muni-year
total_cases <- cases_raw %>%
  filter(CRITERIO %in% c("1", "2","3")) %>%                           # keep confirmed case
  group_by(MUNI, NU_ANO) %>%                            # group by muni, year
  summarise(total_cases = n(), .groups = "drop") %>%    # count cases in each group
  rename(
    muni = MUNI,
    year = NU_ANO
  )


# 3. aggregate cases by sex per muni-year
cases_sex <- cases_raw %>%
  filter(CRITERIO %in% c("1", "2","3")) %>%                           # keep confirmed case
  group_by(MUNI, NU_ANO, CS_SEXO) %>%
  summarise(cases_sex = n(), .groups = "drop")  %>%
  pivot_wider(
    names_from = CS_SEXO,
    values_from = cases_sex,
    names_prefix = "cases_"
  ) %>%
  mutate(
    cases_M = replace_na(cases_M, 0),
    cases_F = replace_na(cases_F, 0),
    cases_I = replace_na(cases_I, 0)
  ) %>%
  rename(
    muni = MUNI,
    year = NU_ANO
  )


# 4. prep population data
pop_total <- pop %>%
  select(municipality_code_age, year, Total, Masculino, Feminino) %>%
  rename(
    muni = municipality_code_age,
    pop_tot = Total,
    pop_m = Masculino,
    pop_f = Feminino
  )

# 5. join cases with population
incidence_data <- total_cases %>%
  left_join(cases_sex, by = c("muni", "year")) %>%
  left_join(pop_total, by = c("muni", "year"))

# 6. calculate incidence per 100k
incidence_data <- incidence_data %>%
  rename(
    cases_tot = total_cases,
    cases_m = cases_M,
    cases_f = cases_F,
    cases_i = cases_I
  ) %>%
  mutate(
    incidence_tot = round(cases_tot / pop_tot * 100000, 3),
    incidence_m   = round(cases_m / pop_m * 100000, 3),
    incidence_f   = round(cases_f / pop_f * 100000, 3)
  )

# 7. save output
write_csv(incidence_data, "/Users/allison/Desktop/Research_BentoLab/SP26/data/dengue/deng_incidence_2016_2025.csv")
