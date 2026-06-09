##################################################################################################
# Brazil Arbovirus (DENV, ZIKA, CHIKV) Transmission (2016-2025)
# Author: 
# Date: 4/27/2026
# Datasets: 
#   - DENV:
#   - ZIKA: 
#   - CHIKV:
#   - pop: Pop_size_estimative_2010_2025_IBGE_collected_in_April_14_2026_clean.csv
# Goal: Merge datasets and calculate incidence (yearly, monthly) for Brazilian municipalities
# Notes:
#   - Population data is yearly -> assume the same values for monthly population
#   - Incidence per 100k: Number of Cases / Municipality Population * 100,000
##################################################################################################

rm(list = ls())

library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

##################################################################################################
# 1. Load datasets 
##################################################################################################

# DENV from ()
# Columns: muni, date, DENV_per_100k 
dengue <- read_csv("data/DENV-serotypes_1996-2025_monthly_mun.csv") %>%
  rename(muni = CD_MUN) %>%
  filter(!is.na(muni)) %>%
  mutate(
    muni = substr(muni, 1, 6),
    date = ymd(date),
    year = year(date),
    month = month(date),
    dengueCases = DENV_total) %>%
  filter(between(year, 2016, 2025))     # removes rows without data
glimpse(dengue)
colnames(dengue)


# ZIKA from https://dadosabertos.saude.gov.br
# Municipality is determined primarily using COMUNINF. Source indicators used when COMUNINF != residence, or COMUNINF is blank
#   - if COMUNINF exists and == residence,  muni = COMUNINF and source_indicator = autochthnonous_residence
#   - if COMUNINF exists and != residence,  muni = COMUNINF and source_indicator = infection_recorded
#   - if COMUNINF is blank,                 muni = ID_MN_RESI and source_indicator = residence_fallback 
# Columns: TP_NOT,DT_NOTIFIC, SEM_NOT, NU_ANO, MUNI, CLASSI_FIN, CRITERIO, CS_SEXO, NU_IDADE_N, source_indicator, UF_CODE

zika <- read_csv("data/combined_zika_2016_2025.csv") %>%
  filter(between(NU_ANO, 2016, 2025)) %>%
  rename(muni = MUNI) %>%
  filter(!is.na(muni)) %>%
  mutate(
    date = ymd(DT_NOTIFIC),
    year = year(date),
    month = month(date)
  ) %>%
  filter(!is.na(year) & !is.na(month))
glimpse(zika) # check
colnames(zika)

# CHIKV from winson smile emoji (???)
chikv <- read_csv("data/bd_CHIKV_positive_2013_2025_english.csv") %>% 
  filter (between(Collection_Year, 2016, 2025)) %>%
  rename(muni = Municipality_Code) %>%
  filter(!is.na(muni)) %>%
  mutate(
    date = lubridate::mdy(Collection_Date),
    year = year(date),
    month = month(date)
  )
glimpse(chikv)
colnames(chikv)


# population
pop <- read_csv("data/Pop_size_estimative_2010_2025_IBGE_collected_in_April_14_2026_clean.csv") %>% 
  filter(Year >= 2016 & Year <= 2025) %>%
  filter(code_muni != 110000) %>% # no population data for state department
  select(code_muni, Year, Pop_estimated, Sigla) %>%
    mutate(
      muni = substr(code_muni, 1, 6),
      year = Year,
      pop_tot = as.numeric(Pop_estimated)
    )
glimpse(pop)
  
##################################################################################################
# 2. Clean datasets to include 6 digit muni, year, month, number of cases/incidence
##################################################################################################

dengueClean <- dengue %>%
  select(muni, year, month, dengueCases)

zikaClean <- zika %>% 
  group_by(muni, year, month) %>%
  filter(!is.na(muni)) %>%
  summarise(zikaCases = n(), .groups = "drop")

chikvClean <- chikv %>%
  group_by(muni, year, month) %>%
  filter(!is.na(muni)) %>%
  summarise(chikvCases = n(), .groups = "drop")


# clean population dataset to use annual data for each month
popClean <- pop %>%
  select(muni, year, pop_tot) %>%
  distinct() %>%                        # one row per muni-year
  expand_grid(month = 1:12) %>%  # create month 1–12 for each muni-year
  arrange(muni, year, month)

# ensure all columns have same type 
dengueClean <- dengueClean %>% 
  mutate(muni = as.character(muni),
        year = as.integer(year),
        month = as.integer(month))

zikaClean <- zikaClean %>%
  mutate(muni = as.character(muni),
         year = as.integer(year),
         month = as.integer(month))

chikvClean <- chikvClean %>%
  mutate(muni = as.character(muni),
         year = as.integer(year),
         month = as.integer(month))

popClean <- popClean %>%
  mutate(muni = as.character(muni),
         year = as.integer(year),
         month = as.integer(month))

##################################################################################################
# 3. Merge datasets and calculate incidence
##################################################################################################

# merge disease datasets to population 
merged <- popClean %>%
  left_join(dengueClean, by = c("muni", "year", "month")) %>%
  left_join(zikaClean,   by = c("muni", "year", "month")) %>%
  left_join(chikvClean,  by = c("muni", "year", "month")) 
glimpse(merged)

merged <- merged %>%
  mutate(
    dengueIncidence = dengueCases / pop_tot * 100000,
    zikaIncidence   = zikaCases / pop_tot * 100000,
    chikvIncidence  = chikvCases / pop_tot * 100000
  ) %>%
  select(
    muni, year, month, pop_tot, dengueCases, zikaCases, chikvCases, dengueIncidence, zikaIncidence, chikvIncidence)

# check: range excluding NA
range(merged$dengueCases, na.rm = TRUE) 
range(merged$zikaCases, na.rm = TRUE)
range(merged$chikvCases, na.rm = TRUE)

# check: how many entries are NA?
mean(is.na(merged$dengueCases)) 
mean(is.na(merged$zikaCases))
mean(is.na(merged$chikvCases))

# to replace NA with 0:
#merged <- merged %>%
#  mutate(
#    dengueCases = replace_na(dengueCases, 0),
#    zikaCases   = replace_na(zikaCases, 0),
#    chikvCases  = replace_na(chikvCases, 0),

    # recalc incidence
#    dengueIncidence = dengueCases / pop_tot * 100000,
#    zikaIncidence   = zikaCases / pop_tot * 100000,
#    chikvIncidence  = chikvCases / pop_tot * 100000,
#  )

merged <- merged |>
  mutate(year_month = as.Date(sprintf("%d-%02d-01", year, month)))

glimpse(merged)
write_csv(merged, "Brazil_arbovirus_monthly_data_2016_2025.csv")


# convert the dataset into long format
mergedLong <- merged %>%
  select(year_month, 
         dengueIncidence,
         zikaIncidence,
         chikvIncidence) %>%
  pivot_longer(
    cols = -year_month,
    names_to = "disease",
    values_to = "incidence"
  ) %>%
  mutate(
    disease = recode(disease,
      dengueIncidence = "Dengue",
      zikaIncidence = "Zika",
      chikvIncidence = "Chikungunya")
  )


mergedLong <- mergedLong %>%
  mutate(disease = factor(disease, levels = c(
    "Dengue", "Zika", "Chikungunya"
  ))) 

mergedLong_muni <- merged %>%
  pivot_longer(
    cols = c(dengueIncidence, zikaIncidence, chikvIncidence),
    names_to = "disease",
    values_to = "incidence"
  ) %>%
  mutate(
    disease = recode(disease,
                     dengueIncidence = "Dengue",
                     zikaIncidence = "Zika",
                     chikvIncidence = "Chikungunya"    )
  )

# aggregate (mean incidence across all municipalities )
munMean <- mergedLong %>%
  group_by(year_month, disease) %>%
  summarise(incidence = mean(incidence, na.rm = T), .groups = "drop")

# national incidence
natIncidence <- merged %>%
  group_by(year_month) %>%
  summarise(
    dengueIncidence = sum(dengueCases, na.rm = TRUE) / sum(pop_tot) * 100000,
    zikaIncidence   = sum(zikaCases, na.rm = TRUE) / sum(pop_tot) * 100000,
    chikvIncidence  = sum(chikvCases, na.rm = TRUE) / sum(pop_tot) * 100000,
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -year_month,
    names_to = "disease",
    values_to = "incidence"
  ) %>%
  mutate(
    disease = recode(disease,
                     dengueIncidence = "Dengue",
                     zikaIncidence = "Zika",
                     chikvIncidence = "Chikungunya"
    )
  )

cols = c("Dengue" = "#F8766D",
             "Zika" = "#7CAE00",
             "Chikungunya" = "#00BFC4")

# PLOT MUNICIPALITY AVG INCIDENCE

ggplot(munMean, aes(x = year_month, y = incidence, color = disease, fill = disease)) + 
  geom_ribbon(aes(ymin = 0, ymax = incidence), alpha = 0.2, color = NA) +
  geom_line(alpha = 0.7) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  guides(fill = "none") +
  labs(
    x = "Year-Month",
    y = "Average Municipality Incidence per 100k",
    color = "Disease"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

natIncidence_log <- natIncidence %>%
  filter(incidence > 0)

# PLOTS FOR NATIONAL INCIDENCE (NOT WORKING)
ggplot(natIncidence_log, aes(x = year_month, y = incidence, color = disease, fill = disease)) + 
  geom_ribbon(aes(ymin = 0, ymax = incidence), alpha = 0.2, color = NA) +
  geom_line(alpha = 0.7) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  guides(fill = "none") +
  labs(
    x = "Year-Month",
    y = "National Incidence per 100k",
    color = "Disease"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# PLOT FOR NATIONAL INCIDENCE (log scale) (NOT WORKING)
ggplot(natIncidence_log, aes(x = year_month, y = incidence, color = disease, fill = disease)) + 
  geom_ribbon(aes(ymin = 0, ymax = incidence), alpha = 0.2, color = NA) +
  geom_line(alpha = 0.7) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  scale_y_log10(oob = scales::squish_infinite) +
  guides(fill = "none") +
  labs(
    x = "Year-Month",
    y = "National Incidence per 100k (Log scale)",
    color = "Disease"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

# HEATMAPS
ggplot(mergedLong_muni,
       aes(x = year_month, y = factor(muni), fill = incidence)) +
  geom_tile(color = NA, linewidth = 0) +
  facet_wrap(~ disease, ncol = 1) +
  scale_fill_gradient(
    low="#FCF0CE",
    high = "#D4180A",
    na.value = "white") + 
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  labs(x = "Time", y = "Municipality", fill = "Incidence") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggplot(mergedLong_muni,
       aes(x = year_month, y = factor(muni), fill = incidence)) +
  geom_tile(color = NA, linewidth = 0) +
  facet_wrap(~ disease, ncol = 1) +
  scale_fill_viridis_c(option = "G",
    na.value = "white",
    trans = "log10") + 
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  labs(x = "Time", y = "Municipality", fill = "Incidence") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

