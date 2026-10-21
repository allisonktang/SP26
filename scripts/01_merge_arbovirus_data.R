##################################################################################################
# Brazil Arbovirus (DENV, ZIKA, CHIKV) Transmission (2016-2025)
# Author: 
# Date: 4/27/2026
# Datasets: 
#   - DENV:   (collected october 13, 2026)
#   - ZIKA:   (collected october 13, 2026)
#   - CHIKV:  (collected october 13, 2026)
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

# DENGUE from https://dadosabertos.saude.gov.br
# Municipality is determined primarily using COMUNINF. Source indicators used when COMUNINF != residence, or COMUNINF is blank
#   - if COMUNINF exists and == residence,  muni = COMUNINF and source_indicator = autochthnonous_residence
#   - if COMUNINF exists and != residence,  muni = COMUNINF and source_indicator = infection_recorded
#   - if COMUNINF is blank,                 muni = ID_MN_RESI and source_indicator = residence_fallback 
# Columns: TP_NOT,DT_NOTIFIC, SEM_NOT, NU_ANO, MUNI, CLASSI_FIN, CRITERIO, CS_SEXO, NU_IDADE_N, source_indicator, UF_CODE
dengue <- read_csv("data/ncombined_dengue_2016_2025.csv") %>%
  filter(between(NU_ANO, 2016, 2025)) %>%
  filter(CLASSI_FIN %in% c(1, 10, 11, 12)) %>%
  mutate(
    muni= as.character(MUNI),
    date = ymd(DT_NOTIFIC),
    year = year(date),
    month = month(date)
  ) %>%
  filter(!is.na(muni), !is.na(year), !is.na(month))
glimpse(dengue) # check
colnames(dengue)

# ZIKA from https://dadosabertos.saude.gov.br
# Municipality is determined primarily using COMUNINF. Source indicators used when COMUNINF != residence, or COMUNINF is blank
#   - if COMUNINF exists and == residence,  muni = COMUNINF and source_indicator = autochthnonous_residence
#   - if COMUNINF exists and != residence,  muni = COMUNINF and source_indicator = infection_recorded
#   - if COMUNINF is blank,                 muni = ID_MN_RESI and source_indicator = residence_fallback 
# Columns: TP_NOT,DT_NOTIFIC, SEM_NOT, NU_ANO, MUNI, CLASSI_FIN, CRITERIO, CS_SEXO, NU_IDADE_N, source_indicator, UF_CODE

zika <- read_csv("data/ncombined_zika_2016_2025.csv") %>%
  filter(between(NU_ANO, 2016, 2025)) %>%
  filter(CLASSI_FIN %in% c(1)) %>%
  mutate(
    muni= as.character(MUNI),
    date = ymd(DT_NOTIFIC),
    year = year(date),
    month = month(date)
  ) %>%
  filter(!is.na(muni), !is.na(year), !is.na(month))
glimpse(zika) # check
colnames(zika)

# CHIKV from https://dadosabertos.saude.gov.br
# Municipality is determined primarily using COMUNINF. Source indicators used when COMUNINF != residence, or COMUNINF is blank
#   - if COMUNINF exists and == residence,  muni = COMUNINF and source_indicator = autochthnonous_residence
#   - if COMUNINF exists and != residence,  muni = COMUNINF and source_indicator = infection_recorded
#   - if COMUNINF is blank,                 muni = ID_MN_RESI and source_indicator = residence_fallback 
# Columns: TP_NOT,DT_NOTIFIC, SEM_NOT, NU_ANO, MUNI, CLASSI_FIN, CRITERIO, CS_SEXO, NU_IDADE_N, source_indicator, UF_CODE

chikv <- read_csv("data/ncombined_chikv_2016_2025.csv") %>%
  filter(between(NU_ANO, 2016, 2025)) %>%
  filter(CLASSI_FIN %in% c(1, 13)) %>%
  mutate(
    muni= as.character(MUNI),
    date = ymd(DT_NOTIFIC),
    year = year(date),
    month = month(date)
  ) %>%
  filter(!is.na(muni), !is.na(year), !is.na(month))
glimpse(chikv) # check
colnames(chikv)

sum(is.na(dengue$date))
sum(is.na(zika$date))
sum(is.na(chikv$date))

#check date range
range(dengue$date)
range(zika$date)
range(chikv$date)


# population
pop <- read_csv("data/Pop_size_estimative_2010_2025_IBGE_collected_in_April_14_2026_clean.csv") %>% 
  filter(between(Year, 2016, 2025)) %>%
  filter(code_muni != 110000) %>% # no population data for state department
  mutate(
    Pop_estimated = trimws(Pop_estimated),
    Pop_estimated = sub("\\s*\\([0-9]+\\)$", "", Pop_estimated),
    Pop_estimated = gsub("\\.", "", Pop_estimated)
  ) %>% # remove footnote-- sometimes formatted with a space, ranges from 1-19 
  transmute(
    muni = substr(as.numeric(code_muni), 1, 6),
    year = Year,
    pop_tot = as.numeric(Pop_estimated)
    )
glimpse(pop)

# check for any duplicate populations
pop %>%
  count(muni, year) %>%
  filter(n > 1)

# check for missing values
pop %>%
  filter(is.na(pop_tot) | pop_tot <= 0)
  
##################################################################################################
# 2. Aggregate confirmed cases by muni, year, month 
##################################################################################################

dengueClean <- dengue %>%
  group_by(muni, year, month) %>%
  summarise(
    dengueCases = n(),
    .groups = "drop"
  )

zikaClean <- zika %>% 
  group_by(muni, year, month) %>%
  summarise(
    zikaCases = n(), 
    .groups = "drop"
  )

chikvClean <- chikv %>%
  group_by(muni, year, month) %>%
  summarise(
    chikvCases = n(), 
    .groups = "drop"
    )

# check for duplicate muni-month cases
dengueClean %>%
  count(muni, year, month) %>%
  filter(n > 1)

zikaClean %>%
  count(muni, year, month) %>%
  filter(n > 1)

chikvClean %>%
  count(muni, year, month) %>%
  filter(n > 1)

# check range
range(dengueClean$dengueCases)
range(zikaClean$zikaCases)
range(chikvClean$chikvCases)

##################################################################################################
# 3. Merge datasets and calculate incidence
##################################################################################################
# expand population dataset to use annual data for each month
popClean <- pop %>%
  distinct(muni, year, pop_tot) %>%                        # one row per muni-year
  expand_grid(month = 1:12) %>%  # create month 1–12 for each muni-year
  mutate(month = as.integer(month)) %>%
  arrange(muni, year, month)

# merge disease datasets to population 
merged <- popClean %>%
  left_join(dengueClean, by = c("muni", "year", "month")) %>%
  left_join(zikaClean,   by = c("muni", "year", "month")) %>%
  left_join(chikvClean,  by = c("muni", "year", "month")) 
glimpse(merged)

merged %>%
  filter(is.na(pop_tot) | pop_tot <= 0)

# calculate incidence per 100k
merged <- merged %>%
  mutate(
    dengueIncidence = dengueCases / pop_tot * 100000,
    zikaIncidence   = zikaCases / pop_tot * 100000,
    chikvIncidence  = chikvCases / pop_tot * 100000
  ) %>%
  select(
    muni, year, month, pop_tot, 
    dengueCases, zikaCases, chikvCases, 
    dengueIncidence, zikaIncidence, chikvIncidence
    )

# check: range excluding NA
range(merged$dengueCases, na.rm = TRUE) 
range(merged$zikaCases, na.rm = TRUE)
range(merged$chikvCases, na.rm = TRUE)

# check: how many entries are NA?
mean(is.na(merged$dengueCases)) 
mean(is.na(merged$zikaCases))
mean(is.na(merged$chikvCases))

# check: duplicates?
merged %>%
  count(muni, year, month) %>%
  filter(n > 1)

merged <- merged %>%
  mutate(year_month = as.Date(sprintf("%d-%02d-01", year, month)))

# sanity check: verify incidence for one muni-year (Sao Paulo City, 2023)
merged %>%
  filter(muni == "355030", year == 2023) %>%
  summarise(
    pop = first(pop_tot),
    
    dengue_cases_sum  = sum(dengueCases, na.rm = T),
    zika_cases_sum    = sum(zikaCases,   na.rm = T),
    chikv_cases_sum   = sum(chikvCases,  na.rm = T),
    
    dengue_manual = dengue_cases_sum / pop * 100000,
    zika_manual   = zika_cases_sum   / pop * 100000,
    chikv_manual  = chikv_cases_sum  / pop * 100000,
    
    dengue_from_col = sum(dengueIncidence, na.rm = T),
    zika_from_col   = sum(zikaIncidence,   na.rm = T),
    chikv_from_col  = sum(chikvIncidence,  na.rm = T)
  ) %>%
  glimpse()


write_csv(merged, "data/nBrazil_arbovirus_monthly_data_2016_2025.csv")


# count cases
summary_counts <- merged %>%
  summarise(
    # 1. num municipalities
    total_municipalities = n_distinct(muni),
    
    # 2. num muni-month observations (total rows in the dataset)
    total_muni_month_obs = n(),
    
    # 3. total case counts (ignores NA)
    total_dengue_cases = sum(dengueCases, na.rm = TRUE),
    total_zika_cases   = sum(zikaCases, na.rm = TRUE),
    total_chikv_cases  = sum(chikvCases, na.rm = TRUE)
  )

print(summary_counts)


# plotting
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



# HEATMAPS (non log)
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


# incidence heatmap
# goal: municipality level incidence heatmap (ungrouped), raw + log10 scale
# input: "Brazil_arbovirus_monthly_data_2016_2025.csv"
muni_heatmap <- ggplot(mergedLong_muni,
       aes(x = year_month, y = factor(muni), fill = incidence)) +
  geom_tile(color = NA, linewidth = 0) +
  facet_wrap(~ disease, ncol = 1) +
  scale_fill_gradient(
    low="#FCF0CE",
    high = "#D4180A",
    na.value = "white", 
    trans = "log10") + 
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  labs(x = "Time", y = "Municipality", fill = "Incidence") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

print(muni_heatmap)

