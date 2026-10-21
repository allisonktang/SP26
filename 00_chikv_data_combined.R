# 00: download and combine chikungunya virus files
# data collected from datasus: october 13, 2026

rm(list = ls())

library(readr)
library(dplyr)
library(stringr)

# condensing municipality information rules
# if COMUNINF exists, use that
# otherwise, use residence
# source indicators will track what municipality ID was used: 
# infection_recorded --> COMUNINF exists, different from residence
# autochthonous_residence --> COMUNINF exists, same as residence
# residence_fallback --> COMUNINF does not exist, used residence instead

years <- 16:25
folder_path <- "/Users/allison/Desktop/Research_BentoLab/SP26/data"

# columns to keep when cleaning
req_col <- c("TP_NOT", "DT_NOTIFIC", "SEM_NOT", "NU_ANO", "ID_MN_RESI", "COMUNINF", "CLASSI_FIN", "CRITERIO", "CS_SEXO", "NU_IDADE_N")

# lists to store dataframes
clean_list <- list() # list with cleaned dataframes

for (year in years) {
  name <- paste0("nCHIKBR",year,".csv")
  path <- file.path(paste0(folder_path,"/chikv"), name)
  
  if (file.exists(path)) {
    df <- read_csv(path, show_col_types = FALSE)
    
    remove <- setdiff(req_col, names(df))
    df[remove] <- NA
    df <- df[req_col]
    clean_list[[length(clean_list)+1]] <- df
    cat("Processed", name, "\n")
  } 
  else {
    cat(name, "not found\n")
  }
}  

# valid IBGE codes
pop_codes <- read_csv("data/Pop_size_estimative_2010_2025_IBGE_collected_in_April_14_2026_clean.csv",
                      show_col_types = FALSE) %>%
  filter(code_muni != 110000) %>%                        
  transmute(muni = substr(as.character(code_muni), 1, 6)) %>%
  distinct() %>%
  pull(muni)

combined <- bind_rows(clean_list)

combined <- combined %>%
  mutate(
    COMUNINF = na_if(as.character(COMUNINF), ""),
    COMUNINF = na_if(as.character(COMUNINF), "0"),
    COMUNINF = ifelse(COMUNINF %in% pop_codes, COMUNINF, NA)   # invalid COMUNINF --> missing
  ) %>%
  mutate(
    source_indicator = case_when(
      !is.na(COMUNINF) & COMUNINF != ID_MN_RESI ~ "infection_recorded",
      !is.na(COMUNINF) & COMUNINF == ID_MN_RESI ~ "autochthonous_residence",
      is.na(COMUNINF) ~ "residence_fallback"
    ),
    MUNI = ifelse(!is.na(COMUNINF), COMUNINF, ID_MN_RESI)
  ) %>%
  select(TP_NOT,DT_NOTIFIC, SEM_NOT, NU_ANO, MUNI, CLASSI_FIN, CRITERIO, CS_SEXO, NU_IDADE_N, source_indicator)


uf_lookup <- tibble(
  UF_code = c(
    "11","12","13","14","15","16","17",
    "21","22","23","24","25","26","27","28","29",
    "31","32","33","35",
    "41","42","43",
    "50","51","52","53"
  ),
  UF = c(
    "RO","AC","AM","RR","PA","AP","TO",
    "MA","PI","CE","RN","PB","PE","AL","SE","BA",
    "MG","ES","RJ","SP",
    "PR","SC","RS",
    "MS","MT","GO","DF"
  )
)



combined <- combined %>%
  mutate(UF_code = substr(MUNI, 1, 2)) %>%
  left_join(uf_lookup, by = "UF_code") %>% 
  filter(NU_ANO <= 2025)

excluded_muni <- combined %>% filter(!MUNI %in% pop_codes)  
nrow(excluded_muni)                                           
excluded_muni %>% filter(CLASSI_FIN %in% c(1,13)) %>% nrow()  

write_csv(excluded_muni, "data/chikv_unmatched_muni_codes.csv") 
# Thirty of 1,238,492 confirmed chikungunya notifications (0.002%) had 
# municipality codes that could not be matched to IBGE municipality codes 
# and were excluded from municipality-level analyses

combined <- combined %>%
  filter(MUNI %in% pop_codes)

write_csv(combined, file.path(folder_path, "ncombined_chikv_2016_2025.csv"))

# CALCULATE INCIDENCE BY MUNI, YEAR
# output columns:
# MUNI  year cases_tot pop_tot incidence_tot 
# Incidence calculated per 100k
library(tidyr)

# 1. load datasets
cases_raw <- combined
cases_raw <- cases_raw %>%
  filter(NU_ANO >= 2016 & NU_ANO <= 2025)
pop <- read_csv("data/Pop_size_estimative_2010_2025_IBGE_collected_in_April_14_2026_clean.csv") %>%
  filter(Year >= 2016 & Year <= 2025) 

# 2. aggregate total cases per muni-date
total_cases <- cases_raw %>%
  filter(CLASSI_FIN %in% c(1, 13)) %>%                           # keep confirmed case: use portalsinan dictionary (5 - discard, 13 = chikungunya)
  transmute(
    muni = as.character(MUNI),
    year = NU_ANO
  ) %>%
  group_by(muni, year) %>%    
  summarise(total_cases = n(), .groups = "drop")

# 3. prep population data
pop_total <- pop %>%
  select(code_muni, Year, Pop_estimated, Sigla) %>%
  mutate(Pop_estimated = trimws(Pop_estimated),
         Pop_estimated = sub("\\s*\\([0-9]+\\)$", "", Pop_estimated),
         Pop_estimated = gsub("\\.", "", Pop_estimated)) %>% # remove footnote-- sometimes formatted with a space, ranges from 1-19 
  transmute(
    muni = substr(as.numeric(code_muni), 1, 6),
    year = Year,
    pop_tot = as.numeric(Pop_estimated)
  ) 
 

# 4. join cases with population
incidence_data <- total_cases %>%
  left_join(pop_total, by = c("muni", "year")
  )

# 5. calculate incidence per 100k
incidence_data <- incidence_data %>%
  rename(
    cases_tot = total_cases
  ) %>%
  mutate(
    incidence = round(cases_tot / pop_tot * 100000, 3)
  )

# 6. checks
# check population data for duplicate muni-year combinations
pop_total %>%
  count(muni, year) %>%
  filter(n > 1)

# check for missing population after joining
incidence_data %>%
  filter(is.na(pop_tot) | pop_tot <= 0)

# 7. save output
incidence_data <- incidence_data %>% 
  filter(!is.na(pop_tot)) 

write_csv(incidence_data, "/Users/allison/Desktop/Research_BentoLab/SP26/nchikv_incidence_2016_2025.csv")

# compare with old data (TEMP)
old_data <- read_csv("data/bd_CHIKV_positive_2013_2025_english.csv") %>% 
  filter(Auxiliary_Virus == "Chikungunya") %>%
  filter (between(Collection_Year, 2016, 2025))

new_data <- read_csv("data/ncombined_chikv_2016_2025.csv")
new_data <- new_data %>%
  filter(CLASSI_FIN %in% c(1, 13))

nrow(old_data)
nrow(new_data)

nrow(new_data) - nrow(old_data)


library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)

# Count old dataset by month
old_monthly <- old_data %>%
  mutate(month_date = floor_date(mdy(Collection_Date), "month")) %>%
  filter(!is.na(month_date)) %>%
  count(month_date, name = "old_cases")

# Count new dataset by month
new_monthly <- new_data %>%
  mutate(month_date = floor_date(as.Date(DT_NOTIFIC), "month")) %>%
  filter(!is.na(month_date)) %>%
  count(month_date, name = "new_cases")

compare_counts <- full_join(
  old_monthly, new_monthly,
  by = "month_date"
) %>%
  replace_na(list(old_cases = 0, new_cases = 0))

# Plot 1: Monthly differences
compare_counts <- compare_counts %>%
  mutate(difference = new_cases - old_cases)

ggplot(compare_counts, aes(x = month_date, y = difference)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Monthly Difference in Chikv Notifications",
    subtitle = "New dataset minus old dataset",
    x = "Month",
    y = "Difference in case counts"
  ) +
  theme_classic()


# Plot 2: Overlay old and new monthly counts
compare_long <- compare_counts %>%
  pivot_longer(
    cols = c(old_cases, new_cases),
    names_to = "dataset",
    values_to = "n"
  ) %>%
  mutate(
    dataset = recode(
      dataset,
      old_cases = "Old",
      new_cases = "New"
    )
  )

ggplot(compare_long,
       aes(x = month_date, y = n, color = dataset, group = dataset)) +
  geom_line() +
  labs(
    title = "Chikv Monthly Notifications: Old vs. New",
    x = "Month",
    y = "Number of notifications",
    color = "Dataset"
  ) +
  theme_classic()
