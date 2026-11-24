# 07: windowed analysis
# goal: one file per window length (3, 6, 12 months), one row per municipality x year-month, with:
#         - cumulative cases for each disease over the window
#         - population
#         - cumulative incidence per 100k for each disease
#         - temperature, precipitation, pressure, humidity: min, mean, max over the window
#         - ENSO: mean over the window
# input:  data/nBrazil_arbovirus_monthly_data_2016_2025.csv
#         data/combined_all_climate_data_2016_2024(7).rds, data/climate_daily_2025.csv
#         data/ocean_climate_oscillations_update_2025.csv
# output: data/analysis_table_3month.csv, analysis_table_6month.csv, analysis_table_12month.csv
# notes:
#   - window = the row's month plus the previous k-1 months (trailing). Windows that reach back
#     before Jan 2016 are NA, so the first k-1 months of each muni are NA for window variables
#   - cases: counts come from individual notification records, so a month with no records counts as 0
#   - cumulative incidence = cumulative cases in the window / population in the row's month * 100,000
#   - climate min/max = lowest/highest daily value in the window
#     climate mean = average of the monthly means
#   - precipitation sum/mean/max use daily precipitation totals (mm/day)
#   - ENSO is same value for every municipality

rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(zoo)

windows    <- c(3, 6, 12)
all_months <- seq(as.Date("2016-01-01"), as.Date("2025-12-01"), by = "month")

state_lookup <- c(
  "11" = "RO", "12" = "AC", "13" = "AM", "14" = "RR", "15" = "PA",
  "16" = "AP", "17" = "TO", "21" = "MA", "22" = "PI", "23" = "CE",
  "24" = "RN", "25" = "PB", "26" = "PE", "27" = "AL", "28" = "SE",
  "29" = "BA", "31" = "MG", "32" = "ES", "33" = "RJ", "35" = "SP",
  "41" = "PR", "42" = "SC", "43" = "RS", "50" = "MS", "51" = "MT",
  "52" = "GO", "53" = "DF"
)

# ================================================================================
#            1. cases + population (muni-month)
# ================================================================================
merged <- read_csv("data/nBrazil_arbovirus_monthly_data_2016_2025.csv", show_col_types = FALSE)

# complete() gives every muni all 120 months
cases <- merged %>%
  mutate(muni = as.character(muni),
         year_month = as.Date(year_month)) %>%
  select(muni, year_month, pop_tot, dengueCases, zikaCases, chikvCases) %>%
  complete(muni, year_month = all_months) %>%
  # no records = 0 cases, but only in months where the municipality exists (has population)
  mutate(across(c(dengueCases, zikaCases, chikvCases),
                ~ if_else(is.na(pop_tot), NA_real_, coalesce(.x, 0)))) %>%
  mutate(dengue_incidence = dengueCases / pop_tot * 100000,
         zika_incidence   = zikaCases   / pop_tot * 100000,
         chikv_incidence  = chikvCases  / pop_tot * 100000)


# ================================================================================
#            2. climate: convert daily observations into muni-month observations (min / mean / max)
# ================================================================================
climate24 <- readRDS("data/combined_all_climate_data_2016_2024(7).rds")
climate25 <- read_csv("data/climate_daily_2025.csv") %>% select(-1, -last_col())

climate <- bind_rows(climate24, climate25) 

# check
range(climate$date) 
climate %>% count(geocodigo) %>% count(n)        # days per municipality: should be 3653 for nearly all

climate_monthly <- suppressWarnings(           # min/max of an all-NA month warns; handled below
  climate %>%
    mutate(date = as.Date(date)) %>%
    filter(date >= as.Date("2016-01-01"), date <= as.Date("2025-12-31")) %>%
    mutate(muni = substr(as.character(geocodigo), 1, 6),
           year_month = as.Date(format(date, "%Y-%m-01"))) %>%
    group_by(muni, year_month) %>%
    summarise(
      temp_min    = min(temp_min,  na.rm = TRUE),
      temp_mean   = mean(temp_med, na.rm = TRUE),
      temp_max    = max(temp_max,  na.rm = TRUE),
      precip_sum  = if (all(is.na(precip_tot))) NA_real_ else sum(precip_tot, na.rm = TRUE),   # monthly total (mm)
      press_min   = min(pressao_min,  na.rm = TRUE),
      press_mean  = mean(pressao_med, na.rm = TRUE),
      press_max   = max(pressao_max,  na.rm = TRUE),
      humid_min   = min(umid_min,  na.rm = TRUE),
      humid_mean  = mean(umid_med, na.rm = TRUE),
      humid_max   = max(umid_max,  na.rm = TRUE),
      .groups = "drop"
    )
) %>%
  # all-NA months come out as Inf / NaN; make them NA
  mutate(across(-c(muni, year_month), ~ if_else(is.finite(.x), .x, NA_real_)))

rm(climate)
rm(climate24)
rm(climate25)

# ================================================================================
#            3. ENSO: weekly -> monthly mean
# ================================================================================

ocean <- read_csv("data/ocean_climate_oscillations_update_2025.csv", show_col_types = FALSE)

enso_monthly <- ocean %>%
  distinct(date, .keep_all = TRUE) %>%
  mutate(year_month = as.Date(format(date, "%Y-%m-01"))) %>%
  filter(year_month >= as.Date("2016-01-01"), year_month <= as.Date("2025-12-01")) %>%
  group_by(year_month) %>%
  summarise(enso = mean(enso, na.rm = TRUE), .groups = "drop") %>%
  complete(year_month = all_months)

# ================================================================================
#            4. combine into one monthly table
# ================================================================================
monthly <- cases %>%
  left_join(climate_monthly, by = c("muni", "year_month")) %>%
  left_join(enso_monthly,    by = "year_month") %>%
  arrange(muni, year_month)

# ================================================================================
#            5. build one table per window
# ================================================================================
build_window_table <- function(k) {
  roll <- function(x, f) zoo::rollapplyr(x, width = k, FUN = f, fill = NA)  # NA if any month in window is NA
  
  monthly %>%
    group_by(muni) %>%
    mutate(
      dengue_cum_cases = roll(dengueCases, sum),
      zika_cum_cases   = roll(zikaCases,   sum),
      chikv_cum_cases  = roll(chikvCases,  sum),
      
      dengue_cum_incidence = dengue_cum_cases / pop_tot * 100000,
      zika_cum_incidence   = zika_cum_cases   / pop_tot * 100000,
      chikv_cum_incidence  = chikv_cum_cases  / pop_tot * 100000,
      
      precip_sum = roll(precip_sum, sum),        # total mm over the window
      
      across(c(temp_min, press_min, humid_min), ~ roll(.x, min)),
      across(c(temp_mean, press_mean, humid_mean), ~ roll(.x, mean)),
      across(c(temp_max, press_max, humid_max), ~ roll(.x, max)),
      enso_mean = roll(enso, mean)
    ) %>%
    ungroup() %>%
    filter(!is.na(pop_tot)) %>%                  # incidence is undefined without population
    mutate(state = unname(state_lookup[substr(muni, 1, 2)])) %>%
    select(year_month, state, muni,
           dengue_cum_cases, zika_cum_cases, chikv_cum_cases,
           pop_tot,
           dengue_cum_incidence, zika_cum_incidence, chikv_cum_incidence,
           temp_min, temp_mean, temp_max,
           precip_sum,
           press_min, press_mean, press_max,
           humid_min, humid_mean, humid_max,
           enso_mean)
}

tables <- list()

for (k in windows) {
  tables[[as.character(k)]] <- build_window_table(k)
  write_csv(tables[[as.character(k)]], sprintf("data/%dmonth_table.csv", k))
  cat("Saved:", k, "month table,", nrow(tables[[as.character(k)]]), "rows\n")
}

# ================================================================================
#            6. checks
# ================================================================================
# municipalities in the table with no climate data in at least one month
monthly %>%
  filter(!is.na(pop_tot), is.na(temp_mean)) %>%
  distinct(muni) %>%
  nrow()

# rows dropped for missing population
sum(is.na(monthly$pop_tot))

# missing values per column (first k-1 months should be NA for window variables, nothing else)
colSums(is.na(tables[["3"]]))

# rolling sum vs. manual sum: Sao Paulo, 3-month window ending June 2023 (Apr + May + Jun)
monthly %>%
  filter(muni == "355030", year_month >= as.Date("2023-04-01"), year_month <= as.Date("2023-06-01")) %>%
  summarise(manual_dengue = sum(dengueCases))

tables[["3"]] %>%
  filter(muni == "355030", year_month == as.Date("2023-06-01")) %>%
  select(dengue_cum_cases)

# cumulative incidence across a year boundary: Sao Paulo, 3-month window ending Feb 2023 (Dec 2022 + Jan + Feb)
# should equal the sum of the three monthly incidences, each using its own year's population
monthly %>%
  filter(muni == "355030", year_month >= as.Date("2022-12-01"), year_month <= as.Date("2023-02-01")) %>%
  summarise(manual_dengue_incidence = sum(dengue_incidence))

tables[["3"]] %>%
  filter(muni == "355030", year_month == as.Date("2023-02-01")) %>%
  select(dengue_cum_incidence)
