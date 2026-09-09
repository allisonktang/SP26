# state aggregation
# goal: recalculate incidence at state level by summing municipality cases and dividing by the sum of muni population for all three arboviruses
#       plot as a heatmap (state x year, log10)
# input: Brazil_arbovirus_monthly_data_2016_2025.csv
# notes: state level incidence is calculated by (sum cases/ sum population), not an average of municipality incidence rates

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

merged <- read_csv("Brazil_arbovirus_monthly_data_2016_2025.csv")

state_lookup <- c(
  "11" = "RO", "12" = "AC", "13" = "AM", "14" = "RR", "15" = "PA",
  "16" = "AP", "17" = "TO", "21" = "MA", "22" = "PI", "23" = "CE",
  "24" = "RN", "25" = "PB", "26" = "PE", "27" = "AL", "28" = "SE",
  "29" = "BA", "31" = "MG", "32" = "ES", "33" = "RJ", "35" = "SP",
  "41" = "PR", "42" = "SC", "43" = "RS", "50" = "MS", "51" = "MT",
  "52" = "GO", "53" = "DF"
)

merged <- merged %>% 
  mutate(state_code = substr(muni,1,2),
         state = state_lookup[state_code])

# CHECK: municipality codes map to a known state, expected val = 0
merged %>% 
  filter(is.na(state)) %>% 
  distinct(muni, state_code) %>% 
  nrow()

# state-year incidence
pop_muni_year <- merged %>% 
  distinct(muni, state, year, pop_tot) # muni monthly population is currently duplicated 12x, remove duplicate

cases_state_year <- merged %>%
  group_by(state,year) %>%
  summarise(
    dengueCases = sum(dengueCases, na.rm = T),
    zikaCases = sum(zikaCases, na.rm = T),
    chikvCases = sum(chikvCases, na.rm = T),
    .groups = "drop"
  )

pop_state_year <- pop_muni_year %>%
  group_by(state, year) %>%
  summarise(pop_tot = sum(pop_tot, na.rm = T), .groups = "drop")

state_year <- cases_state_year %>%
  left_join(pop_state_year, by = c("state", "year")) %>%
  mutate(
    dengueIncidence = dengueCases / pop_tot * 100000,
    zikaIncidence = zikaCases / pop_tot * 100000,
    chikvIncidence = chikvCases / pop_tot * 100000
  ) %>%
  pivot_longer(cols = ends_with("Incidence"), names_to = "disease", values_to = "incidence") %>%
  mutate(disease = recode(disease,
                          dengueIncidence = "Dengue",
                          zikaIncidence = "Zika",
                          chikvIncidence = "Chikungunya"),
         disease = factor(disease, levels = c("Dengue", "Zika", "Chikungunya")))

base_theme <- theme_classic() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.line = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b=10))
  )



state_year_heatmap <- ggplot(state_year, aes(x=factor(year), y=reorder(state, incidence, FUN = max), fill = incidence)) +
  geom_tile(color = NA, linewidth = 0) +
  facet_wrap(~ disease, ncol = 1) +
  scale_fill_gradient(
    low="grey90",
    high = "#D4180A",
    na.value = "white",
    trans = "log",
    breaks = trans_breaks("log10",function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "Year", y = "State", fill = "Incidence\nper 100k", title = "State-Level Yearly Incidence (log scale)") + 
  base_theme 

state_year_heatmap


# create state-month incidence (finer resolution, needed for muni-level grouped by state)

state_month <- merged %>%
  group_by(state, year, month, year_month) %>%
  summarise(
    dengueCases = sum(dengueCases, na.rm = T),
    zikaCases   = sum(zikaCases,   na.rm = T),
    chikvCases  = sum(chikvCases,  na.rm = T),
    pop_tot     = sum(pop_tot,     na.rm = T),
    .groups = "drop"
  ) %>%
  mutate(
    dengueIncidence = dengueCases / pop_tot * 100000,
    zikaIncidence   = zikaCases   / pop_tot * 100000,
    chikvIncidence  = chikvCases  / pop_tot * 100000
  )

write_csv(state_month, "Brazil_arbovirus_state_monthly_2016_2025.csv")
write_csv(state_year,  "Brazil_arbovirus_state_yearly_2016_2025.csv")

