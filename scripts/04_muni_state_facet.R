# muni grouped by state incidence
# goal: plot muni-level incidence grouped by state for all three arboviruses
# input: Brazil_arbovirus_monthly_data_2016_2025.csv, Brazil_arbovirus_state_yearly_2016_2025.csv
# notes: state level incidence is calculated by (sum cases/ sum population), not an average of municipality incidence rates

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidytext)
library(scales)

merged     <- read_csv("Brazil_arbovirus_monthly_data_2016_2025.csv")
state_year <- read_csv("Brazil_arbovirus_state_yearly_2016_2025.csv")

state_lookup <- c(
  "11" = "RO", "12" = "AC", "13" = "AM", "14" = "RR", "15" = "PA",
  "16" = "AP", "17" = "TO", "21" = "MA", "22" = "PI", "23" = "CE",
  "24" = "RN", "25" = "PB", "26" = "PE", "27" = "AL", "28" = "SE",
  "29" = "BA", "31" = "MG", "32" = "ES", "33" = "RJ", "35" = "SP",
  "41" = "PR", "42" = "SC", "43" = "RS", "50" = "MS", "51" = "MT",
  "52" = "GO", "53" = "DF"
)

merged <- merged %>%
  mutate(state_code = substr(muni, 1, 2), state = state_lookup[state_code])

state_order <- state_year %>%
  group_by(state) %>%
  summarise(max_Incidence = max(incidence, na.rm = TRUE)) %>%
  arrange(desc(max_Incidence)) %>%
  pull(state)

muniLong <- merged %>%
  filter(!is.na(state)) %>%
  mutate(state = factor(state, levels = state_order)) %>%
  pivot_longer(
    cols = c(dengueIncidence, zikaIncidence, chikvIncidence),
    names_to = "disease", values_to = "incidence"
  ) %>%
  mutate(disease = recode(disease,
                          dengueIncidence = "Dengue", zikaIncidence = "Zika", chikvIncidence = "Chikungunya"))


plot_muni_by_disease <- function(df, disease_name) {
  d <- df %>% filter(disease == disease_name, !is.na(incidence), incidence > 0)
  
  ggplot(d, aes(
    x = factor(year),
    y = reorder_within(muni, incidence, state),
    fill = incidence
  )) +
    geom_tile(color = "#fcf0ce") +
    facet_wrap(~ state, scales = "free_y") +
    scale_y_reordered() +
    scale_fill_gradient(
      low = "#FCF0CE", high = "#D4180A", trans = "log",
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    labs(x = "Year", y = "Municipality", fill = paste(disease_name, "\nincidence\nper 100k"),
         title = paste(disease_name, "incidence by municipality, grouped by state")) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
    #  strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 10))
    )
}

dengue_plot <- plot_muni_by_disease(muniLong, "Dengue")
zika_plot   <- plot_muni_by_disease(muniLong, "Zika")
chikv_plot  <- plot_muni_by_disease(muniLong, "Chikungunya")

print(dengue_plot)
print(zika_plot)
print(chikv_plot)


