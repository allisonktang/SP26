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

merged <- read_csv("data/Brazil_arbovirus_monthly_data_2016_2025.csv")

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


# edit state order
state_order <- c(
  "AC", "AM", "PA", "RR", "RO", "AP", "TO",
  "PI", "BA", "MA", "PE", "CE", "AL", "SE", "RN", "PB",
  "MT", "GO", "MS","DF",
  "MG", "ES", "RJ", "SP",
  "PR", "SC", "RS"
)

state_year <- state_year %>%
  mutate(state = factor(state, levels = rev(state_order)))


base_theme <- theme_classic() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.line = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b=10))
  )



state_year_heatmap <- ggplot(state_year, aes(x=factor(year), y=state, fill = incidence)) +
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

write_csv(state_month, "data/Brazil_arbovirus_state_monthly_2016_2025.csv")
write_csv(state_year,  "data/Brazil_arbovirus_state_yearly_2016_2025.csv")

# save as png, tiff, and pdf
ggsave(
  "figures/state_year_heatmap.png",
  plot = state_year_heatmap,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/state_year_heatmap.tiff",
  plot = state_year_heatmap,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "figures/state_year_heatmap.pdf",
  plot = state_year_heatmap,
  width = 8,
  height = 10,
  units = "in"
)

###########################################
#         plot geographic map             #
###########################################
library(geobr)
library(sf)

 br_states <- read_state(
   year = 2020,
   simplified = T
 )

names(br_states)

state_labels <- data.frame(
  state = c("RN", "PB", "PE", "AL", "SE", "ES", "RJ", "SC"),
  x = c(-34, -34, -34, -35, -36, -38.5, -40, -47),
  y = c(-5, -7, -8.1, -10, -12, -20, -23, -28)
)

# create plotting function
plot_state_map <- function(data, geography, disease_name, year_value) {
    map_data <- geography %>%
      left_join(
        data %>%
          filter(
            disease == disease_name,
            year == year_value
          ),
        by = c("abbrev_state" = "state")
      )
    
  ggplot(map_data) + 
    geom_sf(
      aes(fill = incidence),
      color = "black",
      linewidth = 0.2
    ) +
    geom_sf_text(
      data = map_data %>%
        filter(!abbrev_state %in% c("RN", "PB", "PE", "AL", "SE", "ES", "RJ", "SC")),
      aes(label = abbrev_state),
      size = 3
    ) +
    geom_text(
      data = state_labels,
      aes(x = x, y = y, label = state),
      size = 3
    )+
    scale_fill_gradient(
      low = "#FCF0CE",
      high = "#D4180A",
      na.value = "white",
      trans = "log",
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    labs(
      fill = "Incidence\nper 100k",
      title = paste0(
        disease_name,
        " Incidence by State, ",
        year_value,
        " (log scale)"
      )
    ) +
    base_theme +
    theme(
      panel.border = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    ) +
    coord_sf(datum = NA)
}


dengue_map <- plot_state_map(
  state_year,
  br_states,
  "Dengue",
  2025
)

zika_map <- plot_state_map(
  state_year,
  br_states,
  "Zika",
  2025
)

chikv_map <- plot_state_map(
  state_year,
  br_states,
  "Chikungunya",
  2025
)

dengue_map
zika_map
chikv_map

ggsave(
  "figures/dengue_geo_heatmap_2025.png",
  plot = dengue_map,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/zika_geo_heatmap_2025.png",
  plot = zika_map,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/chikv_geo_heatmap_2025.png",
  plot = chikv_map,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)
