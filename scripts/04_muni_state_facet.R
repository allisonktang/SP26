# muni grouped by state incidence
# goal: plot muni-level incidence grouped by state for all three arboviruses
# input: Brazil_arbovirus_monthly_data_2016_2025.csv, Brazil_arbovirus_state_yearly_2016_2025.csv
# notes: state level incidence is calculated by (sum cases/ sum population), not an average of municipality incidence rates

rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidytext)
library(scales)

merged     <- read_csv("data/nBrazil_arbovirus_monthly_data_2016_2025.csv")
state_year <- read_csv("data/nBrazil_arbovirus_state_yearly_2016_2025.csv")

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
  summarise(max_incidence = max(incidence, na.rm = T)) %>%
  arrange(desc(max_incidence)) %>%
  pull(state)

muniLong <- merged %>%
  filter(!is.na(state)) %>%
  mutate(state = factor(state, levels = state_order)) %>%
  group_by(muni, state, year) %>%
  summarise(
    dengueCases = sum(dengueCases, na.rm = T),
    zikaCases = sum(zikaCases,   na.rm = T),
    chikvCases  = sum(chikvCases,  na.rm = T),
    pop_tot     = first(pop_tot), #repeat value for each month
    .groups = "drop"
  ) %>%
  mutate(
    dengueIncidence = dengueCases / pop_tot * 100000,
    zikaIncidence   = zikaCases   / pop_tot * 100000,
    chikvIncidence  = chikvCases  / pop_tot * 100000
  ) %>%
  pivot_longer(
    cols = c(dengueIncidence, zikaIncidence, chikvIncidence),
    names_to = "disease", values_to = "incidence"
  ) %>%
  mutate(
    muni = sprintf("%06d", as.integer(muni)),
    disease = recode(
      disease,
      dengueIncidence = "Dengue",
      zikaIncidence = "Zika",
      chikvIncidence = "Chikungunya"
    )
  )

plot_muni_by_disease <- function(df, disease_name) {
  d <- df %>% 
    filter(disease == disease_name, 
           !is.na(incidence),
           incidence > 0)
  
  ggplot(d, aes(
    x = factor(year),
    y = reorder_within(muni, incidence, state),
    fill = incidence
  )) +
    geom_tile() +
    facet_wrap(~ state, scales = "free_y") +
    scale_y_reordered() +
    scale_fill_gradient(
      low = "#FCF0CE", high = "#D4180A", trans = "log10", na.value="white",
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
      strip.text = element_text(face = "bold"),
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

ggsave(
  "figures/nchikv_muni_state_heatmap.png",
  plot = chikv_plot,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/nzika_muni_state_heatmap.png",
  plot = zika_plot,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/ndengue_muni_state_heatmap.png",
  plot = dengue_plot,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)

# plot with a fixed scale
shared_breaks <- c(0, 10, 50, 100, 500, 1000, Inf)
shared_colors <- colorRampPalette(c("#FCF0CE", "#D4180A"))(6)
bin_labels <- c("0-10", "10-50", "50-100", "100-500", "500-1000", "1000+")

plot_muni_by_disease_fixed <- function(df, disease_name, breaks, colors, labels) {
  d <- df %>% 
    filter(disease == disease_name, !is.na(incidence), incidence >= 0) %>%
    mutate(
      incidence_bin = cut(incidence, breaks = breaks, labels = labels,
                          include.lowest = TRUE, right = FALSE)
    )
  
  ggplot(d, aes(
    x = factor(year),
    y = reorder_within(muni, incidence, state),
    fill = incidence_bin
  )) +
    geom_tile() +
    facet_wrap(~ state, scales = "free_y") +
    scale_y_reordered() +
    scale_fill_manual(values = setNames(colors, labels), na.value = "white", drop = FALSE) +
    labs(x = "Year", y = "Municipality", fill = paste(disease_name, "\nincidence\nper 100k"),
         title = paste(disease_name, "incidence by municipality, grouped by state")) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.text = element_text(),
      legend.position = "bottom",
      plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 10))
    )
}

dengue_plot_fixed <- plot_muni_by_disease_fixed(muniLong, "Dengue", shared_breaks, shared_colors, bin_labels)
zika_plot_fixed   <- plot_muni_by_disease_fixed(muniLong, "Zika", shared_breaks, shared_colors, bin_labels)
chikv_plot_fixed  <- plot_muni_by_disease_fixed(muniLong, "Chikungunya", shared_breaks, shared_colors, bin_labels)

print(dengue_plot_fixed)

ggsave(
  "figures/nchikv_fixed_muni_state_heatmap.png",
  plot = chikv_plot_fixed,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/nzika_fixed_muni_state_heatmap.png",
  plot = zika_plot_fixed,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/ndengue_fixed_muni_state_heatmap.png",
  plot = dengue_plot_fixed,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)

# check SC
muniLong %>%
  filter(
    state == "SC",
    disease == "Chikungunya"
  ) %>%
  summarise(
    min_incidence = min(incidence, na.rm = TRUE),
    median_incidence = median(incidence, na.rm = TRUE),
    mean_incidence = mean(incidence, na.rm = TRUE),
    max_incidence = max(incidence, na.rm = TRUE)
  )
muniLong %>%
  filter(
    state == "SC",
    disease == "Chikungunya"
  ) %>%
  arrange(desc(incidence)) %>%
  select(muni, year, incidence) %>%
  head(20)

muniLong %>%
  filter(
    disease == "Chikungunya",
    incidence > 0
  ) %>%
  summarise(
    min = min(incidence),
    q25 = quantile(incidence, 0.25),
    median = median(incidence),
    q75 = quantile(incidence, 0.75),
    max = max(incidence)
  )
muniLong %>%
  filter(
    disease == "Chikungunya",
    state == "SC",
    incidence > 0
  ) %>%
  summarise(
    min = min(incidence),
    q25 = quantile(incidence, 0.25),
    median = median(incidence),
    q75 = quantile(incidence, 0.75),
    max = max(incidence)
  )

sc_plot <- muniLong %>%
  filter(
    state == "SC",
    disease == "Chikungunya",
    incidence > 0
  ) %>%
  ggplot(
    aes(
      x = factor(year),
      y = reorder_within(muni, incidence, state),
      fill = incidence
    )
  ) +
  geom_tile() +
  facet_wrap(~state, scales = "free_y") +
  scale_y_reordered() +
  scale_fill_gradient(
    low = "#FCF0CE",
    high = "#D4180A",
    trans = "log"
  ) + 
  labs(x = "Year", y = "Municipality", fill = paste("Chikungunya", "\nincidence\nper 100k")) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_blank(),    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.line = element_blank(),
    strip.background = element_blank(),
  )

sc_plot

ggsave(
  "figures/sc_plot_test.png",
  plot = sc_plot,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)


# see global scale
muniLong %>%
  filter(
    disease == "Chikungunya",
    incidence > 0
  ) %>%
  summarise(
    min = min(incidence, na.rm = TRUE),
    q25 = quantile(incidence, 0.25, na.rm = TRUE),
    median = median(incidence, na.rm = TRUE),
    q75 = quantile(incidence, 0.75, na.rm = TRUE),
    max = max(incidence, na.rm = TRUE)
  )

muniLong %>%
  filter(
    disease == "Chikungunya",
    incidence > 0
  ) %>%
  arrange(desc(incidence)) %>%
  select(muni, state, year, incidence) %>%
  head(10)

muniLong %>%
  filter(
    state == "MG",
    disease == "Chikungunya",
    incidence > 0
  ) %>%
  summarise(
    min = min(incidence),
    q25 = quantile(incidence, 0.25),
    median = median(incidence),
    q75 = quantile(incidence, 0.75),
    max = max(incidence)
  )

###########################################
#         plot geographic map             #
###########################################
library(geobr)
library(sf)

br_munis <- read_municipality(
  year = 2020,
  simplified = TRUE,
  cache = FALSE
)

names(br_munis)

br_munis <- br_munis %>%
  mutate(
    code_muni = substr(as.character(code_muni),1,6)
  )

brazil_border <- br_munis %>%
  sf::st_union() %>%
  sf::st_cast("MULTIPOLYGON") %>%
  sf::st_cast("POLYGON") %>%
  sf::st_exterior_ring()

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

# create plotting function that facets by year (if you want to do individual year, add year_value to parameters, filter by year, remove facet wrap)
plot_muni_map <- function(data, geography, disease_name) {
  
  map_data <- geography %>%
    left_join(
      data %>%
        filter(
          disease == disease_name),
      by = c("code_muni" = "muni")
    )
  
  ggplot(map_data) + 
    geom_sf(
      aes(fill = incidence),
      color = NA
    ) +
  geom_sf(
    data = brazil_border,
    fill = NA,
    color = "black",
    linewidth = 0.1
    ) +
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
        " Incidence by Municipality (log scale)"
      )
    ) +
    base_theme +
    theme(
      panel.border = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    ) +
    coord_sf(datum = NA) +
    facet_wrap(~year)
}


# create three graphs showing geographic heatmap throughout time
diseases <- unique(muniLong$disease)

for (disease_name in diseases) {
  
  p <- plot_muni_map(
    muniLong,
    br_munis,
    disease_name
  )
  
  ggsave(
    filename = paste0(
      "figures/n",
      tolower(ifelse(disease_name == "Chikungunya", "chikv", disease_name)),
      "_geo_muni.png"    ),
    plot = p,
    width = 12,
    height = 8,
    units = "in",
    dpi = 600
  )
  
  cat("Saved:", disease_name, "\n")
}

# cumulative incidence

muni_cumulative <- muniLong %>%
  group_by(muni, disease) %>%
  summarise(
    cumulative_incidence = sum(incidence, na.rm = TRUE),
    .groups = "drop"
  )

glimpse(muni_cumulative)

plot_cumulative_muni_map <- function(data, geography, disease_name) {
  
  map_data <- geography %>%
    left_join(
      data %>%
        filter(disease == disease_name),
      by = c("code_muni" = "muni")
    )
  
  ggplot(map_data) +
    geom_sf(
      aes(fill = cumulative_incidence),
      color = NA
    ) +
    geom_sf(
      data = brazil_border,
      fill = NA,
      color = "black",
      linewidth = 0.2
    ) +
    scale_fill_gradient(
      low = "#FCF0CE",
      high = "#D4180A",
      na.value = "white",
      trans = "log10",
      breaks = trans_breaks(
        "log10",
        function(x) 10^x
      ),
      labels = trans_format(
        "log10",
        math_format(10^.x)
      )
    ) +
    
    labs(
      fill = "Cumulative incidence\nper 100k",
      title = paste0(
        disease_name,
        " Cumulative Incidence by Municipality, 2016–2025"
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


# create one map for each disease
dengue_cumulative <- plot_cumulative_muni_map(
  muni_cumulative,
  br_munis,
  "Dengue"
)

zika_cumulative <- plot_cumulative_muni_map(
  muni_cumulative,
  br_munis,
  "Zika"
)

chikv_cumulative <- plot_cumulative_muni_map(
  muni_cumulative,
  br_munis,
  "Chikungunya"
)

ggsave(
  "figures/ndengue_cumulative_muni_map.png",
  plot = dengue_cumulative,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/nzika_cumulative_muni_map.png",
  plot = zika_cumulative,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/nchikv_cumulative_muni_map.png",
  plot = chikv_cumulative,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

