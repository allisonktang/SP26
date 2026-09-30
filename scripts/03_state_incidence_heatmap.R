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
    low="#FCF0CE",
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

# state incidence heatmap, grouped by region, per disease
region_lookup <- c(
  "RO" = "North", "AC" = "North", "AM" = "North", "RR" = "North",
  "PA" = "North", "AP" = "North", "TO" = "North",
  "MA" = "Northeast", "PI" = "Northeast", "CE" = "Northeast", "RN" = "Northeast",
  "PB" = "Northeast", "PE" = "Northeast", "AL" = "Northeast", "SE" = "Northeast", "BA" = "Northeast",
  "MG" = "Southeast", "ES" = "Southeast", "RJ" = "Southeast", "SP" = "Southeast",
  "PR" = "South", "SC" = "South", "RS" = "South",
  "MS" = "Central-West", "MT" = "Central-West", "GO" = "Central-West", "DF" = "Central-West"
)

state_month_region <- state_month %>%
  mutate(region = region_lookup[state]) %>%
  filter(!is.na(region)) %>%
  pivot_longer(
    cols = ends_with("Incidence"), names_to = "disease", values_to = "incidence"
  ) %>%
  mutate(
    disease = recode(disease,
                     dengueIncidence = "Dengue", zikaIncidence = "Zika", chikvIncidence = "Chikungunya"),
    region  = factor(region, levels = c("North", "Northeast", "Central-West", "Southeast", "South"))
  )

#check
state_month_region %>% filter(is.na(region)) %>% distinct(state) %>% nrow()

plot_state_by_region <- function(df, disease_name) {
  d <- df %>% filter(disease == disease_name)
  
  ggplot(d, aes(
    x = year_month,
    y = tidytext::reorder_within(state, incidence, region, fun = max, na.rm = T),
    fill = incidence
  )) +
    geom_tile(color = NA) +
    facet_wrap(~ region, ncol = 1, scales = "free_y", strip.position = "right") +
    tidytext::scale_y_reordered() +
    scale_fill_gradient(
      low = "#FCF0CE", high = "#D4180A", na.value = "white",
      trans = "log10",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
    labs(x = "Year-Month", y = "State", fill = "Incidence\n(log10)",
         title = paste(disease_name, "incidence by state, grouped by region")) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text.y.right = element_text(angle = 0),
      plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 10))
    )
}

dengue_region <- plot_state_by_region(state_month_region, "Dengue")
zika_region   <- plot_state_by_region(state_month_region, "Zika")
chikv_region  <- plot_state_by_region(state_month_region, "Chikungunya")

ggsave(
  "figures/dengue_state_year_region_heatmap.png",
  plot = dengue_region,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/zika_state_year_region_heatmap.png",
  plot = zika_region,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/chikv_state_year_region_heatmap.png",
  plot = chikv_region,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)


# combined peak incidence order: rank states by their single highest incidence value
# across ALL THREE diseases (not per-disease), so all plots share one consistent order
incidence_state_order <- state_year %>%
  group_by(state) %>%
  summarise(peak_incidence = max(incidence, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(peak_incidence)) %>%
  pull(state)

plot_state_by_region_ordered <- function(df, disease_name, state_order) {
  d <- df %>% 
    filter(disease == disease_name) %>%
    mutate(state = factor(state, levels = state_order))
  
  ggplot(d, aes(x = year_month, y = state, fill = incidence)) +
    geom_tile(color = NA) +
    facet_wrap(~ region, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_fill_gradient(
      low = "#FCF0CE", high = "#D4180A", na.value = "white",
      trans = "log10",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
    labs(x = "Year-Month", y = "State", fill = "Incidence\n(log10)",
         title = paste(disease_name, "incidence by state, grouped by region")) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text.y.right = element_text(angle = 0),
      plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 10))
    )
}

dengue_region_order <- plot_state_by_region_ordered(state_month_region, "Dengue", incidence_state_order)
zika_region_order   <- plot_state_by_region_ordered(state_month_region, "Zika", incidence_state_order)
chikv_region_order  <- plot_state_by_region_ordered(state_month_region, "Chikungunya", incidence_state_order)

ggsave(
  "figures/dengue_ordered_state_year_region_heatmap.png",
  plot = dengue_region_order,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/zika_ordered_state_year_region_heatmap.png",
  plot = zika_region_order,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/chikv_ordered_state_year_region_heatmap.png",
  plot = chikv_region_order,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
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

# create plotting function that facets by year (if you want to do individual year, add year_value to parameters, filter by year, remove facet wrap)
plot_state_map <- function(data, geography, disease_name) {
    map_data <- geography %>%
      left_join(
        data %>%
          filter(
            disease == disease_name),
        by = c("abbrev_state" = "state")
      )
    
  ggplot(map_data) + 
    geom_sf(
      aes(fill = incidence),
      color = "black",
      linewidth = 0.1
    ) +
    geom_sf_text(
      data = map_data %>%
        filter(!abbrev_state %in% c("RN", "PB", "PE", "AL", "SE", "ES", "RJ", "SC")),
      aes(label = abbrev_state),
      size = 1.5
    ) +
    geom_text(
      data = state_labels,
      aes(x = x, y = y, label = state),
      size = 1.5
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
        " Incidence by State (log scale)"
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
diseases <- unique(state_year$disease)

for (disease_name in diseases) {
  
  p <- plot_state_map(
    state_year,
    br_states,
    disease_name
  )
  
  ggsave(
    filename = paste0(
      "figures/",
      tolower(ifelse(disease_name == "Chikungunya", "chikv", disease_name)),
      "_geo.png"    ),
    plot = p,
    width = 12,
    height = 8,
    units = "in",
    dpi = 600
  )
  
  cat("Saved:", disease_name, "\n")
}

# cumulative incidence

state_cumulative <- state_year %>%
  group_by(state, disease) %>%
  summarise(
    cumulative_incidence = sum(incidence, na.rm = TRUE),
    .groups = "drop"
  )

glimpse(state_cumulative)

plot_cumulative_state_map <- function(data, geography, disease_name) {
  
  map_data <- geography %>%
    left_join(
      data %>%
        filter(disease == disease_name),
      by = c("abbrev_state" = "state")
    )
  
  ggplot(map_data) +
    geom_sf(
      aes(fill = cumulative_incidence),
      color = "black",
      linewidth = 0.1
    ) +
    
    # state labels that fit
    geom_sf_text(
      data = map_data %>%
        filter(
          !abbrev_state %in%
            c("RN", "PB", "PE", "AL", "SE", "ES", "RJ", "SC")
        ),
      aes(label = abbrev_state),
      size = 1.5
    ) +
    
    # manually positioned labels
    geom_text(
      data = state_labels,
      aes(x = x, y = y, label = state),
      size = 1.5
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
        " Cumulative Incidence by State, 2016–2025"
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
dengue_cumulative <- plot_cumulative_state_map(
  state_cumulative,
  br_states,
  "Dengue"
)

zika_cumulative <- plot_cumulative_state_map(
  state_cumulative,
  br_states,
  "Zika"
)

chikv_cumulative <- plot_cumulative_state_map(
  state_cumulative,
  br_states,
  "Chikungunya"
)

ggsave(
  "figures/dengue_cumulative_state_map.png",
  plot = dengue_cumulative,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/zika_cumulative_state_map.png",
  plot = zika_cumulative,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)

ggsave(
  "figures/chikv_cumulative_state_map.png",
  plot = chikv_cumulative,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600
)
