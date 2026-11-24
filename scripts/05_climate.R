# 05: climate
# goal: plot temp and precip as heatmap and time series
#         - min, moving average, max 
#         - national, state, muni (ungrouped & grouped by state)
#         - monthly resolution
# input: data/combined_all_climate_data_2016_2024(7).rds, data/climate_daily_2025.csv

rm(list = ls())

dir.create("figures/climate", recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("cmocean", quietly = TRUE)) install.packages("cmocean")

library(cmocean)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(viridis)
library(zoo)

climate24 <- readRDS("data/combined_all_climate_data_2016_2024(7).rds")
climate25 <- read_csv("data/climate_daily_2025.csv") %>% select(-1, -last_col())

climate <- bind_rows(climate24, climate25) 
  

# ================================================================================
#            1. look ups 
# ================================================================================
base_theme <- theme_classic() +
  theme(
    axis.line = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_text(size = 4.5),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 5),
    plot.title = element_blank()
  )

strip_theme <- theme(
  strip.background = element_rect(fill = "grey90", color = NA),
  strip.text.y.right = element_text(angle = 0, size = 8),
  panel.spacing.y = unit(0.15, "lines")
)

uf_lookup <- data.frame(
  uf_code = c("11","12","13","14","15","16","17",
              "21","22","23","24","25","26","27","28","29",
              "31","32","33","35","41","42","43","50","51","52","53"),
  state   = c("RO","AC","AM","RR","PA","AP","TO",
              "MA","PI","CE","RN","PB","PE","AL","SE","BA",
              "MG","ES","RJ","SP","PR","SC","RS","MS","MT","GO","DF")
)

state_order <- c(
  "AC", "AM", "PA", "RR", "RO", "AP", "TO",
  "PI", "BA", "MA", "PE", "CE", "AL", "SE", "RN", "PB",
  "MT", "GO", "MS", "DF",
  "MG", "ES", "RJ", "SP",
  "PR", "SC", "RS"
)

region_lookup <- data.frame(
  state  = state_order,
  region = c(rep("North", 7), rep("Northeast", 9), rep("Center-West", 4),
             rep("Southeast", 4), rep("South", 3))
)

# ================================================================================
#            2. helper functions 
# ================================================================================

# k = num months for moving average, partial = TRUE shrinks window at edge
ma <- function(x, k) rollapply(x, width = k, FUN = mean, partial = TRUE,
                               align = "center", na.rm = TRUE)

# municipality-month table for one variable
# col = daily column; fun = mean for state vars, sum for totals
build_muni_month <- function(col, fun = mean) {
  climate %>%
    mutate(date = as.Date(date),
           year_month = as.Date(format(date, "%Y-%m-01")),
           uf_code = substr(as.character(geocodigo), 1, 2)) %>%
    filter(date >= as.Date("2016-01-01"), date <= as.Date("2025-12-31")) %>%
    group_by(geocodigo, uf_code, year_month) %>%
    summarise(value = fun(.data[[col]], na.rm = TRUE), .groups = "drop") %>%
    left_join(uf_lookup, by = "uf_code")
}

# national time series: percentile band across municipalities + moving averages
plot_national <- function(df, title, ylab, cols = c("#24476B", "#668EAC", "#D8D3B9", "#D99A65", "#A94F35", "#633B30"), windows = c(3, 6, 12), band = c(0.05, 0.95)) {
  monthly <- df %>%
    group_by(year_month) %>%
    summarise(lo   = quantile(value, band[1], na.rm = TRUE),
              mean = mean(value, na.rm = TRUE),
              hi   = quantile(value, band[2], na.rm = TRUE),
              .groups = "drop") %>%
    arrange(year_month)
  
  ma_long <- bind_rows(lapply(windows, function(k)
    monthly %>% transmute(year_month, window = paste0(k, "-month"), ma = ma(mean, k)))) %>%
    mutate(window = factor(window, levels = paste0(windows, "-month")))
  
  ggplot(monthly, aes(x = year_month)) +
    geom_linerange(
      aes(ymin = lo, ymax = hi),
      color = "grey80",
      linewidth = 0.1
      ) +
    geom_line(
      data = ma_long,
      aes(y = ma, linetype = window, linewidth = window),
      ) +
    scale_linewidth_manual(
      values = c("3-month" = 0.3,
                 "6-month" = 0.3,
                 "12-month" = 0.3)
      ) +
    scale_linetype_manual(
      values = c("3-month" = "solid",
                 "6-month" = "dashed",
                 "12-month" = "twodash")
      ) +
    scale_color_gradientn(
      colors = cols,
      name = ylab
    ) +
    scale_x_date(
      date_breaks = "6 months",
      labels = function(x) {
        ifelse(format(x, "%m") == "01", format(x, "%Y"), "")
        },
      expand = c(0, 0)
      ) +
    labs(title = title, 
         x = NULL, 
         y = NULL,
         linetype = "Moving average",
         ) +
    base_theme +
    
    theme(
      axis.line = element_line(linewidth = 0.2),
      axis.ticks.y = element_line(linewidth = 0.2),
      axis.ticks.x = element_line(linewidth = 0.2),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 4.5),
      legend.spacing.y = unit(0, "mm"),
      legend.key.height = unit(3, "mm"),
      legend.margin = margin(4, 0, 4, 4),
      legend.box.spacing = unit(1, "mm"),
      plot.margin = margin(1, 1, 1, 1, unit = "mm")
    ) +
    
    guides(
      color = guide_colorbar(
        order = 1,
        barheight = unit(16, "mm"),
        barwidth = unit(4, "mm"),
        ticks = TRUE
      ),
      linetype = guide_legend(order = 2),
      linewidth = "none"
    )
}

# state x month heatmap, grouped by region
plot_state_heatmap <- function(df, title, fill_label, fill_scale, anomaly = FALSE) {
  d <- df %>%
    filter(state %in% state_order) %>%
    group_by(state, year_month) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(state) %>%
    mutate(value = if (anomaly) value - mean(value, na.rm = TRUE) else value) %>%
    ungroup() %>%
    left_join(region_lookup, by = "state") %>%
    mutate(region = factor(region, levels = unique(region_lookup$region)),
           state  = factor(state, levels = rev(state_order)))
  
  ggplot(d, aes(x = year_month, y = state, fill = value)) +
    geom_tile() +
    facet_grid(region ~ ., scales = "free_y", space = "free_y") +
    scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m", expand = c(0, 0)) +
    fill_scale +
    labs(title = title, x = "Year-Month", y = "State", fill = fill_label) +
    base_theme +
    theme(strip.background = element_rect(fill = "grey90", color = NA),
          strip.text.y.right = element_text(angle = 0, size = 8),
          panel.spacing.y = unit(0.15, "lines"))
}

# muni x month, ordered by state order and grouped by region
plot_muni_heatmap <- function(df, title, fill_label, fill_scale, standardize = "z") {
  d <- df %>%
    filter(state %in% state_order) %>%
    mutate(month = as.integer(format(year_month, "%m"))) %>%
    group_by(geocodigo, month) %>%
    mutate(anom = if (standardize == "z") (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)
           else value - mean(value, na.rm = TRUE)) %>%
    ungroup() %>%
    left_join(region_lookup, by = "state") %>%
    mutate(region = factor(region, levels = unique(region_lookup$region)))
  
  muni_levels <- d %>%
    distinct(geocodigo, state) %>%
    arrange(match(state, state_order), geocodigo) %>%
    pull(geocodigo)
  
  d <- d %>% mutate(geocodigo = factor(geocodigo, levels = rev(muni_levels)))
  
  lim <- quantile(abs(d$anom), 0.99, na.rm = TRUE)   # clip extremes so colors stay readable
  
  ggplot(d, aes(x = year_month, y = geocodigo, fill = anom)) +
    geom_tile() +
    facet_grid(region ~ ., scales = "free_y", space = "free_y") +
    scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m", expand = c(0, 0)) +
    fill_scale(limits = c(-lim, lim), oob = scales::squish) +
    labs(title = title, x = "Year-Month", y = "Municipality (ordered by state)",
         fill = fill_label) +
    base_theme + 
    strip_theme +
    theme(axis.text.y = element_blank())
}

# ================================================================================
#            3. temp data prep (mean low/high temp instead of percentile)
# ================================================================================

muni_month <- build_muni_month("temp_med") %>%
  rename(temp_avg = value) %>%
  inner_join(build_muni_month("temp_min") %>% select(geocodigo, year_month, tmin = value),
             by = c("geocodigo", "year_month")) %>%
  inner_join(build_muni_month("temp_max") %>% select(geocodigo, year_month, tmax = value),
             by = c("geocodigo", "year_month"))

temp_monthly <- muni_month %>%
  group_by(year_month) %>%
  summarise(lo   = mean(tmin, na.rm = TRUE),
            mean = mean(temp_avg, na.rm = TRUE),
            hi   = mean(tmax, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(year_month) %>%
  mutate(`3-month`  = ma(mean, 3),
         `6-month`  = ma(mean, 6),
         `12-month` = ma(mean, 12))

temp_ma_long <- temp_monthly %>%
  pivot_longer(c(`3-month`, `6-month`, `12-month`),
               names_to = "window", values_to = "ma") %>%
  mutate(window = factor(window, levels = c("3-month", "6-month", "12-month")))

# ================================================================================
#            4. temp plots
# ================================================================================

# A1. national time series, overlapping

temp_cols <- c(
  "#24476B", 
  "#668EAC", 
  "#D8D3B9",
  "#D99A65", 
  "#A94F35", 
  "#633B30"
)

ggplot(temp_monthly, aes(x = year_month)) +
  geom_linerange(
    aes(ymin = lo, ymax = hi),
    color = "grey80",
    linewidth = 0.1
   ) +
  geom_line(
    data = temp_ma_long,
    aes(y = ma, color = ma, linetype = window, linewidth = window),
   ) +
  scale_linewidth_manual(
    values = c("3-month" = 0.3,
               "6-month" = 0.3,
               "12-month" = 0.3)
    ) +
  scale_linetype_manual(
    values = c("3-month" = "solid",
               "6-month" = "dashed",
               "12-month" = "twodash")
    ) +
  scale_color_gradientn(
    colors = temp_cols,
    name = "Mean temperature\n(°C)"
    ) +
  scale_x_date(
    date_breaks = "6 months",
    labels = function(x) {
      ifelse(format(x, "%m") == "01", format(x, "%Y"), "")
    },
    expand = c(0, 0)
    ) +
  
  labs(
    title = "National Monthly Temperature",
    x = NULL,
    y = NULL,
    color = "Mean temperature (°C)",
    linetype = "Moving average"
    ) +
  
  base_theme +
  
  theme(
    axis.line = element_line(linewidth = 0.2),
    axis.ticks.y = element_line(linewidth = 0.2),
    axis.ticks.x = element_line(linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 4.5),
    legend.spacing.y = unit(0, "mm"),
    legend.key.height = unit(3, "mm"),
    legend.margin = margin(4, 0, 4, 4),
    legend.box.spacing = unit(1, "mm"),
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  ) +
  
  guides(
    color = guide_colorbar(
      barheight = unit(16, "mm"),
      barwidth = unit(4, "mm"),
      ticks = TRUE
    ),
    linewidth = "none"
  )

ggsave("figures/climate/temp_national_timeseries.pdf", 
       width = 89, 
       height = 42, 
       units = "mm",
       dpi = 600,
       bg = "transparent")

# A2. national time series, faceted
ggplot(temp_ma_long, aes(x = year_month)) +
  geom_linerange(data = temp_monthly,                     
                 aes(ymin = lo, ymax = hi),
                 color = "#d0e2f2", linewidth = 0.4) +
  geom_line(aes(y = ma), color = "#3c78b0", linewidth = 0.8) +
  facet_wrap(~ window, ncol = 1) +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
  labs(title = "National Monthly Temperature",
       x = NULL, y = "Avg monthly temperature (°C)",
       caption = "Bars: national mean of monthly min to max temperature; line: centered moving average of the national mean") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0))

ggsave("figures/climate/temp_national_timeseries_facets.png", width = 12, height = 8, dpi = 300)

# A3. 6 month ma only
ggplot(filter(temp_ma_long, window == "6-month"), aes(x = year_month)) +
  geom_linerange(data = temp_monthly,
                 aes(ymin = lo, ymax = hi),
                 color = "#d0e2f2", linewidth = 0.4) +
  geom_line(aes(y = ma), color = "#3c78b0", linewidth = 0.9) +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
  labs(title = "National Monthly Temperature",
       x = NULL, y = "Avg monthly temperature (°C)",
       caption = "min to max temperature; line: 6-month centered moving average of the national mean") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.minor = element_blank())

ggsave("figures/climate/temp_national_ma6.png", width = 12, height = 4, dpi = 300)

# B. seasonality: each line is a year
temp_monthly %>%
  mutate(year  = factor(format(year_month, "%Y")),
         month = as.integer(format(year_month, "%m"))) %>%
  ggplot(aes(month, mean, group = year, color = year)) +
  geom_line() +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(title = "Seasonal Cycle by Year",
       x = NULL, 
       y = "Mean temperature (°C)", 
       color = NULL) +
  theme_minimal()

ggsave("figures/climate/temp_seasonal_by_year.png", width = 8, height = 4, dpi = 300)


# C. state-month heatmap
plot_state_heatmap(muni_month %>% mutate(value = temp_avg),
                   "Monthly Average Temperature by State, Grouped by Region",
                   "Avg temp (°C)",
                   scale_fill_distiller(palette = "RdBu", direction = -1))
                   #scale_fill_distiller(palette = "RdYlBu", direction = -1))
                   #scale_fill_cmocean(name = "thermal"))
                   # scale_fill_cmocean(name = "balance"))

ggsave("figures/climate/temp_state_heatmap_thermal.png", width = 12, height = 6, dpi = 300)

# D. muni-month heatmap anomaly
plot_muni_heatmap(muni_month %>% mutate(value = temp_avg),
                   "Monthly Temperature Anomaly by Municipality, Grouped by Region",
                   "z-score",
                  function(...) scale_fill_distiller(palette = "RdBu", direction = -1))
                  # scale_fill_cmocean(name = "thermal"))
                  # scale_fill_cmocean(name = "balance"))

ggsave("figures/climate/temp_muni_heatmap.png", width = 12, height = 8, dpi = 300)
                   

# ================================================================================
#            5. precipitation (monthly total, mm)
# ================================================================================

fun_safe <- function(x, f) if (all(is.na(x))) NA_real_ else f(x, na.rm = TRUE)
summarise(value = fun_safe(.data[[col]], fun), ...) 

precip_muni <- build_muni_month("precip_tot", fun = sum)


precip_cols <- c("#DEF5E5FF",
                 "#60CEACFF",
                 "#3497A9FF",
                 "#395D9CFF",
                 "#382A54FF",
                 "#0B0405FF"
)

# A. National time series
plot_national(precip_muni, "National Monthly Precipitation",
              "Total monthly\nprecipitation (mm)", color = "grey20", precip_cols)

ggsave("figures/climate/precip_national_timeseries.pdf", 
       width = 89, 
       height = 42, 
       units = "mm",
       dpi = 600,
       bg = "transparent")

# B. state-month heatmap
plot_state_heatmap(precip_muni, "Monthly Precipitation by State, Grouped by Region",
                   "Precip (mm)",
                   scale_fill_distiller(palette = "YlGnBu", direction = 1,
                                        transform = "sqrt"))
                   #scale_fill_cmocean(name = "deep", transform = "sqrt"))

ggsave("figures/climate/precip_state_heatmap.png", width = 12, height = 6, dpi = 300)

# C. muni-month heatmap anomaly
plot_muni_heatmap(precip_muni, "Monthly Precipitation Anomaly by Municipality, Grouped by Region",
                  "z-score",
                  function(...) scale_fill_distiller(palette = "BrBG", direction = 1, ...))

ggsave("figures/climate/precip_muni_heatmap.png", width = 12, height = 8, dpi = 300)


# ================================================================================
#            6. humidity (monthly mean of daily median)
# ================================================================================
umid_muni <- build_muni_month("umid_med")
umid_cols <- c("#edf8fb",
               "#bfd3e6",
               "#9ebcda",
               "#8c96c6",
               "#8856a7",
               "#810f7c")

# A. National time series
plot_national(umid_muni, "National Monthly Relative Humidity",
              "Relative humidity (%)", color = "#6a3d9a")

ggsave("figures/climate/umid_national_timeseries.png", width = 12, height = 4, dpi = 300)


plot_national(umid_muni, "National Monthly Relative Humidity",
              "Relative humidity (%)", color = "grey20", umid_cols)

ggsave("figures/climate/umid_national_timeseries.pdf", 
       width = 89, 
       height = 42, 
       units = "mm",
       dpi = 600,
       bg = "transparent")

# B. state-month heatmap
# need to add custom scale, since data concentrated near max
library(scales)
pow2_trans <- trans_new(
  name = "pow2",
  transform = function(x) x^2,
  inverse = function(x) sqrt(x)
)
plot_state_heatmap(umid_muni, "Monthly Humidity by State, Grouped by Region",
                   "Humidity (%)",
                   scale_fill_distiller(palette = "BuPu", direction = 1, trans = pow2_trans))
                   # scale_fill_cmocean(name = "dense") 
ggsave("figures/climate/umid_state_heatmap.png", width = 12, height = 6, dpi = 300)

# C. muni-month heatmap anomaly
plot_muni_heatmap(umid_muni, "Monthly Humidity Anomaly by Municipality, Grouped by Region",
                  "z-score",
                  function(...) scale_fill_distiller(palette = "PuOr", direction = 1, ...))
ggsave("figures/climate/umid_muni_heatmap.png", width = 12, height = 8, dpi = 300)

# ================================================================================
#            7. pressure (state heatmap as anomaly from state mean)
# ================================================================================
pressao_muni <- build_muni_month("pressao_med")
pressao_cols <- c("#FBEAE3", 
                  "#F4B59B", 
                  "#E2688A", 
                  "#B03A8E", 
                  "#6A2A7A", 
                  "#2E1A47")

# A. National time series
plot_national(pressao_muni, "National Monthly Pressure (IQR)", "Pressure (atm)", color = "grey20", pressao_cols, band = c(0.25,0.75))

ggsave("figures/climate/pressao_national_timeseries.pdf", 
       width = 89, 
       height = 42, 
       units = "mm",
       dpi = 600,
       bg = "transparent")

# B. state-month heatmap
plot_state_heatmap(pressao_muni, "Monthly Pressure Anomaly by State, Grouped by Region",
                   "Anomaly",
                   scale_fill_distiller(palette = "RdBu", direction = -1),
                   anomaly = TRUE)
                   # scale_fill_cmocean(name = "balance") +
ggsave("figures/climate/pressao_state_heatmap.png", width = 12, height = 6, dpi = 300)

# C. muni-month heatmap anomaly
plot_muni_heatmap(pressao_muni, "Monthly Pressure Anomaly by Municipality, Grouped by Region",
                  "z-score",
                  function(...) scale_fill_distiller(palette = "RdBu", direction = -1, ...))
ggsave("figures/climate/pressao_muni_heatmap.png", width = 12, height = 8, dpi = 300)

# ================================================================================
#            8. lagged analysis table (municipality x month)
# ================================================================================
# one row per municipality-month with lag 0-3 of each variable
# lag1 on the March row = February's value
all_months <- seq(as.Date("2016-01-01"), as.Date("2025-12-01"), by = "month")

climate_monthly_lagged <- muni_month %>%
  select(geocodigo, year_month, temp_avg) %>%
  left_join(precip_muni  %>% select(geocodigo, year_month, precip_tot   = value),
            by = c("geocodigo", "year_month")) %>%
  left_join(umid_muni    %>% select(geocodigo, year_month, umid_med     = value),
            by = c("geocodigo", "year_month")) %>%
  left_join(pressao_muni %>% select(geocodigo, year_month, pressao_med  = value),
            by = c("geocodigo", "year_month")) %>%
  complete(geocodigo, year_month = all_months) %>%   # insert missing months so lag() is truly 1 month
  arrange(geocodigo, year_month) %>%
  group_by(geocodigo) %>%
  mutate(across(c(temp_avg, precip_tot, umid_med, pressao_med),
                list(lag1 = ~ lag(.x, 1), lag2 = ~ lag(.x, 2), lag3 = ~ lag(.x, 3)))) %>%
  ungroup() %>%
  mutate(uf_code = substr(as.character(geocodigo), 1, 2)) %>%
  left_join(uf_lookup, by = "uf_code") %>%
  select(geocodigo, state, year_month, everything(), -uf_code)

saveRDS(climate_monthly_lagged, "data/climate_monthly_lagged.rds")

# ================================================================================
#            9. geo state maps by year (annual temperature anomaly)
# ================================================================================
library(sf)
library(geobr)   # needs internet the first time

states_sf <- geobr::read_state(year = 2020, showProgress = FALSE,  cache = FALSE) %>%
  select(abbrev_state)

temp_state_year <- muni_month %>%
  filter(state %in% state_order) %>%
  mutate(year = as.integer(format(year_month, "%Y"))) %>%
  group_by(state, year) %>%
  summarise(temp = mean(temp_avg, na.rm = TRUE), .groups = "drop") %>%
  group_by(state) %>%
  mutate(anom = temp - mean(temp)) %>%
  ungroup()

lim <- max(abs(temp_state_year$anom))

states_sf %>%
  left_join(temp_state_year, by = c("abbrev_state" = "state")) %>%
  ggplot() +
  geom_sf(aes(fill = anom), color = "white", linewidth = 0.1) +
  facet_wrap(~ year, nrow = 2) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-lim, lim)) +
  labs(title = "Annual Temperature Anomaly by State (vs. 2016-2025 mean)", fill = "°C") +
  theme_void() +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5))

ggsave("figures/climate/temp_state_map_by_year.png", width = 14, height = 6, dpi = 300)


# ================================================================================
#            10. climate with incidence on top
# ================================================================================
library(patchwork)
library(readr)

x_year <- scale_x_date(breaks = seq(as.Date("2016-01-01"), as.Date("2025-01-01"), by = "1 year"),
                       date_labels = "%Y",
                       limits = c(as.Date("2015-12-01"), as.Date("2026-01-15")))

x_heat <- scale_x_date(breaks = seq(as.Date("2016-01-01"), as.Date("2025-07-01"), by = "6 months"),
                       date_labels = "%Y-%m",
                       limits = c(as.Date("2015-12-17"), as.Date("2025-12-16")),
                       expand = c(0, 0))

# national monthly incidence from 01
inc <- read_csv("data/nBrazil_arbovirus_monthly_data_2016_2025.csv", show_col_types = FALSE)

inc_nat <- inc %>%
  group_by(year_month) %>%
  summarise(Dengue      = sum(dengueCases, na.rm = TRUE) / sum(pop_tot) * 1e5,
            Zika        = sum(zikaCases,   na.rm = TRUE) / sum(pop_tot) * 1e5,
            Chikungunya = sum(chikvCases,  na.rm = TRUE) / sum(pop_tot) * 1e5,
            .groups = "drop") %>%
  pivot_longer(-year_month, names_to = "disease", values_to = "incidence") %>%
  mutate(disease = factor(disease, levels = c("Dengue", "Zika", "Chikungunya")))

inc_cols <- c("Dengue" = "#F8766D", "Zika" = "#7CAE00", "Chikungunya" = "#00BFC4")

# top panel; x_scale must match the plot underneath (x_year or x_heat)
plot_incidence_top <- function(x_scale, diseases = c("Dengue", "Zika", "Chikungunya")) {
  ggplot(filter(inc_nat, disease %in% diseases),
         aes(year_month, incidence, color = disease)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = inc_cols) +
    x_scale +
    scale_y_log10() +
    labs(x = NULL, y = "Incidence\nper 100k", color = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top")
}

# stack incidence over any climate plot
# kind = "ts" for line plots (x_year), "heat" for heatmaps (x_heat)
with_incidence <- function(climate_plot, kind = "ts", diseases = c("Dengue", "Zika", "Chikungunya"),
                           heights = NULL) {
  x_scale <- if (kind == "ts") x_year else x_heat
  if (is.null(heights)) heights <- if (kind == "ts") c(1, 2.5) else c(1, 4)
  plot_incidence_top(x_scale, diseases) / (climate_plot + x_scale) +
    plot_layout(heights = heights)
}

# ---- time series + incidence
p_temp <- ggplot(filter(temp_ma_long, window == "6-month"), aes(x = year_month)) +
  geom_linerange(data = temp_monthly, aes(ymin = lo, ymax = hi),
                 color = "#d0e2f2", linewidth = 0.4) +
  geom_line(aes(y = ma), color = "#3c78b0", linewidth = 0.9) +
  labs(title = "National Monthly Temperature", x = NULL, y = "Avg monthly temperature (°C)",
       caption = "Bars: national mean of monthly min to max temperature; line: 6-month moving average") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

ggsave("figures/climate/temp_with_incidence.png",
       with_incidence(p_temp), width = 12, height = 6, dpi = 300)

ggsave("figures/climate/precip_with_incidence.png",
       with_incidence(plot_national(precip_muni, "National Monthly Precipitation",
                                    "Total monthly precipitation (mm)", color = "#1b7837")),
       width = 12, height = 6, dpi = 300)

ggsave("figures/climate/umid_with_incidence.png",
       with_incidence(plot_national(umid_muni, "National Monthly Relative Humidity",
                                    "Relative humidity (%)", color = "#6a3d9a")),
       width = 12, height = 6, dpi = 300)

# ---- state heatmaps + incidence
ggsave("figures/climate/temp_state_heatmap_with_incidence.png",
       with_incidence(plot_state_heatmap(muni_month %>% mutate(value = temp_avg),
                                         "Monthly Average Temperature by State, Grouped by Region",
                                         "Avg temp (°C)",
                                         scale_fill_distiller(palette = "RdBu", direction = -1)),
                      kind = "heat"),
       width = 12, height = 8, dpi = 300)

ggsave("figures/climate/precip_state_heatmap_with_incidence.png",
       with_incidence(plot_state_heatmap(precip_muni, "Monthly Precipitation by State, Grouped by Region",
                                         "Precip (mm)",
                                         scale_fill_distiller(palette = "YlGnBu", direction = 1,
                                                              transform = "sqrt")),
                      kind = "heat"),
       width = 12, height = 8, dpi = 300)

ggsave("figures/climate/umid_state_heatmap_with_incidence.png",
       with_incidence(plot_state_heatmap(umid_muni, "Monthly Humidity by State, Grouped by Region",
                                         "Humidity (%)",
                                         scale_fill_distiller(palette = "BuPu", direction = 1,
                                                              transform = pow2_trans)),
                      kind = "heat"),
       width = 12, height = 8, dpi = 300)


# ================================================================================
#            10. climate with incidence overlaid
# ================================================================================

library(patchwork)
library(readr)

# ---- x-axis
x_year <- scale_x_date(
  breaks = seq(as.Date("2016-01-01"), as.Date("2025-01-01"), by = "1 year"),
  date_labels = "%Y",
  limits = c(as.Date("2015-12-01"), as.Date("2026-01-15"))
)

x_heat <- scale_x_date(
  breaks = seq(as.Date("2016-01-01"), as.Date("2025-07-01"), by = "6 months"),
  date_labels = "%Y-%m",
  limits = c(as.Date("2015-12-17"), as.Date("2025-12-16")),
  expand = c(0, 0)
)


# ================================================================================
#            10A. national monthly incidence
# ================================================================================

inc <- read_csv(
  "data/nBrazil_arbovirus_monthly_data_2016_2025.csv",
  show_col_types = FALSE
) %>%
  mutate(year_month = as.Date(year_month))

# Calculate national incidence:
# total cases / total population * 100,000
inc_nat <- inc %>%
  group_by(year_month) %>%
  summarise(
    Dengue = sum(dengueCases, na.rm = TRUE) /
      sum(pop_tot, na.rm = TRUE) * 1e5,
    
    Zika = sum(zikaCases, na.rm = TRUE) /
      sum(pop_tot, na.rm = TRUE) * 1e5,
    
    Chikungunya = sum(chikvCases, na.rm = TRUE) /
      sum(pop_tot, na.rm = TRUE) * 1e5,
    
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Dengue, Zika, Chikungunya),
    names_to = "disease",
    values_to = "incidence"
  ) %>%
  mutate(
    disease = factor(
      disease,
      levels = c("Dengue", "Zika", "Chikungunya")
    )
  )

inc_cols <- c(
  "Dengue" = "#F8766D",
  "Zika" = "#7CAE00",
  "Chikungunya" = "#00BFC4"
)


# ================================================================================
#            10B. prepare national climate time series
# ================================================================================

# Precipitation and humidity use the same approach as plot_national():
# monthly mean across municipalities, with 5th-95th percentile range.

make_climate_monthly <- function(df) {
  
  monthly <- df %>%
    group_by(year_month) %>%
    summarise(
      lo = quantile(value, 0.05, na.rm = TRUE),
      mean = mean(value, na.rm = TRUE),
      hi = quantile(value, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(year_month) %>%
    mutate(
      ma6 = ma(mean, 6)
    )
  
  return(monthly)
}

precip_monthly <- make_climate_monthly(precip_muni)

umid_monthly <- make_climate_monthly(umid_muni)


# ================================================================================
#            10C. function to overlay incidence on climate time series
# ================================================================================

plot_climate_incidence <- function(
    climate_monthly,
    title,
    ylab,
    climate_color,
    diseases = c("Dengue", "Zika", "Chikungunya")
) {
  
  # Keep only the selected diseases
  incidence_data <- inc_nat %>%
    filter(disease %in% diseases)
  
  # Range of the climate variable
  climate_range <- range(
    climate_monthly$lo,
    climate_monthly$hi,
    na.rm = TRUE
  )
  
  # Range of incidence across the selected diseases
  incidence_range <- range(
    incidence_data$incidence,
    na.rm = TRUE
  )
  
  # Convert incidence values to the climate variable's scale
  # This allows both variables to be drawn in the same plotting panel.
  scale_factor <- diff(climate_range) / diff(incidence_range)
  
  incidence_data <- incidence_data %>%
    mutate(
      incidence_scaled =
        climate_range[1] +
        (incidence - incidence_range[1]) * scale_factor
    )
  
  ggplot() +
    
    # Climate range
    geom_linerange(
      data = climate_monthly,
      aes(
        x = year_month,
        ymin = lo,
        ymax = hi
      ),
      color = "grey80",
      linewidth = 0.4
    ) +
    
    # Climate 6-month moving average
    geom_line(
      data = climate_monthly,
      aes(
        x = year_month,
        y = ma6
      ),
      color = climate_color,
      linewidth = 0.9
    ) +
    
    # Incidence curves
    geom_line(
      data = incidence_data,
      aes(
        x = year_month,
        y = incidence_scaled,
        color = disease
      ),
      linewidth = 0.7
    ) +
    
    # Left axis: climate
    # Right axis: original incidence values
    scale_y_continuous(
      name = ylab,
      sec.axis = sec_axis(
        ~ (. - climate_range[1]) / scale_factor + incidence_range[1],
        name = "Incidence per 100,000"
      )
    ) +
    
    scale_color_manual(values = inc_cols) +
    
    x_year +
    
    labs(
      title = title,
      x = NULL,
      color = NULL
    ) +
    
    theme_minimal() +
    
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      axis.title.y.left = element_text(color = climate_color),
      axis.title.y.right = element_text(color = "black")
    )
}


# ================================================================================
#            10D. temperature + incidence
# ================================================================================

# Temperature uses the existing temp_monthly data:
# lo = national mean monthly minimum temperature
# mean = national mean monthly average temperature
# hi = national mean monthly maximum temperature
# 6-month = centered 6-month moving average

temp_monthly_inc <- temp_monthly %>%
  mutate(ma6 = `6-month`)

p_temp_inc <- plot_climate_incidence(
  climate_monthly = temp_monthly_inc,
  title = "National Monthly Temperature and Arbovirus Incidence",
  ylab = "Temperature (°C)",
  climate_color = "#3c78b0"
)

ggsave(
  "figures/climate/temp_with_incidence.png",
  p_temp_inc,
  width = 12,
  height = 5,
  dpi = 300
)


# ================================================================================
#            10E. precipitation + incidence
# ================================================================================

p_precip_inc <- plot_climate_incidence(
  climate_monthly = precip_monthly,
  title = "National Monthly Precipitation and Arbovirus Incidence",
  ylab = "Precipitation (mm)",
  climate_color = "#1b7837"
)

ggsave(
  "figures/climate/precip_with_incidence.png",
  p_precip_inc,
  width = 12,
  height = 5,
  dpi = 300
)


# ================================================================================
#            10F. humidity + incidence
# ================================================================================

p_umid_inc <- plot_climate_incidence(
  climate_monthly = umid_monthly,
  title = "National Monthly Relative Humidity and Arbovirus Incidence",
  ylab = "Relative humidity (%)",
  climate_color = "#6a3d9a"
)

ggsave(
  "figures/climate/umid_with_incidence.png",
  p_umid_inc,
  width = 12,
  height = 5,
  dpi = 300
)


# ================================================================================
#            10G. state heatmaps
# ================================================================================

# Incidence is not overlaid on these heatmaps because their y-axis
# represents state rather than a continuous climate measurement.
# Keep the climate heatmaps as separate figures.

# Temperature
p_temp_heat <- plot_state_heatmap(
  muni_month %>% mutate(value = temp_avg),
  "Monthly Average Temperature by State, Grouped by Region",
  "Avg temp (°C)",
  scale_fill_distiller(palette = "RdBu", direction = -1)
)

ggsave(
  "figures/climate/temp_state_heatmap.png",
  p_temp_heat,
  width = 12,
  height = 6,
  dpi = 300
)


# Precipitation
p_precip_heat <- plot_state_heatmap(
  precip_muni,
  "Monthly Precipitation by State, Grouped by Region",
  "Precip (mm)",
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    transform = "sqrt"
  )
)

ggsave(
  "figures/climate/precip_state_heatmap.png",
  p_precip_heat,
  width = 12,
  height = 6,
  dpi = 300
)


# Humidity
p_umid_heat <- plot_state_heatmap(
  umid_muni,
  "Monthly Humidity by State, Grouped by Region",
  "Humidity (%)",
  scale_fill_distiller(
    palette = "BuPu",
    direction = 1,
    transform = pow2_trans
  )
)

ggsave(
  "figures/climate/umid_state_heatmap.png",
  p_umid_heat,
  width = 12,
  height = 6,
  dpi = 300
)

