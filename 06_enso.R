# 06: ocean (ENSO)
# goal: plot ENSO as monthly bars colored by sign and strength
#         - monthly mean of weekly values
#         - event thresholds at +/-0.5
#         - 180 x 72 mm
# input: data/ocean_climate_oscillations_update_2025.csv (columns: date, enso, iod, pdo; starts 1993)
# output: figures/ocean/enso_timeseries.pdf

rm(list = ls())

dir.create("figures/ocean", recursive = TRUE, showWarnings = FALSE)

library(readr)
library(dplyr)
library(ggplot2)

ocean <- read_csv("data/ocean_climate_oscillations_update_2025.csv", show_col_types = FALSE)

range(ocean$date)
colSums(is.na(ocean))

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

# diverging, blue = La Nina, red = El Nino
enso_cols <- c("#2B5C8A", "#6FA3C7", "#D6E7F0", "#F7D9C8", "#E08B6A", "#A63D2F")

# conventional ENSO event thresholds
enso_thresholds <- c(-0.5, 0.5)

# keep one month before the study period so lagged effects can be computed later
plot_start <- as.Date("2015-12-01")
plot_end   <- as.Date("2025-12-01")

# ================================================================================
#            2. helper functions
# ================================================================================

# monthly ENSO bars colored on a diverging scale, thresholds dashed
plot_enso <- function(monthly, title, ylab, cols = enso_cols) {
  
  lim <- ceiling(max(abs(monthly$value), na.rm = TRUE) * 2) / 2   # symmetric, rounded to 0.5
  
  ggplot(monthly, aes(x = year_month)) +
    geom_col(
      aes(y = value, fill = value),
      width = 25
    ) +
    geom_hline(
      yintercept = 0,
      color = "grey30",
      linewidth = 0.2
    ) +
    geom_hline(
      yintercept = enso_thresholds,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.2
    ) +
    scale_fill_gradientn(
      colors = cols,
      limits = c(-lim, lim),
      name = ylab
    ) +
    scale_x_date(
      breaks = seq(as.Date("2016-01-01"), max(monthly$year_month), by = "6 months"),
      labels = function(x) {
        ifelse(format(x, "%m") == "01", format(x, "%Y"), "")
      },
      limits = c(min(monthly$year_month) - 15, max(monthly$year_month) + 16),
      expand = c(0, 0)
    ) +
    labs(title = title,
         x = NULL,
         y = NULL) +
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
      legend.box.spacing = unit(1, "mm"),
      plot.margin = margin(1, 1, 1, 1, unit = "mm")
    ) +
    
    guides(
      fill = guide_colorbar(
        barheight = unit(16, "mm"),
        barwidth = unit(4, "mm"),
        ticks = TRUE
      )
    )
}

# ================================================================================
#            3. ENSO data prep (weekly -> monthly mean)
# ================================================================================

# full series kept (1993 onward) so lags can be computed later
enso_full <- ocean %>%
  mutate(date = as.Date(date),
         year_month = as.Date(format(date, "%Y-%m-01"))) %>%
  group_by(year_month) %>%
  summarise(value = mean(enso, na.rm = TRUE), .groups = "drop") %>%
  arrange(year_month)

# display window
enso_monthly <- enso_full %>%
  filter(year_month >= plot_start, year_month <= plot_end)

# sanity check: Dec 2015 up to the last month in the file (115 rows if it ends Jun 2025)
nrow(enso_monthly)
range(enso_monthly$year_month)

# ================================================================================
#            4. ENSO plots
# ================================================================================

# A. monthly bars
plot_enso(enso_monthly, "Monthly ENSO Index", "ENSO index")

ggsave("figures/ocean/enso_timeseries.pdf",
       width = 175,
       height = 66,
       units = "mm",
       dpi = 600,
       bg = "transparent")

