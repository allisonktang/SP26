# incidence heatmap
# goal: municipality level incidence heatmap (ungrouped), raw + log10 scale
# input: "Brazil_arbovirus_monthly_data_2016_2025.csv"
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

merged <- read_csv("Brazil_arbovirus_monthly_data_2016_2025.csv")

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
                     chikvIncidence = "Chikungunya"),
    disease = factor(disease, levels = c("Dengue", "Zika", "Chikungunya"))
  )

base_theme <- theme_classic() + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    axis.text.x = element_text(angle = 90, hjust = 1) + 
    plot.title = element_text(
      size = 16,
      hjust = 0.5,          
      margin = margin(b = 10))
  )

# raw/linear scale

raw_plot <- ggplot(mergedLong_muni, aes(x = year_month, y = factor(muni), fill = incidence)) +
  geom_tile(color = NA, linewidth = 0) +
  facet_wrap(~ disease, ncol = 1) +
  scale_fill_gradient(
    low="grey90",
    high = "#D4180A",
    na.value = "white") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  labs(x = "Time", y = "Municipality", fill = "Incidence", title = "Municipality Incidence (linear scale)") + 
  base_theme 


# log10 scale
log_plot <-  ggplot(mergedLong_muni, aes(x = year_month, y = factor(muni), fill = incidence)) +
   geom_tile(color = NA, linewidth = 0) +
   facet_wrap(~ disease, ncol = 1) +
   scale_fill_gradient(
     low="grey90",
     high = "#D4180A",
     na.value = "white",
     trans = "log10",
     breaks = scales::trans_breaks("log10", function(x) 10^x),
     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
   scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
   labs(x = "Time", y = "Municipality", fill = "Incidence", title = "Municipality Incidence (log scale)") + 
   base_theme 