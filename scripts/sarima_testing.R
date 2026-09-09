rm(list=ls())

library(readr)
library(forecast)
library(dplyr)

combined <- read_csv("data/Brazil_arbovirus_monthly_data_2016_2025.csv")

sub <- combined %>%
  filter(muni == 355030) %>%
  arrange(year, month)
  
nrow(sub)
sum(is.na(sub$dengueIncidence))
summary(sub$dengueIncidence)
var(sub$dengueIncidence, na.rm = TRUE)
  
ts_m <- ts(sub$dengueIncidence,
           start = c(min(sub$year), min(sub$month)),
           frequency = 12) # turn into monthly time series

sum(is.na(ts_m)) # count NA months that we want to fill

fit <- auto.arima(ts_m, seasonal = TRUE) # apply arima model

summary(fit)

checkresiduals(fit) # check if residuals look like white noise... not too sure how to interpret

# forecast
fc <- forecast(fit, h=36)

plot(fc)

ts_filled <- na.interp(ts_m) # interpolate missing data
sum(is.na(ts_filled)) # number of NA months --> 0

data.frame(
  original = as.numeric(ts_m),
  interpolated = as.numeric(ts_filled)
)

plot(ts_m,
     main = "Dengue incidence with missing values for muni=355030")

lines(ts_filled,
      lty = 2)
