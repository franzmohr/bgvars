# Downloaded January, 2024.
# https://www.mohaddes.org/gvar
# https://data.mendeley.com/datasets/kfp5fhgkvf/1

rm(list = ls())

library(dplyr)
library(lubridate)
library(readxl)
#library(seasonal)
library(tidyr)
library(zoo)

# Country table ----
countries <- as.data.frame(read_xls("data-raw/gvar2023/Country Codes.xls", col_types = "text")) %>%
  mutate(iso = c("AR", "AU", "AT", "BE", "BR", "CA", "CN", "CL", "FI", "FR", "DE",
                 "IN", "ID", "IT", "JP", "KR", "MY", "MX", "NL", "NO", "NZ", "PE",
                 "PH", "ZA", "SA", "SG", "ES", "SE", "CH", "TH", "TR", "GB", "US"))

#  PPP Data ----
PPP <- as.data.frame(read_xls("data-raw/gvar2023/PPP-GDP WDI (1990-2018).xls", sheet = "WDI"))
PPP <- na.omit(PPP)[, -which(names(PPP) %in% c("Country Code", "Indicator Name", "Indicator Code"))]
nam <- PPP[, "Country Name"]
time <- names(PPP)[-1]
PPP <- data.frame(time, t(PPP[,-1]))
nam[which(nam == "Korea, Rep.")] <- "Korea"
nam[which(nam == "United Kingdom")] <- "UK"
nam[which(nam == "United States")] <- "USA"
names(PPP) <- c("Year", nam)
rownames(PPP) <- NULL
dates <- as.character(PPP[, "Year"])
PPP <- ts(PPP[, -1], start = 1990, frequency = 1)
dimnames(PPP)[[1]] <- dates
dimnames(PPP)[[2]] <- countries$iso

# Country data ----
variables <- c("y", "Dp", "eq", "ep", "r", "lr")

country_data <- NULL
nam <- c()

for (i in variables) {
  temp <- read_xls("data-raw/gvar2023/GVAR_2023Q3.xls", sheet = i) %>%
    pivot_longer(cols = -c("date"), names_to = "ctry") %>%
    mutate(var = i) %>%
    group_by(ctry) %>%
    mutate(cond = !all(value == 123456789)) %>%
    ungroup() %>%
    filter(cond) %>%
    select(-cond)
  country_data <- bind_rows(country_data, temp)
}

country_data <- country_data %>%
  left_join(countries, by = c("ctry" = "Country Short Name")) %>%
  select(date, ctry = iso, var, value)

ctry_names <- unique(pull(country_data, "ctry"))
ctry_names <- ctry_names[order(ctry_names)]

final_data <- NULL
for (i in ctry_names) {
  temp <- country_data %>%
    filter(ctry == i) %>%
    pivot_wider(names_from = "var", values_from = "value") %>%
    arrange(date)
  dates <- as.character(pull(temp, "date"))
  temp <- ts(temp[, which(names(temp) %in% variables)], start = as.numeric(zoo::as.yearqtr(pull(temp, "date"))[1]), frequency = 4)
  dimnames(temp)[[1]] <- dates
  
  final_data <- c(final_data, list(temp))
}

names(final_data) <- ctry_names
country_data <- final_data; rm(final_data)


# Global data ----
global_variables <- c("poil", "pmat", "pmetal")
global_data <- NULL
for (i in global_variables) {
  temp <- read_xls("data-raw/gvar2023/GVAR_2023Q3.xls", sheet = i) %>%
    pivot_longer(cols = -c("date"), names_to = "ctry") %>%
    mutate(var = i)
  global_data <- bind_rows(global_data, temp)
}
global_data <- global_data %>%
  select(date, var, value) %>%
  mutate(date = zoo::as.yearqtr(date)) %>%
  pivot_wider(names_from = "var", values_from = "value") %>%
  arrange(date)
dates <- pull(global_data, "date")

global_data <- ts(global_data[, which(names(global_data) %in% global_variables)],
                  start = c(1979, 2), frequency = 4)
dimnames(global_data)[[1]] <- dates

#### Weight matrix ####
#weights <- array(NA, dim = c(length(country.data), length(country.data), 37))
#dimnames(weights) <- list(countries[, "Country Code"], countries[, "Country Code"], as.character(1980:2016))
weights <- c()
w.names <- c()
for (i in countries[, "Country Code"]) {
  temp <- as.data.frame(read_excel("data-raw/gvar2023/flows.xls",
                                   sheet = i, na = "NaN"))
  # Extract years
  dates <- temp[, 1]
  # Drop year and home country information
  temp <- temp[, -(1:2)]
  # Reorder columns
  temp <- temp[ , countries[, "Country Code"]]
  # Set home country values to zero
  temp[, i] <- 0
  # Replace country code with ISO code
  dimnames(temp)[[2]] <- countries[, "iso"]
  # Generate ts object
  temp <- ts(temp, start = dates[1], frequency = 1)
  dimnames(temp)[[1]] <- as.character(dates)
  weights <- c(weights, list(temp))
  w.names <- c(w.names, countries[which(countries["Country Code"] == i), "iso"])
  #weights[i , ,] <- t(as.matrix(temp))
}
names(weights) <- w.names
weight_data <- weights

submodel_data <- NULL
for (i in ctry_names) {
  submodel_data[[i]] <- list("endogen" = country_data[[i]],
                             "weights" = weights[[i]])
}
class(submodel_data) <- list("submodeldata", "list")

# Save result ----
gvar2023 <- list("submodel_data" = submodel_data,
                 "global_data" = global_data,
                 "region_weights" = PPP)

save(gvar2023, file = "data/gvar2023.rda", version = 2)
