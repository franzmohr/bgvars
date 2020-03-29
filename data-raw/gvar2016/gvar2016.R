# Downloaded Feb 19, 2018.
# https://sites.google.com/site/gvarmodelling/gvar-toolbox/download

rm(list = ls())

library(dplyr)
library(lubridate)
library(readxl)
#library(seasonal)
library(tidyr)
library(zoo)

####  PPP Data ####
PPP <- as.data.frame(read_xls("data-raw/gvar2016/PPP-GDP WDI (1990-2016).xls", sheet = "WDI"))
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

#### Country data ####
variables <- c("y", "Dp", "eq", "ep", "r", "lr")
countries <- as.data.frame(read_xls("data-raw/gvar2016/Country Codes.xls", col_types = "text"))

data <- c()
nam <- c()
for (i in countries[, "Country Name"]) {
  temp <- as.data.frame(read_xls("data-raw/gvar2016/Country Data (1979Q2-2016Q4).xls", sheet = i))
  dates <- as.character(zoo::as.yearqtr(temp[, "date"]))
  temp <- ts(temp[, which(names(temp)%in%variables)], start = as.numeric(zoo::as.yearqtr(temp[, "date"])[1]), frequency = 4)
  dimnames(temp)[[1]] <- dates
  #temp <- zoo::zoo(temp[, which(names(temp)%in%variables)], order.by = zoo::as.yearqtr(temp[, "date"]))
  data <- c(data, list(temp))
  nam <- c(nam, i)
  rm(temp)
}

capcase <- function(x) { # from ?tolower
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

nam <- unlist(lapply(tolower(nam), capcase))
countries <- cbind(countries, "Name" = nam)
countries$Name <- as.character(countries$Name)
countries$Name[which(countries$Name == "United Kingdom")] <- "UK"
countries$Name[which(countries$Name == "Usa")] <- "USA"
names(data) <- countries$Name
country.data <- data

#### Add credit data to country data ####
#download.file("https://www.bis.org/statistics/full_bis_total_credit_csv.zip", destfile = "data-raw/gvar2016/bis.zip")
#download.file("https://www.bis.org/statistics/full_webstats_long_cpi_dataflow_csv.zip", destfile = "data-raw/gvar2016/cpi.zip")
#unzip("data-raw/gvar2016/bis.zip", exdir = "data-raw/gvar2016")
#unzip("data-raw/gvar2016/cpi.zip", exdir = "data-raw/gvar2016")

#bis <- read.csv("data-raw/gvar2016/WEBSTATS_TOTAL_CREDIT_DATAFLOW_csv_col.csv", stringsAsFactors = FALSE) %>%
#  filter(FREQ == "Q",
#         TC_BORROWERS == "P",
#         TC_LENDERS == "A",
#         UNIT_TYPE == "XDC",
#         TC_ADJUST == "A",
#         VALUATION == "M") %>%
#  gather(key = "date", value = "value", -(FREQ:Time.Period)) %>%
#  mutate(date = as.yearqtr(date, format = "X%Y.Q%q")) %>%
#  rename(cty = Borrowers..country) %>%
#  select(date, cty, value) %>%
#  filter(!is.na(value)) %>%
#  rename(credit = value)

#cpi <- read.csv("data-raw/gvar2016/WEBSTATS_LONG_CPI_DATAFLOW_csv_col.csv", stringsAsFactors = FALSE) %>%
#  filter(FREQ == "M",
#         UNIT_MEASURE == 628) %>%
#  gather(key = "date", value = "value", -(FREQ:Time.Period)) %>%
#  rename(cty = Reference.area) %>%
#  select(date, cty, value) %>%
#  mutate(month = as.yearmon(date, format = "X%Y.%m"),
#         date = as.yearqtr(month)) %>%
#  filter(!is.na(value)) %>%
#  group_by(date, cty) %>%
#  filter(n() == 3) %>%
#  summarise(cpi = mean(value)) %>%
#  ungroup()


#### SA Adjust ####
#setwd("~/Dropbox/R/bgvars/data-raw/gvar2016")

# credit
#temp <- bis %>%
#  mutate(credit = log(credit)) %>%
#  filter(!is.na(credit),
#         abs(credit) != Inf) %>%
#  select(date, cty, credit) %>%
#  spread(cty, credit) %>%
#  filter(date >= "1950 Q1")

#dates <- unique(temp$date)
#dates_credit <- dates[order(dates)]

#avail <- apply(temp, 2, function(x){length(na.omit(x))}) >= 151

#temp <- temp[, avail] %>%
#  rename(UK = `United Kingdom`,
#         USA = `United States`)

#temp <- temp[,dimnames(temp)[[2]] %in% names(country.data)]
#temp <- ts(temp, start = 1950, frequency = 4)

#for (i in 1:ncol(temp)) {
#  temp_i <- seas(temp[, i], transform.function = "none", na.action = na.exclude)
#  temp[, i] <- final(temp_i)
  #temp_i <- x12(na.omit(temp[, i]))
  #temp[(length(temp[, i]) - length(temp_i@d11) + 1):length(temp[, i]), i] <-  temp_i@d11
#}

#credit <- exp(temp)

# cpi
#temp <- cpi %>%
#  select(date, cty, cpi) %>%
#  mutate(cpi = log(cpi)) %>%
#  filter(abs(cpi) != Inf) %>%
#  spread(cty, cpi) %>%
#  filter(date >= "1970 Q1")

#dates <- unique(temp$date)
#dates_cpi <- dates[order(dates)]

#avail <- apply(temp, 2, function(x){length(na.omit(x))}) >= 151

#temp <- temp[, avail] %>%
#  rename(UK = `United Kingdom`,
#         USA = `United States`)

#temp <- temp[,dimnames(temp)[[2]] %in% names(country.data)]
#temp <- ts(temp, start = 1950, frequency = 4)

#for (i in 1:ncol(temp)) {
#  temp_i <- seas(temp[, i], transform.function = "none", na.action = na.exclude)
#  temp[, i] <- final(temp_i)
#  #temp_i <- x12(na.omit(temp[, i]))
#  #temp[(length(temp[, i]) - length(temp_i@d11) + 1):length(temp[, i]), i] <-  temp_i@d11
#}

#cpi <- exp(temp)

#setwd("/home/franz/Dropbox/R/bgvars")

#credit <- as.data.frame(credit) %>%
#  mutate(date = dates_credit) %>%
#  gather(key = "cty", value = "credit", -date)

#cpi <- as.data.frame(cpi) %>%
#  mutate(date = dates_cpi) %>%
#  gather(key = "cty", value = "cpi", -date)

#### Add credit series ####s
#bis <- full_join(credit, cpi, by = c("date", "cty")) %>%
#  filter(date >= "1979 Q2",
#         date <= "2016 Q4") %>%
#  mutate(rcredit = credit / cpi) %>%
#  filter(!is.na(rcredit),
#         abs(rcredit) != Inf) %>%
#  arrange(date) %>%
#  group_by(cty) %>%
#  mutate(index = mean(rcredit[substring(date, 1, 4) == "2010"])) %>%
#  ungroup() %>%
#  mutate(credit = log(rcredit / index)) %>%
#  select(date, cty, credit) %>%
#  spread(cty, credit)
  
#avail <- apply(bis, 2, function(x){length(na.omit(x))}) == nrow(bis)

#bis <- bis[, avail]
#bis <- ts(bis[, -1], start = 1979.25, frequency = 4)

#for (i in dimnames(bis)[[2]]) {
#  if (i %in% names(country.data)) {
#    dim_nam <- dimnames(country.data[[i]])
#    dim_nam[[2]] <- c(dim_nam[[2]], "cred")
#    country.data[[i]] <- cbind(country.data[[i]], bis[,i])
#    dimnames(country.data[[i]]) <- dim_nam
#  } else {
#    next
#  }
#}

#### Global data ####
global.variables <- c("poil", "pmat", "pmetal")
global.data <- as.data.frame(read_xls("data-raw/gvar2016/Country Data (1979Q2-2016Q4).xls", sheet = 1))
dates <- as.character(zoo::as.yearqtr(global.data[, "date"]))
global.data <- ts(global.data[, which(names(global.data)%in%global.variables)], start = as.numeric(zoo::as.yearqtr(global.data[, "date"])[1]), frequency = 4)
dimnames(global.data)[[1]] <- dates

#global.data <- zoo::zoo(global.data[, which(names(global.data)%in%global.variables)], order.by = zoo::as.yearqtr(global.data[, "date"]))

#### Weight matrix ####
#weights <- array(NA, dim = c(length(country.data), length(country.data), 37))
#dimnames(weights) <- list(countries[, "Country Code"], countries[, "Country Code"], as.character(1980:2016))
weights <- c()
w.names <- c()
for (i in countries[, "Country Code"]) {
  temp <- as.data.frame(read_xlsx("data-raw/gvar2016/Trade Flows (1980-2016).xlsx",
                                  sheet = i, na = "NaN"))
  #dimnames(temp)[[1]] <- temp[, 1]
  dates <- temp[, 1]
  temp <- temp[, -(1:2)]
  # Reorder columns
  temp <- temp[ , countries[, "Country Code"]]
  # Set home country to zero
  temp[, i] <- 0
  temp <- ts(temp, start = dates[1], frequency = 1)
  dimnames(temp) <- list(as.character(dates), countries[, "Name"])
  weights <- c(weights, list(temp))
  w.names <- c(w.names, countries[which(countries["Country Code"] == i), "Name"])
  #weights[i , ,] <- t(as.matrix(temp))
}
names(weights) <- w.names
#for (i in 1:dim(weights)[3]) {
#  weights[,, i] <- t(apply(weights[,, i], 1, function(x){s <- sum(x); return(x/s)}))
#}
#dimnames(weights)[[1]] <- countries$Name
#dimnames(weights)[[2]] <- countries$Name
weight.data <- weights

#### Rename countries ####
ctr <- c("AR", "AU", "AT", "BE", "BR", "CA", "CN", "CL", "FI", "FR", "DE",
         "IN", "ID", "IT", "JP", "KR", "MY", "MX", "NL", "NO", "NZ", "PE",
         "PH", "ZA", "SA", "SG", "ES", "SE", "CH", "TH", "TR", "GB", "US")

names(country.data) <- ctr
dimnames(PPP)[[2]] <- ctr

names(weight.data) <- ctr
for (i in ctr) {
  dimnames(weight.data[[i]])[[2]] <- ctr
}

#### Save result ####
gvar2016 <- list("country_data" = country.data,
                 "global_data" = global.data,
                 "region_weights" = PPP,
                 "weight_data" = weight.data)
save(gvar2016, file = "data/gvar2016.rda", version = 2)

# country_data <- country.data
# global_data <- global.data
# region_weights <- PPP
# weight_data <- weight.data
# 
# save(country_data, file = "data/country_data.rda")
# save(global_data, file = "data/global_data.rda")
# save(region_weights, file = "data/region_weights.rda")
# save(weight_data, file = "data/weight_data.rda")

