# Downloaded March 1, 2018.
# http://qed.econ.queensu.ca/jae/2007-v22.1/dees-di_mauro-pesaran-smith/

rm(list = ls())

library(readxl)

####  PPP Data ####
PPP <- as.data.frame(read_xls("data-raw/dees2007/pppgdp9901.XLS", sheet = "Sheet1"))
PPP <- na.omit(PPP)
names(PPP)[1] <- "Country Name"
nam <- PPP[, "Country Name"]
time <- names(PPP)[-1]
PPP <- data.frame(time, t(PPP[,-1]))
nam[which(nam == "Korea, Rep.")] <- "Korea"
nam[which(nam == "United Kingdom")] <- "UK"
nam[which(nam == "United States")] <- "USA"
names(PPP) <- c("Year", nam)
rownames(PPP) <- NULL
dates <- substring(as.character(PPP[, "Year"]), 2, 5)
PPP <- ts(PPP[, -1], start = 1999, frequency = 1)
dimnames(PPP)[[1]] <- dates

# Country data
global.variables <- "poil"
countries <- as.data.frame(read_xls("data-raw/dees2007/countrycodes.xls", sheet = "Sheet1", col_types = "text"))
capcase <- function(x) { # from ?tolower
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}
countries[,1] <- unlist(lapply(tolower(countries[,1]), capcase))
countries[which(countries[,1] == "South Africa"), 1] <- "SAfrica"
countries[which(countries[,1] == "United Kingdom"), 1] <- "UK"
countries[which(countries[,1] == "Usa"), 1] <- "USA"

vars <- c("y", "p", "rs", "rl", "q", "ep")

data <- c()
nam <- c()
for (i in countries[, "Country Names"]) {
  temp <- as.data.frame(read_xls("data-raw/dees2007/data33.xls", sheet = i))
  names(temp) <- tolower(names(temp))
  # Check for constant series of ones (for USA model)
  temp <- temp[, !apply(temp, 2, function(x){length(unique(x)) == 1})]
  # Omit global vars
  if (any(names(temp) %in% global.variables)) {
    temp <- temp[, -which(names(temp)%in%global.variables)]
  }
  # Rename if ex-rate is an index
  if (any(names(temp) %in% "eindex")) {names(temp)[which(names(temp) == "eindex")] <- "e"}
  
  ts_start <- temp[1, "date"]
  ts_start <- as.numeric(c(substring(ts_start, 1, 4), substring(ts_start, 6, 6)))
  temp <- ts(temp[, -which(names(temp) == "date")], start = ts_start, frequency = 4)
  
  temp_raw <- temp
  names_raw <- dimnames(temp)[[2]]
  
  temp_dees <- NULL
  names_dees <- NULL
  if ("y" %in% names_raw) {
    temp_dees <- cbind(temp_dees, log(temp_raw[, "y"]))
    names_dees <- c(names_dees, "y")
  }
  if ("cpi" %in% names_raw) {
    temp_dees <- cbind(temp_dees, log(temp_raw[, "cpi"]))
    names_dees <- c(names_dees, "p")
  }
  if ("eq" %in% names_raw & "cpi" %in% names_raw) {
    temp_dees <- cbind(temp_dees, log(temp_raw[, "eq"] / temp_raw[, "cpi"]))
    names_dees <- c(names_dees, "q")
  }
  if ("e" %in% names_raw & "cpi" %in% names_raw) {
    temp_dees <- cbind(temp_dees, log(temp_raw[, "e"]) - log(temp_raw[, "cpi"]))
    names_dees <- c(names_dees, "ep")
  }
  if ("rs" %in% names_raw) {
    temp_dees <- cbind(temp_dees, .25 * log(1 + temp_raw[, "rs"] / 100))
    names_dees <- c(names_dees, "rs")
  }
  if ("rl" %in% names_raw) {
    temp_dees <- cbind(temp_dees, .25 * log(1 + temp_raw[, "rl"] / 100))
    names_dees <- c(names_dees, "rl")
  }
  
  names_temp <- c(names_dees, paste(names_raw, "raw", sep = "_"))
  temp <- cbind(temp_dees, temp_raw)
  dimnames(temp) <- list(NULL, names_temp)
  
  pos_i <- names_temp[which(names_temp %in% vars)]
  pos_i <- vars[which(vars %in% pos_i)]
  temp <- temp[, pos_i]
  
  data <- c(data, list(temp))
  nam <- c(nam, countries[i, 1])
  rm(temp)
}

names(data) <- countries[, 1]
country.data <- data

# Global data
global.variables <- "POIL"
global.data <- as.data.frame(read_xls("data-raw/dees2007/data33.xls", sheet = 1))[-1,]
global.data <- ts(matrix(global.data[, which(names(global.data) %in% global.variables)]), start = c(1979, 2), frequency = 4)
dimnames(global.data) <- list(NULL, "poil")
global.data <- log(global.data)


# Weight matrix
#weights <- array(NA, dim = c(length(country.data), length(country.data), 37))
#dimnames(weights) <- list(countries[, "Country Code"], countries[, "Country Code"], as.character(1980:2016))
weights <- c()
w.names <- c()
w_data <- read_xls("data-raw/dees2007/tradematrix8003.xls", sheet = "tradematrix8003_ne")
for (i in as.numeric(countries[, "Country Codes"])) {
  temp <- w_data[w_data$country_primary == i,]
  #dimnames(temp)[[1]] <- temp[, 1]
  dates <- temp[, 1]
  temp <- temp[, -(1:2)]
  # Reorder columns
  temp <- temp[ , paste("_", countries[, "Country Codes"], sep = "")]
  # Set home country to zero
  temp[, paste("_", i, sep = "")] <- 0
  temp <- ts(temp, start = dates$year[1], frequency = 1)
  dimnames(temp) <- list(as.character(dates$year), countries[, "Country Names"])
  weights <- c(weights, list(temp))
  w.names <- c(w.names, countries[which(countries["Country Codes"] == i), "Country Names"])
}
names(weights) <- w.names
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
dees2007 <- list("country_data" = country.data,
                 "global_data" = global.data,
                 "region_weights" = PPP,
                 "weight_data" = weight.data)

rm(list = ls()[-which(ls() == "dees2007")])

devtools::load_all(".")

country_data <- dees2007$country_data
global_data <- dees2007$global_data
region_weights <- dees2007$region_weights
weight_data <- dees2007$weight_data

# Generate EA area region with 2 year (2 periods), rolling window weights
ea <- c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")

# Adjust EA exchange rate
ea_fx <- exp(read.csv("data-raw/dees2007/Euroarea26.csv", stringsAsFactors = FALSE)[, "ep"])
ea_fx <- ts(ea_fx, start = c(1979, 2), frequency = 4)
ea_fx <- mean(ea_fx[which(floor(time(ea_fx)) == 2000)])
for (i in ea) {
  denum <- mean(exp(country_data[[i]][which(floor(time(country_data[[i]])) == 2000), "ep"]))
  country_data[[i]][, "ep"] <- log(ea_fx * exp(country_data[[i]][, "ep"]) / denum)
}

temp <- create_regions(country_data = country_data,
                       regions = list("EA" = ea),
                       period = 1999:2001,
                       region_weights = region_weights,
                       weight_data = weight_data)

country_data <- temp$country_data
weight_data <- temp$weight_data


  # Generate weight matrices as 3 year, rolling window averages
gvar_weights <- create_weights(weight_data = weight_data, period = 1999:2001,
                               country_data = country_data)


#### Specs ####
country_data <- diff_variables(country_data, variables = c("p"))

# Take values directly form 
country_data$EA[, "ep"] <- read.csv("data-raw/dees2007/Euroarea26.csv", stringsAsFactors = FALSE)[, "ep"]

# Create an object with country model specifications
model_specs <- create_specifications(country_data = country_data,
                                     global_data = global_data,
                                     domestic = list(variables = c("y", "p", "rs", "rl", "q", "ep"), lags = 1),
                                     foreign = list(variables = c("y", "p", "rs", "rl", "q", "ep", "poil"), lags = 1),                                     
                                     deterministic = list(const = "unrestricted", trend = "restricted"),
                                     countries = NULL,
                                     iterations = 10000,
                                     burnin = 5000,
                                     r = 1,
                                     type = "VEC")

# Argentina
model_specs$AR$domestic$variables <- c("y", "p", "rs", "q", "ep")
model_specs$AR$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$AR$domestic$lags <- 2
model_specs$AR$rank <- 2

# Australia
model_specs$AU$domestic$variables <- c("y", "p", "rs", "rl", "q", "ep")
model_specs$AU$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$AU$rank <- 4

# Brazil
model_specs$BR$domestic$variables <- c("y", "p", "rs", "ep")
model_specs$BR$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$BR$domestic$lags <- 2
model_specs$BR$rank <- 1

# Canada
model_specs$CA$domestic$variables <- c("y", "p", "rs", "rl", "q", "ep")
model_specs$CA$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$CA$rank <- 4

# China
model_specs$CN$domestic$variables <- c("y", "p", "rs", "ep")
model_specs$CN$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$CN$domestic$lags <- 2
model_specs$CN$rank <- 1

# Chile
model_specs$CL$domestic$variables <- c("y", "p", "rs", "q", "ep")
model_specs$CL$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$CL$domestic$lags <- 2
model_specs$CL$rank <- 2

# India
model_specs$IN$domestic$variables <- c("y", "p", "rs", "q", "ep")
model_specs$IN$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$IN$domestic$lags <- 2
model_specs$IN$rank <- 2

# Indonesia
model_specs$ID$domestic$variables <- c("y", "p", "rs", "ep")
model_specs$ID$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$ID$domestic$lags <- 2
model_specs$ID$rank <- 3

# Japan
model_specs$JP$domestic$variables <- c("y", "p", "rs", "rl", "q", "ep")
model_specs$JP$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$JP$rank <- 4

# South Korea
model_specs$KR$domestic$variables <- c("y", "p", "rs", "rl", "q", "ep")
model_specs$KR$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$KR$domestic$lags <- 2
model_specs$KR$rank <- 4

# Malaysia
model_specs$MY$domestic$variables <- c("y", "p", "rs", "q", "ep")
model_specs$MY$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$MY$domestic$lags <- 2

# Mexico
model_specs$MX$domestic$variables <- c("y", "p", "rs", "ep")
model_specs$MX$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$MX$rank <- 3

# Norway
model_specs$NO$domestic$variables <- c("y", "p", "rs", "rl", "q", "ep")
model_specs$NO$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$NO$domestic$lags <- 2
model_specs$NO$rank <- 2

# New Zealand
model_specs$NZ$domestic$variables <- c("y", "p", "rs", "rl", "q", "ep")
model_specs$NZ$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$NZ$domestic$lags <- 2
model_specs$NZ$rank <- 3

# Peru
model_specs$PE$domestic$variables <- c("y", "p", "rs", "ep")
model_specs$PE$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$PE$domestic$lags <- 2
model_specs$PE$rank <- 3

# Phillipines
model_specs$PH$domestic$variables <- c("y", "p", "rs", "q", "ep")
model_specs$PH$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$PH$domestic$lags <- 2
model_specs$PH$rank <- 2

# South Africa
model_specs$ZA$domestic$variables <- c("y", "p", "rs", "rl", "q", "ep")
model_specs$ZA$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$ZA$domestic$lags <- 2

# Saudia Arabia
model_specs$SA$domestic$variables <- c("y", "p", "ep")
model_specs$SA$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$SA$domestic$lags <- 2

# Singapore
model_specs$SG$domestic$variables <- c("y", "p", "rs", "q", "ep")
model_specs$SG$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$SG$domestic$lags <- 2
model_specs$SG$rank <- 3

# Sweden
model_specs$SE$domestic$variables <- c("y", "p", "rs", "rl", "q", "ep")
model_specs$SE$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$SE$domestic$lags <- 2
model_specs$SE$rank <- 3

# Switzerland
model_specs$CH$domestic$variables <- c("y", "p", "rs", "rl", "q", "ep")
model_specs$CH$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$CH$rank <- 3

# Thailand
model_specs$TH$domestic$variables <- c("y", "p", "rs", "q", "ep")
model_specs$TH$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$TH$rank <- 3

# Turkey
model_specs$TR$domestic$variables <- c("y", "p", "rs", "ep")
model_specs$TR$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$TR$domestic$lags <- 2

# United Kingdom
model_specs$GB$domestic$variables <- c("y", "p", "rs", "rl", "q", "ep")
model_specs$GB$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$GB$domestic$lags <- 2
model_specs$GB$rank <- 3

# Euro area
model_specs$EA$domestic$variables <- c("y", "p", "rs", "rl", "q", "ep")
model_specs$EA$foreign$variables <- c("y", "p", "q", "rs", "rl")
model_specs$EA$domestic$lags <- 2
model_specs$EA$foreign$lags <- 2
model_specs$EA$rank <- 2

# United States
model_specs$US$domestic$variables <- c("y", "p", "rs", "rl", "q", "poil")
model_specs$US$foreign$variables <- c("y", "p", "ep")
model_specs$US$domestic$lags <- 2
model_specs$US$foreign$lags <- 2
model_specs$US$rank <- 2


dees2007 <- list("country_data" = country_data,
                 "global_data" = global_data,
                 #"region_weights" = region_weights,
                 "weight_data" = gvar_weights,
                 "model_specs" = model_specs)

usethis::use_data(dees2007, overwrite = TRUE, version = 2)
