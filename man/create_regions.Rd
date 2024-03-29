% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_regions.R
\name{create_regions}
\alias{create_regions}
\title{Create Regional Series}
\usage{
create_regions(country_data, weight_data, region_weights, regions, period)
}
\arguments{
\item{country_data}{a named list of time-series objects containing country data.}

\item{weight_data}{a named list of named time-series objects. See 'Details'.}

\item{region_weights}{a multivariate time-series object containing data used
to weight the observations in \code{country_data}.}

\item{regions}{a named list of character vectors containing specifications for countries
in \code{country_data}, which should be combined to a region. The name of the
respective list element will become the name of the region.}

\item{period}{either a single integer or a numeric vector specifiying the periods in
\code{region_weights}, which should be used to construct weights.}
}
\value{
A list containing the following elements:
\item{country_data}{a named list of updated time-series objects with country data.}
\item{weight_data}{a named list of named time-series objects for the
generation of country-specific weight matrices.}
}
\description{
Combines multiple country-specific time series to regional series.
}
\details{
If a numeric vector is provided for argument \code{period}, the function weights
country-specific observations based on the sums over the specified periods.
If an integer is proved, the country-specific observations are weighted
according to rolling sums over the last \code{period} periods. If country
data start earlier than the series in \code{region_weights}, the sums over
the first \code{period} observations of \code{region_weights} are used
until the periods match.

The elements of \code{weight_data} should be time-series objects,
where colums are named after the trading partners of the respective country
and the names of the elements should match the country names in
\code{country_data}. For each country list element the time-series object should
contain zeros in the column of the respective country. For examle, the column
"US" in the "US"-element of the list should only contain zeros.
}
\examples{
# Load data
data("gvar2019")

# Create regions
temp <- create_regions(country_data = gvar2019$country_data,
             weight_data = gvar2019$weight_data,
             region_weights = gvar2019$region_weights,
             regions = list(EA =  c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")),
             period = 3)

# New data sets
country_data <- temp$country_data
weight_data <- temp$weight_data 

}
