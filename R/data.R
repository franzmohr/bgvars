#' Global VAR Database (1979Q2-2016Q4)
#'
#' The data set of Mohaddes and Raissi (2018) contains economic time series for
#' 33 countries and 3 commodities from 1979Q2 to 2016Q4.
#' 
#' @usage data("gvar2016")
#' 
#' @format A named list containing the following elements:
#' \describe{
#'   \item{\strong{country_data}}{A named list of 33 time-series objects with 151 quarterly
#'   observations from 1979Q2 to 2016Q4 of the following variables:
#'   
#'   \itemize{
#'   \item \strong{y:} log real GDP
#'   \item \strong{Dp:} rate of inflation
#'   \item \strong{r:} short-term interest rate
#'   \item \strong{lr:} long-term interest rate
#'   \item \strong{ep:} log deflated exchange rate
#'   \item \strong{eq:} log real equity prices
#'   }
#'   
#'   The variables are available for the following countries
#' \tabular{lccccccc}{
#' \strong{Country} \tab \strong{Code}  \tab \strong{\code{y}} \tab \strong{\code{Dp}}  \tab \strong{\code{r}}
#'  \tab \strong{\code{lr}} \tab \strong{\code{ep}} \tab \strong{\code{eq}} \cr
#' Argentina \tab AR \tab x \tab x \tab x \tab   \tab x \tab x \cr
#' Australia \tab AU \tab x \tab x \tab x \tab x \tab x \tab x \cr
#' Austria \tab AT \tab x \tab x \tab x \tab x \tab x \tab x \cr
#' Belgium \tab BE \tab x \tab x \tab x \tab x \tab x \tab x  \cr
#' Brazil \tab BR \tab x \tab x \tab x \tab  \tab x \tab   \cr
#' Canada \tab CA  \tab x \tab x \tab x \tab x \tab x \tab x \cr
#' China \tab CN  \tab x \tab x \tab x \tab  \tab x \tab   \cr
#' Chile \tab CL  \tab x \tab x \tab x \tab  \tab x \tab x \cr
#' Finland \tab FI  \tab x \tab x \tab x \tab  \tab x \tab x  \cr
#' France \tab FR  \tab x \tab x \tab x \tab x \tab x \tab x  \cr
#' Germany \tab DE \tab x \tab x \tab x \tab x \tab x \tab x  \cr
#' India \tab IN \tab x \tab x \tab x \tab x \tab x \tab x  \cr
#' Indonesia \tab ID \tab x \tab x \tab x \tab  \tab x \tab    \cr
#' Italy \tab IT \tab x \tab x \tab x \tab x \tab x \tab x  \cr
#' Japan \tab JP \tab x \tab x \tab x \tab x \tab x \tab x  \cr
#' Korea \tab KR \tab x \tab x \tab x \tab x \tab x \tab x  \cr
#' Malaysia \tab MY  \tab x \tab x \tab x \tab  \tab x \tab x  \cr
#' Mexico \tab MX  \tab x \tab x \tab x \tab  \tab x \tab  \cr
#' Netherlands \tab NL \tab x \tab x \tab x \tab x \tab x \tab x \cr
#' Norway \tab NO  \tab x \tab x \tab x \tab x \tab x \tab x  \cr
#' New Zealand \tab NZ  \tab x \tab x \tab x \tab x \tab x \tab x  \cr
#' Peru \tab PE  \tab x \tab x \tab x \tab  \tab x \tab  \cr
#' Philippines \tab PH \tab x \tab x \tab x \tab  \tab x \tab x \cr
#' South Africa \tab ZA  \tab x \tab x \tab x \tab x \tab x \tab x  \cr
#' Saudi Arabia \tab SA  \tab x \tab x \tab  \tab  \tab x \tab   \cr
#' Singapore \tab SG \tab x \tab x \tab x \tab  \tab x \tab x  \cr
#' Spain \tab ES \tab x \tab x \tab x \tab x \tab x \tab x  \cr
#' Sweden \tab SE  \tab x \tab x \tab x \tab x \tab x \tab x  \cr
#' Switzerland \tab CH \tab x \tab x \tab x \tab x \tab x \tab x  \cr
#' Thailand \tab TH  \tab x \tab x \tab x \tab  \tab x \tab x \cr
#' Turkey \tab TR  \tab x \tab x \tab x \tab  \tab x \tab \cr
#' United Kingdom \tab  GB \tab x \tab x \tab x \tab x \tab x \tab x  \cr
#' USA \tab US  \tab x \tab x \tab x \tab x \tab  \tab x  \cr
#'}
#'   }
#'   \item{\strong{global_data}}{A time-series object containing 151 quarterly observations
#'   from 1979Q2 to 2016Q4 of the following variables
#'   \itemize{
#'   \item \strong{poil:} log of oil prices
#'   \item \strong{pmat:} log of agricultural raw material prices
#'   \item \strong{pmetal:} log of metals prices
#'   }}
#'   
#'   \item{\strong{region_weights}}{A time-series object containing annual data on PPP-GDP
#'   from 1990 to 2016 for the countries listed above.}
#'   
#'   \item{\strong{weight_data}}{A named list of time-series objects containing annual data on
#'   trade flows from 1980 to 2016 for the countries listed above.}
#' }
#' 
#' @details The variables in \code{country_data} were constructed as
#'   \itemize{
#'   \item \strong{y:} \eqn{y_{it} = ln(GDP_{it})}
#'   \item \strong{Dp:} \eqn{Dp_{it} = p_{it} - p_{it-1}} with \eqn{p_{it} = ln(CPI_{it})}
#'   \item \strong{r:} \eqn{r_{it} = 0.25 ln(1 + R_{it}^{S} / 100)}
#'   \item \strong{lr:} \eqn{r_{it} = 0.25 ln(1 + R_{it}^{L} / 100)}
#'   \item \strong{ep:} \eqn{ep_{it} = ln(E_{it} / CPI_{it})}
#'   \item \strong{eq:} \eqn{eq_{it} = ln(EQ_{it} / CPI_{it})}
#   \item \strong{cred:} \eqn{credit_{it} = ln(CREDIT_{it} / CPI_{it})}
#'   }
#'  where \eqn{GDP_{it}} is real gross domestic product at time \eqn{t} for country \eqn{i}, \eqn{CPI_{it}}
#'  is the consumer price index, \eqn{E_{it}} is the nominal exchange rate in terms of US dollar, \eqn{EQ_{it}}
#'  is the nominal equity price index, and \eqn{R_{it}^{S}} and \eqn{R_{it}^{L}} are short-term and long-term
#'  interest rates, respectively.
#' 
#' See Mohaddes and Raissi (2018) for further details on the sources and compilation of all variables.
#' 
#' @references
#' 
#' Mohaddes, K., & Raissi, M. (2018). \emph{Compilation, revision and updating of the global VAR (GVAR) database, 1979Q2--2016Q4} (mimeo). University of Cambridge: Faculty of Economics.
#' The paper and dataset can be downloaded from \href{https://www.mohaddes.org/gvar}{https://www.mohaddes.org/gvar}.
#' 
"gvar2016"

#' Dees, di Mauro, Pesaran \& Smith data set
#'
#' The data set of Dees et al. (2007) contains economic time series for
#' 26 economies and oil prices from 1979Q2 to 2003Q4.
#' 
#' @usage data("dees2007")
#' 
#' @format A named list containing the following elements:
#' \describe{
#'   \item{\strong{country_data}}{A named list of 26 time-series objects with 99 quarterly
#'   observations from 1979Q2 to 2003Q4 of the following variables:
#'   
#'   \itemize{
#'   \item \strong{y:} log real GDP
#'   \item \strong{p:} rate of inflation
#'   \item \strong{rs:} short-term interest rate
#'   \item \strong{rl:} long-term interest rate
#'   \item \strong{q:} log real equity prices
#'   \item \strong{ep:} log real exchange rate
#'   }
#'   }
#'   
#'   \item{\strong{global_data}}{A time-series object containing 99 quarterly observations
#'   from 1979Q2 to 2003Q4 of log of oil prices.}   
#'   \item{\strong{region_weights}}{A time-series object containing annual data on PPP-GDP
#'   from 1999 to 2001 for 26 economies.}
#'   \item{\strong{weight_data}}{A named list of time-series objects containing annual data on
#'   trade flows from 1980 to 2003 for 26 economies.}
#'   \item{\strong{dees_specs}}{A named list of country model specifications as used in
#'   Dees et al. (2007).}
#' }
#' 
#' See the supplementary material to Dees et al. (2007) for further details on the
#' sources and compilation of the data (\url{http://qed.econ.queensu.ca/jae/datasets/dees001/}).
#' 
#' @references
#' 
#' Dees, S., di Mauro, F., Pesaran, M. H., & Smith, V. (2007). Exploring the
#' international linkages of the euro area: A global VAR analysis.
#' \emph{Journal of Applied Econometrics, 22}(1), 1--38.
#' \url{https://doi.org/10.1002/jae.932}
#' 
"dees2007"

