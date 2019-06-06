#' Summarize climatic series to higher temporal resolution
#'
#' Function to convert climatic series (provided as \code{RasterStack}) into a
#' coarser time frequency series for a period of interest. This function transforms the \code{RasterStack}
#' into an \code{xts} time series object to extract the values for the period of interest and
#' apply some summary function. It is mainly a wrapper from the \code{apply.} function family
#' in the package xts (Ryan and Ulrich 2017).
#'
#' @usage sumSeries(r, p, yr0, l, fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months",
#' freqout = "years")
#'
#' @param r \code{RasterStack} containing the time series of the climatic variable.
#' @param p \code{character string} defining the period to extract for the calculation
#' of the series (see examples).
#' @param yr0 \code{character string} specifying the first (yr0) year in the series (see examples).
#' @param l \code{integer} length of the input time series.
#' @param fun \code{logical} summary function to be computed. Summary functions need to be applied by cell (columns)
#' so should have the structure 'function(x) apply(x, 2, function(y){...})'. For convinience, sumSeries imports
#' colMaxs, and colMins from package ‘matrixStats’ (Bengtsson 2018) so they can be called in directly.
#' @param freqin \code{character string} specifying the original time frequency of the series.
#' @param freqout \code{character string} specifying the desired time frequency of the new series.
#' Must be one of the following: "weeks", "months", "quarters", "years", "other". Argument "other"
#' allows for user-defined functions to be applied on the 'xts' time series object over the period of interest (see examples).
#'
#' @return A \code{RasterStack} with the new series.
#'
#' @references \href{https://CRAN.R-project.org/package=xts}{Ray and Ulrich. 2017}. xts: eXtensible Time Series. R package version 0.10-1. \cr
#' \href{https://CRAN.R-project.org/package=matrixStats}{Bengtsson 2018}. matrixStats: Functions that Apply to Rows and Columns
#' of Matrices (and to Vectors). R package version 0.53.1.
#'
#' @importFrom matrixStats colMaxs colMins
#' @importFrom xts xts apply.weekly apply.monthly apply.quarterly apply.yearly
#' @export
#' @author Jorge Garcia Molinos
#' @examples
#' # Monthly mean SST (HadISST) data for Europe Jan-1950 to Dec-2010
#'
#' ?HSST
#'
#' # Calculate mean annual monthly SST
#'
#' yrSST <- sumSeries(HSST, p = "1969-01/2009-12", yr0 = "1955-01-01", l = nlayers(HSST),
#' fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years")
#'
#' # Extract Jul Aug mean SST each year (xts months are indexed from 0 to 11)
#'
#' myf = function(x, m = c(7,8)){
#' x[xts::.indexmon(x) %in% (m-1)]
#' }
#'
#' JlAugSST <- sumSeries(HSST, p = "1969-01/2009-12", yr0 = "1950-01-01", l = raster::nlayers(HSST),
#' fun = myf, freqin = "months", freqout = "other")
#'
#' # Same but calculating the annual variance of the two months
#'
#' myf = function(x, m = c(7,8)){
#' x1 <- x[xts::.indexmon(x) %in% (m-1)]
#' xts::apply.yearly(x1, function(y) apply(y, 2, function(y){var(y, na.rm = TRUE)}))
#' }
#'
#' meanJASST <- sumSeries(HSST, p = "1969-01/2009-12", yr0 = "1950-01-01", l = raster::nlayers(HSST),
#' fun = myf, freqin = "months", freqout = "other")
#'
#' @rdname sumSeries

sumSeries <- function(r, p, yr0, l = nlayers(r), fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years"){
# construct xts object
m <- t(getValues(r))
dates <- seq(as.Date(yr0), length = l, by = freqin)
ts1 <- xts(m, order.by = dates)
# subset for the period of interest
x <- ts1[p]

# calculate the annual series
if(freqout == "weeks"){s <- apply.weekly(x, fun)}
if(freqout == "months"){s <- apply.monthly(x, fun)}
if(freqout == "quarters"){s <- apply.quarterly(x, fun)}    # Jan-Mar (Q1); Apr-Jn (Q2); Jl-Sep(Q3); Oct-Dec (Q4)
if(freqout == "years"){s <- apply.yearly(x, fun)}
if(freqout == "other"){s <- fun(x)}

# create raster stack
r1 <- stack()
for(i in 1:nrow(s)){
r2 <- raster(r[[1]])
r2[] <-  as.numeric(s[i,])
r1 <- addLayer(r1, r2)
}
if(freqout != "other"){names(r1) <- seq(start(x), length = nlayers(r1), by = freqout)}
return(r1)
}



