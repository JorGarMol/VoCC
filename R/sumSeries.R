#' Summarize climatic series to higher temporal resolutions
#'
#' Function to convert climatic series (provided as \code{RasterStack}) into a
#' coarser time frequency series for a period of interest. This function is a wrapper
#' of the apply. function family from the package xts.
#'
#' @param r \code{RasterStack} containing the time series of the climatic variable.
#' @param p \code{character string} defining the period to extract for the calculation
#' of the series (see examples).
#' @param yr0 \code{character string} specifying the first (yr0) year in the series (see examples).
#' @param l \code{integer} length of the input time series.
#' @param fun \code{logical} summary function to be computed. Summary functions need to be applied by cell (columns)
#' so should have the structure 'function(x) apply(x, 2, function(y){...})'. For convinience, sumSeries imports colMeans,
#' colMaxs, and colMins from package ‘matrixStats’ so they can be called in directly.
#' @param freqin \code{character string} specifying the original time frequency of the series.
#' @param freqout \code{character string} specifying the desired time frequency of the new series.
#' Must be one of the following: "weeks", "months", "quarters", "years". Alternatively,
#' an \code{integer} representing a month (1-12) can be provided as argument where an
#' annual series is required for a specific month.
#'
#' @return A \code{RasterStack} with the new series.
#' @references Jeffrey A. Ryan and Joshua M. Ulrich (2017). xts: eXtensible Time Series. R package version 0.10-1.
#'  \url{https://CRAN.R-project.org/package=xts}
#'  Henrik Bengtsson (2018). matrixStats: Functions that Apply to Rows and Columns of Matrices (and to Vectors). R package version 0.53.1.
#'  \url{https://CRAN.R-project.org/package=matrixStats}
#' @import raster xts
#' @importFrom matrixStats colMaxs colMins colMeans
#' @export
#' @author Jorge Garcia Molinos
#' @rdname sumSeries

sumSeries <- function(r, p = "1969-01/2009-12", yr0 = "1870-01-01", l = nlayers(r), fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years"){
t0 = as.Date(yr0)
# construct xts object
m <- t(getValues(r))
dates <- seq(as.Date(t0), length = l, by = freqin)
ts1 <- xts(m, order.by = dates)
# subset for the period of interest
ts2 <- ts1[p]

# calculate the annual series
if(freqout == "weeks"){s <- apply.weekly(ts2, fun)}
if(freqout == "months"){s <- apply.monthly(ts2, fun)}
if(freqout == "quarters"){s <- apply.quarterly(ts2, fun)}    # Jan-Mar (Q1); Apr-Jn (Q2); Jl-Sep(Q3); Oct-Dec (Q4)
if(freqout == "years"){s <- apply.yearly(ts2, fun)}
if(is.numeric(freqout)){s <- ts2[.indexmon(ts2) %in% (freqout-1)]}             # xts months are indexed from 0 to 11

# create raster stack
r1 <- stack()
for(i in 1:nrow(s)){
r2 <- raster(r[[1]])
r2[] <-  as.numeric(s[i,])
r1 <- addLayer(r1, r2)
}

return(r1)
}



