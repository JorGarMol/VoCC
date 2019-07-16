#' Shift in timing of seasonal climatology
#'
#' Function to calculate the seasonal shift in the arriving of typical seasonal climates
#' for a given period of interest as per Burrows et al. (2011).
#'
#' @usage shiftTime(r, yr1, yr2, yr0, th, m)
#'
#' @param r \code{stack} with monthly values of the climatic
#' variable for the period of interest.
#' @param yr0 \code{integer} specifying the first year in the series.
#' @param yr1 \code{integer} specifying the initial year for the period of interest.
#' @param yr2 \code{integer} specifying the end year for the period of interest.
#' @param th \code{integer} minimum number of non NAs in the series needed to
#' calculate the trend (default 3).
#' @param m \code{integer} number (1-12) of the month for which the shift is to be calculated
#'
#' @return a \code{stack} with the long-term monthly trend (C/year for temperature in degrees; "mTrend"),
#' seasonal rate of change (C/month; "seaRate"), and seasonal shift (day/decade; "seaShift").
#'
#' @references \href{http://science.sciencemag.org/content/334/6056/652}{Burrows et al. 2011}. The pace of
#' shifting climate in marine and terrestrial ecosystems. Science, 334, 652-655.
#'
#' @export
#' @author Jorge Garcia Molinos and Michael T. Burrows
#' @examples
#' Apr <- shiftTime(HSST, yr1 = 1960, yr2 = 2009, yr0 = 1955, th = 10, m = 4)
#'
#' plot(Apr)
#'
#' @rdname shiftTime


shiftTime <- function(r, yr1, yr2, yr0, th, m){
# 1. Long term trends in monthly values (e.g. deg/year if temperature)
m1 <- ((yr1-yr0)*12)+ m
m2 <- ((yr2-yr0)*12)+ m
r1 <- r[[seq(m1, m2, by = 12)]]
trend <- tempTrend(r1, th)[[1]]
# 2. seasonal rate of shift in temperature centered on each month (deg/month) = difference in
# temperature between  preceding and following months divided by 2 months (slope) preceding month
b <- ifelse((m-1) == 0, 12, (m-1))
m1 <- ((yr1-yr0)*12)+b
m2 <- ((yr2-yr0)*12)+b
x2 <- r[[seq(m1, m2, by = 12)]]
# following month
b <- ifelse((m+1) == 13, 1, (m+1))
m1 <- ((yr1-yr0)*12)+ b
m2 <- ((yr2-yr0)*12)+ b
x3 <- r[[seq(m1, m2, by = 12)]]
# slope
x4 <- raster::mean((x3-x2)/2, na.rm = TRUE)

# 3. seasonal shifts (month/year) converted to days per decade by multiplying by 10 years, 365.25 days per year and dividing by 12 months
sShift <- (trend/x4)*(3652.5/12)
sShift[sShift == Inf | sShift == -Inf] <- NA
r2 <- raster::stack(trend, x4, sShift)
names(r2) <- c("mTrend", "seaRate", "seaShift")
return(r2)
}
