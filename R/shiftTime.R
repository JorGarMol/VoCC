#' Shift in timing of seasonal climatology
#'
#' Function to calculate the seasonal shift in the arriving of seasonal climatology
#' for a given period of interest as per Burrows et al. 2011.
#'
#' @usage shiftTime(r, yr1, yr2, yr0, th, month)
#'
#' @param r \code{stack} with monthly values of the climatic
#' variable for the period of interest.
#' @param yr0 \code{integer} specifying the first year in the series.
#' @param yr1 \code{integer} specifying the initial year for the period of interest.
#' @param yr2 \code{integer} specifying the end year for the period of interest.
#' @param th \code{integer} max number of non NAs in the series needed to
#' calculate the trend (default 3).
#' @param month \code{integer} number (1-12) of the month for which the shift is to be calculated
#'
#' @return a \code{stack} with the long-term monthly trend (C/year; "mTrend"),
#' seasonal rate of change (C/month; "seaRate"), and seasonal shift (day/decade; "seaShift").
#'
#' @references \href{http://science.sciencemag.org/content/334/6056/652}{Burrows et al. 2011}. The pace of shifting climate in marine and terrestrial ecosystems. Science, 334, 652-655.
#'
#' @export
#' @author Jorge Garcia Molinos and Michael T. Burrows
#' @examples
#' data(HSST)
#' Apr <- shiftTime(HSST, yr1 = 1960, yr2 = 2009, yr0 = 1955, th = 10, month = 4)
#' plot(Apr)
#'
#' @rdname shiftTime


shiftTime <- function(r, yr1, yr2, yr0, th, month){
# 1. Long term trends in monthly temperatures(deg/year)
m1 <- ((yr1-yr0)*12)+ month
m2 <- ((yr2-yr0)*12)+ month
r1 <- r[[seq(m1, m2, by = 12)]]
trend <- tempTrend(r1, th)[[1]]

# 2. seasonal rate of shift in temperature centered on each month (deg/month) = difference in temperature between  preceding and following months divided by 2 months (slope)
# preceding month
b <- ifelse((month-1) == 0, 12, (month-1))
m1 <- ((yr1-yr0)*12)+b
m2 <- ((yr2-yr0)*12)+b
x2 <- r[[seq(m1, m2, by = 12)]]
# following month
b <- ifelse((month+1) == 13, 1, (month+1))
m1 <- ((yr1-yr0)*12)+ b
m2 <- ((yr2-yr0)*12)+ b
x3 <- r[[seq(m1, m2, by = 12)]]
# slope
x4 <- mean((x3-x2)/2, na.rm = TRUE)

# 3. seasonal shifts (month/year) converted to days per decade by multiplying by 10 years, 365.25 days per year and dividing by 12 months
sShift <- (trend/x4)*(3652.5/12)
# replace Inf by NA (spurious values resulting from sea ice cells having 0 x4)
sShift[sShift == Inf | sShift == -Inf] <- NA
r2 <- raster::stack(trend, x4, sShift)
names(r2) <- c("mTrend", "seaRate", "seaShift")
return(r2)
}
