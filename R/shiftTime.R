#' Shift in timing of seasonal temperatures
#'
#' Function to calculate the seasonal shift in the arriving of mean monthly temperatures 
#' for the study period as per Burrows et al. 2011.
#' 
#' @param rasStack \code{stack} with monthly values of the climatic 
#' variable for the period of interest.
#' @param yr0,1,2 \code{integer} specifying the first (yr0) year in the series, 
#' and the initial (yr1) and end (yr2) years for the period of interest.
#' @param th \code{integer} max number of non NAs in the series needed to 
#' calculate the trend (default 3).  
#' @param month \code{integer} number (1-12) of the month for which the shift is to be calculated
#'
#' @return a \code{stack} with the long-term monthly trend (C/year; "mTrend"), 
#' seasonal rate of change (C/month; "seaRate"), and seasonal shift (day/month; "seaShift").
#' 
#' @import raster
#' @export
#' @author Jorge Garcia Molinos and Michael T. Burrows
#' @rdname shiftTime


shiftTime <- function(rasStack, yr1, yr2, yr0, th, month){
# 1. Long term trends in monthly temperatures(deg/year)
m1 <- ((yr1-yr0)*12)+ month    
m2 <- ((yr2-yr0)*12)+ month
r <- rasStack[[seq(m1, m2, by = 12)]]
trend <- tempTrend(r, th)[[1]]
                    
# 2. seasonal rate of shift in temperature centered on each month (deg/month) = difference in temperature between  preceding and following months divided by 2 months (slope)
# preceding month
b <- ifelse((month-1) == 0, 12, (month-1))
m1 <- ((yr1-yr0)*12)+b    
m2 <- ((yr2-yr0)*12)+b
x2 <- rasStack[[seq(m1, m2, by = 12)]]
# following month
b <- ifelse((month+1) == 13, 1, (month+1))
m1 <- ((yr1-yr0)*12)+ b    
m2 <- ((yr2-yr0)*12)+ b
x3 <- rasStack[[seq(m1, m2, by = 12)]]
# slope
x4 <- mean((x3-x2)/2, na.rm = TRUE)

# 3. seasonal shifts (month/year) converted to days per decade by multiplying by 10 years, 365.25 days per year and dividing by 12 months
sShift <- (trend/x4)*(3652.5/12)
# replace Inf by NA (spurious values resulting from sea ice cells having 0 x4)
sShift[sShift == Inf | sShift == -Inf] <- NA
r <- stack(trend, x4, sShift)
names(r) <- c("mTrend", "seaRate", "seaShift")
return(r)
}