#' Summary annual values for the climatic series.
#'
#' Function to calculate the mean/max/min monthly annual series from a monthly  
#' series for a period of interest.
#' 
#' @param r \code{RasterStack} with mean monthly values of the climatic variable
#'  for the period of interest. Only series with complete years (no missing data)
#'  should be used.
#' @param yr0,1,2 \code{Integers} specifying the first (yr0) year in the series, 
#' the initial (yr1) and the last (yr2) years for the period of interest.
#' @param comp \code{logical} specifying the summary statistic to be computed. 
#' Must be one of mean (default), max or min.
#'
#' @return A \code{RasterStack} with the annual mean, max, or min values.
#' @import raster
#' @export
#' @author Jorge Garcia Molinos
#' @rdname anSeries

anSeries <- function(r, yr1, yr2, yr0, comp = "mean"){

r1 <- stack()
m1 <- (yr1-yr0)*12+1   # position starting month 
m2 <- ((yr2-yr0)*12+1)+11            # position ending month

x <- getValues(r[[m1:m2]])
ocean <- which(rowSums(is.na(x))!= ncol(x))    # remove land cells  (whole row NAs)
x <- t(x[ocean,])
f <- rep(1:(yr2-yr1+1), each = 12)

for(i in 1:(yr2-yr1+1)){  # iterate from the beginning at 12-month steps
rr <- raster(r[[1]])
if(comp == "mean"){rr[ocean] <- colMeans(x[which(f == i),], na.rm = T)}
if(comp == "max"){rr[ocean] <- do.call(pmax, c(as.data.frame(t(x[which(f == i),])), na.rm = TRUE))}
if(comp == "min"){rr[ocean] <- do.call(pmin, c(as.data.frame(t(x[which(f == i),])), na.rm = TRUE))}
r1 <- addLayer(r1, rr)
}

return(r1)
}
