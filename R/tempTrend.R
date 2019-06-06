#' Long-term local climatic trends
#'
#' Function to calculate temporal trend from a raster series
#' of a climatic variable. This trend is to be used for the calculation of the
#' gradient-based climate velocity using gVoCC.
#'
#' @usage tempTrend(r, th)
#'
#' @param r \code{RasterStack} containing a time series of (annual, seasonal, monthly...) values of
#' the climatic variable for the period of interest.
#' @param th \code{Integer} minimum number of observations in the series needed to
#' calculate the trend at each cell.
#'
#' @return A \code{RasterStack} containing the cell-specific temporal trends
#' extracted from simple linear regressions of the climatic variable against time
#' ("slpTrends" in degree Celsius per year), together with their standard
#' errors ("seTrends") and statistical significance ("sigTrends").
#'
#' @seealso{\code{\link{spatGrad}}, \code{\link{gVoCC}}}
#'
#' @export
#' @author Jorge Garcia Molinos and Christopher J. Brown
#' @examples
#'
#' yrSST <- sumSeries(HSST, p = "1969-01/2009-12", yr0 = "1955-01-01", l = nlayers(HSST),
#' fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years")
#'
#' # Mean annual SST trend (minimum threshold of 10 years of data), with SE and p-values.
#'
#' tr <- tempTrend(yrSST, th = 10)
#'
#' plot(tr)
#'
#' @rdname tempTrend

tempTrend <- function(r, th) {
y <- getValues(r)
ocean <- which(rowSums(is.na(y))!= ncol(y))    # remove land cells
y <- t(y[ocean, ])
N <- apply(y, 2, function(x) sum(!is.na(x)))
ind <- which(N >= th)
y <- y[,ind]  # drop cells with less than th observations
N <- apply(y, 2, function(x) sum(!is.na(x)))
x <- matrix(nrow = nlayers(r), ncol = ncol(y))
x[] <- 1:nlayers(r)
# put NA values into the x values so they correspond with y
x1 <- y
x1[!is.na(x1)] <- 1
x <- x*x1
# calculate the sum terms
sx <- apply(x, 2, sum, na.rm = T)
sy <- apply(y, 2, sum, na.rm = T)
sxx <- apply(x, 2, function(x) sum(x^2, na.rm = T))
syy <- apply(y, 2, function(x) sum(x^2, na.rm = T))
xy <- x*y
sxy <- apply(xy, 2, sum, na.rm = T)
# Estimate slope coefficients and associated standard errors and p-values
slope <- (N*sxy-(sx*sy))/(N*sxx-sx^2)
sres <- (N*syy-sy^2-slope^2*(N*sxx-sx^2))/(N*(N-2))
SE <- suppressWarnings(sqrt((N*sres)/(N*sxx-sx^2)))
Test <- slope/SE
p <- mapply(function(x,y) (2*pt(abs(x), df = y-2, lower.tail = FALSE)), x = Test, y = N)

slpTrends <- sigTrends <- seTrends <- raster(r[[1]])
slpTrends[ocean[ind]] <- slope
seTrends[ocean[ind]] <- SE
sigTrends[ocean[ind]] <- p
output <- stack(slpTrends,seTrends,sigTrends)
names(output) <- c("slpTrends", "seTrends", "sigTrends")
return(output)
}
