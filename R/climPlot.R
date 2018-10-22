#' Binned scatter plot for 2-dimensional climate space
#'
#' Function to create a binned scatter plot of two climate variables.
#'
#' @usage   climPlot(xy, x.binSize, y.binSize)
#'
#' @param xy \code{data.frame} with cells as rows and 4 columns representing the present and future local values for the two variables (V1p, V1f, V2p, V2f).
#' @param x.binSize \code{numeric} the bin size for the first variable.
#' @param y.binSize \code{numeric} the bin size for the second variable.
#'
#' @return a series of \code{plot} objects displaying the (i) present and (ii) future
#' cell frequency for each combination of local climates,
#' and (iii) the location of remnant, novel and disappearing climates between both periods.
#'
#' @seealso{\code{\link{climAna}}, \code{\link{climPCA}}}
#'
#' @import data.table
#' @importFrom RColorBrewer brewer.pal
#' @importFrom fields image.plot
#' @export
#' @author Jorge Garcia Molinos
#' @examples
#'
#' data(JapTC)
#' # Plot climate space for the two first variables(annual precipitation and maximum temperature)
#' xy <- na.omit(data.frame(getValues(JapTC[[1]]), getValues(JapTC[[2]]),
#' getValues(JapTC[[3]]), getValues(JapTC[[4]])))
#' climPlot(xy, x.binSize = 5, y.binSize = 0.2)
#'
#' @rdname climPlot


climPlot <- function(xy, x.binSize, y.binSize){
xp = xy[,1]
yp = xy[,3]
xf = xy[,2]
yf = xy[,4]
# bins per axis
x.nbins = floor((abs(range(xp, xf)[2]-range(xp, xf)[1]))/x.binSize)
y.nbins = floor((abs(range(yp, yf)[2]-range(yp, yf)[1]))/y.binSize)
x.bin <- seq(floor(min(cbind(xp, xf))), ceiling(max(cbind(xp, xf))), length = x.nbins)
y.bin <- seq(floor(min(cbind(yp, yf))), ceiling(max(cbind(yp, yf))), length = y.nbins)

# define palette
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(64)

# present
freq <-  as.data.frame(table(findInterval(xp, x.bin), findInterval(yp, y.bin)))
freq[,1] <- as.numeric(freq[,1])
freq[,2] <- as.numeric(freq[,2])
freq2Dp <- diag(nrow = x.nbins, ncol = y.nbins)
freq2Dp[cbind(freq[,1], freq[,2])] <- freq[,3]
freq2Dp[freq2Dp[] == 0] <- NA

# future
freq <-  as.data.frame(table(findInterval(xf, x.bin), findInterval(yf, y.bin)))
freq[,1] <- as.numeric(freq[,1])
freq[,2] <- as.numeric(freq[,2])
freq2Df <- diag(nrow = x.nbins, ncol = y.nbins)
freq2Df[cbind(freq[,1], freq[,2])] <- freq[,3]
freq2Df[freq2Df[] == 0] <- NA

# Adjust values above the upper limit (minimum of the two maxima) to the upper limit so image.plot plot all on the same color scale
UL <- min(max(freq2Dp, na.rm=TRUE), max(freq2Df, na.rm=TRUE))
freq2Dp[freq2Dp > UL] <- UL
freq2Df[freq2Df > UL] <- UL
# plot climate space
fields::image.plot(x.bin, y.bin, freq2Dp, col=r, xlab = "Variable 1", ylab = "Variable 2", main = 'Frequency of local climates present conditions')
fields::image.plot(x.bin, y.bin, freq2Df, col=r, xlab = "Variable 1", ylab = "Variable 2", main = 'Frequency of local climates future conditions')

# novel (in future but not present, 2), remnant (in both, 1), and dissapearing (in present but not future, 3) climates
freq2D <- diag(nrow = x.nbins, ncol = y.nbins)
freq2D[] <- NA
for(i in 1:x.nbins){
for(j in 1:y.nbins){
freq2D[i,j] <- ifelse(is.na(freq2Dp[i,j]) & !is.na(freq2Df[i,j]), 1, ifelse(!is.na(freq2Dp[i,j]) & is.na(freq2Df[i,j]), 2, ifelse(is.na(freq2Dp[i,j]) & is.na(freq2Df[i,j]), NA, 0)))
}}
image(x.bin, y.bin, freq2D, col = c("#56B4E9", "#009E73", "#D55E00"), xlab = "Variable 1", ylab = "Variable 2", main = 'Remnant (blue), novel (green), and disappearing (red) climates')
}




