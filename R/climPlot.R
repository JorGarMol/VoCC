#' Binned scatter plot for 2-dimensional climate space
#'
#' Function to create a binned scatter plot of two climate variables.
#'
#' @usage   climPlot(xy, x.binSize, y.binSize, x.name="V1", y.name="V2")
#'
#' @param xy \code{data.frame} with cells as rows and 4 columns representing the present and future local values for the two variables (V1p, V1f, V2p, V2f).
#' @param x.binSize \code{numeric} the bin size for the first variable.
#' @param y.binSize \code{numeric} the bin size for the second variable.
#' @param x.name \code{character} the variable name for the first variable. Used to label the plot.
#' @param y.name \code{character} the variable name for the second variable. Used to label the plot.
#'
#' @return A series of \code{plot} objects displaying the (i) present and (ii) future
#' cell frequency for each combination of local climates,
#' and (iii) the location of remnant, novel and disappearing climates between both periods.
#'
#' @seealso{\code{\link{dVoCC}}, \code{\link{climPCA}}}
#'
#' @import data.table
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom cowplot plot_grid
#' @export
#' @author Jorge Garcia Molinos and Naoki H. Kumagai
#' @examples
#'
#' # Plot climate space for the two first variables(annual precipitation and maximum temperature)
#' xy <- na.omit(data.frame(getValues(JapTC[[1]]), getValues(JapTC[[2]]),
#' getValues(JapTC[[3]]), getValues(JapTC[[4]])))
#'
#' out <- climPlot(xy, x.binSize = 5, y.binSize = 0.2, x.name="Precipitation (mm)",
#' y.name="Temperature max (°C)")
#'
#' # output plots can be saved as:
#' ggplot2::ggsave(plot=out, filename=paste0(getwd(), "/example_plot.pdf"), width=17, height=17, unit="cm")
#'
#'
#' @rdname climPlot


climPlot <- function(xy, x.binSize, y.binSize, x.name="V1", y.name="V2"){
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
  rf <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,'Spectral')))
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

  # novel (in future but not present, 2), remnant (in both, 1), and dissapearing (in present but not future, 3) climates
  freq2D <- diag(nrow = x.nbins, ncol = y.nbins)
  freq2D[] <- NA
  for(i in 1:x.nbins){
    for(j in 1:y.nbins){
      freq2D[i,j] <- ifelse(is.na(freq2Dp[i,j]) & !is.na(freq2Df[i,j]), 1, ifelse(!is.na(freq2Dp[i,j]) & is.na(freq2Df[i,j]), 2, ifelse(is.na(freq2Dp[i,j]) & is.na(freq2Df[i,j]), NA, 0)))
    }}

  # plot climate space
  Freq2Dpf <- rbind( data.frame(x=x.bin, y=rep(y.bin, each=length(x.bin)), freq=c(freq2Dp)),
                     data.frame(x=x.bin, y=rep(y.bin, each=length(x.bin)), freq=c(freq2Df)) )
  Freq2Dpf$clim <- factor(rep(c("Historical","Present"), each=nrow(Freq2Dpf)/2), levels=c("Historical", "Present"))
  Freq2Dpf <- Freq2Dpf[!is.na(Freq2Dpf$freq),]

  Freq2D <- factor(c("Novel climate", "Remnant climate", "Disappearing climate")[c(freq2D)+1],
                   levels=c("Novel climate", "Remnant climate", "Disappearing climate"))
  Freq2D <- data.frame(x=x.bin, y=rep(y.bin, each=length(x.bin)), freq=Freq2D)
  Freq2D <- Freq2D[!is.na(Freq2D$freq),]

  panelAB <- ggplot(Freq2Dpf, aes(x=x, y=y, fill=freq)) + geom_raster() + scale_fill_gradientn(
    colors = r, name="Cell count", breaks=seq(0,UL,20), guide=guide_colorbar(ticks.linewidth=1.5)) + facet_wrap(~clim, scale="free_y",
                                                                                                                ncol=2) + scale_x_continuous(limits=c(min(x.bin), max(x.bin))) + scale_y_continuous(limits=c(min(y.bin), max(y.bin))) + labs(
                                                                                                                  x=x.name, y=y.name) + theme(legend.position="bottom",legend.justification=c(1,0), legend.key.width=unit(2, "cm"), strip.text=element_blank())

  panelC <- ggplot(Freq2D, aes(x=x, y=y, fill=freq)) + geom_raster() + scale_fill_manual(values=c("#56B4E9", "#009E73", "#D55E00"),
                                                                                         name="Climate type") + labs(x=x.name, y=y.name)

  panels <- cowplot::plot_grid(panelAB, panelC, nrow=2, rel_heights=c(1.3, 1.0))
  return(panels)
}




