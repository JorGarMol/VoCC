#' Climatic residence time of a polygon
#'
#' Function to calculate VoCC-based residence time of isotherms within a polygon after Loaire et al. (2009)
#'
#' @usage resTime(pg, vel, areapg = NA)
#'
#' @param pg \code{SpatialPolygon} or a \code{SpatialPolygonsDataFrame} containing the polygons for which
#' the residence time is to be calculated. The polygons must be on the same coordinate system as vel.
#' @param vel \code{raster} with climate velocity (km/year) for the period of interest.
#' @param areapg \code{vector} with the area (in km2) of the polygons. Use NA (default) to calculate internally if field not avilable.
#'
#' @return a \code{data.frame} containing for each polygon its ID, mean velocity (km/yr),
#' diameter of the equivalent circle (km), and residence time (years) as the ratio D/vel.
#'
#' @references \href{https://www.nature.com/articles/nature08649}{Loarie et al. 2009}. The velocity of climate change. Nature, 462, 1052-1055.
#'
#' @seealso{\code{\link{gVoCC}}}
#'
#' @import data.table sp
#' @export
#' @author Jorge Garcia Molinos
#' @examples
#'
#' # Load example Exclusive Economic Zone polygon
#'
#' ?EEZ
#' ?HSST
#' yrSST <- sumSeries(HSST, p = "1969-01/2009-12", yr0 = "1955-01-01", l = nlayers(HSST),
#' fun = function(x) colMeans(x, na.rm = TRUE),
#' freqin = "months", freqout = "years")
#' tr <- tempTrend(yrSST, th = 10)
#' sg <- spatGrad(yrSST, th = 0.0001, projected = FALSE)
#' v <- gVoCC(tr,sg)
#' vel <- v[[1]]
#' # Calculating area internally
#' a1 <- resTime(EEZ, vel, areapg = NA)
#' a1
#' # Using the area field from the polygon data table
#' a2 <- resTime(EEZ, vel, areapg = as.numeric(as.numeric(levels(EEZ$Area_km2))[EEZ$Area_km2]))
#' a2
#' # Using a user defined polygon
#' x_coord <- c(-28, -20, -20.3, -25.5)
#' y_coord <- c(60, 61, 63, 62)
#' p <- Polygon(cbind(x_coord, y_coord))
#' sps <- SpatialPolygons(list(Polygons(list(p),1)))
#' a3 <- resTime(sps, vel, areapg = NA)
#'
#' plot(vel)
#' plot(EEZ, add = TRUE)
#' plot(sps, add = TRUE)
#'
#' @rdname resTime


resTime <- function(pg, vel, areapg = NA){
RT <- suppressWarnings(data.table(ID = sp::getSpPPolygonsIDSlots(pg)))
RT[, v := raster::extract(vel, pg, small=TRUE, fun=mean, na.rm=TRUE)]
# If not provided, calculate the area of the polygon
if(all(is.na(areapg))){
areapg <- geosphere::areaPolygon(pg)/1000000              # area polygon in km2
}
RT[, d := 2*sqrt((areapg/pi))]
RT[, resTim := abs(d/v)]
return(RT[])
}





