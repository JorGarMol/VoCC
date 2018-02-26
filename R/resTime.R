#' Climatic residence time of a polygon
#'
#' Function to calculate vocc residence time of a polygon after Loaire et al. (2009)
#'
#' @param pg \code{SpatialPolygon} or a \code{SpatialPolygonsDataFrame} containing the polygons for which
#' the residence time is to be calculated. The polygons must be on the same coordinate system as vel.
#' @param vel \code{raster} with local climate velocity (km/year) for the period of interest.
#' @param areapg \code{vector} with the area (in km2) of the polygons. Use NA (default) to calculate internally if field not avilable.
#'
#' @return a \code{data.frame} containing for each polygon its ID, mean velocity (km/yr),
#' diameter of the equivalent circle (km), and residence time (years) as the ratio D/vel.
#'
#' @import data.table
#' @export
#' @author Jorge Garcia Molinos
#' @examples
#'
#' # Load example Exclusive Economic Zone polygons
#'
#' data(EEZs)
#' # Calculating area internally
#' a1 <- resTime(EEZs, vel, areapg = NA)
#' a1
#' # Using the area field from the polygon data table
#' a2 <- resTime(EEZs, vel, areapg = EEZs$Area_km2)
#' a2
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





