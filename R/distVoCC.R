#' Calculate velocity of climate change using the nearest distance method
#'
#' @usage distVoCC(pre, post, tdiff, tol, denom = 1,
#'   distfun = st_distance, ...)
#'
#' @param pre \code{raster} giving raster of initial climate conditions
#' @param post \code{raster} giving raster of future climate conditions
#' @param tdiff \code{numeric} giving time difference between pre and post
#' @param tol \code{numeric} giving tolerance for classification of
#' analogous climate conditions.
#' @param denom \code{numeric} giving time giving denominator for
#' conversion of distance units (default is metres)
#' @param distfun \code{function} that will calculate distances between
#' analogous climate conditions.
#' @param ... \code{function} other parameters passed to \code{distfun}
#'
#' @return \code{raster} with the velocity of climate change.
#' Default units are metres per unit time (in units of \code{tdiff}).
#'
#' @details This algorithm is based on the one described in
#' Hamann et al. 2015 Velocity of climate change algorithms for
#' guiding conservation and management. Global Change Biology 21: 997-1004
#' The primary difference is that we use spatial data structures from
#' \code{raster} and \code{sf} in the calculations, we use the \code{cut}
#' function to identify analogous climate conditions, whereas Hamann et al.
#' used rounding, and we allow flexible specification of the calculation
#' of distances.
#' As a consequence you must specify the projection of input raster layers.
#' \code{pre} and \code{post} must be matching rasters and
#' must contain CRS details (i.e. \code{!is.na(proj4string(pre))}).
#'
#' The default distance is the Euclidean distance or, in the case of
#'  unprojected coordinates the great circle distance calculated using the
#' \code{st_distance()} function.
#' User defined distance functions should take \code{sf} points objects for
#' the historical and future conditions as the first two arguments
#' respectively.
#'
#' @author Christopher J. Brown
#' @examples
#' data(HSST)
#'
#' tol <- 0.1
#' pre <- raster(HSST,1)
#' post <- raster(HSST, 50*12)
#' tdiff <- 50
#' units <- 1000 #convert to km
#' rvocc <- distvocc(pre, post, tdiff, tol, denom = units)
#' plot(rvocc)
#' @rdname distVoCC
#' @export

distVoCC <- function(pre, post, tdiff, tol, denom = 1, distfun = st_distance, ...){

	names(pre) <- "prevar"
	names(post) <- "postvar"

	datpre <- raster::rasterToPoints(pre, spatial = T) %>%
		sf::st_as_sf()
	datpost <- raster::rasterToPoints(post, spatial = T) %>%
	sf::st_as_sf()

	dat <- sf::st_join(datpre, datpost)
	ymin <- min(c(datpre$prevar, datpre$postvar))
	ymax <- max(c(datpre$prevar, datpre$postvar))
	breaks <- seq(ymin-tol, ymax+tol, by = tol)

	dat$prefact <- cut(dat$pre, breaks = breaks, labels = FALSE)
	dat$postfact <- cut(dat$post, breaks = breaks, labels = FALSE)

	unq <- unique(dat$prefact)
	xout <- purrr::map(unq, ~.getmindist(.x, dat, distfun = st_distance, ...))

	dvocc <- do.call("rbind", xout)
	dvocc$vocc <- (dvocc$dist/denom)/tdiff
	spvocc <- as(dvocc, "Spatial")
	r <- raster::rasterize(spvocc, pre, field = "vocc")
	return(r)
	}
