#' @rdname .getmindist
.getmindist <- function(inum, dat, distfun = sf::st_distance, ...){
	dtemppre <- subset(dat, prefact == inum)
	dtemppost <- subset(dat, postfact == inum)

	if (nrow(dtemppost) == 0){
		dtemppre$dist <- Inf
	} else {
		dmat <- distfun(dtemppre, dtemppost, ...)
		dtemppre$dist <- apply(dmat, 1, min)
	}
	dtemppre
}
