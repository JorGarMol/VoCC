#' Internal. Angle associated to the spatial gradient
#' @param dx \code{numeric} giving the longitudinal gradient component
#' @param dy \code{numeric} giving the latitudinal gradient component
#' @author Jorge Garcia Molinos and David S. Schoeman
#' angulo()

angulo <- function(dx, dy){
		d <- cbind(dx, dy)
			angline <- function(rw){
				angle <- ifelse(rw[2] < 0, 180 + CircStats::deg(atan(rw[1]/rw[2])),
			 	ifelse(rw[1] < 0, 360 + CircStats::deg(atan(rw[1]/rw[2])), CircStats::deg(atan(rw[1]/rw[2]))))
				return(angle)
				}
		return(apply(d, 1, angline))
		}















