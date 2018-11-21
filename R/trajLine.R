#' Climate velocity trajectory spatial lines
#'
#' Create a spatial line data frame object from trajectory points.
#'
#' @usage trajLine(x)
#'
#' @param x \code{data.frame} containing the coordinates (x, y) of the constituent
#' points and identification number (trajIDs) for each trajectory as returned by VoCCTraj.
#'
#' @return A \code{SpatialLinesDataFrame} with one line per trajectory as specified in x.
#' To avoid artifacts, trajectories crossing the date line need to be splitted into two segments.
#' Where the trajectory on one side of the date line is only composed of a single point,
#' the trajectory won't be displayed (no line object created). The function assumes
#' a -180 to 180 longitudinal arrengment.
#'
#' @seealso{\code{\link{voccTraj}}}
#'
#' @export
#' @author Jorge Garcia Molinos
#' @examples
#'
#' data(HSST)
#' yrSST <- sumSeries(HSST, p = "1969-01/2009-12", yr0 = "1955-01-01", l = nlayers(HSST),
#' fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years")
#' tr <- tempTrend(yrSST, th = 10)
#' sg <- spatGrad(yrSST, th = 0.0001, projected = FALSE)
#' v <- gVoCC(tr,sg)
#' vel <- v[[1]]
#' ang <- v[[2]]
#'
#' # calculate the annual SST mean over the period
#' mn <- mean(yrSST, na.rm = T)
#'
#' # get the set of starting cells for the trajectories
#' lonlat <- na.omit(data.frame(xyFromCell(vel, 1:ncell(vel)), vel[], ang[], mn[]))[,1:2]
#'
#' # Calculate trajectories.
#' traj <- voccTraj(lonlat, vel, ang, mn, tyr = 50)
#'
#' # create a spatial line data frame from traj
#' lns <- trajLine(traj)
#' plot(mn)
#' plot(lns, add = TRUE)
#'
#' # Export as ESRI shape file
#'
#' rgdal::writeOGR(lns, ".", "velTraj", "ESRI Shapefile")
#'
#' @rdname trajLine

trajLine <- function (x){

  spl <- split(x, x$trajIDs)
  # remove traj consisiting of a single point
  i <- sapply(spl, function(x) {nrow(x) == 1})
  spl <- spl[!i]

  lns <- vector("list", length(spl))
  for (i in 1:length(spl)) {
        s <- which(abs(diff(coordinates(spl[[i]][,1:2]))) > 180)
        if(length(s) > 0){
        SPL <- split(as.data.frame(coordinates(spl[[i]][,1:2])), cumsum(1:nrow(coordinates(spl[[i]][,1:2])) %in% (s+1)))
        lns[[i]] <- Lines(lapply(SPL, function(x) Line(coordinates(x))), ID = i)
        }else{
        lns[[i]] <- Lines(list(Line(coordinates(spl[[i]][,1:2]))), ID = i)
     }}
  SpatialLinesDataFrame(SpatialLines(lns, proj4string = CRS(proj4string(vel))), data = data.frame(trajIDs = unique(x$trajIDs)))
}



