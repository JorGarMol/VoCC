#' Climate velocity trajectories
#'
#' Function to calculate vocc trajectories after Burrows et al (2014). Trajectories
#' are calculated by propagating climatic isopleths using the magnitude and direction of
#' local (cell) velocities. This is a slightly modified version of the original Burrows et al. (2014)
#' approach in that iterations of a trajectory are based on cumulative time travelled instead of using fixed time steps.
#'
#' @usage voccTraj(lonlat, vel, ang, mn, tyr, trajID = 1:nrow(lonlat), correct = FALSE)
#'
#' @param lonlat \code{data.frame} with the longitude and latitude (in decimal degrees)
#' of the points to project.
#' @param vel \code{raster} with the magnitude of gradient-based climate velocity.
#' @param ang \code{raster} with velocity angles in degrees.
#' @param mn \code{raster} with the overall mean climatic value over the period of interest.
#' @param tyr \code{integer} temporal length of the period of interest.
#' @param trajID \code{integer} specifying the identifiers for the trajectories.
#' @param correct \code{logical} does the input raster need to be corrected to account for cropped margins?
#' Unless the raster extent is global, calculation of trajectories will throw an error at the margins
#' as the trajectories go beyond the raster extent (no input values). To avoid this, an option is given for
#' expanding the extent by the resolution of the raster (1 column/row) with NAs. Note that those trajectories
#' reaching the extent limits will be artificially bounced back so should be discarded at that point.
#' Alternatively, users may choose to crop to a larger extent to the domain of interest (appropriately
#' defined by lonlat), so the extra extent buffer for those trajectories getting to the border
#' of the raster.
#'
#' @return a \code{data.frame} containing the coordinates ("x", "y") of the constituent
#' points and identification number ("trajIDs") for each trajectory.
#'
#' @references \href{https://www.nature.com/articles/nature12976}{Burrows et al. 2014}. Geographical limits to species-range shifts are suggested by climate velocity. Nature, 507, 492-495.
#'
#' @seealso{\code{\link{gVoCC}}, \code{\link{trajClas}}}
#' @importFrom geosphere destPoint distGeo
#' @export
#' @author Jorge Garcia Molinos, David S. Schoeman and Michael T. Burrows
#' @examples
#'
#' yrSST <- sumSeries(HSST, p = "1969-01/2009-12", yr0 = "1955-01-01", l = nlayers(HSST),
#' fun = function(x) colMeans(x, na.rm = TRUE),
#' freqin = "months", freqout = "years")
#' tr <- tempTrend(yrSST, th = 10)
#' sg <- spatGrad(yrSST, th = 0.0001, projected = FALSE)
#' v <- gVoCC(tr,sg)
#' vel <- v[[1]]
#' ang <- v[[2]]
#'
#' # Calculate the annual SST mean over the period
#' mn <- mean(yrSST, na.rm = T)
#'
#' # Get the set of starting cells for the trajectories
#' lonlat <- na.omit(data.frame(xyFromCell(vel, 1:ncell(vel)), vel[], ang[], mn[]))[,1:2]
#'
#' # Calculate trajectories
#' # The following throws an error due to the trajectories moving beyond the raster extent
#' traj <- voccTraj(lonlat, vel, ang, mn, tyr = 50)
#'
#' # This accounts for the extent issue
#' traj <- voccTraj(lonlat, vel, ang, mn, tyr = 50, correct = TRUE)
#'
#' @rdname voccTraj


voccTraj <- function(lonlat, vel, ang, mn, tyr, trajID = 1:nrow(lonlat), correct = FALSE){

  if(correct == TRUE){
    e <- extent(vel)+(rep(res(vel), each = 2)*c(-1,1,-1,1))
    vel <- extend(vel, e, value = NA)
    ang <- extend(ang, e, value = NA)
    mn <- extend(mn, e, value = NA)
  }

# make sure all raster has consistent NAs
  vel[is.na(ang)] <- NA
  vel[is.na(mn)] <- NA
  ang[is.na(vel)] <- NA
  mn[is.na(vel)] <- NA
# Set up variables to catch results, allocating the right amount of memory
  nc <- nrow(lonlat)
  remaining <- rep(tyr, nc) # String containing the time remaining for each trajectory
  llon <- rep(NA, (nc * tyr) + nc)  # Starting lons, plus one more set for each iteration
  llat <- rep(NA, (nc * tyr) + nc)
  rem <- rep(NA, (nc * tyr) + nc) # this is to keep track of trajectories trapped in internal sinks
# populate the first n slots with starting points
	llon[1:nc] <- lonlat[,1]
	llat[1:nc] <- lonlat[,2]
  rem[1:nc] <- remaining
  i <- 0      # set the iteration counter
  land <- "N"  # set to not hit land
  bounce <- "N"  # set to not bounced
# Calculate the trajectories
  pb <- utils::txtProgressBar(min = 0, max = 100, style = 3)
  while(sum(remaining <= 0) != nc){  # while there is at least one trajectory active
    utils::setTxtProgressBar(pb, 100*(sum(remaining <= 0)/nc))

    llold <- lonlat # Take a copy of lonlat (starting cell xy)
    resto <- which(remaining > 0)    # index with remaining active trajectories
    fcells <- raster::cellFromXY(vel, llold)     # focal cells
    # Extract lon and lat of landing point for the remaining active trajectories
    # limit the displacement to 2 cell lengths to reduce later the number of intermediate points
    dth <- max(res(vel))*222000
    dis <- ifelse(dth < (abs(vel[fcells[resto]])*remaining[resto]*1000), dth,(abs(vel[fcells[resto]])*remaining[resto]*1000))
    lonlat[resto,] <- as.data.frame(destPoint(llold[resto,], ang[fcells[resto]], dis)) # distance input in meters. The function adjusts internally for -180-180 and pole crossings
    tcells <- raster::cellFromXY(vel, lonlat)

# Step 1. where the trajectory is still in the same cell by tyr, it has terminated.
      # Flag those trajectories by resetting the reminding time to 0 to get them out of the next iteration.
      remaining[fcells == tcells] <- 0
      if(sum(remaining == 0) == nc){break}   # to avoid error when only 1 trajectory is left and finishes inside a cell

# Step 2. For the rest, get the last point in the starting cell and the first point in a new cell
      resto <- which(remaining > 0)     # update resto
      d <- round((distGeo(llold[resto,], lonlat[resto,])/1000), 0)
      d[d == 0] <- 1
      Trajxy = splitLine(A = llold[resto,], B = lonlat[resto,], n = d)
      # now get a list with the cells for the points in each trajectory
      if(is.list(Trajxy)){
      Trajcells <- lapply(Trajxy, cellFromXY, object = vel)
      # the first new cell in the trajectory
      index <- lapply(Trajcells, function(x) which(x != x[1])[1])
      newcell <- mapply(function(X,Y) {X[Y]}, X = Trajcells, Y = index)
      # the coordinates for that first point out of the focal cell
      newxy <- as.data.frame(t(mapply(function(X,Y) {X[Y,]}, X=Trajxy, Y=index)))
      # the coordinates for the last focal cell point
      oldxy <- as.data.frame(t(mapply(function(X,Y) {X[Y-1,]}, X=Trajxy, Y=index)))
      }else{       # if 1 trajectory the output from splitLine is a single matrix instead of a list
      Trajcells <- apply(Trajxy, 1, cellFromXY, object = vel)
      newcell <- Trajcells[which(Trajcells != Trajcells[1])[1]]
      newxy <- data.frame(x = Trajxy[which(Trajcells != Trajcells[1])[1],1], y = Trajxy[which(Trajcells != Trajcells[1])[1],2])
      oldxy <- data.frame(x = Trajxy[(which(Trajcells != Trajcells[1])[1])-1,1], y = Trajxy[(which(Trajcells != Trajcells[1])[1])-1,2])
      }
      # starting and destination cell ids
      oldcell <- fcells[resto]
      # Get the new velocity at the new cells
      velend <- raster::extract(vel, newcell)
      # set remaining time for each of the running trajectories
      remaining[resto][is.na(velend)] <- remaining[resto][is.na(velend)]-((distGeo(llold[resto,][is.na(velend),], oldxy[is.na(velend),])/1000)/abs(vel[fcells[resto]][is.na(velend)]))
      remaining[resto][!is.na(velend)] <- remaining[resto][!is.na(velend)]-((distGeo(llold[resto,][!is.na(velend),], newxy[!is.na(velend),])/1000)/abs(vel[fcells[resto]][!is.na(velend)]))
      # For those ending in marine cells update lonlat info to the new point
      lonlat[resto,][!is.na(velend),] <- newxy[!is.na(velend),]
      # For those ending in land cells update lonlat info to the last marine point
      lonlat[resto,][is.na(velend),] <- oldxy[is.na(velend),]

# Step 3. From step 3 onwards check if the trajectory has bounced back to the origin cell. If so, redirect along cell border.
     if(i >= 1){
      current <- newcell[!is.na(velend)]   # current cell
      last <- oldcell[!is.na(velend)] # last cell
      if(land == "Y" & bounce == "Y"){ # Given the counter increases each time the trajectory table is updated,
      # where some traj hit land AND bounced in the previous iteration the values were repeated twice for all cells. Hence, need to go 3 steps back instead of 1
      last2 <- raster::cellFromXY(vel, cbind(llon[(((i-3) * nc) + 1):(((i-3) * nc) + nc)][resto][!is.na(velend)], llat[(((i-3) * nc) + 1):(((i-3) * nc) + nc)][resto][!is.na(velend)]))  # last but one cell
      land <- "N"  # set back to no hit land
      bounce <- "N"
      }else if(xor(land == "Y", bounce == "Y")){  # if one of the two conditions then need to go two steps back
      last2 <- raster::cellFromXY(vel, cbind(llon[(((i-2) * nc) + 1):(((i-2) * nc) + nc)][resto][!is.na(velend)], llat[(((i-2) * nc) + 1):(((i-2) * nc) + nc)][resto][!is.na(velend)]))
      land <- "N"
      bounce <- "N"
      }else{    # if non of the two, then just one step back
      last2 <- raster::cellFromXY(vel, cbind(llon[(((i-1) * nc) + 1):(((i-1) * nc) + nc)][resto][!is.na(velend)], llat[(((i-1) * nc) + 1):(((i-1) * nc) + nc)][resto][!is.na(velend)]))
      }
      # Identify the bouncing trajectories
      ind <- which(current == last2 & current != last)
      if(length(ind) > 0){
      bounce <- "Y"
      # take the mean velocity from the appropriate lat/lon vel components.
      vb <- mapply(function(X,Y) {ifelse(abs(X-Y) > 1, mean(c((vel[X]*sin(pi*ang[X]/180)),(vel[Y]*sin(pi*ang[Y]/180)))),mean(c((vel[X]*cos(pi*ang[X]/180)),(vel[Y]*cos(pi*ang[Y]/180)))))},
      X = current[ind], Y = last[ind])
      # take the corresponding angle (0/180 / 90/270 for lat / lon movements)
      ab <- mapply(function(X,Y,Z) {ifelse(abs(X-Y) > 1 & Z > 0, 90, ifelse(abs(X-Y) > 1 & Z < 0, 270, ifelse(abs(X-Y) == 1 & Z > 0, 0, 180)))}, X = current[ind], Y = last[ind], Z = vb)
      if(is.matrix(ab)){ab <- ab[,1]}
      # Extract lon and lat of point previous to bounce and new destination point for the remaining active trajectories
      p <- oldxy[!is.na(velend),][ind,]
      dis <- ifelse(dth < (abs(vb)*remaining[resto][!is.na(velend)][ind]*1000), dth, (abs(vb)*remaining[resto][!is.na(velend)][ind]*1000))
      destp <- as.data.frame(destPoint(p, ab, dis))
      # where the trajectory is still in the same cell by tyr, it has terminated. Flag those trajectories by resetting the reminding time to 0 to get them out of the next iteration.
      same <- which(cellFromXY(vel, p) == cellFromXY(vel, destp))
      remaining[resto][!is.na(velend)][ind][same] <- 0
      # For the rest, get the correct destination point
      rest <- which(cellFromXY(vel, p) != cellFromXY(vel, destp))     # update rest
      if(length(rest)>0){
      d <- round((distGeo(p[rest,], destp[rest,])/1000), 0)
      d[d == 0] <- 1
      Trajxy = splitLine(A = p[rest,], B = destp[rest,], n = d)
      # now get a list with the cells for the points in each trajectory
      if(is.list(Trajxy)){
      Trajcells <- lapply(Trajxy, cellFromXY, object = vel)
      # the first new cell in the trajectory
      index <- lapply(Trajcells, function(x) which(x != x[1])[1])
      # update the coordinates for that first point out of the focal cell
      newxy[!is.na(velend),][ind,][rest,] <- as.data.frame(t(mapply(function(X,Y) {X[Y,]}, X=Trajxy, Y=index)))
      # update points info for the old position in case the trajectory hit land and need to be repositioned
      i <- i+1
      llon[((i * nc) + 1):((i * nc) + nc)] <- lonlat[,1]
      llat[((i * nc) + 1):((i * nc) + nc)] <- lonlat[,2]
      llon[((i * nc) + 1):((i * nc) + nc)][resto][!is.na(velend)][ind][rest] <- oldxy[!is.na(velend),][ind,][rest,1]
      llat[((i * nc) + 1):((i * nc) + nc)][resto][!is.na(velend)][ind][rest] <- oldxy[!is.na(velend),][ind,][rest,2]
      # update the coordinates for the last focal cell point
      oldxy[!is.na(velend),][ind,][rest,] <- as.data.frame(t(mapply(function(X,Y) {X[Y-1,]}, X=Trajxy, Y=index)))
      }else{       # if all the trajectories have same number of points the output from gcIntermediate is a matrix instead of a list
      Trajcells <- apply(Trajxy, 1, cellFromXY, object = vel)
      newxy[!is.na(velend),][ind,][rest,] <- data.frame(x = Trajxy[which(Trajcells != Trajcells[1])[1], 1], y = Trajxy[which(Trajcells != Trajcells[1])[1],2])
      i <- i+1
      llon[((i * nc) + 1):((i * nc) + nc)] <- lonlat[,1]  # necessary to not leave the other cells as NAs
      llat[((i * nc) + 1):((i * nc) + nc)] <- lonlat[,2]
      llon[((i * nc) + 1):((i * nc) + nc)][resto][!is.na(velend)][ind][rest] <- oldxy[!is.na(velend),][ind,][rest,1]
      llat[((i * nc) + 1):((i * nc) + nc)][resto][!is.na(velend)][ind][rest] <- oldxy[!is.na(velend),][ind,][rest,2]
      oldxy[!is.na(velend),][ind,][rest,] <- data.frame(x = Trajxy[(which(Trajcells != Trajcells[1])[1])-1,1], y = Trajxy[(which(Trajcells != Trajcells[1])[1])-1,2])
      }
      # Update main traj info
      # get remaining time for each of the running trajectories
      remaining[resto][!is.na(velend)][ind][rest] <- remaining[resto][!is.na(velend)][ind][rest]-((distGeo(newxy[!is.na(velend),][ind,][rest,], p[rest,])/1000)/abs(vb[rest]))
      rem[((i * nc) + 1):((i * nc) + nc)] <- remaining
      lonlat[resto,][!is.na(velend),][ind,][rest,] <- newxy[!is.na(velend),][ind,][rest,]
      # update the velocity at new cells (some of the redirected bouncing trajectories might have ended in a land cell)
      newcell[!is.na(velend)][ind][rest] <- cellFromXY(vel, newxy[!is.na(velend),][ind,][rest,])
      velend[!is.na(velend)][ind][rest] <- raster::extract(vel, cellFromXY(vel, newxy[!is.na(velend),][ind,][rest,]))
      }}
     }

# Step 4. For those traj ending on land redirect the trajectories if possible
     if(sum(is.na(velend)) > 0){
      land <- "Y"      # to know if land was hit in the next loop when looking at bouncing trajectories
      onland <- which(is.na(velend))  # Identify which rows of velend are on land
      # break and produce error if trajectories go beyond raster extent
      if(anyNA(newcell[onland])) {stop('Trajectories beyond raster extent - set correct = TRUE or limit lonlat to a sufficiently smaller domain relative to current raster extent')}
      # For each cell that ends on land, look for a suitable target cell...
      # Make list of candidate cell IDs: overlap between cells adjacent to "from" (fcells[is.na(sflags)][onland]) AND "to" (newcell[onland])cells
      # that are in the direction of movement (any other cells would mean "un-natural" movements).
      ccells <- suppressWarnings(mapply(function(X, Y) {
      Z <- as.numeric(intersect(adjacent(vel, X, directions = 8), adjacent(vel, Y, directions = 8)))
      Z[!is.na(vel[Z])]
      }, X = oldcell[onland], Y = newcell[onland], SIMPLIFY = FALSE))

      # Now check if a suitable cell, other than the focal cell, is available along the line of movement
      # given departure direction. Diagonals included.
      # First, calculate the velocity associated to each focal cell with which to move along the cell border given direction
      # if transfer is vertical (upper/lower cell) use horizontal velocity component, else (horizontal) use vertical component
      # negative values indicate E-W and N-S movements
      v <- mapply(function(X,Y) {ifelse(abs(X-Y) > 1, abs(vel[Y])*sin(pi*ang[Y]/180), abs(vel[Y])*cos(pi*ang[Y]/180))}, X = newcell[onland], Y = oldcell[onland])
      a <- rep(NA, length(v))  # to store angles for border movement

      for(k in 1:length(v)){
      target <- newcell[onland][k]
      focal <- oldcell[onland][k]
      if((target-focal) > 1 & v[k] > 0){    # vertical transfer, horziontal movement
      ccells[[k]] <- subset(ccells[[k]], ccells[[k]] %in% c(focal+1,focal+ncol(vel)+1))
      a[k] <- 90
      }else if((target-focal) > 1 & v[k] < 0){
      ccells[[k]] <- subset(ccells[[k]], ccells[[k]] %in% c(focal-1,focal+ncol(vel)-1))
      a[k] <- 270
      }else if((target-focal) < -1 & v[k] > 0){
      ccells[[k]] <- subset(ccells[[k]], ccells[[k]] %in% c(focal+1,focal-ncol(vel)+1))
      a[k] <- 90
      }else if((target-focal) < -1 & v[k] < 0){
      ccells[[k]] <- subset(ccells[[k]], ccells[[k]] %in% c(focal-1,focal-ncol(vel)-1))
      a[k] <- 270
      } else if((target-focal) == 1 & v[k] > 0){    # horizontal transfer, vertical movement
      ccells[[k]] <- subset(ccells[[k]], ccells[[k]] %in% c(focal-ncol(vel),focal-ncol(vel)+1))
      a[k] <- 0
      }else if((target-focal) == 1 & v[k] < 0){
      ccells[[k]] <- subset(ccells[[k]], ccells[[k]] %in% c(focal+ncol(vel),focal+ncol(vel)+1))
      a[k] <- 180
      }else if((target-focal) == -1 & v[k] > 0){
      ccells[[k]] <- subset(ccells[[k]], ccells[[k]] %in% c(focal-ncol(vel),focal-ncol(vel)-1))
      a[k] <- 0
      }else if((target-focal) == -1 & v[k] < 0){
      ccells[[k]] <- subset(ccells[[k]], ccells[[k]] %in% c(focal+ncol(vel),focal+ncol(vel)-1))
      a[k] <- 180
      }}

      # Flag out the cells for which no suitable potential neighbours are available (all potential neighbours are NA; the trajectory has nowhere to go)
      empty <- as.logical(lapply(ccells, function(x) identical(x, numeric(0))))
      remaining[resto][onland][empty] <- 0
      lonlat[resto,][onland,][empty,] <- oldxy[onland,][empty,]
      # if there are cells with available neighbours
      if(sum(!empty) > 0){
      # select those cells
      leftcell <- ccells[!empty]
      focal <- oldcell[onland][!empty]
      target <- newcell[onland][!empty]
      nxy <- newxy[onland,][!empty,]
      oxy <- oldxy[onland,][!empty,]
      vleft <- v[!empty]
      aleft <- a[!empty]
      # Check if there is an actual suitable neighbour
      sst <- lapply(leftcell, function(x) mn[x])   # sst at focal and neighbouring cells
      # if warming (cooling), need a suitable neighbour with cooler (warmer) local temperatures
      for(k in 1:length(leftcell)){
      leftcell[k] <- ifelse(vel[focal[k]]>0 & sum(sst[[k]] < mn[focal[k]]) > 0, leftcell[[k]][which.min(sst[[k]])],
      ifelse(vel[focal[k]]<0 & sum(sst[[k]] > mn[focal[k]]) > 0, leftcell[[k]][which.max(sst[[k]])], NA))
      }

      # Flag cells for which no suitable neighbours are available (no suitable sst)
      remaining[resto][onland][!empty][is.na(as.numeric(leftcell))] <- 0
      lonlat[resto,][onland,][!empty,][is.na(as.numeric(leftcell)),] <- oldxy[onland,][!empty,][is.na(as.numeric(leftcell)),]

      # For those traj with a suitable neighbour, relocate the trajectory to the nearest newcell point and calculate the time needed to reach it
      if(length(target[!is.na(as.numeric(leftcell))])>0){
      dudcent <- xyFromCell(vel, target[!is.na(as.numeric(leftcell))]) # Coordinates of dud target cell on land
		  newcent <- xyFromCell(vel, as.numeric(leftcell[!is.na(as.numeric(leftcell))])) # Coordinates of new target cell
		  oldcent <- xyFromCell(vel, focal[!is.na(as.numeric(leftcell))]) # Coordinates of departure cell
		  loncent <- data.frame(oldx = oldcent[,1], NAx = dudcent[,1], newx = newcent[,1]) # An object containing the longitudes of the cells involved - needed for fixing dateline

      for(k in 1:nrow(loncent)){
      if(abs(max(loncent[k,]) - min(loncent[k,])) > 180){loncent[k,][loncent[k,] < 0] <- loncent[k,][loncent[k,] < 0] + 360}  # Remove sign change for dateline, if needed
      loncent$lonline[k] <- ifelse(abs(loncent[k,2] - loncent[k,3]) > abs(loncent[k,2] - loncent[k,1]),    # if lon difference between NA cell and new cell is larger than NA to old cell
			mean(c(loncent[k,1], loncent[k,3])), mean(c(loncent[k,1], loncent[k,2]))) # Figure out position of dividing longitude
      loncent$lonline[k] <- loncent$lonline[k] - (360 * floor((loncent$lonline[k] + 180) / 360)) # Return to -180o to 180o format
      loncent$lonline[k] <- ifelse(loncent$lonline[k] == 180 && newcent[k,1] < 0, -180, loncent$lonline[k])
      loncent$lonline[k] <- ifelse(loncent$lonline[k] == -180 && newcent[k,1] > 0, 180, loncent$lonline[k]) # Arrange lonline to accommodate dateline
      loncent$latline[k] <- ifelse(abs(dudcent[k,2] - newcent[k,2]) > abs(dudcent[k,2] - oldcent[k,2]), mean(c(oldcent[k,2], newcent[k,2])), mean(c(oldcent[k,2], dudcent[k,2])))
      loncent$dirlon[k] <- ifelse(newcent[k,1] > loncent$lonline[k], 1, -1) # Adjust direction of movement relative to lonline
	    loncent$dirlat[k] <- ifelse(newcent[k,2] > loncent$latline[k], 1, -1) # Same for latline
      loncent$lonnew[k] <- loncent$lonline[k] + (loncent$dirlon[k] * 0.0001) # add a small offset to put the traj on the new cell
      loncent$latnew[k] <- loncent$latline[k] + (loncent$dirlat[k] * 0.0001)
      loncent$latnew[k] <- ifelse(loncent$latnew[k] > 90, 90, ifelse(loncent$latnew[k] < -90, -90, loncent$latnew[k]))
      loncent$lonnew[k] <- loncent$lonnew[k] - (360 * floor((loncent$lonnew[k] + 180) / 360))
      }

      loncent$oldcell <- focal[!is.na(as.numeric(leftcell))]
      loncent$lonold  <- oxy[!is.na(as.numeric(leftcell)),1]
      loncent$latold  <- oxy[!is.na(as.numeric(leftcell)),2]
      loncent$dis <- (distGeo(oxy[!is.na(as.numeric(leftcell)),], loncent[,8:9])/1000)        # distance from dold to new point
      loncent$dur <- loncent$dis/abs(vleft[!is.na(as.numeric(leftcell))])                 # time taken to get there
      loncent$remaining <- remaining[resto][onland][!empty][!is.na(as.numeric(leftcell))]-loncent$dur

     # where remaining is < 0, the trajectory terminates before reaching new destination. Flag those traj out and place them in the corresponding final coordinates.
     if(sum(loncent$remaining < 0) > 0){
     loncent[loncent$remaining < 0, 8:9] <- as.matrix(destPoint(oxy[!is.na(as.numeric(leftcell)),][loncent$remaining < 0,], aleft[!is.na(as.numeric(leftcell))][loncent$remaining < 0],
     (abs(vleft[!is.na(as.numeric(leftcell))][loncent$remaining < 0])*(remaining[resto][onland][!empty][!is.na(as.numeric(leftcell))][loncent$remaining < 0])*1000)))
     loncent$remaining[loncent$remaining < 0] <- 0
     }

     # Finally, update lonlat
     i <- i+1 # increase the counter by 1
     # first input the point before hitting land
     lonlat[resto,][onland,][!empty,][!is.na(as.numeric(leftcell)),] <- loncent[,11:12]
     llon[((i * nc) + 1): ((i * nc) + nc)] <- lonlat[,1]  # Add final lon to the list
	   llat[((i * nc) + 1): ((i * nc) + nc)] <- lonlat[,2] # Add final lat to the list
     # then the new point to carry the trajectory on with
     remaining[resto][onland][!empty][!is.na(as.numeric(leftcell))] <- loncent$remaining
     rem[((i * nc) + 1):((i * nc) + nc)] <- remaining
     lonlat[resto,][onland,][!empty,][!is.na(as.numeric(leftcell)),] <- loncent[,8:9]
     }}
    }

# Step 5. Update register before moving to the next projection step
    i <- i+1 # increase the counter by 1
    llon[((i * nc) + 1):((i * nc) + nc)] <- lonlat[,1]  # Add final lon to the list
    llat[((i * nc) + 1):((i * nc) + nc)] <- lonlat[,2] # Add final lat to the list
    rem[((i * nc) + 1):((i * nc) + nc)] <- remaining

# Step 6. Check for trajectories trapped in internal sinks and terminate them
    if(i >= 8){
    cs <- which(remaining != 0)
    # index the trajectories trapped in a sink (those for which distance travelled over the last 4 iterations is less than 0.5 km)
    if(length(cs) > 0){
    ind <- which(unlist(lapply(cs, function(x) all(tail(abs(diff(rem[seq(x, by = nc, length.out = i+1)], 1))[abs(diff(rem[seq(x, by = nc, length.out = i+1)], 1)) != 0], 4) < 0.5))))
    # Terminate trapped trajectories
    if(length(ind)>0){
    remaining[cs][ind] <- 0
    }}}

  }
  close(pb)

# prepare output
trajIDs = rep(trajID,(length(llon)/nc)) # rep(initialcells,(length(llon)/nc))
traj <- na.omit(data.frame(x = llon, y = llat, trajIDs = trajIDs))
# remove duplicated points (to keep track of the trajectorie's ID all trajectories are repeated over iterations even if they had terminated)
traj <- traj[!duplicated(traj),]

return(traj)
}
