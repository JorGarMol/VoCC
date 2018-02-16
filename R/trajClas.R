#' Spatial classification based on climate velocity trajectories
#'
#' Function for the spatial classification of cells based on vocc trajectories after Burrows et al (2014).
#'
#' @param traj \code{data.frame} as retuned by voccTraj containing the coordinates 
#' and identification number for each trajectory.
#' @param vel \code{raster} with the magnitude of local climate velocity.
#' @param ang \code{raster} with velocity angles.
#' @param trajSt \code{integer} number of trajectories starting from each cell or spatial unit.
#' @param tyr \code{integer} number of years comprising the period of interest.
#' @param nmL \code{numeric} upper threshold (distance units as per vel object) up to which 
#' a trajectory is considered to have traveled a negiglibe distance over the study period (non-moving).
#' @param smL \code{numeric} upper threshold up to which a trajectory is considered to have traveled a small
#' distance over the study period (slow-moving).
#' @param Nend \code{numeric} the percentage of trajectories ending to be used as threshold in the classification (see Burrows et al. 2014).
#' @param Nst \code{numeric} the percentage of trajectories starting to be used as threshold in the classification (see Burrows et al. 2014).
#' @param NFT \code{numeric} the percentage of trajectories flowing through to be used as threshold in the classification (see Burrows et al. 2014).
#'
#' @return a \code{raster.stack} containing the trajectory classification ("TrajClas"), 
#' as well as those based on trajectory length ("ClassL"; 1 non-moving, 2 slow-moving, 3 fast-moving cells),
#' boundrary ("BounS") and internal sinks ("IntS"), and the proportion of trajectories ending("PropEnd"), 
#' flowing through ("PropFT") and starting ("PropSt"). The trajectory classes ("TrajClas") are (1) non-moving, 
#' (2) slow-moving, (3) internal sinks, (4) boundary sinks, (5) sources, (6) relative sinks, (7) corridors, 
#' (8) divergence and (9) convergence.
#'
#' @examples
#' \dontrun{
#' data(traj25)
#' clas <- trajClas(traj25, vel, ang, trajSt = 16, tyr = 50, nmL = 20, smL = 100, Nend = 45, Nst = 15, NFT = 70)
#' my_col = c('gainsboro','darkseagreen1','coral4','firebrick2','mediumblue','darkorange1','magenta1','cadetblue1','yellow1')
#' plot(clas[[7]], legend = FALSE, col = my_col)
#' legend(x='bottomleft', legend = c("N-M", "S-M", "IS", "BS", "Srce", "RS", "Cor", "Div", "Con"), fill = my_col,horiz = TRUE, cex = 0.7)
#' }
#'
#' @import raster data.table
#' @export
#' @author Jorge Garcia Molinos
#' @rdname trajClas

trajClas <- function(traj, vel, ang, trajSt, tyr, nmL, smL , Nend, Nst, NFT){

TrajEnd <- TrajFT <- TrajSt <- IntS <- BounS <- TrajClas <- raster(ang)

# add cell ID to the data frame                                
traj <- data.table(traj)
traj$cid <- cellFromXY(ang, traj[,1:2])

# A. Number of traj starting from each cell
TrajSt[!is.na(ang[])] <- trajSt

# B. Number of traj ending in each cell
tr <- traj[, .SD[.N], by = trajIDs]  # subset last point of each trajectory
enTraj <- tr[,.N, by = cid]
TrajEnd[!is.na(vel)] <- 0
TrajEnd[enTraj$cid] <- enTraj$N

# C. Number of traj flowing through each cell
cxtrj <- unique(traj, by = c("trajIDs", "cid")) 
TotTraj <- cxtrj[,.N, by = cid]       # total number of trajectories per cell
TrajFT[!is.na(vel)] <- 0
TrajFT[TotTraj$cid] <- TotTraj$N
TrajFT <- TrajFT - TrajEnd - TrajSt   # subtract traj starting and ending to get those actually transversing the cell
TrajFT[TrajFT[] < 0] <- 0   # to avoid negative values in ice covered cells

# C. Identify cell location for internal sinks (groups of 4 endorheic cells with angles pointing inwards)
ll <- data.table(xyFromCell(ang, 1:ncell(ang)))
ll[,1:2] <- ll[,1:2] + 0.1   # add small offset to the centroid
a <- fourCellsFromXY(ang, as.matrix(ll[,1:2]))
a <- t(apply(a, 1, sort))
# correct sequences on date line
a[seq(360, by = 360, length = nrow(ang)),] <- t(apply(a[seq(360, by = 360, length = nrow(ang)),], 1, function(x) {x[c(2,1,4,3)]}))
b <- matrix(extract(ang, as.vector(a)), nrow = ncell(ang), ncol = 4, byrow = FALSE)        # extract the angles for each cell
ll[, c("c1", "c2", "c3", "c4", "ang1", "ang2", "ang3", "ang4") := data.frame(a, b)]
# now look if the 4 angles point inwards (internal sink)
ll[, c("d1", "d2", "d3", "d4") := .(((ang1 - 180)  *  (90 - ang1)), ((ang2 - 270)  *  (180 - ang2)), ((ang3 - 90)  *  (0 - ang3)), ((ang4 - 360)  *  (270 - ang4)))]
ll[, isink := 0L]
ll[d1 > 0 & d2 > 0 & d3 > 0 & d4 > 0, isink := 1L]
# get the cids for the cells contained in the sinks
celdas <- ll[isink == 1, 3:6]
IntS[!is.na(vel)] <- 0
IntS[c(celdas[[1]], celdas[[2]], celdas[[3]], celdas[[4]])] <- 1

# D. Identify cell location for boundary sinks (coastal cells which are disconected from cooler climates under warming or warmer climates under cooling)
# detect coastal cells
coast <- suppressWarnings(boundaries(vel, type = 'inner', asNA = TRUE))       # to avoid warning for coercing NAs via asNA = TRUE
# make a list of vel values and SST values for each coastal cells and their marine neighbours
cc <- na.omit(data.table(cid = 1:ncell(vel), coast = coast[]))
ad <- adjacent(vel, cc$cid, 8, sorted = TRUE, include = TRUE) # matrix with adjacent cells
ad <- data.table(coastal = ad[,1], adjacent = ad[,2], cvel = vel[ad[,1]], ctemp = mn[ad[,1]], atemp = mn[ad[,2]])
# locate the sinks
ad <- na.omit(ad[ad$cvel != 0,])      # remove cells with 0 velocity (ice) and with NA (land neighbours)
j <- ad[, ifelse(.SD$cvel > 0, all(.SD$ctemp <= .SD$atemp), all(.SD$ctemp >= .SD$atemp)), by = coastal]
setkey(j)
j <- unique(j)
BounS[!is.na(vel)] <- 0
BounS[unique(subset(j$coastal, j$V == TRUE))] <- 1

# Total number of trajectories per cell and proportions per cell
TrajTotal <- calc(stack(TrajSt, TrajFT, TrajEnd), sum, na.rm = TRUE)
TrajTotal[is.na(ang[])] <- NA
PropTrajEnd <- (TrajEnd/TrajTotal)*100                          
PropTrajFT <- (TrajFT/TrajTotal)*100
PropTrajSt <- (TrajSt/TrajTotal)*100

# reclassify by traj length
rclM <- matrix(c(0, (nmL/tyr), 1, (nmL/tyr), (smL/tyr), 2, (smL/tyr), Inf, 3), ncol=3, byrow=TRUE) 
v <- raster(vel)
v[] <- abs(vel[])
ClassMov <- reclassify(v, rclM)

# Classify the cells
TrajClas[!is.na(vel)] <- 0
# capture non-moving (1)
TrajClas[ClassMov[] == 1] <- 1
# capture slow-moving (2)
TrajClas[ClassMov[] == 2] <- 2
# capture internal (3) and (4) boundary sinks
TrajClas[IntS[] == 1] <- 3
TrajClas[BounS[] == 1] <- 4
# Classify remaining cells into sources(5), rel sinks(6), corridors(7), divergence(8) and convergence(9)
d <- na.omit(data.table(cid = 1:ncell(TrajClas), val = TrajClas[]))
d <- d[val == 0, 1]
d[, Nend := PropTrajEnd[d$cid]]
d[, Nst := PropTrajSt[d$cid]]
d[, NFT := PropTrajFT[d$cid]]
d$clas <- ifelse(d$Nend == 0, 5, ifelse(d$Nend > Nend & d$Nst < Nst, 6, ifelse(d$NFT > NFT, 7,ifelse(d$Nend < d$Nst, 8, 9))))           
TrajClas[d$cid] <- d$clas

# return raster
s <- stack(PropTrajEnd, PropTrajFT, PropTrajSt, ClassMov, IntS, BounS, TrajClas)
names(s) <- c("PropEnd", "PropFT", "PropSt", "ClassL", "IntS", "BounS", "TrajClas")
return(s)
}







 
 
 
 
 
 
 
 
 