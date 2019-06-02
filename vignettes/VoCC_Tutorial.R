## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
# If not installed already get the VoCC package first
# devtools::install_github("JorGarMol/VoCC")
library(VoCC)
library(rasterVis)
library(gridExtra)
library(doParallel)
library(foreach)
library(scales)
library(data.table)
library(mapplots)

## ------------------------------------------------------------------------
#?marshift

## ------------------------------------------------------------------------
# info on the sea surface temperature data set
#?HSST
# monthly to annual averages
r <- sumSeries(HSST, p = "1960-01/2009-12", yr0 = "1955-01-01", l = nlayers(HSST), fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years")
# temporal trend
vt <- tempTrend(r, th = 10)
# spatial gradient
vg <- spatGrad(r, th = 0.0001, projected = FALSE)
# climate velocity
gv <- gVoCC(vt, vg)

# Now the distance-based velocities
# Take 1960-1970 as base period against 2000-2009
r2 <- stack(mean(r[[1:10]], na.rm = TRUE), mean(r[[41:50]], na.rm = TRUE))
# prepare the data frame with the necessary variables
clim <- na.omit(data.frame(getValues(r2), cid = 1:ncell(r)))
clim[,c("x","y")] <- xyFromCell(r, clim$cid)
# climate velocity; 1965-2004 (40 yr), 500 km search radius
v <- dVoCC(clim, n = 1, tdiff = 40, method = "Single", climTol = 0.1, geoTol = 500, distfun = "GreatCircle", trans = NA, lonlat = TRUE)

## ------------------------------------------------------------------------
# Change sign as needed and create the distance-based velocity raster
ind <-  which(r2[[1]][v$focal] > r2[[2]][v$target])
v$velBis <- v$vel
v$velBis[ind] <- v$vel[ind] * -1
# put output in raster format
dv <- raster(gv)
dv[v$focal] <- v$velBis

data(marshift)
marshift$GV <- with(marshift, raster::extract(abs(gv[[1]]), cbind(long, lat), buffer = (Shift* (timespan/10)*1000), fun = mean, na.rm = TRUE, small = TRUE))
marshift$DV <- with(marshift, raster::extract(abs(dv), cbind(long, lat), buffer = (Shift* (timespan/10)*1000), fun = mean, na.rm = TRUE, small = TRUE))
# retrieve the closest marine cell for those centroids falling on land
marshift$GV <- with(marshift, ifelse(is.na(GV), gv[[1]][which.min(replace(distanceFromPoints(gv[[1]], cbind(long,lat)), is.na(gv[[1]]), NA))], GV))
marshift$DV <- with(marshift, ifelse(is.na(DV), dv[which.min(replace(distanceFromPoints(dv, cbind(long, lat)), is.na(dv), NA))],DV))
# fit the regression models
Mgv <- lm(Shift^(1/4) ~ I((GV*10)^(1/4)), data = marshift, weights=years_data)
Mdv <- lm(Shift^(1/4) ~ I((DV*10)^(1/4)), data = marshift, weights=years_data)
summary(Mgv) 
summary(Mdv)

