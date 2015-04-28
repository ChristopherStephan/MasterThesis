library(raster)
library(rgdal)
library(bfast)
library(classInt)
library(RColorBrewer)
library(rts)
library(stringr)
library(M3)
require(animation)
require(date)
require(igraph)
library(maptools)
library(dichromat)
library(gridExtra)
library(rasterVis)

##############################
### Create satellite image time series
##############################

# location of files
pathFilt = "/media/christopher/backup_2/boku/c84a8757358a86a8f675f0591ec1f7e5/filt"
pathRaw = "/media/christopher/backup_2/boku/c84a8757358a86a8f675f0591ec1f7e5/raw"

# list of raster files
lst = list.files(path=pathFilt, pattern="EVI", full.names=TRUE)

# extracting julian date of file name (?<=\/filt\/)(.*)(?=.MOD)
dates = lapply(lst, function(x) str_extract(x, perl("(?<=\\/filt\\/).*?(?=.MOD|$)")))

# julian date to date
dates = as.Date(decipher.M3.date(as.numeric(unlist(dates))))

# creating a RasterStack object
rasterStack <- stack(lst)

writeRaster(rasterStack, filename="manicore_filt.grd", format="raster", overwrite=TRUE, options=c("COMPRESS=NONE"))

##############################
### Process raster, scale EVI, trim time series
##############################

rasterStack = "/home/christopher/Documents/Master Thesis/Spatial Analysis/Site3/manicore_filt.gri"

rasterStack = stack(rasterStack)

prodesPolygon = readOGR("/home/christopher/Documents/Master Thesis/Spatial Analysis/Site3", layer="prodesPolygon2011_reprj")

# crop scene
spatialExtent = readOGR("/home/christopher/Documents/Master Thesis/PRODES/Manicore/newExtent", layer="ext")
rasterStack = crop(rasterStack, spatialExtent)
plot(rasterStack,1)

# trim satellite image time series at the end of PRODES year 2011
# trim from 1st of August 2011 - Julian day 213
# http://landweb.nascom.nasa.gov/browse/calendar.html
# names(rasterStack)
rasterStack = subset(rasterStack, c(1:264))

# scale the EVI values
rasterStack = calc(rasterStack, function(x) x/10000)

plot(rasterStack, 1)
plot(prodesPolygon, add=T, col="red")

writeRaster(rasterStack, filename="site3_filt.grd", format="raster", overwrite=TRUE, options=c("COMPRESS=NONE"))

##############################
### APPLY BFASTMONITOR ON RASTER STACK
##############################

rasterStack = "/home/christopher/Documents/Master Thesis/Spatial Analysis/Site3/site3_filt.gri"
rasterStack = "/home/christopher/Documents/Master Thesis/Spatial Analysis/Site3/Raster files/site3_raw.gri"

rasterStack = stack(rasterStack)

# helper function
xbfastmonitor <- function(x, dates) {
  evits <- bfastts(x, dates, type = c("16-day"))
  # call of bfastmonitor
  bfm <- bfastmonitor(data = evits, start=c(2010, 12), formula = response ~ harmon, history=c(2006,1))  
  
  return(cbind(bfm$breakpoint, bfm$magnitude, bfm$history[1], bfm$history[2], bfm$monitor[1], bfm$monitor[2]))
}

ptm = proc.time()

timeofbreak = calc(rasterStack, fun=function(x){
  res = t(apply(x, 1, xbfastmonitor, dates))
  return(res)
})

runtimeInMinutes = (proc.time() - ptm) / 60

print(runtimeInMinutes)

# set meaningful layernames
names(timeofbreak) = c("breakpoint", "magnitude", "historyStart", "historyEnd", "monitorStart", "monitorEnd")

# save the resulting raster object
writeRaster(timeofbreak, filename="site3Bfm_filt.grd", format="raster", overwrite=TRUE, options=c("COMPRESS=NONE"))


##############################
### Forest mask 
##############################

site3Bfm_filt = "/home/christopher/Documents/Master Thesis/Spatial Analysis/Site3/Raster files/site3Bfm_filt.gri"
site3BfmStack_filt = stack(site3Bfm_filt)
plot(site3BfmStack_filt, 2)
plot(prodesPolygon, add=TRUE)

site3Bfm_raw = "/home/christopher/Documents/Master Thesis/Spatial Analysis/Site3/Raster files/site3Bfm_raw.gri"
site3BfmStack_raw = stack(site3Bfm_raw)
plot(site3BfmStack_raw, 2)
plot(prodesPolygon, add=TRUE)

# create mask
maskPolygon = readOGR("/home/christopher/Documents/Master Thesis/Spatial Analysis/Site3", layer="mask_reprj")
plot(maskPolygon, add=TRUE)

mm = mask(raster(site3BfmStack_raw, 2), maskPolygon)
mm = reclassify(mm, c(-1, 1, 1))
plot(mm, col="lightsalmon", legend=F)
levelplot(mm, margin=F, col.regions="lightsalmon", colorkey=F, main="Rasterized forest mask and clear-cut deforestation polygons from PRODES year 2011") + layer(sp.lines(prodesPolygonShp, lwd=0.8))

site3BfmStack_raw = mask(site3BfmStack_raw, maskPolygon, inverse=TRUE)
site3BfmStack_filt = mask(site3BfmStack_filt, maskPolygon, inverse=TRUE)

##############################
### Plotting the change magnitude
##############################

site3BfmRaster_filt_brks = raster(site3BfmStack_filt, 1)
plot(site3BfmRaster_filt_brks)
site3BfmRaster_filt_brks_rcl = reclassify(site3BfmRaster_filt_brks, c(0, 2013, 1))
plot(site3BfmRaster_filt_brks_rcl)

site3BfmRaster_filt_chmg = raster(site3BfmStack_filt, 2)
site3BfmRaster_filt = mask(site3BfmRaster_filt_chmg, site3BfmRaster_filt_brks, updatevalue=0)
plot(site3BfmRaster_filt)

site3BfmRaster_raw_brks = raster(site3BfmStack_raw, 1)
plot(site3BfmRaster_raw_brks)
site3BfmRaster_raw_brks_rcl = reclassify(site3BfmRaster_raw_brks, c(0, 2013, 1))
plot(site3BfmRaster_raw_brks_rcl)

site3BfmRaster_raw_chmg = raster(site3BfmStack_raw, 2)
site3BfmRaster_raw = mask(site3BfmRaster_raw_chmg, site3BfmRaster_raw_brks, updatevalue=0)
plot(site3BfmRaster_raw)

magVals_filt = na.omit(getValues(site3BfmRaster_filt))
magVals_raw = na.omit(getValues(site3BfmRaster_raw))

# negative magnitudes
redPurples = c("#f1eef6", "#d7b5d8", "#df65b0", "#dd1c77", "#980043")

# positive magnitudes
yellowGreens = c("#addd8e", "#78c679", "#41ab5d", "#238443", "#005a32")

# estimate class breaks for negative change magnitudes
fj5NegMag_filt = classIntervals(magVals_filt[magVals_filt < 0], n=5, style="fisher", precision = 3)
fj5NegMag_raw = classIntervals(magVals_raw[magVals_raw < 0], n=5, style="fisher", precision = 3)

# estimate class breaks for positve change magnitudes
fj5PosMag_filt = classIntervals(magVals_filt[magVals_filt > 0], n=5, style="fisher", precision = 3)
fj5PosMag_raw = classIntervals(magVals_raw[magVals_raw > 0], n=5, style="fisher", precision = 3)

# subsetting the change magnitude map
xdel = abs(61.80164 - 61.16771)/2
ydel = abs(8.107933 - 7.628022)/2

left = extent(c(-61.80164, -61.16771 - xdel, -8.107933, -7.628022))
right = extent(c(-61.80164 + xdel, -61.16771, -8.107933, -7.628022))

site3BfmRaster_filt_right = crop(site3BfmRaster_filt, left)
site3BfmRaster_raw_right = crop(site3BfmRaster_raw, left)
plot(site3BfmRaster_filt_right)
plot(site3BfmRaster_raw_right)

prodesPolygonShp = readShapeLines("/home/christopher/Documents/Master Thesis/Spatial Analysis/Site3/prodesPolygon2011_reprj.shp")
maskPolygonShp = readShapeLines("/home/christopher/Documents/Master Thesis/Spatial Analysis/Site3/mask_reprj.shp")

#right_1
brksFilt = round(c(fj5NegMag_filt$brks, fj5PosMag_filt$brks[-1]), 3)
brksRaw = round(c(fj5NegMag_raw$brks, fj5PosMag_raw$brks[-1]), 3)

redsAndGreens = c(rev(redPurples), yellowGreens)

myColorkeyFilt = list(at=brksFilt, labels=list(at=brksFilt))
filt = levelplot(site3BfmRaster_filt_right, col.regions=redsAndGreens, at=brksFilt,
                 colorkey=myColorkeyFilt, pretty=TRUE, margin=F, main="Smoothed MODIS EVI")

filt = filt + layer(sp.lines(prodesPolygonShp, lwd=0.8))

myColorkeyRaw = list(at=brksRaw, labels=list(at=brksRaw))
raw = levelplot(site3BfmRaster_raw_right, col.regions=redsAndGreens, at=brksRaw,
                colorkey=myColorkeyRaw, pretty=TRUE, margin=F, main="Raw MODIS EVI")
raw = raw + layer(sp.lines(prodesPolygonShp, lwd=0.8)) # + layer(sp.lines(maskPolygonShp, fill='darkgray'))

print(filt, split=c(1, 1, 2, 1) , more=TRUE)
print(raw, split=c(2, 1, 2, 1))


#right_2
brksFilt = round(c(fj5NegMag_filt$brks), 3)[-(5:6)]
brksRaw = round(c(fj5NegMag_raw$brks), 3)[-(5:6)]

redsAndGreens = c(rev(redPurples[-1]))

myColorkeyFilt = list(at=brksFilt, labels=list(at=brksFilt))
filt = levelplot(site3BfmRaster_filt_right, col.regions=rev(redPurples)[-(5:6)], at=brksFilt,
                 colorkey=myColorkeyFilt, pretty=TRUE, margin=F, main="Smoothed MODIS EVI", col.na="red")
filt = filt + layer(sp.lines(prodesPolygonShp, lwd=0.8))

myColorkeyRaw = list(at=brksRaw, labels=list(at=brksRaw))
raw = levelplot(site3BfmRaster_raw_right, col.regions=rev(redPurples)[-(5:6)], at=brksRaw,
                colorkey=myColorkeyRaw, pretty=TRUE, margin=F, main="Raw MODIS EVI")
raw = raw + layer(sp.lines(prodesPolygonShp, lwd=0.8)) # + layer(sp.lines(maskPolygonShp, fill='darkgray'))

print(filt, split=c(1, 1, 2, 1) , more=TRUE)
print(raw, split=c(2, 1, 2, 1))

#right_3
brksFilt = round(c(fj5NegMag_filt$brks), 3)[-(4:6)]
brksRaw = round(c(fj5NegMag_raw$brks), 3)[-(4:6)]

redsAndGreens = c(rev(redPurples[-1]))

myColorkeyFilt = list(at=brksFilt, labels=list(at=brksFilt))
filt = levelplot(site3BfmRaster_filt_right, col.regions=rev(redPurples)[-(3:6)], at=brksFilt,
                 colorkey=myColorkeyFilt, pretty=TRUE, margin=F, main="Smoothed MODIS EVI", col.na="red")
filt = filt + layer(sp.lines(prodesPolygonShp, lwd=0.8))

myColorkeyRaw = list(at=brksRaw, labels=list(at=brksRaw))
raw = levelplot(site3BfmRaster_raw_right, col.regions=rev(redPurples)[-(3:6)], at=brksRaw,
                colorkey=myColorkeyRaw, pretty=TRUE, margin=F, main="Raw MODIS EVI")
raw = raw + layer(sp.lines(prodesPolygonShp, lwd=0.8)) # + layer(sp.lines(maskPolygonShp, fill='darkgray'))

print(filt, split=c(1, 1, 2, 1) , more=TRUE)
print(raw, split=c(2, 1, 2, 1))

#right_4
brksFilt = round(c(fj5NegMag_filt$brks), 3)[-(3:6)]
brksRaw = round(c(fj5NegMag_raw$brks), 3)[-(3:6)]

redsAndGreens = c(rev(redPurples[-1]))

myTheme <- BTCTheme()
myTheme$panel.background$col = 'grey'

myColorkeyFilt = list(at=brksFilt, labels=list(at=brksFilt))
filt = levelplot(site3BfmRaster_filt_right,  par.settings=myTheme, col.regions=c("#dd1c77"), at=brksFilt,
                 colorkey=myColorkeyFilt, pretty=TRUE, margin=F, main="Smoothed MODIS EVI")
filt = filt + layer(sp.lines(prodesPolygonShp, lwd=0.8))

myColorkeyRaw = list(at=brksRaw, labels=list(at=brksRaw))
raw = levelplot(site3BfmRaster_raw_right, col.regions=rev(redPurples)[-(3:6)], at=brksRaw,
                colorkey=myColorkeyRaw, pretty=TRUE, margin=F, main="Raw MODIS EVI")
raw = raw + layer(sp.lines(prodesPolygonShp, lwd=0.8)) # + layer(sp.lines(maskPolygonShp, fill='darkgray'))

print(filt, split=c(1, 1, 2, 1) , more=TRUE)
print(raw, split=c(2, 1, 2, 1))



##############################
### Error matrix
##############################

# edit vector for building threshold range types
brksFilt = round(c(fj5NegMag_filt$brks), 3)#[-(5:6)]
brksRaw = round(c(fj5NegMag_raw$brks), 3)[-(5:6)]

# edit first parameter and the brks to generate the corresponding error matrices for smoothed and raw data
site3BfmRasterRecl = reclassify(site3BfmRaster_filt, rcl=c(first(brksFilt)-0.1, last(brksFilt), 1, last(brksFilt), 1, 0))

plot(site3BfmRasterRecl, colNA="red")
values(site3BfmRasterRecl)

site3BfmRasterRecl_filt_right = reclassify(site3BfmRaster_filt_right, rcl=c(first(brksFilt)-0.1, last(brksFilt), 1, last(brksFilt), 1, 0))
site3BfmRasterRecl_raw_right = reclassify(site3BfmRaster_raw_right, rcl=c(first(brksRaw)-0.1, last(brksRaw), 1, last(brksRaw), 1, 0))

# Binary map showing all negative magnitudes as green and all positve magnitueds as red
plot(site3BfmRasterRecl, main="Red: BFAST Monitor identified change", legend=F)
plot(prodesPolygon, add=T, lwd=0.1)

rc1 = levelplot(site3BfmRasterRecl_filt_right, margin=F, main="Smoothed MODIS EVI", colorkey=F,
           par.settings = myTheme, col.regions=c("lightgray", "red"))
rc1 = rc1 + layer(sp.lines(prodesPolygonShp, lwd=0.8))


rc2 = levelplot(site3BfmRasterRecl_raw_right, margin=F, main="Raw MODIS EVI", colorkey=F,
                par.settings = myTheme, col.regions=c("lightgray", "red"))
rc2 = rc2 + layer(sp.lines(prodesPolygonShp, lwd=0.8))

print(rc1, split=c(1, 1, 2, 1) , more=TRUE)
print(rc2, split=c(2, 1, 2, 1))


#create a dummy raster and reset the raster values
dummyRaster = reclassify(site3BfmRaster_raw, rcl=c(-1, 1, 0))
# test if the dummyRaster was fully resetted, TRUE is expected
all(getValues(dummyRaster) == 0, na.rm=TRUE)

# rasterize the prodes polygon where deforestation was monitored in the PRODES year 2011
# update the value of the cells that overlap with the prodes polygon with the value 1
# prodesPolygonRaster = rasterize(prodesPolygon, dummyRaster, update=TRUE, field=1)
prodesPolygonRaster = rasterize(as(prodesPolygon, "SpatialLines"), dummyRaster, update=TRUE, field=1)
prodesPolygonRaster = rasterize(prodesPolygon, prodesPolygonRaster, update=TRUE, field=1)
plot(prodesPolygonRaster, main="Gray: rasterized PRODES polygon", breaks=c(0.5, 1), col="lightgray", legend=F)
plot(prodesPolygon, add=T, lwd=1)

# rasterize the prodes polygon where deforestation was monitored before 2013
prodesPolygon = readShapeLines("/home/christopher/Documents/Master Thesis/Spatial Analysis/Site3/prodesPolygonsSmaller2013.shp")
prodesPolygonRaster = rasterize(as(prodesPolygon, "SpatialLines"), dummyRaster, update=TRUE, field=1)
prodesPolygonRaster = rasterize(prodesPolygon, prodesPolygonRaster, update=TRUE, field=1)
plot(prodesPolygonRaster, main="Gray: rasterized PRODES polygon", breaks=c(0.5, 1), col="lightgray", legend=F)
plot(prodesPolygon, add=T, lwd=1)

valuesProdes = getValues(prodesPolygonRaster)
samplesProdesDeforested = length(valuesProdes[valuesProdes == 1])
samplesProdesForested = length(valuesProdes[valuesProdes == 0])

prodes = getValues(prodesPolygonRaster) # PRODES
bfm = getValues(site3BfmRasterRecl) # change map 

confusionMatrix = aggregate(x=prodes, by=list(prodes, bfm), FUN = length)

# group 1: PRODES, group 2: BFAST Monitor
confusionMatrix

confusionMatrix$x[1]+confusionMatrix$x[2] 
confusionMatrix$x[3]+confusionMatrix$x[4]

confusionMatrix$x[1]+confusionMatrix$x[3] 
confusionMatrix$x[2]+confusionMatrix$x[4]

# functions to compute confustion matrix based statistics
# overall accuracy
# the relative observed agreement among PRODES and BFAST Monitor
getOverallAccuracy = function(confusionMatrix){
  return((confusionMatrix$x[1] + confusionMatrix$x[4]) /
           (confusionMatrix$x[1] + confusionMatrix$x[2] +
              confusionMatrix$x[3] + confusionMatrix$x[4])) 
}

getUsersAccuracies = function(confusionMatrix, class){
  # forest / no change
  if(class == 0)
    userAcc = (confusionMatrix$x[1]) / (confusionMatrix$x[1] + confusionMatrix$x[2])
  # deforested / change
  else if(class == 1)
    userAcc = (confusionMatrix$x[4]) / (confusionMatrix$x[4] + confusionMatrix$x[3])
  
  return(userAcc)
}

getProducersAccuracies = function(confusionMatrix, class){
  # forest / no change
  if(class == 0)
    prodAcc = (confusionMatrix$x[1]) / (confusionMatrix$x[1] + confusionMatrix$x[3])
  # deforested / change
  else if(class == 1)
    prodAcc = (confusionMatrix$x[4]) / (confusionMatrix$x[2] + confusionMatrix$x[4])
  
  return(prodAcc)
}

getOverallAccuracy(confusionMatrix)
getUsersAccuracies(confusionMatrix, 1)
getUsersAccuracies(confusionMatrix, 0)
getProducersAccuracies(confusionMatrix, 1)
getProducersAccuracies(confusionMatrix, 0)

##############################
### Grid distances for false positves of BAST Monitor
##############################
prodesPolygonExt = readShapeLines("/home/christopher/Documents/Master Thesis/Spatial Analysis/Site3/prodesPolygon2011_reprj_ext.shp")

gd =  gridDistance(prodesPolygonRaster, origin=1)
breaks = c(0, 1, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000)
gdCols = rev(c("#ffffcc", "#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#b10026", "white"))
plot(gd, colNA="black", col=gdCols, main="Distances from deforested areas in m", breaks=breaks)

gdm = mask(gd, maskPolygon, updatevalue=F, inverse=T)
plot(prodesPolygonExt, add=TRUE)
plot(gdm)

myColorkey = list(at=breaks, labels=list(at=breaks))
l = levelplot(gdm, col.regions=gdCols, at=breaks, colorkey=myColorkey, margin=F, main="Distances from deforested areas in m")
l = l + layer(sp.lines(prodesPolygonShp))
l



gdStack =  stack(gdm, site3BfmRasterReclFilt)
gdMask = mask(gdm, site3BfmRasterReclFilt,  updatevalue=F, maskvalue=0)
#plot(gdMask, col=gdCols, breaks=breaks)
l = levelplot(gdMask, col.regions=gdCols, at=breaks, colorkey=myColorkey, margin=F,
              main="Distances from false positives to deforested areas in m")
l = l + layer(sp.lines(prodesPolygonExt))
l


# Cumulative relative frequency graph
distances = getValues(gdMask)
# remove null values
distances = distances[distances != 0]
plot(ecdf(distances), verticals = T, col.points = "blue", col.hor = "red",
     main="Cumulative relative frequency graph for the distances",
     ylab="Cumulative distance proportion", xlab="Distances to PRODES polygon [m]")

# histogram for distance classes
hist(distances)

##############################
### Plotting the time of change point  
##############################

minisampleBfmBrkpt = raster(minisampleBfmStack, 1)

plot(minisampleBfmBrkpt,  col = heat.colors(8))
plot(prodesPolygon, add=TRUE)

plot(minisampleBfmBrkpt,  col = heat.colors(5), zlim=c(2004.4347, 2007))
plot(prodesPolygon, add=TRUE)

##############################
### Filtering the BFAST Monitor Change map
##############################

# 3x3 filter that filters out pixels that have no neighbor euqals 1 
filter = function(x){
  seedHasNeighbor = any(x[c(1:4, 6:9)] > 0, na.rm=TRUE)
  if((!is.na(x[5])) && (x[5] > 0) && (seedHasNeighbor == FALSE)) {
    return(0)
  }
  return(x[5])
}

site3BfmRasterReclFilt = focal(site3BfmRasterRecl, w=matrix(1/9, nrow=3, ncol=3), fun=filter, pad=TRUE, padValue=0) 

# re-establishing comparable values for binary classes
vals = getValues(site3BfmRasterReclFilt)
vals = replace(vals, vals == 0, 0)
vals = replace(vals, vals >= 0.1, 1)
site3BfmRasterReclFilt = setValues(site3BfmRasterReclFilt, vals)

# compare the filtered maps
plot(site3BfmRasterRecl)
plot(site3BfmRasterReclFilt)

# confusion matrix
prodes = getValues(prodesPolygonRaster) # PRODES
bfm = getValues(site3BfmRasterReclFilt) # change map 

confusionMatrix = aggregate(x=prodes, by=list(prodes, bfm), FUN = length)
confusionMatrix

getOverallAccuracy(confusionMatrix)
getUsersAccuracies(confusionMatrix, 1)
getUsersAccuracies(confusionMatrix, 0)
getProducersAccuracies(confusionMatrix, 1)
getProducersAccuracies(confusionMatrix, 0)


##############################
### Plotting the forest mask
##############################
plot(dummyRaster, col="white", legend=F)
plot(maskPolygon, add=T, col="lightsalmon", lwd=0.5)

# number of pixeld covered by the forest mask
165256 - 113724 + 1029

##############################
### Create Raster Time Series
##############################

# get layer names from raster stack
layerNames = names(rasterStack)

# extracting julian date of layer name 
dates = lapply(layerNames, function(x) str_extract(x, perl("(?<=X).*?(?=.MOD|$)")))

# transform julian date to date
dates = as.Date(decipher.M3.date(as.numeric(unlist(dates))))

rasterStackTs = rts(rasterStack, dates)

plot(site3BfmRasterRecl)
plot(prodesPolygon, add=T)

##############################
### Selecting BFAST Monitor's true and error classifications
##############################

click(site3BfmRasterRecl, cell=T, xy=T)

lines(cl$x, cl$y, type = "p", col="blue")
### raw
# false negative:
# 25591 A
# 10593 B -61.61302 -7.711728
# 23121 C -61.54159 -7.809942
# 43929 D -61.37195 -7.972888
# E -61.3965 -7.754138 16086
# F -61.22239 -7.876906 31784
# G -61.57284 -7.941638 39863

# false positives
# -61.35409 -7.660388 4177 -0.04547702 H
# -61.52373 -7.69387  8361 -0.0380592 I
# -61.47239 -7.825567 25140 -0.05334639 J 
# 61.59963 -7.986281 45531 -0.04111024 K
# -61.63088 -7.676013 6041 -0.107695 L
# -61.76704 -8.093424 58803 -0.02548427 M

# true positives
# -61.45231 -7.972888 43893 -0.08494597 N
# -61.41436 -7.910388 35958 -0.06306232 O 
# -61.52373 -7.937174 39317 -0.05896039 P
# -61.68891 -7.872442 31007 -0.07542304 Q
# -61.49248 -7.698335 8943 -0.0466025 R
# -61.23579 -7.946103 40582 -0.1078322 S

### filt
# false negative:
# -61.72686 -7.841192 26731 0
# -61.61302 -7.711728 10593 0
# -61.60409 -7.827799 25365 0 597

# false positives
# -61.47239 -7.825567 25140 1
# -61.4032 -8.086728 58399  1
# -61.5215 -7.69387 8362    1

# true positives
# -61.45231 -7.975121 44178 1
# -61.72909 -7.825567 25025 1
# -61.46347 -7.832263 25996 1


# Time Series of pixels where change was reported by PRODES 2006
timeSeries = ts(as.numeric(rasterStackTs[25996,]), start=c(2000, 4), frequency=23)
plot(timeSeries)
# type = c("OLS-CUSUM", "OLS-MOSUM", "RE", "ME", "fluctuation")

# original call of bfastmonitor with the same parameter set
bfm =  bfastmonitor(timeSeries, start=c(2010, 12), formula = response ~ harmon, history=c(2006, 1))
plot(bfm)
plot(bfm, main="I: Break detected at: 2011(6)")


##############################
### Plot locations of true and error classifications
##############################

prodesPolygon2012 = readShapeLines("/home/christopher/Documents/Master Thesis/Spatial Analysis/Site3/prodesPolygon2012_reprj")
myTheme <- BTCTheme()
myTheme$panel.background$col = 'white'
#main="Red: BFAST Monitor's reclassified change magnitudes of smoothed MODIS EVI with threshold type 3"
rc1 = levelplot(site3BfmRasterRecl, margin=F,colorkey=F,
                par.settings = myTheme, col.regions=c("lightgray", "red"))
rc1 = rc1 + layer(sp.lines(prodesPolygonShp, lwd=0.8)) + layer(sp.lines(prodesPolygon2012, lwd=0.8, col="blue"))
rc1 = rc1 + layer(sp.points(pointsFalseNegs, col = 'blue', cex=2))
rc1 = rc1 + layer(sp.text(c(-61.73132 + 0.01, -7.830031 + 0.01), txt="A", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.61302 + 0.01, -7.711728 + 0.01), txt="B", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.54159 + 0.01, -7.809942 + 0.01), txt="C", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.37195 + 0.01, -7.972888 + 0.01), txt="D", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.3965 + 0.01, -7.754138 + 0.01), txt="E", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.22239 + 0.01, -7.876906 + 0.01), txt="F", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.57284 + 0.01, -7.941638 + 0.01), txt="G", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.35409 + 0.01, -7.660388 + 0.01), txt="H", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.52373 + 0.01, -7.69387 + 0.01), txt="I", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.47239 + 0.01, -7.825567 + 0.01), txt="J", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.59963 + 0.01, -7.986281 + 0.01), txt="K", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.63088 + 0.01, -7.67601 + 0.01), txt="L", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.76704 + 0.01, -8.093424  + 0.01), txt="M", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.45231 + 0.01, -7.972888+ 0.01), txt="N", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.41436 + 0.01, -7.910388 + 0.01), txt="O", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.52373 + 0.01, -7.937174 + 0.01), txt="P", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.68891 + 0.01, -7.872442 + 0.01), txt="Q", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.49248 + 0.01, -7.698335 + 0.01), txt="R", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.23579 + 0.01, -7.946103 + 0.01), txt="S", cex=1.5))
rc1

x = c(-61.73132, -61.61302, -61.54159, -61.37195, -61.3965, -61.22239, -61.57284, -61.35409, -61.52373, -61.47239,- 61.59963, -61.63088, -61.76704, -61.45231, -61.41436, -61.52373, -61.68891, -61.49248, -61.23579) 
y = c(-7.830031, -7.711728, -7.809942, -7.972888, -7.754138, -7.876906, -7.941638, -7.660388, -7.69387, -7.825567, -7.986281, -7.67601, -8.093424, -7.972888, -7.910388, -7.937174, -7.872442, -7.698335, -7.946103)
xy = cbind(x,y)
v1 = c("A", "B", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s")
v2 = c("A", "B", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s")
attributes <- as.data.frame(cbind(v1,v2))
pointsFalseNegs = SpatialPointsDataFrame(xy, attributes) #, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))


rc1 = levelplot(site3BfmRasterRecl, margin=F,colorkey=F,
                par.settings = myTheme, col.regions=c("lightgray", "red"))
rc1 = rc1 + layer(sp.lines(prodesPolygonShp, lwd=0.8)) + layer(sp.lines(prodesPolygon2012, lwd=0.8, col="blue2"))
rc1 = rc1 + layer(sp.points(pointsFalseNegs, col = 'blue', cex=2)) 
rc1 = rc1 + layer(sp.text(c(-61.72686  + 0.01, -7.841192 + 0.01), txt="A", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.61302 + 0.01, -7.711728 + 0.01), txt="B", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.60409 + 0.01, -7.827799 + 0.01), txt="C", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.47239 + 0.01, -7.825567 + 0.01), txt="D", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.4032 + 0.01, -8.086728 + 0.01), txt="E", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.5215 + 0.01, -7.69387 + 0.01), txt="F", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.45231 + 0.01, -7.975121 + 0.01), txt="G", cex=1.5))
rc1 = rc1 + layer(sp.text(c(-61.72909 + 0.01, -7.825567 + 0.01), txt="H", cex=1.5))
rc1 = rc1 + layer(sp.text(c( -61.46347 + 0.01, -7.832263 + 0.01), txt="I", cex=1.5))
rc1 

x = c(-61.72686, -61.61302, -61.60409, -61.47239, -61.4032, -61.5215, -61.45231, -61.72909, -61.46347) 
y = c(-7.841192, -7.711728, -7.827799, -7.825567, -8.086728, -7.69387, -7.975121, -7.825567, -7.832263)
xy = cbind(x,y)
v1 = c("A", "B", "c", "d", "e", "f", "g", "h", "i")
v2 = c("A", "B", "c", "d", "e", "f", "g", "h", "i")
attributes <- as.data.frame(cbind(v1,v2))
pointsFalseNegs = SpatialPointsDataFrame(xy, attributes) #, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))


##############################
### Stratified sampling
##############################

userAccuracyNoChangeVec = c()
userAccuracyChangeVec = c()
prodAccuracyNoChangeVec = c()
prodAccuracyChangeVec = c()


for(i in 1:100) {
  sample = sampleStratified(site3BfmRasterReclFilt, size=50) #prodesPolygonRaster
  sampleCell = sample[,1]
  
  prodesSampleCellValVector = c()
  for(i in c(1:length(sampleCell))){
    sampleCellVal = sampleCell[i]
    prodesSampleCellVal = prodesPolygonRaster[sampleCellVal]
    prodesSampleCellValVector = append(prodesSampleCellValVector, prodesSampleCellVal)
  }
  
  prodes = prodesSampleCellValVector
  bfm = sample[,2]
  
  confusionMatrix = aggregate(x=prodes, by=list(prodes, bfm), FUN = length)
  # print(confusionMatrix)
  
  userAccuracyNoChangeVec = append(userAccuracyNoChangeVec, getUsersAccuracies(confusionMatrix, 0))
  userAccuracyChangeVec = append(userAccuracyChangeVec, getUsersAccuracies(confusionMatrix, 1))
  
  prodAccuracyNoChangeVec = append(prodAccuracyNoChangeVec, getProducersAccuracies(confusionMatrix, 0))
  prodAccuracyChangeVec = append(prodAccuracyChangeVec, getProducersAccuracies(confusionMatrix, 1))
  
}

##############################
### Create a movie of satellite image time series
##############################

names(rasterStack) = dates

saveGIF({
  for(i in c(211:length(names(rasterStack)))){
    plot(rasterStack, 264, col=(heat.colors(20)), zlim=c(0, 1))
    plot(prodesPolygon, add=TRUE)    
  }
}, interval=0.2, movie.name="movie_filt_evi_site_1.gif")

plot(rasterStack, 264, col=(heat.colors(10)), zlim=c(0, 1))
plot(prodesPolygon, add=TRUE)

# plotting the scene with visual change
plot(rasterStack, 264, col=(heat.colors(20)), zlim=c(0, 1))
plot(prodesPolygon, add=TRUE)

plot(rasterStack, c(259:264), col=(heat.colors(10)), zlim=c(0, 1))
plot(prodesPolygon, add=TRUE)
