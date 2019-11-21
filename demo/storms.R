require(XML); 
library(sp);
library(spacetime);
library(trajectories)
extract.track=function(year = 2012, p = TRUE) {
# based on # ARTHUR CHARPENTIER # http://freakonometrics.hypotheses.org/17113
	loc <- paste("http://weather.unisys.com/hurricane/atlantic/",year,"/index.php",sep="")
	tabs <- readHTMLTable(htmlParse(loc)) 
	storms <- unlist(strsplit(as.character(tabs[[1]]$Name),split=" "))
	index <- storms %in% c("Tropical","Storm", paste("Hurricane-",1:6,sep=""),
		"Depression","Subtropical","Extratropical","Low",
		paste("Storm-",1:6,sep=""), "Xxx")
	nstorms  <- storms[!index]
	tracks = list()
	for(i in length(nstorms):1) {
		loc=paste("http://weather.unisys.com/hurricane/atlantic/",
			year, "/", nstorms[i], "/track.dat", sep="")
		track=read.fwf(loc, skip=3, widths = c(4,6,8,12,4,6,20))
		names(track)=c("ADV", "LAT", "LON", "TIME", "WIND", "PR", "STAT")
		track$LAT=as.numeric(as.character(track$LAT))
		track$LON=as.numeric(as.character(track$LON))
		track$WIND=as.numeric(as.character(track$WIND))
		track$PR=as.numeric(as.character(track$PR))
		track$year=year
		pts = SpatialPoints(cbind(track$LON, track$LAT), 
			CRS("+proj=longlat +datum=WGS84 +datum=WGS84"))
		time = as.POSIXct(paste(year, track$TIME, sep="/"), format="%Y/ %m/%d/%HZ  ",tz="UTC")
		tr = Track(STIDF(pts, time, track))
		tracks[nstorms[i]] = tr
		if (p==TRUE)
			cat(year,i,nstorms[i],nrow(track),"\n")
	}
	return(Tracks(tracks))
}
#m=extract.track(2012)
#m=extract.track(2011:2012)
TOTTRACK=list()
#for(y in 2012:1851)
for(y in 2012:2009)
	#TOTTRACK=rbind(TOTTRACK, extract.track(y))
	# http://robertgrantstats.wordpress.com/2014/10/01/transparent-hurricane-paths-in-r/
	if (!inherits(try(x <- extract.track(y)), "try-error"))
		TOTTRACK[as.character(y)] = x

storms = TracksCollection(TOTTRACK)

library(maps)
map("world",xlim=c(-80,-40),ylim=c(10,50),col="light yellow",fill=TRUE)
plot(storms, col = sp::bpy.colors(4, alpha = .25), lwd = 8, add = TRUE)

plot(storms)
x = approxTracksCollection(storms, by = "30 min", FUN = spline)
plot(x, col = 'red', add = TRUE)

TOTTRACK = as(storms, "data.frame")
library(ks)
U=TOTTRACK[,c("LON","LAT")]
U=U[!is.na(U$LON),]
H=diag(c(.2,.2))
# note that this might be not meaningful, as coords are longlat:
fat=kde(U,H,xmin=c(min(U[,1]),min(U[,2])),xmax=c(max(U[,1]),max(U[,2])))
z=fat$estimate
long = fat$eval.points[[1]]
lat = fat$eval.points[[2]]
image(long, lat, z)
plot(storms, add=TRUE)
map("world",add=TRUE)

plotLLSlice = function(obj, xlim, ylim, newCRS, maps = TRUE, gridlines = TRUE,..., mar = par('mar') - 2, f = .05) {
	stopifnot(require(rgdal))
	stopifnot(require(rgeos))
	l.out = 100 # length.out
	# draw outer box:
	p = rbind(cbind(xlim[1], seq(ylim[1],ylim[2],length.out = l.out)), 
          cbind(seq(xlim[1],xlim[2],length.out = l.out),ylim[2]), 
	      cbind(xlim[2],seq(ylim[2],ylim[1],length.out = l.out)), 
          cbind(seq(xlim[2],xlim[1],length.out = l.out),ylim[1]))
	LL = CRS("+init=epsg:4326")
	bb = SpatialPolygons(list(Polygons(list(Polygon(list(p))),"bb")), proj4string = LL)

	expBboxLL = function(x, f = f) { x@bbox[,1] = bbox(x)[,1] - f * apply(bbox(x), 1, diff); x }
	plot(expBboxLL(spTransform(bb, newCRS), f), mar = mar)

	if (gridlines)
		plot(spTransform(gridlines(bb), newCRS), add = TRUE)

	if (maps) {
		stopifnot(require(maps))
		m = map(xlim = xlim, ylim = ylim, plot = FALSE, fill = TRUE)
		IDs <- sapply(strsplit(m$names, ":"), function(x) x[1])
		library(maptools)
		m.sp <- map2SpatialPolygons(m, IDs=IDs, proj4string = LL)
		m = gIntersection(m.sp, bb) # cut map slice in WGS84
		plot(spTransform(m, newCRS), add = TRUE, col = grey(0.8))
	}

	if (is(obj, "TracksCollection"))
		obj = as(obj, "SpatialLinesDataFrame")

	obj = spTransform(obj, CRS(proj4string(bb))) # to WGS84, might be obsolete
	obj = gIntersection(obj, bb) # cut Slice
	plot(spTransform(obj, newCRS), add = TRUE, ...)
	text(labels(gridlines(bb), newCRS))
}

laea = CRS("+proj=laea +lat_0=30 +lon_0=-80")
data(storms)
plotLLSlice(storms, xlim = c(-100.01,-19.99), ylim = c(10, 55), laea, col = 'orange', lwd = 2)
plotLLSlice(storms, xlim = c(-90,-80), ylim = c(20, 30), laea, col = 'orange', lwd = 2)
