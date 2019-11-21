require(RCurl)
require(rgdal)
require(rjson)
require(sp)
require(spacetime)
importEnviroCar = function(trackID, url = "https://envirocar.org/api/stable/tracks/") {
	url = getURL(paste(url, trackID, sep = ""), 
		.opts = list(ssl.verifypeer = FALSE)) # .opts needed for Windows
	# Read data into spatial object.
	spdf = readOGR(dsn = url, layer = "OGRGeoJSON", verbose = FALSE)
	# Convert time from factor to POSIXct.
	#time = as.POSIXct(spdf$time, format = "%Y-%m-%dT%H:%M:%SZ")
	time = as.POSIXct(paste0(as.character(spdf$time),"00"), format = "%Y/%m/%d %H:%M:%S%z")
	# Convert phenomena from JSON to data frame.
	phenomena = lapply(as.character(spdf$phenomenons), fromJSON)
	values = lapply(phenomena, function(x) as.data.frame(lapply(x, function(y) y$value)))
	# Get a list of all phenomena for which values exist.
	names = vector()
	for(i in values)
		names = union(names, names(i))
	# Make sure that each data frame has the same number of columns.
	values = lapply(values, function(x) {
		xNames = names(x)
		# Get the symmetric difference.
		diff = setdiff(union(names, xNames), intersect(names, xNames))
		if(length(diff) > 0)
			x[diff] = NA
		x
	})
	# Bind values together.
	data = do.call(rbind, values)
	sp = SpatialPoints(coords = coordinates(spdf), 
		proj4string = CRS("+proj=longlat +datum=WGS84"))
	stidf = STIDF(sp = sp, time = time, data = data)
	Track(track = stidf)
}
A3 = importEnviroCar("528cf1a3e4b0a727145df093")
plot(A3)
