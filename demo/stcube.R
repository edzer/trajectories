library(spacetime)
if (!requireNamespace("enviroCaR", quietly = TRUE)) {
# Import enviroCar track.
  importEnviroCar = function(trackID, url = "https://envirocar.org/api/stable/tracks/") {
  	require(RCurl)
  	require(rgdal)
  	require(rjson)
  	url = getURL(paste(url, trackID, sep = ""), .opts = list(ssl.verifypeer = FALSE))
  	# Read data into spatial object.
  	spdf = readOGR(dsn = url, layer = "OGRGeoJSON", verbose = FALSE)
  	# Convert time from factor to POSIXct.
  	time = as.POSIXct(spdf$time)
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

  A1 = importEnviroCar("528cf1a3e4b0a727145df093")
  A2 = importEnviroCar("528cf19ae4b0a727145deb40")
  A3 = importEnviroCar("528cf194e4b0a727145de63d")
  
  A = Tracks(tracks = list(A1, A2))
  B = Tracks(tracks = list(A3))
  Tr = TracksCollection(tracksCollection = list(A, B))
  
} else {
  A1 = enviroCaR::importEnviroCar("https://envirocar.org/api/stable/", "528cf1a3e4b0a727145df093")
  A2 = enviroCaR::importEnviroCar("https://envirocar.org/api/stable/", "528cf19ae4b0a727145deb40")
  A3 = enviroCaR::importEnviroCar("https://envirocar.org/api/stable/", "528cf194e4b0a727145de63d")
  
  A = Tracks(tracks = list(A1@tracksCollection$Tracks1@tracks$Track1, A2@tracksCollection$Tracks1@tracks$Track1))
  B = Tracks(tracks = list(A3@tracksCollection$Tracks1@tracks$Track1))
  Tr = TracksCollection(tracksCollection = list(A, B))
}
# Plot tracks in a space-time cube. As fetching the basemap takes its time,
# "showMap" defaults to FALSE. See ?stcube for detailed information.

stcube(A, showMap = FALSE)
stcube(B, showMap = TRUE)
stcube(Tr, showMap = TRUE)
