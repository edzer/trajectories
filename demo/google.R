# read a file obtained by "Export to kml" from https://maps.google.com/locationhistory/
filename = system.file("history-06-09-2015.kml", package="trajectories")

library(XML)
kml <- xmlToList(filename)

tr = kml$Document$Placemark$Track
cc = which(names(tr) == "coord")
coord = t(sapply(kml$Document$Placemark$Track[cc], function(x) scan(text = x, quiet = TRUE)))[,1:2]
when = which(names(tr) == "when")
# convert the "-07:00" into " -0700" with sub:
#time = strptime(sub("-08:00$", " -0800", unlist(kml$Document$Placemark$Track[when])),
time = strptime(sub("([+\\-])(\\d\\d):(\\d\\d)$", " \\1\\2\\3", 
	unlist(kml$Document$Placemark$Track[when])), "%Y-%m-%dT%H:%M:%OS %z")

library(sp)
library(spacetime)
library(trajectories)
track = Track(STI(SpatialPoints(coord, CRS("+proj=longlat +ellps=WGS84")), time))
summary(track)
head(as(track, "data.frame"))
plot(track, axes = TRUE)
