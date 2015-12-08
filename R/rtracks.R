rTrack = function(n = 100, origin = c(0,0), start = as.POSIXct("1970-01-01"), ar = .8, 
		step = 60, sd0 = 1, ...) { 
	if (length(ar) == 1 && ar == 0)
		xy = cbind(cumsum(rnorm(n, sd = sd0)) + origin[1], cumsum(rnorm(n, sd = sd0)) + origin[2])
	else 
		xy = cbind(origin[1] + cumsum(as.vector(arima.sim(list(ar = ar), n, sd = sd0, ...))),
			       origin[2] + cumsum(as.vector(arima.sim(list(ar = ar), n, sd = sd0, ...))))
	T = start + 0:(n-1) * step
	sti = STI(SpatialPoints(xy), T)
	Track(sti)
}

rTracks = function(m = 20, start = as.POSIXct("1970-01-01"), delta = 7200, sd1 = 0, 
		origin = c(0,0), ...)
	Tracks(lapply(0:(m-1) * delta, function(x) rTrack(start = start + x, 
		origin = origin + rnorm(2, sd = sd1), ...)))

rTracksCollection = function(p = 10, sd2 = 0, ...)
	TracksCollection(lapply(1:p, function(x) rTracks(origin = rnorm(2, sd = sd2) ,...)))
