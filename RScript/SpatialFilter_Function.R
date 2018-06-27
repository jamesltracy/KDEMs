SpatialFilter <- function(xy, dist, mapUnits = F) {
## NOTE: Probably should always work with data in geographic projection with WGS84 datum for this function
    #mapUnits=T
    #xy=monrst.spdf
    #dist=1
    ## Code by Pascal Title, Univ. Michigan, Ecology and Evol. Biology
    ## From: http://stackoverflow.com/questions/22051141/spatial-filtering-by-proximity-in-r
    #xy can be either a SpatialPoints or SPDF object, or a matrix
    # calculate desired buffer distance around presence points
    #dist is in km if mapUnits=F, in mapUnits otherwise
    if (!mapUnits) {
        d <- spDists(xy,longlat=T)
    }
    if (mapUnits) {
        d <- spDists(xy,longlat=F)
    }
    diag(d) <- NA
    close <- (d <= dist)
    diag(close) <- NA
    closePts <- which(close,arr.ind=T)
    discard <- matrix(nrow=2,ncol=2)
    if (nrow(closePts) > 0) {
            while (nrow(closePts) > 0) {
                if ((!paste(closePts[1,1],closePts[1,2],sep='_') %in% paste(discard[,1],discard[,2],sep='_')) & (!paste(closePts[1,2],closePts[1,1],sep='_') %in% paste(discard[,1],discard[,2],sep='_'))) {
                discard <- rbind(discard, closePts[1,])
                closePts <- closePts[-union(which(closePts[,1] == closePts[1,1]), which(closePts[,2] == closePts[1,1])),]
                }
            }
        discard <- discard[complete.cases(discard),]
        return(xy[-discard[,1],])
    }
    if (nrow(closePts) == 0) {
        return(xy)
    }
}

###########################################################################
#require(rgeos)
#require(sp)
#pts <- readWKT("MULTIPOINT ((3.5 2), (1 1), (2 2), (4.5 3), (4.5 4.5), (5 5), (1 5))")
#
#pts2 <- filterByProximity(pts,dist=2, mapUnits=T)
#
#plot(pts, axes=TRUE)
#
#
#apply(as.data.frame(pts),1,function(x) plot(gBuffer(SpatialPoints(coords=matrix(c(x[1],x[2]),nrow=1)),width=2),add=T))
#plot(pts2,add=T,col='blue',pch=20,cex=2)