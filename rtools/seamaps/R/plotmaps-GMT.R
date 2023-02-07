#----------------------------------------------------------------
# FUNCTION using GMT (modified Anders Nielsen's function 'plotmap'
# allowing overplotting on the map (add=TRUE), showing the coastline and grid)
#----------------------------------------------------------------
  plot.map<-function (x1, x2, y1, y2, resolution = 3, grid = TRUE, add = FALSE,
    save = FALSE, landcolor = rgb(0.85,0.82,0.66), seacolor = "white",
    zoom = FALSE, coastline=FALSE, coastcolor=rgb(0.37,0.35,0.3))
  {
    #comment: GMT functions need always incresing or decreasing vector of coordinates
    if ((x1<0) & (x2<0)) {x1<-x1+360; x2<-x2+360}
    if ((x1>=180) & (x2>360)) {x1<-x1-360; x2<-x2-360}

    GMT <- function(x1, x2, y1, y2, resolution = 3) {
        read.ps.line <- function(txt) {
            txt.split <- strsplit(txt, split = " ")[[1]]
            ret <- c(NA, NA)
            if (length(txt.split) == 3) {
                if (txt.split[3] %in% c("M", "moveto", "D")) {
                  ret <- as.numeric(txt.split[1:2])
                }
            }
            return(ret)
        }
        if (resolution < 1 || resolution > 5)
            stop("resolution from 1 (full) to 5(crude)")
        res <- c("f", "h", "i", "l", "c")[resolution]
        filen <- tempfile("gmtmap")
        on.exit(unlink(c(filen, ".gmtcommands4")))
        cmd <- paste("gmt pscoast -Rd",x1,"/",x2,"/",y1,"/",y2, " -Jx0.5id -P -G -D", 
	             res, " -X0 -Y0 > ", filen,sep = "")
        system(cmd)
        
        txt <- readLines(filen); 
        mat <- matrix(unlist(lapply(txt, read.ps.line)), ncol = 2,
            byrow = TRUE)
        for (i in 2:nrow(mat)) {
            if (!is.na(mat[i, 1]) & !is.na(mat[i - 1, 1]))
                mat[i, ] <- mat[i, ] + mat[i - 1, ]
        }
        maxx <- max(mat[, 1], na.rm = TRUE)
        maxy <- max(mat[, 2], na.rm = TRUE)
        mat[, 1] <- mat[, 1]/600 + x1
        mat[, 2] <- mat[, 2]/600 + y1
        return(mat)
    }
    junk <- GMT(x1, x2, y1, y2, resolution = resolution)
    if (!add) {
        plot(c(x1, x2), c(y1, y2), type = "n", ylab = "", xlab = "",
            xaxs = "i", yaxs = "i",xaxt="n",yaxt="n")
        rect(x1, y1, x2, y2, col = seacolor)
    }
    if (grid){
            axis(1, tck = 1, labels=FALSE, lwd=0.5, lty=2, col=rgb(0.2,0.2,0.2)) # no ugly labels
            axis(2, tck = 1, labels=FALSE, lwd=0.5, lty=2, col=rgb(0.2,0.2,0.2))
    }
    polygon(junk, border = landcolor, col = landcolor)
    if (zoom) {
        ret <- locator(2)
        if (length(ret$x) < 2) {
            zoom <- FALSE
        }
        else {
            x1 <- min(ret$x)
            x2 <- max(ret$x)
            y1 <- min(ret$y)
            y2 <- max(ret$y)
        }
        plot.map(x1, x2, y1, y2, resolution, grid, add, save,landcolor, seacolor, zoom)
    }
    if (save) {
        dimnames(junk)[[2]] <- c("longitude", "latitude")
        return(junk)
    }
   # add coastline
   if (coastline){# & !(x1<180 & x2>360)){
     filen <- tempfile("gmtmap") ;  on.exit(unlink(c(filen, ".gmtcommands4")))
     res <- c("f", "h", "i", "l", "c")[resolution]
     cmd <- paste("gmt pscoast -Rd",x1,"/",x2,"/",y1,"/",y2," -Jx2id -W -M  -D",res," > ",filen,sep="")
     system(cmd)
     dat <- readLines(filen);
     ff<-function(str)if(regexpr("#",str)>0){
       c(NA, NA)
     }else{
       xy<-as.numeric(strsplit(str, split="\t")[[1]])
       xy[1]<-ifelse(xy[1]<0,xy[1]+360,xy[1])
       return(xy)      
     }
     lines(t(sapply(dat, ff)),col=coastcolor,lwd=0.5)
    }
  }

  #' Adding coast contour to existing map
  #' @param reso is the coastline resolution, from 1 (full resolution) to  5 (crude). Recommended value for basin wide maps is 4.
  #' @examples 
  #' #TBD 
  #' @export
  addcoast<-function(reso=3, coastcolor=rgb(0.37,0.35,0.3),...) {
    # ----------------------------------------------------------------------------
    # "THE BEER-WARE LICENSE":
    # <anders.nielsen@hawaii.edu> wrote this function. As long as you retain this notice you
    # can do whatever you want with this stuff. If we meet some day, and you think
    # this stuff is worth it, you can buy me a beer in return. Anders Nielsen
    # ----------------------------------------------------------------------------

    usr<-par("usr"); x1<-usr[1]; x2<-usr[2]; y1<-usr[3]; y2<-usr[4];
    filen <- tempfile("gmtmap") ;  on.exit(unlink(c(filen, ".gmtcommands4")))
    res <- c("f", "h", "i", "l", "c")[reso]
    cmd <- paste("gmt pscoast -R",x1,"/",x2,"/",y1,"/",y2," -Jx2id -W -M -D",res," > ",filen,sep="")
    system(cmd)
    dat <- readLines(filen);
    ff<-function(str)if(regexpr("#",str)>0){
      c(NA, NA)
    }else{
      xy<-as.numeric(strsplit(str, split="\t")[[1]])
     if (x1<180 & x2>180) xy[1]<-ifelse(xy[1]<0,xy[1]+360,xy[1])
      return(xy)      
    }
    lines(t(sapply(dat, ff)),col=coastcolor, lwd=1)
  }


  Wpos<-function(x){
    return(ifelse(x<0,x+360,x))
  }

  Wneg<-function(x){
    return(ifelse(x>180,x-360,x))
  }


 #' Plotting map with nice labels and formatting
 #' @param X,Y are the vectors of longitude and latitude, can be of length two, with map borders.
 #' @param reso is the map resolution, from 1 (full resolution) to  5 (crude). Recommended value for basin wide maps is 4.
 #' @param grid is boolean. If TRUE, then the grid will be plotted at intervals determined by function 'axis'.
 #' @param col is the color of the land. 
 #' @param col.sea is the color of the sea. 
 #' @param mar are the user margins passed to function 'par'
 #' @param cex defines the size of the axis labels passed as cex.axis in function 'par'
 #' @examples 
 #' #Pacific-wide map with default settings 
 #' nice.map(120:290,-50:50)
 #' nice.map(c(120,290),c(-50,50)) # same output
 #' @export
 nice.map<-function(X,Y,reso=4,grid=TRUE,col="lightgrey",col.sea="white",mar=c(4,6,6,2),cex=1.0){

   par(mar=mar,cex.axis=cex,tck=-0.01,mgp=c(0.25,0.8,0))
   nx<-length(X)
   if (X[nx]>360 | (X[1]<0 & X[nx]>0))
     plot.map(X[1],X[nx],Y[1],Y[length(Y)],res=reso,grid,landcolor=dark.grey,coastline=FALSE,seacolor=col.sea)
   else plot.map(X[1],X[nx],Y[1],Y[length(Y)],res=reso,grid,landcolor=col,coastline=TRUE,seacolor=col.sea)
  
   box(lwd=2.0)
   #----------------------------------------------
   # MAP LABELS
   #----------------------------------------------
   usr=par("usr"); x1<-usr[1]; x2<-usr[2];
   X<-x1:x2
   xlabels<-pretty(X); 	ylabels<-pretty(Y)
   xlabnames<-ifelse((Wneg(xlabels)>=0 & xlabels!=180),
 		    parse(text=paste(Wneg(xlabels),"^o * E",sep="")),
 		    ifelse(xlabels==180,parse(text=paste(xlabels,"^o",sep="")),
 			   parse(text=paste((-Wneg(xlabels)),"^o * W",sep=""))))
   ylabnames<-ifelse(ylabels<0,parse(text=paste(abs(ylabels),"^o * S",sep="")),
 		    parse(text=paste(ylabels,"^o * N",sep=""))); 
   ylabnames<-ifelse(ylabels==0,parse(text=paste(abs(ylabels),"^o",sep="")),ylabnames)
   tlen<-0.0
   axis(1,at=xlabels,lab=xlabnames,tck=tlen); 
   axis(2,at=ylabels,lab=ylabnames,las=1,tck=tlen); 
   return(list(xlab=xlabels,ylab=ylabels))
 }

 #' Plotting map with nice labels and formatting on the existing plot
 #' @param X,Y are the vectors of longitude and latitude, can be of length two, with map borders.
 #' @param reso is the map resolution, from 1 (full resolution) to  5 (crude). Recommended value for basin wide maps is 4.
 #' @param grid is boolean. If TRUE, then the grid will be plotted at intervals determined by function 'axis'.
 #' @param col is the color of the land. 
 #' @param col.sea is the color of the sea. 
 #' @param mar are the user margins passed to function 'par'
 #' @param cex defines the size of the axis labels passed as cex.axis in function 'par'
 #' @examples 
 #' #Pacific-wide map with default settings 
 #' nice.map(120:290,-50:50)
 #' nice.map(c(120,290),c(-50,50)) # same output
 #' @export
 nice.map.over<-function(X,Y,reso=5,grid=TRUE,col=light.grey,cex=1.0){

   par(cex.axis=cex,tck=-0.01,mgp=c(1,0.5,0))
   nx<-length(X);
   if (X[nx]>360 | (X[1]<0 & X[nx]>0))
     plot.map(X[1],X[nx],Y[1],Y[length(Y)],res=reso,grid=FALSE,landcolor=col,coastline=FALSE,add=TRUE)
   else plot.map(X[1],X[nx],Y[1],Y[length(Y)],res=reso,grid=FALSE,coastline=TRUE,landcolor=col,add=TRUE)
 
   box(lwd=2.0)
   #----------------------------------------------
   # MAP LABELS
   #----------------------------------------------
   xlabels<-pretty(X); 
 
   ylabels<-pretty(Y)
   ny<-length(ylabels)
   if (ny>5) ylabels<-ylabels[seq(1,ny,2)]
   xlabnames<-ifelse((Wneg(xlabels)>=0 & xlabels!=180),
 		    parse(text=paste(Wneg(xlabels),"^o * E",sep="")),
 		    ifelse(xlabels==180,parse(text=paste(xlabels,"^o",sep="")),
 			   parse(text=paste((-Wneg(xlabels)),"^o * W",sep=""))))
   ylabnames<-ifelse(ylabels<0,parse(text=paste(abs(ylabels),"^o * S",sep="")),
 		    parse(text=paste(ylabels,"^o * N",sep=""))); 
   ylabnames<-ifelse(ylabels==0,parse(text=paste(abs(ylabels),"^o",sep="")),ylabnames)
   axis(2,at=ylabels,lab=ylabnames,las=1);
   par(mgp=c(1,0.85,0))
   axis(1,at=xlabels,lab=xlabnames)
   if (grid){
      axis(1, tck = 1, at=xlabels, labels=FALSE, lwd=0.5, lty=2, col=rgb(0.5,0.5,0.5))
      axis(2, tck = 1, at=ylabels, labels=FALSE, lwd=0.5, lty=2, col=rgb(0.5,0.5,0.5))
   }
 }


