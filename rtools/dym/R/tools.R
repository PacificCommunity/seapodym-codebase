# A set of functions to manipulate DYM files: do arithmetic operations on multiple files, subset files, aggregate through space (to time series) or through time (to 2D distribution).

is.decimal <- function(date){
	if (is.numeric(date)) return(TRUE)
	else return(FALSE)
}

is.date <- function(date){
	if (class(date)=="Date") return(TRUE)
	else return(FALSE)
}

#' Converts DYM's decimal date to the Date format
#'
#' @param dec.date is the date or an array of dates in decimal format as used in DYM files, i.e., year+(month-1)/12+(day_of_month-1)/366. 
#' @return corresponding date(s) in the Date format.
#' @examples 
#' dec.2date(1975.0); #returns "1974-12-31"
#' dec.2date(1975.01); #returns "1975-01-03"
#' @export
dec.2date <- function(dec.date){
	year  <- trunc(dec.date)
	days <- trunc((dec.date - year)*365);
	date<-as.Date(paste(year,1,1,sep="-"))+days-1
	return(date)
}

#' Converts a vector of monthly dates to a vector of DYM's decimal dates. Useful when working with 360-days calendar to generate a tvect to be written into the DYM file.
#'
#' @param dates a vector with monthly dates as used for the 360-days calendar. 
#' @return corresponding vector with its elements in decimal format.
#' @examples 
#' dates <- gen.monthly.dates(as.Date("1970-01-15"),as.Date("1970-05-15"))
#' tvect <- mdate.2dec(dates);
#' @export
mdate.2dec <- function(dates) as.integer(format(dates,"%Y"))+(as.integer(format(dates,"%m"))-.5)/12


#' Converts date(s) in Date format to date(s) in DYM's decimal dates. Useful for generating a tvect to be written into the DYM file.
#'
#' @param dates a date or a vector of dates in Date format. 
#' @return corresponding vector with its elements in decimal format.
#' @examples 
#' tvect <- date.2dec(as.Date("1975-01-15")); #returns 1975.038
#' @export
date.2dec <- function(dates){
	tvect <- as.integer(format(dates,"%Y"))+
			(as.integer(format(dates,"%m"))-1)/12+
			(as.integer(format(dates,"%d"))-1)/366
	return(tvect)
}

#for file names the same is done by basename(filename), but maybe useful for more general cases
get.filenames <- function(full.names) sapply(lapply(sapply(dyms.in,function(x)strsplit(x,"/")),rev),"[[",1)

cell.surface.area <- function(lat,dx,dy){
	#returns the area of a cell on a sphere in sq.km
	#between 'lat' and 'lat+dy'. Sizes dx and dy are
	#to be entered in minutes, angle 'lat' in degrees.
	R = 6378.1;
	Phi1 = lat*pi/180.0;
	Phi2 = (lat+dy/60.0)*pi/180.0;
	dx_radian = (dx/60.0)*pi/180;
	S = R*R*dx_radian*(sin(Phi2)-sin(Phi1));

	return(S)
}
   
check.dimensions <- function(filename,ts,xs,ys){

	txy <- read.txy.dym(filename)

	flag <- FALSE
	check.range <- function(d,d.dym,type){
		if (d[1] < d.dym[1]){
			d[1] <- d.dym[1]
			flag <<- TRUE
			message("Changed first ",type," to ",d[1],", which is the first DYM ",type,".")
		}
		if (d[1] > d.dym[length(d.dym)] | d[2] < d.dym[1]){
			message("ERROR: ",type," range is outside of DYM's ",type," range. Exit now.")
			return(NULL)
		}
		if (d[2] > d.dym[length(d.dym)]){
			d[2] <- d.dym[length(d.dym)]
			flag <<- TRUE
			message("Changed last ",type," to ",d[2],", which is the last DYM ",type,".")
		} 
		return(d)
	}
	
	if (!is.null(ts)){   
		t.dym <- dec.2date(txy$t)
		ts <- check.range(ts,t.dym,"date")
		if (is.null(ts)) stop("Date problem!")
	}
	x.dym  <- txy$x 
	if (!is.null(xs)){
		xs <- check.range(xs,x.dym,"longitude")
		if (is.null(xs)) stop("Longitude problem!")				
	} else 
		xs <- range(x.dym)

	y.dym  <- txy$y 
	if (!is.null(ys)){
		ys <- check.range(ys,y.dym,"latitude")
		if (is.null(ys)) stop("Latitude problem!")				
	} else 
		ys <- range(y.dym)	

	return(list(t=ts,x=xs,y=ys,changed=flag))		
}

#' Subsetting one or multiple DYM files
#'
#' Function subsets multiple DYM files with identical dimensions. It will subset each file and write the new set of DYMs with the names given in dyms.out. If the latter is not provided, the files with the same names will be saved to the folder 'subset.name-year1-year2_x1-x2_y1-y2'.
#' @param dyms.in the vector of strings containing the full names (with full path or relative to home directory) of the input DYM files.
#' @param dyms.out the name of the output DYM file. If not specified, then the same names will be used for the subsetted files, and saved into a separate folder.
#' @param ts the set of two dates in the Date format c(t1,t2) for subsetting. If left undefined, the original time period will be used.
#' @param xs the set of longitudes c(x1,x2) for subsetting. If left undefined, the original longitudinal range will be used.
#' @param ys the set of latitudes in format c(y1,y2) for subsetting. If left undefined, the original latitudinal range will be used.
#' @examples 
#' #Example 1. Subset the global forcing directory 'forcing_dir_global' for the Western and Central Pacific Ocean and 30-years time period. SEAPODYM environment is needed to run this example.
#' t.range <- c(as.Date("1991-1-1"),as.Date("2020-12-31"))
#' x.range <- c(100,210)
#' y.range <- c(-50,50)
#' dyms.in <- list.files(paste0(Sys.getenv("SEA_FORCING_HOME"),forcing_dir_global,"/"),full.names=TRUE,pattern=".dym")
#' subset.dyms(dyms.in,ts=t.range,xs=x.range,ys=y.range,subset.name="wcpo");
#' @export subset.dyms  
subset.dyms <- function(dyms.in,dyms.out=NULL,ts=NULL,xs=NULL,ys=NULL,subset.name="subset"){
	  
	if (is.null(ts) & is.null(xs) & is.null(ys)){
		stop("Subsetting dimensions not specified, nothing to do - stopped.")
	}
	#check dimensions and conformity: if dimension subsetting defined 
	#and incorrect, then exit without subsetting and issue warning.
	out <- check.dimensions(dyms.in[1],ts,xs,ys)
	ts <- out$t; xs <- out$x; ys <- out$y
	
	t0 <- NULL; tfin <- NULL
	if (!is.null(ts)){
		dates <- sapply(ts,function(x)format(x,c('%Y','%m','%d')))
		t0   <- as.numeric(dates[,1])
		tfin <- as.numeric(dates[,2])
	}

	#Subsetting		
	if (is.null(dyms.out)){
		dir.out <- paste0(subset.name,"-",				  
						  t0[1],"-",tfin[1],
						  "_",xs[1],"-",xs[2],
						  "_",ys[1],"-",ys[2])
		if (!dir.exists(dir.out))
			dir.create(dir.out)

		dyms.out <- paste0(dir.out,"/",basename(dyms.in))
	}
	for (f in 1:length(dyms.in)){
		dym<-read.var.dym(dyms.in[f],t0,tfin,region=c(xs,ys),t.format="dym")
		var<-dym$var; tc<-dym$t; x<-dym$x; y<-dym$y; mask<-dym$landmask
		#if (is.date(tc)) tc <- date.2dec(tc)
		write.var.dym(dyms.out[f],tc,x,y,mask,var)
	}
}

#' Computing the product of two or more DYM files
#'
#' Sum up multiple DYM files with identical dimensions, writing the result to new DYM file and/or returning it.
#' @param dyms.in the vector of strings containing the full names (with full path or relative to home directory) of the input DYM files.
#' @param dym.out the name of the output DYM file. If not specified, then the result will be returned (see out parameter).
#' @param out the logical value (FALSE by default) telling if the result should be returned by the function. It will be automatically turned to TRUE if dym.out is not provided.
#' @param average a logical indicating whether the result should be an average.
#' @returns Optionally, returns a list with dimensional variables, t, x and y, as well as z, a 3d variable.
#' @examples 
#' # Example 1. Get the skipjack denstiy that is accessible to the fishing gear that can fish only in waters with euphotic depth < 100 m. SEAPODYM environment is needed to run this example.
#' dir.in  <- paste0(Sys.getenv("SEA_FORCING_HOME"),"/run-interim_2x30d_po/")
#' dym.1 <- list.files(dir.in,pattern="zeu",full.names=TRUE)
#' dym.2 <- system.file("extdata", "skj.dym", package = "dym")
#' txy <- read.txy.dym(dym.2)
#  # First make dym.1 with identical to dym.2 dimensions
#' dym.1s <- paste0("my-",basename(dym.1))
#' subset.dyms(dym.1,dym.1s,range(dec.2date(txy$t)),range(txy$x),range(txy$y))
#' qfun <- function(x,slope,depth) 1/(1+exp((slope)*(x-depth)))
#' dym.1sq <- paste0("my-q-",basename(dym.1s))
#' eapply.dyms(dym.1s,dym.1sq,FUN="qfun",slope=0.5,depth=100)
#' dym.out <- paste0("accessible-",basename(dym.2))
#' mult.dyms(c(dym.1sq,dym.2),dym.out); # writing the result to dym.out file
#'
#' @export mult.dyms
mult.dyms <- function(dyms.in,dym.out=NULL,out=FALSE){ 

	for (f in 1:length(dyms.in)){
		dym<-read.var.dym(dyms.in[f],apply.mask=TRUE,verbose=TRUE)
		var<-dym$var; tc<-dym$t; x<-dym$x; y<-dym$y; mask<-dym$landmask
			
		if (f==1){
			var.out<-var
			tc1 <- tc; x1 <- x; y1 <- y; mask1 <- mask
		}
		if (f>1){
			var.out<-var.out * var
			if (any(tc != tc1) | any(x != x1) | any(y != y1))
				stop("Different dimensions, exit.")
			if (any(mask != mask1)){
				message("Proceed with caution, masks are different! See file mask-diff.txt ")
				write(mask-mask1,"mask-diff.txt",ncolumn=length(x))
			}
		}
	}

	if (!is.null(dym.out)){
		write.var.dym(dym.out,date.2dec(tc),x,y,mask,var.out)
	} else 
		out <- TRUE
		if (out) return(list(t=tc,x=x,y=y,z=var.out))
}

#' Computing the sum of multiple DYM files
#'
#' Sum up multiple DYM files with identical dimensions, writing the result to new DYM file and/or returning it.
#' @param dyms.in the vector of strings containing the full names (with full path or relative to home directory) of the input DYM files.
#' @param dym.out the name of the output DYM file. If not specified, then the result will be returned (see out parameter).
#' @param out the logical value (FALSE by default) telling if the result should be returned by the function. It will be automatically turned to TRUE if dym.out is not provided.
#' @param average a logical indicating whether the result should be an average.
#' @returns Optionally, returns a list with dimensional variables, t, x and y, as well as z, a 3d variable.
#' @examples 
#' # Example 1. Get the epipelagic forage at night
#' dir.in <- paste0(Sys.getenv("SEA_FORCING_HOME"),"/run-era5np_1x30d-po/forage/")
#' dyms.in <- paste0(dir.in,"Fbiom_",c("epi","mmeso","hmbathy"),".dym")
#' dym.out <- gsub("epi","epi_night",dyms.in[1])
#' sum.dyms(dyms.in,dym.out); # writing the result to DYM file
#'
#' # Example 2. Get the density of recruits to be compared with MFCL class
#' setwd(paste0(Sys.getenv("SKJ_HOME"),"/REF/output/"))
#' dyms.in <- paste0("skj_age",4:6,".dym")
#' dym.out <- "skj_recruits.dym"
#' sum.dyms(dyms.in,dym.out); # writing the result to the DYM file
#' out <- sum.dyms(dyms.in,dym.out,out=TRUE); # writing DYM file and output the result
#' out <- ave.dyms(dyms.in); # output the average density in selected age classes
#' @export sum.dyms
sum.dyms <- function(dyms.in,dym.out=NULL,out=FALSE,average=FALSE){ 

	for (f in 1:length(dyms.in)){
		dym<-read.var.dym(dyms.in[f],t.format="dym",apply.mask=TRUE,verbose=TRUE)
		var<-dym$var; tc<-dym$t; x<-dym$x; y<-dym$y; mask<-dym$landmask
			
		if (f==1){
			var.out<-var
			tc1 <- tc; x1 <- x; y1 <- y; mask1 <- mask
		}
		if (f>1){
			var.out<-var.out + var
			if (any(tc != tc1) | any(x != x1) | any(y != y1))
				stop("Different dimensions, exit.")
			if (any(mask != mask1)){
				message("Proceed with caution, masks are different! See file mask-diff.txt ")
				write(mask-mask1,"mask-diff.txt",ncolumn=length(x))
			}
		}
	}

	if (average==TRUE)
		var.out <- var.out/length(dyms.in)		

		if (!is.null(dym.out)){
			write.var.dym(dym.out,tc,x,y,mask,var.out)
		} else 
			out <- TRUE

	if (out) return(list(t=tc,x=x,y=y,z=var.out))
}


#' @describeIn sum.dyms Average multiple DYM files with identical dimensions, writing the results to new DYM file and/or returning it.
#' @export
ave.dyms <- function(dyms.in,dym.out=NULL,out=FALSE){ 

	res <- sum.dyms(dyms.in,dym.out,out=out,average=TRUE)
	if (out) return(res)		
}
	
#' Computing the difference between two DYM files
#'
#' Function computes the difference between two DYM files 'dym1-dym2' and either writes it to another DYM file or returns the result.
#' @param dym.1 a string containing the full name (with full path or relative to home directory) of the first DYM file.
#' @param dym.2 a string containing the full name (with full path or relative to home directory) of the second DYM file.
#' @param dym.out a string giving the name of output DYM file.
#' @returns Optionally, returns a list with dimensional variables, t, x and y, as well as z, a 3d variable.
#' @examples 
#' # Example 1. Compare seasonal dynamics between two years 
#' dym.in <- system.file("extdata", "skj.dym", package = "dym")
#' txy <- read.txy.dym(dym.in)
#' yy <- unique(floor(txy$t))
#' seasons <- c(1,4,7,10) # first months to split year to seasons
#' fyy <- paste0(yy,".dym")
#' for (f in fyy){
#' 	y <- unlist(strsplit(basename(f),".",fixed=T))[1]
#' 	ts <- c(as.Date(paste0(y,"-01-01")),as.Date(paste0(y,"-12-31")))
#' 	subset.dyms(dym.in,f,ts)
#' 	make.climatology(f,imonths=seasons)
#' } # two files with seasonal means named [f]_clm_quarterly.dym created
#' ff.seasonal <- gsub(".dym","_clm_quarterly.dym",fyy)
#' dat <- diff.dyms(ff.seasonal[1],ff.seasonal[2],out=TRUE)
#' # Now plot the results
#' x <- dat$x; y <- dat$y
#' par(mfrow=c(2,2),mar=c(2,2,2,2),mgp=c(3,0.3,0),tck=0.02)
#' for (n in 1:4)
#' 	image(x,y,dat$z[n,,],zlim=range(dat$z,na.rm=T),
#' 			col=hcl.colors(hcl.pals("diverging")[3],n=15),
#'				main=paste("Quarter",n))
#' @export diff.dyms
diff.dyms <- function(dym1,dym2,dym.out=NULL,out=FALSE){ 

	dym<-read.var.dym(dym1,t.format="dym",apply.mask=TRUE,verbose=TRUE)
	var1<-dym$var; tc<-dym$t; x<-dym$x; y<-dym$y; mask<-dym$landmask
	dym<-read.var.dym(dym2,t.format="dym",apply.mask=TRUE,verbose=TRUE)
	var2<-dym$var;

	if (any(tc != dym$t) | any(x != dym$x) | any(y != dym$y))
		stop("Different dimensions, exit.")
	if (any(mask != dym$landmask)){
		message("Proceed with caution, masks are different! See file mask-diff.txt ")
		write(mask-dym$landmask,"mask-diff.txt",ncolumn=length(x))
	}

	var.out <- var1 - var2

	if (!is.null(dym.out)){
		write.var.dym(dym.out,tc,x,y,mask,var.out)
	} else out <- TRUE

	if (out) return(list(t=tc,x=x,y=y,z=var.out))
}

#' Element-wise operation on DYM file(s) variables
#'
#' Function applies a function to every element of DYM file(s) and writes the result to dyms.out with the same dimensions.
#' @param dyms.in a string or a vector of strings containing the full names (with full path or relative to home directory) of the input DYM files.
#' @param dyms.out a string or a vector of strings giving the full names of output DYM file(s). Must be provided, otherwise will do nothing.
#' @param FUN a function to be applied to the DYM variable on an element level, should be defined elsewhere or directly as an argument.
#' @examples 
#' # Example 1. Multiply all elements of micronekton input data by 0.5. Run this example once the SEAPODYM environment is fully set up and below variables are known
#' dir.in  <- paste0(Sys.getenv("SEA_FORCING_HOME"),my_forcing_dir,"/")
#; mult <- 0.5
#' dir.out <- paste0(dir.in,"forage_x",mult,"/")
#' if (!dir.exists) dir.create(dir.out)
#' dyms.in  <- list.files(paste0(dir.in,"forage/"),pattern="Fbiom",full.names=TRUE)
#' dyms.out <- paste0(dir.out,basename(dyms.in))
#' eapply.dyms(dyms.in,dyms.out,function(x)mult*x);
#'
#' #Example 2. compute the gaussian function of sea surface temperature
#' dym.in  <- list.files(dir.in,pattern="sst",full.names=TRUE)
#' dym.out <- gsub("sst","gauss-sst",dym.in)
#' gauss <- function(x,a,b) exp(-((x-a)^2)/(2*b*b))
#' eapply.dyms(dym.in,dym.out,"gauss",a=29,b=2)   
#' @export
eapply.dyms <- function(dyms.in,dyms.out,FUN,...){ 

	if (!exists("dyms.out") | is.null(dyms.out) | any(is.na(dyms.out)) | length(dyms.out)!=length(dyms.in)) 
	stop("Names of output files not specified. Exit.")
	
	func <- get(FUN)
	#func <- match.fun(FUN)
	for (f in 1:length(dyms.in)){
		dym<-read.var.dym(dyms.in[f],t.format="dym",verbose=TRUE)
		var<-dym$var; tc<-dym$t; x<-dym$x; y<-dym$y; mask<-dym$landmask
		var.out <- func(var,...)
		write.var.dym(dyms.out[f],tc,x,y,mask,var.out)
	}
}

#' Aggregating DYM variable through space and returning time series
#'
#' Function reads DYM file, computes the sum (default) or the average over selected region and returns the time series. The sum always implies conversion to the mass content and sensible only for biomass densities. For averaging, it can return a simple average without modifying the units of the DYM data, or convert the data to the mass content (e.g., for biomass densities) by multiplying by the cell area before averaging, or a weighted average with weights being the grid cell areas (e.g., for temperature). 
#' @param dym.in a string containing the full name (with full path or relative to home directory) of the DYM file.
#' @param ts the set of two dates in Date format c(t1,t2) for subsetting. If left undefined, the original time period will be used.
#' @param xs the set of longitudes c(x1,x2) for subsetting. If lest undefined, the original longitudinal range will be used.
#' @param ys the set of latitudes in format c(y1,y2) for subsetting. If left undefined, the original latitudinal range will be used.
#' @param FUN a name of a function for aggregation, can take values "sum", "mean", "wmean", "mmean" for a sum, a simple mean, weighted mean and mass mean respectively.
#' @param ivalue invalid value to be used to ignore data.
#' @returns Returns a list with two array variables, t being the time vector and y being the aggregated variable.
#' @examples 
#' dym.in <- system.file("extdata", "skj.dym", package = "dym")
#' txy <- read.txy.dym(dym.in)
#' txy
#' x.range <- c(100,180)
#' y.range <- c(-25,25)
#' dat <- get.ts.dym(dym.in,xs=x.range,ys=y.range); 
#; plot(dat$t,dat$y,type="l",main="Total skipjack biomass in selected region",xlab="",ylab="metric tons")
#' @export
get.ts.dym <- function(dym.in,ts=NULL,xs=NULL,ys=NULL,FUN="sum",ivalue=NULL){

	#first deal with dimensions for subsetting		
	dims <- check.dimensions(dym.in,ts,xs,ys)
	ts <- dims$t; xs <- dims$x; ys <- dims$y
	
	t0 <- NULL; tfin <- NULL
	if (!is.null(ts)){
		dates <- sapply(ts,function(x)format(x,c('%Y','%m','%d')))
		t0   <- as.numeric(dates[,1])
		tfin <- as.numeric(dates[,2])
	}

	#by applying mask, getting NA on land		
	dym<-read.var.dym(dym.in,t0,tfin,c(xs,ys),verbose=FALSE,apply.mask=TRUE)
	var<-dym$var; tc<-dym$t; x<-dym$x; y<-dym$y;

	if (!is.null(ivalue)) 
		var <- ifelse(var==ivalue,NA,var)
   
	lat.area<-cell.surface.area(y,60*abs(diff(x)[1]),60*abs(diff(y)[1]))
	
	func <- get("mean")
	if (FUN!="mmean") 
		func <- get("sum")

	corlat <- function(x,...) func(lat.area*t(x),...)

	if (FUN!="mean")
		res <- apply(var,1,corlat,na.rm=T)
	else res <- apply(var,1,mean,na.rm=T)
	
	if (FUN=="wmean") res <- res/(length(x)*sum(lat.area))
		
	#if (is.decimal(tc)) tc <- dec.2date(tc)
	return(list(t=tc,y=res))
}

#' Aggregating DYM variable through time and returning 2d variable
#'
#' Function reads DYM file, computes the average (default) or the sum over selected time period and returns the result as a 2d array and/or writes it to a DYM1 file. Note, the only purpose of writing DYM1 file is to visualize it in SeapodymView.
#' @param dym.in a string containing the full name (with full path or relative to home directory) of the DYM file.
#' @param dym.out a string giving the name of output DYM file if need to be written. In this case, it will be written in the DYM1 format for visualization in SeapodymView.
#' @param ts the set of two dates in Date format c(t1,t2) for subsetting. If left undefined, the original time period will be used.
#' @param xs the set of longitudes c(x1,x2) for subsetting. If lest undefined, the original longitudinal range will be used.
#' @param ys the set of latitudes in format c(y1,y2) for subsetting. If left undefined, the original latitudinal range will be used.
#' @param FUN a function to be used for aggregation, e.g., "sum" or "mean" (default).
#' ivalue invalid value to be used to ignore data.
#' @param out the logical value (FALSE by default) telling if the output should be returned by the function. It will be automatically turned to TRUE if dym.out is not provided.
#' @returns Optionally, returns a list with x, y coordinates and z, a matrix of values of the aggregated variable.
#' @examples 
#' # Example 1. Subsetting and computing average distribution, returning the result.
#' dym.in <- system.file("extdata", "skj.dym", package = "dym");
#' dat <- get.2d.dym(dym.in); 
#' image(dat$x,dat$y,dat$z,zlim=range(dat$z,na.rm=T),col=hcl.colors(hcl.pals("sequential")[11],n=15))
#' # Example 2. Get decadal average of primary production, removing cells with invalid values. SEAPODYM environment is needed for this example.
#' setwd(paste0(Sys.getenv("SEA_FORCING_HOME"),"/run-era5np_1x30d-po/"))
#' dym.in <- list.files(pattern="_pp_")
#' t.range <- c(as.Date("1981-1-1"),as.Date("1990-12-31"))
#' dat <- get.2d.dym(dym.in,ts=t.range,ivalue=-999);
#' image(dat$x,dat$y,dat$z,zlim=range(dat$z,na.rm=T))
#' # Plot with seamaps
#' require(seamaps)
#' image(dat$x,dat$y,dat$z,xlab='',ylab='',xaxt='n',yaxt='n',main="Primary production, 1981-1990")
#' par(mgp=c(3,0.3,0),tck=0.02)
#' nice.map.over(range(dat$x)+c(-.5,.5),range(dat$y)+c(-.5,.5))
#' @export
get.2d.dym <- function(dym.in,dym.out=NULL,ts=NULL,xs=NULL,ys=NULL,FUN="mean",ivalue=NULL,out=FALSE){

	dims <- check.dimensions(dym.in,ts,xs,ys)
	ts <- dims$t; xs <- dims$x; ys <- dims$y
	
		t0 <- NULL; tfin <- NULL
		if (!is.null(ts)){
			dates <- sapply(ts,function(x)format(x,c('%Y','%m','%d')))
		t0   <- as.numeric(dates[,1])
		tfin <- as.numeric(dates[,2])
	}

	#by applying mask, getting NA on land		
	dym<-read.var.dym(dym.in,t0,tfin,c(xs,ys),t.format="dym",verbose=FALSE,apply.mask=TRUE)
	var<-dym$var; tc<-dym$t; x<-dym$x; y<-dym$y; mask<-dym$landmask
	
	if (!is.null(ivalue)) 
		var <- ifelse(var==ivalue,NA,var)
		
	func <- get(FUN)
	res  <- t(apply(var, 2, function(x)apply(x,2,func,na.rm=T)))

	if (!is.null(dym.out)){
		#if (is.date(tc)) tc <- date.2dec(tc)		
		write.var.dym1(dym.out,tc[1],x,y,mask,res)
	} else out <- TRUE
	
	if (out) return(list(x=x,y=y,z=res))
}

#' Making climatology of a DYM variable
#'
#' Function reads DYM file and generates the climatology (average) over selected time period and months. The output DYM file will be written in the same directory and name of the input DYM file with added suffix 'clm_' and either 'monthly' (default), 'quarterly', 'seasonal' or 'other', depending on the selection of months.  
#' @param dym.in a string containing the full name (with full path or relative to home directory) of the DYM file.
#' @param imonths the set of first months indices for averaging. For example, 1:12 corresponds to monthly climatology, c(1,3,7,10) to quarterly. Note, the last month of each climatological interval is always the month preceding the first month of the next, or the first element in a set.
#' @param ts the set of two dates in Date format c(t1,t2) for subsetting. If left undefined, the original time period will be used.
#' @param out the logical value (FALSE by default) telling if the output should be returned by the function, in addition to be writted to the DYM file. 
#' @returns Optionally, returns a list with dimensional variables, x and y, as well as z, a climatological 3d variable.
#' @examples 
#' Example 1. Make monthly climatology of skipjack density 
#' dym.in <- system.file("extdata", "skj.dym", package = "dym")
#' make.climatology(dym.in) 
#'
#' Example 2. Make climatology over calendar seasons (dec-feb, mar-may, jun-aug, sep-nov) of primary production. SEAPODYM environment is needed to run this example.
#' setwd(paste0(Sys.getenv("SEA_FORCING_HOME"),"/run-era5np_1x30d-po/"))
#' dym.in  <- list.files(pattern="_pp_")
#' dym.in  <- unique(gsub("_clm_seasonal","",dym.in))
#' seasons <- c(12,3,6,9)
#' dat <- make.climatology(dym.in,seasons,out=TRUE); # writing to the DYM file [dym.in]_epi_clm_seasonal.dym and returning result
#' # Plot seasonal distributions showing values > 75 mmolC/m2/day in black
#' x <- dat$x; y <- dat$y
#' par(mfrow=c(2,2),mar=c(2,2,2,2),mgp=c(3,0.3,0),tck=0.02)
#' for (n in 1:4)
#' 	image(x,y,dat$z[n,,],zlim=c(0,75),
#' 			col=c(hcl.colors(hcl.pals("sequential")[24],n=15),1),
#'				main=paste(months(as.Date(paste0("2009-",seasons[n],"-1"))),"-",
#'							months(as.Date(paste0("2009-",seasons[c(2:4,1)][n]-1,"-1")))))
#' @export   
make.climatology <- function(dym.in,imonths=1:12,ts=NULL,out=FALSE){

	t0 <- NULL; tfin <- NULL
	if (!is.null(ts)){
		dates <- sapply(ts,function(x)format(x,c('%Y','%m','%d')))
		t0   <- as.numeric(dates[,1])
		tfin <- as.numeric(dates[,2])
	}

	#by applying mask, getting NA on land		
	dym<-read.var.dym(dym.in,t0,tfin,apply.mask=TRUE,verbose=TRUE)
	var<-dym$var; tc<-dym$t; x<-dym$x; y<-dym$y; mask<-dym$landmask
	#if (is.decimal(tc)) tc <- dec.2date(tc)
	
	nbm <- length(imonths)
	if (nbm>=length(tc)) stop("Bad time period selected, stop execution.")
	var.out <- var[1:nbm,,] #clone the first nbm matrices

	#get the intervals for climatology
	two.yrs <- rep(1:12,2)
	ms <- two.yrs[which(two.yrs==imonths[1])[1]:length(two.yrs)][1:12]
	months <- split(ms, cumsum(ms %in% imonths))
	tc.months <- as.numeric(format(tc,"%m"))
	for (i in 1:nbm){
		message("Averaging over months ",paste(months[[i]],collapse=" "))
		ind <- which(!is.na(match(tc.months,months[[i]])))
		var.out[i,,] <- t(apply(var[ind,,], 2, function(x)apply(x,2,mean,na.rm=T)))
	}

	#writing
	clm.type <- ifelse(nbm==12,"monthly",
					   ifelse(all(imonths==seq(1,12,3)),"quarterly",
							  ifelse(all(imonths==c(12,3,6,9)),"seasonal","other")))
	dym.out <- gsub(".dym",paste0("_clm_",clm.type,".dym"),dym.in)
	
	print(dym.out)
	#generate climatological time vector:
	tclm <- date.2dec(as.Date(paste0("1900-",seq(1,12,12/nbm),"-15")))
	write.var.dym(dym.out,tclm,x,y,mask,var.out)

	if (out) return(list(x=x,y=y,z=var.out))
}

