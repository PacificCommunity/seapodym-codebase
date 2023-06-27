#The first version of this function was written by Mary Borderies and later modified and updated by Inna Senina
# **********************************************************************************************************

    #' Writing DYM file
    #'
    #' Writes data to DYM file
    #' @param file.out is the name of the DYM file.
    #' @param tvect is the date vector in decimal format.
    #' @param x is the vector of longitude.
    #' @param y is the vector of latitude.
    #' @param mask is the land mask - the matrix of integers 0 - 2, with 0 for land, 1 for the sea with only one layer (nearshore), 2 for the sea with two layers and 3 for the sea with three pelagic layers though the depth column.
    #' @param data is the 3d array of data to be written in DYM file.
    #' @param verbose flag (default value is TRUE) controlling the prompt for the nominal function execution.
    #' @examples 
    #' write.var.dym("sst.dym",tt,x,y,mask,var);
    #' @export
    write.var.dym<-function(file.out,tvect,x,y,mask,data,verbose=TRUE){

	if (verbose)
	  message("Writing file ",file.out,"...")
	#initialization
	Idfunc<-0;
	nlon<-length(x)
	nlat<-length(y)
	nlevel<-length(tvect)
	xlon<-matrix(x,nrow=nlat,ncol=nlon,byrow=T)
	ylat<-matrix(y,nrow=nlat,ncol=nlon,byrow=F)
	ylat<- apply(t(ylat), 1, rev)

	#writing
	con<-file(file.out,"wb") 
	writeChar("DYM2",con,eos=NULL)
	writeBin(Idfunc,con,size=4)
	writeBin(min(data,na.rm=T),con,size=4)
	writeBin(max(data,na.rm=T),con,size=4)	
	writeBin(nlon,con,size=4)
	writeBin(nlat,con,size=4)
	writeBin(length(tvect),con,size=4)
	writeBin(tvect[1],con,size=4)
	writeBin(tvect[length(tvect)],con,size=4)

	for(i in 1:nlat){
		writeBin(xlon[i,],con,size=4)
	}
	for(i in 1:nlat){
		writeBin(ylat[i,],con,size=4)
	}
	writeBin(tvect,con,size=4) 

	for(i in 1:nlat){
		writeBin(as.integer(mask[i,]),con,size=4)#ATTN: not recognized by C and Java. Need to find a solution!!!
	}
			
	#control in case on nans
	data<-ifelse(is.na(data),-999,data)
	if (length(dim(data))==3 & dim(data)[1]!=nlevel){
          message("WARNING: number of matrices is not equal to nlevel!")		
  	  data[,,]<-data[1:nlevel,,]
	}
	for(ti in 1:nlevel){
		if (length(dim(data))==3)
		  dat<-apply(t(data[ti,,]),2,rev) 
		if (length(dim(data))==2)
		  dat<-apply(t(data),1,rev) 
		for(i in 1:nlat){
			writeBin(dat[i,],con,size=4) 
		}
	}		
	close(con)
    }

    #' Writes DYM file with population density initial conditions
    #' @param file.out is the name of the DYM file.
    #' @param x is the vector of longitude.
    #' @param y is the vector of latitude.
    #' @param mask is the land mask - the matrix of integers 0 - 2, with 0 for land, 1 for the sea with only one layer (nearshore), 2 for the sea with two layers and 3 for the sea with three pelagic layers though the depth column.
    #' @param data is the 3d array of data to be written in DYM file with the structure data[age,x,y]. The mean age of each age class is not written currently in this file, although it can be passed to zlevel. In the current version, zlevel is filled with age indices.
    #' @examples 
    #' write.restart.dym("skj_cohorts.dym",x,y,mask,var);
    #' @export
    write.restart.dym<-function(file.out,x,y,mask,data){

	#initialization
	Idfunc<-0;
	nlon<-length(x)
	nlat<-length(y)
	nlevel<-dim(data)[1]
	zlevel<-seq(1,nlevel,1.0)

	xlon<-matrix(x,nrow=nlat,ncol=nlon,byrow=T)
	ylat<-matrix(y,nrow=nlat,ncol=nlon,byrow=F)
	ylat<- apply(t(ylat), 1, rev)
	
	#writing
	con<-file(file.out,"wb") 
	writeChar("DYM2",con,eos=NULL)
	writeBin(Idfunc,con,size=4)
	writeBin(min(data,na.rm=T),con,size=4)
	writeBin(max(data,na.rm=T),con,size=4)	
	writeBin(nlon,con,size=4)
	writeBin(nlat,con,size=4)
	writeBin(nlevel,con,size=4)
	writeBin(1,con,size=4)
	writeBin(zlevel[nlevel],con,size=4)

	for(i in 1:nlat){
		writeBin(xlon[i,],con,size=4)
	}
	for(i in 1:nlat){
		writeBin(ylat[i,],con,size=4)
	}
	writeBin(zlevel,con,size=4) 

	for(i in 1:nlat){
		writeBin(mask[i,],con,size=4)
	}
			
	#control in case on nans
	data<-ifelse(is.na(data),-999,data)
	if (dim(data)[1]!=nlevel){
          message("WARNING: number of matrices is not equal to nlevel!")		
  	  data[,,]<-data[1:nlevel,,]
	}

	for(ti in 1:nlevel){
		dat<-apply(t(data[ti,,]),2,rev) 
		for(i in 1:nlat){
			writeBin(dat[i,],con,size=4) 
		}
	}		
	close(con)
    }

