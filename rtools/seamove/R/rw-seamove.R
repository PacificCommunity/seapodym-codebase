check.vars<-function(dir)
{
  files<-list.files(paste0(dir,"/"),pattern="FluxesRegion")	
  if (!file.exists(paste0(dir,files[1]))){
    message("Error:: Can't find fluxes file in directory ",dir)
    stop("")
  }
  sp.file<<-unlist(strsplit(files[1],"_"))[1]
  if (!exists("sp")){ 
    sp<<-sp.file
    message("Resetting species name to ",sp.file)
  }
  if (!exists("figs.dir")) figs.dir<<-"./figs/"; 
  if (!dir.exists(figs.dir)) dir.create(figs.dir)
  return(paste0(dir,"/",files))
}

get.files<-function(dir)
{
  files<-list.files(paste0(dir,"/"),pattern="FluxesRegion")	
  if (!file.exists(paste0(dir,files[1]))){
    message("Error:: Can't find fluxes file in directory ",dir)
    stop("")
  }
  return(paste0(dir,"/",files))
}

get.spname<-function(filename)
{
  line1<-scan(filename,what="#",quiet=T,n=5)
  sp.name<-line1[5]
  if (!exists("sp") | is.na(sp) | sp!=sp.name) sp<<-sp.name
}

#' Reading the header of the SEAPODYM output file with movement matrices.
#' @param Filename is the full name (with path) of the file to read.
#' @return Returns the list with species name, mean age and length, as well as number of regions and dates written in the 'fluxes' files.
#' @export
read.data.fluxes.header<-function(filename)
{
 dat.all<-scan(filename,what="#",quiet=T)
 sp<<-dat.all[5]
 A<<-as.numeric(dat.all[8])
 L<<-as.numeric(dat.all[11])
 W<<-as.numeric(dat.all[14])
 nbr<-as.numeric(dat.all[16])
 dates<-as.Date(dat.all[which(dat.all=="#")+1])

 return(list(sp=sp,age=A,len=L,wei=W,nreg=nbr,t=dates))
}

get.years<-function(dates) as.numeric(format(dates,"%Y"))

#' Reading the SEAPODYM output file with movement matrices.
#' @param filename is the full name (with path) of the file to read.
#' @param years can be given in case the movement matrices should be extracted for a given sub-set of years. The default value is NULL. 
#' @return Returns the matrices with SEAPODYM movement fluxes.
#' @examples
#' filename<-"output/skj_FluxesRegion_age31.txt"
#' dat<-read.data.fluxes(filename,2001:2010)
#' @export
read.data.fluxes<-function(filename,years=NULL)
{
 hea<-read.data.fluxes.header(filename)
 nbr<-hea$nreg
 dates<-hea$t

 lskip<-6; lread<-1e5
 if (!is.null(years)){
   ind<-which(!is.na(match(get.years(dates),years)))
   dates<-dates[ind]
   lskip<-lskip+(ind[1]-1)*(nbr+1)
   lread<-length(dates)*nbr
 }
 dat<-read.table(filename,skip=lskip,nrows=lread)
 mat<-split(as.matrix(dat),gl(length(dates),nbr,labels=dates))

 return(mat)
}

read.mfcl.move.file<-function(filename,nbr,nba)
{
  dat<-read.table(filename,skip=1,nrows=1e5)
  file.labs<-paste("Movement period ",rep(1:4,each=nba),rep(paste(" age class ",1:nba),4))
  mat<-split(as.matrix(dat),gl(length(file.labs),nbr,labels=file.labs))	
  return(mat)
}

#' Reading the SEAPODYM output file with movement matrices and reduces their dimensions to the MFCL regional structure, given the provided mapping between MFCL and SEAPODYM regions, where SEAPODYM regions follow: nbr.sea > nbr.mfcl, total.area(regs.sea) = total.area(regs.mfcl), and reg.mfcl.i = sum_j(regs.sea.j). The function will write the SEAPODYM output files with reduced regional structure in the new output directory.
#' @param dir.in is the SEAPODYM output directory containing files [spname]_FluxesRegion_age[a].txt, where a is the age class indices from 0 to A+, with the original regional structure defined in simulations.
#' @param dir.out is the new output directory containing files [spname]_FluxesRegion_age[a].txt with new regional structure.
#' @param age is the vector of age class indices, determines the window of age classes to be extracted from the outputs. Note, the age indices must be a sub-set of age classes written in the outputs.
#' @param reg.mfcl is the list with MFCL regions, where each element of the list contains the indices of the seapodym regions.
#' @examples
#' # Example #1:
#' dir.in<-"./output"
#' dir.out<-"./output/8regs/"
#' reg.mfcl<-list(r1=c(1,2),r2=c(3,4),r3=c(5,6),r4=7,r5=c(8,9),r6=c(10,11),r7=c(12,13,14),r8=15)
#' reg.structure.mapping(dir.in,dir.out,1:48,reg.mfcl)
#' 
#' # Example #2:
#' dir.out<-"output/5regs/"
#' reg.mfcl<-list(r1=1:7,r2=c(8,9),r3=c(10,11),r4=c(12,13,14),r5=15)
#' reg.structure.mapping(dir.in,dir.out,1:48,reg.mfcl)
#' @export
reg.structure.mapping<-function(dir.in,dir.out,age,reg.mfcl)
{
    nbr<-length(reg.mfcl)
    files.in<-get.files(dir.in)
    files.out<-gsub(dir.in,dir.out,files.in)
    
    message("Regional mapping... ")
    pb<-txtProgressBar(style=3)
    for (a in age){
        setTxtProgressBar(pb,(a-1)/length(age))
        ind<-grep(paste0("_age",a,".txt"),files.in)
        if (length(ind)==0) {dev.off(); stop(paste("Error:: No file for age",a))}
        dat<-read.data.fluxes(files.in[ind])
        hea<-read.data.fluxes.header(files.in[ind])

        dat.age<-hea$age
        dat.len<-hea$len
        dat.wei<-hea$wei

        #writing the same information to the header with updated nb. of regions:
        fileout<-files.out[ind]
        write(paste("Regional fluxes predicted for",hea$sp),fileout)
        write(paste("Mean age(months)",dat.age),fileout,append=T)
        write(paste("Mean length(cm)",dat.len),fileout,append=T)
        write(paste("Mean weight(kg)",dat.wei),fileout,append=T)
        write(paste("Nb.regions",nbr,"\n"),fileout,append=T)

        nbr.sea<-sqrt(length(dat[[1]]))
 
        nbt<-length(dat)
        dates<-labels(dat)
        for (t in 1:nbt){
            dat.t<-matrix(dat[[t]],nbr.sea,nbr.sea,byrow=FALSE)
            ij.rm<-NULL;
            for (reg in 1:nbr){
                regs.sea<-reg.mfcl[[reg]]
                if (t==1 & a==1) message("Aggregating seapodym regions ", paste(regs.sea,collapse=","))
  
                #sum-up the rows of seapodym matrix and put zeros to other sub-regions to keep the total sum
                r.cur<-regs.sea[1]; r.0<-regs.sea[-1]; ij.rm<-c(ij.rm,r.0)
  		          if (length(r.0)>0){
  		              dat.t[r.cur,]<-apply(dat.t[regs.sea,],2,sum); dat.t[r.0,]<-0
  		              #same for columns
  		              dat.t[,r.cur]<-apply(dat.t[,regs.sea],1,sum); dat.t[,r.0]<-0
  		          } else {
  		              #if only one seapodym region for a mfcl region
  		              dat.t[r.cur,]<-dat.t[regs.sea,]
  		              dat.t[,r.cur]<-dat.t[,regs.sea]
  		          }
            }
            #now reduce the matrix dimension to remove zero lines and columns
            for (i in rev(ij.rm)) { dat.t<-dat.t[-i,]; dat.t<-dat.t[,-i]}
    
	    write(paste("#",dates[t]),fileout,append=T)
	    write.table(dat.t,fileout,append=T,quote=FALSE,row.names=FALSE,col.names=FALSE)
	}
    }  
    close(pb)
    message("done.")
}

#' Aggregating the age structure of the SEAPODYM movement fluxes to a coarser age structure used in the Multifan-CL model configuration. 
#' @param dir.in is the SEAPODYM output directory containing files [spname]_FluxesRegion_age[a].txt, where a is the age class indices from 0 to A+. 
#' @param dir.out is the new output directory containing files [spname]_FluxesRegion_age[a].txt with the age structure as defined in Multifan-CL. 
#' @param a0, am are the first and last index of age classes to be extracted from the outputs in dir.in and aggregated to coarser age structure. The vector given by a0:am should be a sub-set of age classes in the outputs.
#' @param a.plus has the logical value, indicating whether the last age class in age.in is the A+. If yes, it will be used as is in the aggregation to Multifan-CL structure. Otherwise, it will be included in the aggregation.
#' @param aggregate.by is the size of the Multifan-CL age class in SEAPODYM age class units. Note, the Multifan-CL age indices start with 1 and not 0 as original SEAPODYM files.  
#' @return Returns the number of aggregated age classes.
#' @examples
#' # Example of skipjack configuration with 49 age classes, aggregating
#' # classes older than 3 months of age, i.e., skipping the first three 
#' # monthly classes:
#' dir.in<-"SKJ-DIR/output"
#' dir.out<-"SKJ-DIR/output-mfcl/"
#' nba<-aggregate.age.classes(dir.in,dir.out,a0,am=48,age.plus=TRUE,aggregate.by=3)
#' @export
aggregate.age.classes<-function(dir.in,dir.out,a0,am,a.plus=TRUE,aggregate.by=3)
{
  if ((!a.plus & (am-a0+1)%%aggregate.by!=0)| (a.plus & (am-a0)%%aggregate.by!=0))
    message("WARNING: the number of age classes is not divisible of ",aggregate.by,
    	     ". The last age class (or preceding A+) will be smaller!")

  if (!dir.exists(dir.out)) dir.create(dir.out)
  
  ages.new<-c(seq(a0,am-1,aggregate.by),am)
  if (!a.plus) ages.new<-seq(a0,am,aggregate.by)
  nba<-length(ages.new)

  message("Age structure aggregation... ")
  
  #need to know the species
  if (!exists("sp")){
    ff<-check.vars(dir.in)
    message("Found ",length(ff),"fluxes files in the dir.in for ",sp)
  }
    
  for (age in 1:nba){
    a1<-ages.new[age]
    a2<-min(a1+aggregate.by-1,am)
    #A+ can be aggregated with one or two previous age classes in this routine
    message("New age class ",age, " aggregating ages ", a1,":",a2)
    da<-a2-a1+1
  
    counter<-0;
    for (a in a1:a2){
      counter<-counter+1
      filename<-paste(dir.in,"/",sp,"_FluxesRegion_age",a,".txt",sep="")
      dat<-read.data.fluxes(filename)
      hea<-read.data.fluxes.header(filename)
      if (counter==1){
        dat.sum<-dat
        out.age<-hea$age
        out.len<-hea$len
        out.wei<-hea$wei
      }
      if (counter>1){
	      out.age<-out.age+hea$age
	      out.len<-out.len+hea$len
	      out.wei<-out.wei+hea$wei
	      nbt<-length(dat)
        for (ti in 1:nbt)
          dat.sum[[ti]]<-dat.sum[[ti]]+dat[[ti]]
      }
    }
    if (age<nba){
      out.age<-round(out.age/da,2)
      out.len<-round(out.len/da,2)
      out.wei<-round(out.wei/da,2)
    }

    #now writing
    fileout<-paste(dir.out,"/",sp,"_FluxesRegion_age",age,".txt",sep="")
    write(paste("Regional fluxes predicted for",hea$sp),fileout)
    write(paste("Mean age(months)",out.age),fileout,append=T)
    write(paste("Mean length(cm)",out.len),fileout,append=T)
    write(paste("Mean weight(kg)",out.wei),fileout,append=T)
    write(paste("Nb.regions",hea$nreg,"\n"),fileout,append=T)

    nbt<-length(dat)
    dates<-labels(dat.sum)
    ncol<-hea$nreg
    for (ti in 1:nbt){
      write(paste("#",dates[ti]),fileout,append=T)
      write(c(t(matrix(dat.sum[[ti]],ncol,ncol))),fileout,append=T,ncolumns=ncol)
    }
  }
  message("done.")
  return(nba)
}

#' Computing four seasonal vectors of movement probabilities by age for ri-->rj directly from SEAPODYM output files [spname]_FluxesRegion_age[a].txt. See also function get.movement.prob.by.age.season(), which reads the movement probabilities previously exported to the Multifan-CL format.
#' @param dir is the SEAPODYM output directory containing files [spname]_FluxesRegion_age[a].txt, where a is the age class indices from 0 to A+.
#' @param ri and rj are two regions, so that the vector of probabilities at age will be derived for ri->rj movement.
#' @param age is the vector of selecled age class indices, must be a sub-set of age classes written in the outputs.
#' @param print.out is the logical parameter, can be useful in case of large files to follow progression. 
#' @return Returns the vector of probabilities of movement from region ri to rj for the age classes 'age'.
#' @examples 
#' i<-1; j<-2; age<-0:47
#' move_ij<-comp.movement.prob.by.age.season(dir.out,i,j,age)
#' @export
comp.movement.prob.by.age.season<-function(dir,ri,rj,age=0:50,print.out=FALSE)
{
  #compute movement probability by quarter and age class	
  move_ij<-array(0,c(4,length(age)))	
  if (print.out)
      message("Getting movement probabilities at age ",appendLF=FALSE)  
  for (a in age){
    if (print.out) message(a,"..",appendLF=FALSE)
    filename<-paste0(dir,sp,"_FluxesRegion_age",a,".txt")	  
    dat<-read.data.fluxes(filename)
    if (a==age[1]) {
	    hea<-read.data.fluxes.header(filename)
	    nbr<-hea$nreg
    }
    dates<-as.Date(labels(dat))
    qtrs <-as.numeric(gsub("Q","",quarters(dates)))
    for (qtr in 1:4){
      ind<-which(qtrs==qtr)
      dat.sum<-dat[[ind[1]]]
      for(n in ind[-1])
        dat.sum<-dat.sum+dat[[n]]
        
      mat<-matrix(dat.sum,nbr,nbr)
      pmat<-mat/rowSums(mat)
      move_ij[qtr,a]<-pmat[ri,rj]
    }
  }
  if (print.out) message(".done.")
  return(move_ij)              
}

#' Function reads the Multifan-CL file with movement probabilities with seasonal and age structure.
#' @param filename is the name of the Multifan-CL file with movement probabilities.
#' @param ri and rj are two regions, so that the vector of probabilities at age describe ri->rj movement.
#' @param age is the vector of Multifan-CL age class indices, beginning with 1.
#' @return Returns the vector of probabilities of movement from region ri to rj for the age classes 'age'.
#' @examples 
#' i<-1; j<-2; age<-0:47
#' move_ij<-get.movement.prob.by.age.season(mfcl.file,i,j,age)
#' @export
get.movement.prob.by.age.season<-function(filename,ri,rj,nbr,age)
{ #Note, in MFCL-formatted file, the ri->rj is written in column 'ri' and row 'rj', 
  #while in the SEAPODYM files, the ri->rj is written in row rj and column rj.
  
  nba<-length(age)
  move_ij<-array(0,c(4,nba))	
  
  mat<-read.mfcl.move.file(filename,nbr,nba)
  file.labs<-labels(mat)
  for (qtr in 1:4){
    for (a in age){
      ind<-which(file.labs==paste("Movement period ",qtr," age class ",a))
      pmat<-matrix(mat[[ind]],nbr)
      move_ij[qtr,a]<-pmat[rj,ri] 
    }
  }
  return(move_ij) 
}

#' Converting the SEAPODYM movement fluxes to the Multifan-CL movement probabilities. 
#' @param dir.in is the SEAPODYM output directory containing files [spname]_FluxesRegion_age[a].txt, where a is the age class indices from 0 to A+. It can either be the original output with the original regional and age structures defined in simulations, or the aggregated outputs.
#' @param dir.out is the new output directory containing files [spname]_FluxesRegion_age[a].txt with the regional and age structures defined in Multifan-CL. The file with movement probabilities will be written in this directory as well.
#' @param years can be given in case the movement matrices should be extracted for a given sub-set of years. The default value is NULL. 
#' @param age.in is the vector of age class indices, determines the full or a sub-set of age classes to be extracted from the outputs written in dir.in. 
#' @param age.plus has the logical value, indicating whether the last age class in age.in is the A+. If yes, it will be used as is in the aggregation to Multifan-CL structure. Otherwise, it will be included in the aggregation.
#' @param aggregate.age is the size of the Multifan-CL age class in SEAPODYM age class units. If '0', then no aggregation will be done. 
#' @param nbr.in is the number of regions in the SEAPODYM outputs written in dir.in.
#' @param aggregate.reg has the logical value indicating whether the SEAPODYM regional structure should be mapped to Multifan-CL. See reg.mfcl and function reg.structure.mapping() for more details.
#' @param reg.mfcl is the list with MFCL regions, where each element of the list contains the indices of the seapodym regions. Default value is NULL.
#' @param plot.out activates plotting of resulting movement probabilities for Multifan-CL regions and age structure. Defalt value is TRUE. Note, the PNG files with plots will be written in figs.dir in the working directory.
#' @examples
#' # Example #1: skipjack model with 48 monthly age classes and 15 regions, 
#' # extracting movement rates for age classes older than 3 months of age,
#' # i.e., skipping first three monthly classes, the 48th being the A+:
#' # First, define the parent SKJ.DIR, where the outputs are stored.
#' dir.in<-paste0(SKJ.DIR,"/output/")
#' dir.out<-paste0(SKJ.DIR,"/output-mfcl/")
#' # Second, provide the regional mapping between SEAPODYM and Multifan-CL regions
#' reg.mfcl<-list(r1=c(1,2),r2=c(3,4),r3=c(5,6),r4=7,r5=c(8,9),r6=c(10,11),r7=c(12,13,14),r8=15)
#' # Third, call the function:
#' sea2mfcl.movement(dir.in,dir.out,1983:2010,age.in=2:47,age.plus=TRUE,aggregate.age=3,nbr.in=15,aggregate.reg=TRUE,reg.mfcl)
#'
#' #Example #2: albacore model configuration
#' dir.in<-paste0(ALB.DIR,"/output/")
#' dir.out<-paste0(ALB.DIR,"/output-mfcl/")
#' figs.dir<-paste0(ALB.DIR,"figs/figs-libtest/")
#' sea2mfcl.movement(dir.in,dir.out,1979:2010,age.in=5:148,age.plus=FALSE,aggregate.age=3,nbr.in=4,aggregate.reg=FALSE) 
#' @export
sea2mfcl.movement<-function(dir.in,dir.out,years=NULL,
			    age.in,age.plus=TRUE,aggregate.age=3,
			    nbr.in,aggregate.reg=FALSE,reg.mfcl=NULL,
			    plot.out=TRUE)
{
  nbr<-nbr.in
  if (aggregate.reg){
    nbr<-length(reg.mfcl)
    dir.int<-paste0(dir.in,"/output_",nbr,"regs/")
    if (!dir.exists(dir.int)) dir.create(dir.int)
    reg.structure.mapping(dir.in,dir.int,age.in,reg.mfcl)
    dir.in<-dir.int
  }

  nba<-length(age.in)
  age.out<-age.in
  files<-check.vars(dir.in)

  if (aggregate.age>0)
    nba.out<-aggregate.age.classes(dir.in,dir.out,age.in[1],age.in[nba],age.plus,aggregate.age)
  else if (aggregate.age==0){
    nba.out<-nba
    dir.out<-dir.in
  }
  age.out<-1:nba.out #aggregated classes indices begin with 1
  files<-get.files(dir.out)

  mfcl.file<-paste0(dir.out,"movement-matrices-for-MFCL.txt")
  cat("# movement matrices", "\n", file = mfcl.file)

  dat<-read.data.fluxes(files[1],years)
  dates<-as.Date(labels(dat)); 
  nby<-length(years)
  qtrs <-as.numeric(gsub("Q","",quarters(dates)))
  for (qtr in 1:4){
    ind<-which(qtrs==qtr)
    for (age in age.out){
      f.ind<-grep(paste0("_age",age,".txt"),files)
      dat<-read.data.fluxes(files[f.ind],years)
      #sum-up through all years for quarter 'qtr'
      dat.out<-dat[[ind[1]]]; for (i in ind[2:nby]) dat.out<-dat.out + dat[[i]]
      mat.sea<-matrix(dat.out,nbr,nbr)
      #now transpose to write in MFCL format: rows - in, columns - out.
      p<-t(mat.sea/rowSums(mat.sea)); 
      
      cat("# Movement period ",qtr," age class ",age, "\n", file = mfcl.file, append = T)
      write.table(round(p,10),file=mfcl.file,sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
    }
  }
  message("Written the movement probabilities to ",mfcl.file)
  
  if (plot.out){
    for (i in 1:nbr){
      for (j in 1:nbr){
        move.ij<-get.movement.prob.by.age.season(mfcl.file,i,j,nbr,age.out)	
        plot.movement.prob.by.age.season(list(a=age.out,p=move.ij),i,j)
      }
    }
    message("Done.")
  }
}

