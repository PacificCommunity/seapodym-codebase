library(RColorBrewer)

generic.name<-function(name){
  sp.acronyms<-c("skj","bet","alb","yft","swo")
  sp.names<-c("Skipjack","Bigeye","Albacore","Yellowfin","Swordfish")
  return(sp.names[grep(name,sp.acronyms)])
}

modcols<-function(pal,coef){
  pal<-t(col2rgb(pal))*coef
  out<-rgb(ifelse(pal>255,255,pal),maxColorValue=255)
  return(out)
}

plot.monthly.mean.fluxes<-function(dat,reg,add.B=TRUE,age.title="")
{#IS 20220802: see modifications in plot.qtr.mean.fluxes. Need to integrate them here
  dates<-as.Date(labels(dat))
  mo.names<-substr(month.name,1,3)
  nb.reg<-sqrt(length(dat[[1]]))
  regions<-(1:nb.reg)[-reg]
  names.in<-paste("reg",regions,"in",sep=".")
  names.out<-paste("reg",regions,"out",sep=".")
  fluxes.in <-matrix(rep(0,12*(nb.reg-1)),nb.reg-1,12,byrow=T,
	      dimnames=list(names.in,mo.names))
  fluxes.out<-fluxes.in; rownames(fluxes.out)<-names.out

  B.reg<-array(0,12)
  for (m in 1:12){
    ind<-grep(paste("-",m,"-",sep=""),names(dat)) 
    mat.mean<-array(0,c(nb.reg,nb.reg))
    for (i in ind){
      mat<-matrix(dat[[i]],nb.reg,nb.reg,byrow=F)
      mat.mean<-mat.mean+mat/length(ind)
      #Biomass in the region before the movement:
      B.reg[m]<-sum(mat.mean[reg,])
    }
    fluxes.out[,m]<- -mat.mean[reg,regions]
    fluxes.in[,m]<- mat.mean[regions,reg]
  }
  net.fluxes<-colSums(fluxes.in+fluxes.out) # fluxes.out are negative
  message("In");  print(fluxes.in)
  message("Out"); print(fluxes.out)
  message("Net"); print(net.fluxes)
  max.sum1<-max(colSums(fluxes.in))
  max.sum2<- max(-colSums(fluxes.out))
  print(c(max.sum1,max.sum2))
  max.sum<-1.25*max(max.sum1,max.sum2)
  pal<-brewer.pal(12,"Set3") 
  mycol<-c(pal,modcols(pal,0.7),modcols(pal,1.3)) #36 colors max
  mycol<-mycol[regions]
  
  mylim<-c(-max.sum,max.sum)
  par(las=1,mar=c(3,3,5,3))
  barplot(fluxes.in,ylim=mylim,col=mycol); axis(2,col="white",lwd=3)
  barplot(fluxes.out,add=TRUE,col=mycol)
  abline(h=0,lwd=2)
  legend("topleft",paste("reg",regions,sep=" "),bty="n",fill=mycol,ncol=6)
  par(new=TRUE)
  plot(1:12,net.fluxes,type="b",lwd=1.5,xaxt="n",yaxt="n",xlab='',ylab='',ylim=mylim);
  if (add.B){
    par(new=TRUE)
    print(B.reg)
    plot(1:12,B.reg,type="b",lwd=2,xaxt="n",yaxt="n",xlab='',ylab='');axis(4)
  }
  title(paste0(generic.name(sp)," biomass fluxes for region ",reg," and ",age.title)) 	    
}

plot.qtr.mean.fluxes<-function(dat,reg,add.B=TRUE,age.title="")
{
	
  dates<-as.Date(labels(dat))
  qtrs <-as.numeric(gsub("Q","",quarters(dates)))
  qtr.names<-c("Q1","Q2","Q3","Q4")
  nb.reg<-sqrt(length(dat[[1]]))
  regions<-(1:nb.reg)[-reg]
  names.in<-paste("reg",regions,"in",sep=".")
  names.out<-paste("reg",regions,"out",sep=".")

  print(4*(nb.reg-1))

  fluxes.in <-matrix(rep(0,4*(nb.reg-1)),nb.reg-1,4,byrow=T,
		     dimnames=list(names.in,qtr.names))
  fluxes.out<-fluxes.in; rownames(fluxes.out)<-names.out
  B.reg<-array(0,4)

  for (qtr in 1:4){
    ind<-which(qtrs==qtr)
    mat.mean<-array(0,c(nb.reg,nb.reg))
    for (i in ind){
      mat<-matrix(dat[[i]],nb.reg,nb.reg,byrow=F)
      mat.mean<-mat.mean + mat/length(ind)
      #Biomass in the region before the movement:
      #it should be shifter by 3 months back to correspond 
      #to the model reference biomass (see function get.B.ts)
      #do not modify weights, so the total B corresponds 
      #to the model solution after time and age integration
      qtr.ref<-qtr-1; if (qtr.ref==0) qtr.ref<-4
      B.reg[qtr.ref]<-sum(mat.mean[reg,])
    }
    fluxes.out[,qtr]<- -mat.mean[reg,regions]
    fluxes.in[,qtr]<- mat.mean[regions,reg]
  }
  funit<-"metric tons/quarter"
  bunit<-"metric tons"
  if (mean(B.reg)>1000){
    B.reg<-0.001*B.reg
    fluxes.in<-0.001*fluxes.in
    fluxes.out<-0.001*fluxes.out
    funit<-"thous. metric tons/quarter"
    bunit<-"thous. metric tons"
  }

  net.fluxes<-colSums(fluxes.in+fluxes.out) # fluxes.out are negative
  message("In");  print(fluxes.in)
  message("Out"); print(fluxes.out)
  message("Net"); print(net.fluxes)
  max.sum1<-max(colSums(fluxes.in))
  max.sum2<- max(-colSums(fluxes.out))
  print(max.sum1); print(max.sum2)
  max.sum<-1.25*max(max.sum1,max.sum2)
  pal<-brewer.pal(12,"Set3")[c(2,1,3:12)]
  mycol<-c(pal,modcols(pal,0.7),modcols(pal,1.3)) #36 colors max
  mycol<-mycol[regions]

  mylim<-c(-max.sum,max.sum)
  if (add.B) mylim<-c(-1.2*max(B.reg),1.2*max(B.reg))
  par(las=1,mar=c(3,5,5,4),tck=-0.02,mgp=c(3,0.5,0))
  bp<-barplot(fluxes.in,ylim=mylim,col=mycol,space=0.5,border=mycol); 
  axis(2,col="white",lwd=3);  
  barplot(fluxes.out,add=TRUE,col=mycol,space=0.5,border=mycol)#"gray20")
  abline(h=0)
  legend("topleft",paste("reg",regions,sep=" "),bty="n",fill=mycol,border=mycol,ncol=6)
  lines(bp,net.fluxes,type="b",lwd=1.5);
  if (add.B){
    print(B.reg)
    col.B<-"gray40"
    par(new=TRUE,col.axis=col.B)
    plot(bp,B.reg,type="b",pch=20,cex=0.5,lty=3,lwd=2,xaxt="n",yaxt="n",xlab='',ylab='',
    	ylim=c(0.8*min(B.reg),1.2*max(B.reg)),xlim=c(bp[1]-0.5,bp[4]+0.5),col=col.B);
	axis(4,col=col.B)
    par(col.axis=1)	
  }
  par(las=0); 
  mtext(funit,2,line=3,cex=1.25)
  mtext(bunit,4,line=3,cex=1.25)
  title(paste0(generic.name(sp)," biomass flow rates for region ",reg," and ",age.title))
}

#' Plotting seasonal movement probabilities by age between two regions, either read from 'movement-matrices-for-MFCL.txt' file or computed from SEAPODYM outputs.   
#' @param move_ij is the list containing the age vector and the seasonal movement probabilities at age from region 'i' to region 'j', stored as a 2d array with dimensions 4 (seasons) and the number of age classes. The age structure can be as defined in the simulation, or aggregated into quarterly classes, which is typical resolution of Multifan-CL.     
#' @param ri and rj are regional indices, which should be in the range from 1 to 'nreg', with 'nreg' being the total number of regions in the regional structure defined in seapodym_fluxes simulations. 
#' @examples 
#' Example with MFCL file
#' file.name<-"./movement-matrices-for-MFCL.txt"
#' i<-1; j<-2; 
#' move_ij<-get.movement.prob.by.age.season(file.name,i,j,nbr,nba)	
#' plot.movement.prob.by.age.season(move_ij,i,j)
#' @export plot.movement.prob.by.age.season
plot.movement.prob.by.age.season<-function(move_ij,ri,rj)
{		
  mycol<-brewer.pal(8,"Set2")[c(6,3,2,8)]; col.fg<-"grey40"
  png.file<-paste0(figs.dir,"movement_probability_r",ri,"-to-r",rj,".png")
  message("Plotting movement probabilities at age to ", png.file)
  png(png.file,600,400)
  par(las=1,xaxs="i",yaxs="i",mar=c(3.5,3.75,3.5,1),mgp=c(2.75,0.3,0),
      fg=col.fg,cex=1.25,col.axis=col.fg)
  age <-move_ij$a
  prob<-move_ij$p
  y.pretty<-pretty(prob)
  ymax <-max(y.pretty)
  plot(age,prob[1,],type="n",ylim=c(0,ymax),xlab="",ylab="",lwd=4,tck=1)
  for (qtr in 1:4) lines(age,prob[qtr,],col=mycol[qtr],lwd=4)
  mtext("Age class (quarterly)",1,1.5,col=col.fg,cex=1.25)
  mtext("Probability (quarterly)",2,2.5,col=col.fg,cex=1.25,las=0)
  mtext(paste0("Movement probability R",ri,"-->R",rj),3,line=2,cex=1.5,col="grey30")
  legend(floor(max(age)/5),1.12*ymax,legend=paste0("Q",1:4),lty=rep(1,4),
	 lwd=rep(4,4),col=mycol,ncol=4,bty="n",cex=1,xpd=TRUE)
  dev.off()
}


#' Plotting SEAPODYM fluxes for a given region in PDF file. 
#' 
#' The monthly or quatrerly mean values averaged over the time series written in the output files, or the monthly or quarterly fluxes for a selected year. This function uses global environment variables 'sp' and 'figs.dir'. If 'sp' was not set, it will be reset from the file names, which are conventional. Same for figs.dir, it will be created in the working directory if it was not defined earlier. 
#' @param dir is the SEAPODYM output directory containing files [spname]_FluxesRegion_age[a].txt, where a is the age class indices from 0 to A+.
#' @param reg is the region for which the fluxes, i.e., biomass flow per deltaT from and to 'reg', will be plotted.
#' @param age is the vector of age class indices, selected for plotting. Must be a sub-set of age classes written in the outputs.
#' @param year by default is set to NULL, but it can be set to a value within the temporal range in the outputs.
#' @param treso is the temporal resolution of fluxes that will be plotted, can be either "monthly" or "quarterly". If the outputs are written as quarterly, both options will provide the same quarterly visualization. If outputs are monthly and treso is quarterly, then the monthly fluxes will be aggregated into quarterly. 
#' @param add.B is logical parameter, if set to TRUE, then the total biomass in the region 'reg' before the movement will be added to the plot on the right-hand y axis.
#' @examples 
#' plot.fluxes.region("./output/",1,3:10) # plot fluxes from and into region 1 for only immature adults (if skipjack)
#' @export plot.fluxes.region
plot.fluxes.region<-function(dir,reg,age=0:18,year=NULL,treso="quarterly",add.B=FALSE)
{
  files<-check.vars(dir)
  pdf.name<-paste(figs.dir,"/",sp,"_fluxes_reg",reg,".pdf",sep="")
  if (!is.null(year)) 
    pdf.name<-paste(figs.dir,"/",sp,"_fluxes_year",year,"_reg",reg,".pdf",sep="")
  pdf(pdf.name,7,5)
  for (a in age){
    ind<-grep(paste0("_age",a,".txt"),files)
    if (length(ind)==0) {dev.off(); stop(paste("Error:: No file for age",a))}
  
    dat<-read.data.fluxes(files[ind],year)
    age.title<-paste0("age ",A," months","\n (mean length ",L,"cm, mean weight ",W,"kg)")
    if (treso=="monthly")
        plot.monthly.mean.fluxes(dat,reg,add.B,age.title)
    if (treso=="quarterly")
        plot.qtr.mean.fluxes(dat,reg,add.B,age.title)
  }
  dev.off()
}

#' Plotting SEAPODYM fluxes for a given region in PDF file. 
#' 
#' The monthly or quatrerly mean values averaged over the time series written in the output files, or the monthly or quarterly fluxes for a selected year. This function uses global environment variables 'sp' and 'figs.dir'. If 'sp' was not set, it will be reset from the file names, which are conventional. Same for figs.dir, it will be created in the working directory if it was not defined earlier. 
#' @param dir is the SEAPODYM output directory containing files [spname]_FluxesRegion_age[a].txt, where a is the age class indices from 0 to A+.
#' @param reg is the region for which the fluxes, i.e., biomass flow per deltaT from and to 'reg', will be plotted.
#' @param age is the vector of age class indices, selected as constituting the life stage. Must be a sub-set of age classes written in the outputs. It the responsibility of the user to select the ages that correspond to the life stage. Selected range will be included in the plot title
#' @param year by default is set to NULL, but it can be set to a value within the temporal range in the outputs.
#' @param areso is the age resolution of fluxes that will be plotted, can be either "monthly" or "quarterly". 
#' @param treso is the temporal resolution of fluxes that will be plotted, can be either "monthly" or "quarterly". If the outputs are written as quarterly, both options will provide the same quarterly visualization. If outputs are monthly and treso is quarterly, then the monthly fluxes will be aggregated into quarterly. 
#' @param add.B is logical parameter, if set to TRUE, then the total biomass in the region 'reg' before the movement will be added to the plot on the right-hand y axis.
#' @examples 
#' plot.fluxes.region("./output/",1,3:10) # plot fluxes from and into region 1 for only immature adults (if skipjack)
#' @export plot.fluxes.region.lstage
plot.fluxes.region.lstage<-function(dir,reg,age=0:18,lstage.name="all",
					year=NULL,areso="monthly",treso="quarterly",add.B=FALSE)
{ 
  files<-check.vars(dir)
  if (areso=="monthly")
    age.title<-paste0(lstage.name,"\n(",paste0(min(age)," - ",max(age)+1)," months of age)")
  if (areso=="quarterly")
    age.title<-paste0(lstage.name,"\n(",paste0(range(age),collapse=" - ")," quarters of age)")
  lstage.name<-gsub(" ","-",lstage.name)
  pdf.name<-paste(figs.dir,"/",sp,"_fluxes_reg",reg,"_",lstage.name,".pdf",sep="")
  if (!is.null(year)) 
    pdf.name<-paste(figs.dir,"/",sp,"_fluxes_year",year,"_reg",reg,"_",lstage.name,".pdf",sep="")

  pdf(pdf.name,7,5)
  
  counter<-0
  for (a in age){
    ind<-grep(paste0("_age",a,".txt"),files)
    if (length(ind)==0) {dev.off(); stop(paste("Error:: No file for age",a))}
    counter<-counter+1
  
    dat<-read.data.fluxes(files[ind],year)
    dates<-as.Date(labels(dat))
    nbt<-length(dates)
    if (counter==1)
      dat.sum<-dat
    if (counter>1){
      for (ti in 1:nbt)
        dat.sum[[ti]]<-dat.sum[[ti]]+dat[[ti]]
    } 
  }
  if (treso=="monthly")
      plot.monthly.mean.fluxes(dat.sum,reg,add.B,age.title)
  if (treso=="quarterly")
      plot.qtr.mean.fluxes(dat.sum,reg,add.B,age.title)

  dev.off()
}


#To be described later: this function to be used once plotting the movements with uncertainties
#@examples 
#to read and plot the output files with movement probabilities at age
#move.est<-scan(file.out,skip=2,nlines=1)
#move.min<-scan(file.out,skip=4,nlines=1)
#move.max<-scan(file.out,skip=6,nlines=1)
#polygon.movement.prob.by.age(move.est,move.min,move.max,1,2,3,qtr.col=mycol[1])

polygon.movement.prob.by.age<-function(move_ij_mean,move_ij_min,move_ij_max,qtr,ri,rj,age=1:50,qtr.col=1)
{		
  png.file<-paste0(figs.dir,"movement_probability_r",ri,"-to-r",rj,"_qtr",qtr,"_with-errs.png")
  png(png.file,600,300)
  par(las=1,xaxs="i",yaxs="i",mar=c(3.5,3.75,3.5,1),mgp=c(2.75,0.3,0),
      fg=col.fg,cex=1.25,col.axis=col.fg)
  nba<-length(age)
  y.pretty<-pretty(move_ij_max)
  ymax <-max(y.pretty)
  plot(age,move_ij_mean,type="n",ylim=c(0,ymax),xlab="",ylab="",tck=1)
  #polygon(c(age,rev(age)),c(move_ij_min,move_ij_max[nba:1]),col=rgb(0.9,0.9,0.9),border="gray80")
  #lines(age,move_ij_mean,col=qtr.col,lwd=1)
  polygon(c(age,rev(age)),c(move_ij_min,move_ij_max[nba:1]),col=qtr.col,border="gray80")
  lines(age,move_ij_mean,lwd=2)
  mtext("Age class (quarterly)",1,2,col=col.fg,cex=1.25)
  mtext("Probability (quarterly)",2,2.5,col=col.fg,cex=1.25,las=0)
  mtext(paste0("Movement probability R",ri,"-->R",rj),3,line=2,cex=1.5,col="grey30")
  legend(floor(max(age)/2.5),1.2*ymax,legend=paste0("Q",qtr),
	 lty=1,lwd=4,col=qtr.col,ncol=4,bty="n",cex=1,xpd=TRUE)
  dev.off()
}

