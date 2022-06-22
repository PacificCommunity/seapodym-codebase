get.month<-function(num){
 names<-c('January','February','March','April','May','June','July',
              'August','September','October','November','December')
 return(names[num])
}

sp_name<-function(sp, FULL=TRUE){
   if (sp=="cjm") return("jack mackerel")
   if (FULL)
     species<-c(paste(c("skipjack","bigeye","yellowfin","albacore"),"tuna"),"swordfish")
   if (!FULL)
     species<-c("skipjack","bigeye","yellowfin","albacore","swordfish")
   for (i in 1:length(species)){
       if (!is.na(all(pmatch(substring(sp,1:3,1:3),unlist(strsplit(species[i],NULL)))))){ 
	 name<-species[i]; break;
       }
   }
   return(name)
}

xtoi<-function(x,xmin,dx) {
   return(round((x-xmin)/dx,digits=0)+1)
}

ytoj<-function(y,ymin,dy){
   return(round((ymin-y)/dy,digits=0)+1)
}

xtoi2<-function(xmin,x,dx)
{
  return(round((x-xmin)/dx,digits=0)+1)
}

ytoj2<-function(ymin,y,dy)
{
  return(round((ymin-y)/dy,digits=0)+1)
}


sum2d<-function(vars)
{#use this function for the catch files
  nt<-dim(vars)[1]
  SUM<-array(NA,nt)
  for (i in 1:nt){
      SUM[i]<-sum(vars[i,,],na.rm=TRUE)
  }
  SUM<-ifelse(SUM==0,NA,SUM)
  return(SUM)
}

mean2d<-function(vars)
{
  nt<-dim(vars)[1]
  MEAN<-array(NA,nt)
  for (i in 1:nt){
      MEAN[i]<-mean(vars[i,,],na.rm=TRUE)
  }
  MEAN<-ifelse(MEAN==0,NA,MEAN)
  return(MEAN)
}


sum1d<-function(vars)
{
  nx<-dim(vars)[2]
  ny<-dim(vars)[3]
  SUM<-array(NA,c(nx,ny))
  for (i in 1:nx){
    for (j in 1:ny){
      SUM[i,j]<-sum(vars[,i,j],na.rm=TRUE)
    }
  }
  SUM<-ifelse(SUM==0,NA,SUM)
  return(SUM)
}

cell.surface.area<-function(lat,dx,dy)
{#returns the area of a cell on a sphere in sq.km
     R = 6378.1;
     Phi1 = lat*pi/180.0;
     Phi2 = (lat+dy/60.0)*pi/180.0;
     dx_radian = (dx/60.0)*pi/180;
     S = R*R*dx_radian*(sin(Phi2)-sin(Phi1));

     return(S)
}

cell.surface.area.2<-function(lat,dx,dy)
{
     g = (lat * pi) / 180.0;# transform lat(deg) in lat(radian)
	
     S<-dx*dy*1.852^2*cos(g)

     return(S)
}

sum.nb<-function(x,y,var)
{#use this function for the files with the density in nb/km2 units 
 #the output will be in mln nb
  nx<-length(x)
  dx<-60*(x[2]-x[1])
  dy<-60*(y[2]-y[1])
  area<-cell.surface.area(y,dx,dy)
  SUM<-1e-6*sum(t(var)*area,na.rm=TRUE)
  SUM<-ifelse(SUM==0,NA,SUM)
  return(SUM)
}


sum1d.nb<-function(x,y,vars)
{#use this function for the files with the density in nb/km2 units 
 #the output will be in mln nb
  nt<-dim(vars)[1]
  nx<-length(x)
  dx<-60*(x[2]-x[1])
  dy<-60*(y[2]-y[1])
  SUM<-array(NA,nt)
  area<-cell.surface.area(c(y),dx,dy)
  for (i in 1:nt){
    SUM[i]<-1e-6*sum(t(vars[i,,])*area,na.rm=TRUE)
  }
  SUM<-ifelse(SUM==0,NA,SUM)
  return(SUM)
}

sum1d.d2b<-function(x,y,vars)
{#from 2d density in vars[nt,nx,ny] to time series of total biomass (result in mt)
  nt<-dim(vars)[1]
  nx<-length(x)
  dx<-60*(x[2]-x[1])
  dy<-60*(y[2]-y[1])

  SUM<-array(NA,nt)
  area<-cell.surface.area(y,dx,dy)
  print(area)
  for (i in 1:nt){
    SUM[i]<-sum(t(vars[i,,])*area,na.rm=TRUE)
  }
  SUM<-ifelse(SUM==0,NA,SUM)
  return(SUM)
}

mean1d.d2b<-function(x,y,vars)
{#from 2d density in vars[nt,nx,ny] to time series of average biomass (result in mt)
 #Attn, vector 'y' should have the same order as vars!	
  nt<-dim(vars)[1]
  nx<-length(x)
  dx<-60*(x[2]-x[1])
  dy<-60*(y[2]-y[1])
  MEAN<-array(NA,nt)
  area<-cell.surface.area(y,dx,dy)
  for (i in 1:nt){
    MEAN[i]<-mean(t(vars[i,,])*area,na.rm=TRUE)
  }
  MEAN<-ifelse(MEAN==0,NA,MEAN)
  return(MEAN)
}

gen.monthly.dates<-function(t0,tfin){

 years<-seq(t0[1],tfin[1],1)
 for (y in years){
   dyear<-as.Date(paste(y,1:12,15,sep="-"))
   is.out<-(dyear<as.Date(paste(t0[1],t0[2],15,sep="-"))|
	    dyear>as.Date(paste(tfin[1],tfin[2],15,sep="-")))
   if (any(is.out))
     dyear<-dyear[!is.out]
   if (y==years[1]) dates<-dyear
   if (y!=years[1]) dates<-c(dates,dyear)
 }
 return(dates)
}

dates.to.monthly<-function(dates){
  years<-sort(unique(format(dates,"%Y")))
  for (y in years){
    dyear<-as.Date(paste(y,1:12,15,sep="-"))
    is.out<-(dyear<dates[1] | dyear>dates[length(dates)])
    if (any(is.out))
      dyear<-dyear[!is.out]
    if (y==years[1]) dates.out<-dyear
    if (y!=years[1]) dates.out<-c(dates.out,dyear)
  }
  return(dates.out)
}

#Aggregating MONTHLY data to YEARLY time step:
get.yearly<-function(t.mo,var){		
  t.yr<-as.numeric(format(t.mo,"%Y"))		
  Y<-sort(unique(t.yr))
  var.out<-array(NA,length(Y))
  for (i in 1:length(Y)){
    ind<-which(t.yr==Y[i])
    var.out[i]<-sum(var[ind])
  }
  return(list(t=Y,V=var.out))
}


get.ts.catch<-function(file,t0=c(1990,1),tfin=c(1990,12),coords=c(lonmin,lonmax,latmin,latmax))
{#in Mt
  data<-read.var.reg.dym(file,t0,tfin,coords)
  if (is.null(data)) return(NULL)
  x<-data$xc; y<-data$yc; y<-y[length(y):1]; c.time<-data$t
  Var.sum<-1e-6*sum2d(data$var) #simple sum over spatial dimensions, returns time series
  return(list(t=c.time,V=Var.sum))
}

get.ts.lifestage<-function(file,t0=c(1990,1),tfin=c(1990,12),coords=c(lonmin,lonmax,latmin,latmax))
{
  data<-read.var.reg.dym(file,t0,tfin,coords)
  if (is.null(data)) return(NULL)
  x<-data$xc; y<-data$yc; y<-y[length(y):1]; c.time<-data$t
  Var.stage<-sum1d.nb(x,y,data$var)# this function doesn't change units, i.e. 
  #if in the file the unit is nb/km2 it will give the total number, if mt/km2 then the result is mt.
  #print(Var.stage)
  return(list(t=c.time,V=Var.stage))
}

  DL<-function(jday,lat){
    rum <- (jday-80)/365.25;
    delta <- sin(rum*pi*2)*sin(pi*23.5/180);
    codel = asin(delta);
    phi   = lat*pi/180;
    argu  = tan(codel)*tan(phi);
    argu  = pmin(1.,argu);
    argu  = pmax(-1.,argu);
    dl    = 24.0-2.0*acos(argu)*180.0/pi/15 ;
    dl    = pmax(dl,0.0)
    return(dl)
  }


