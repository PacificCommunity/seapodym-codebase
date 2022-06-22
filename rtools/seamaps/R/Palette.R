#R function to add palette to the color map

library(RColorBrewer)

get.bins<-function(Ncols){
 #choice of nb of bins and step:
 #select between 5 and 15 bins
 n<-0; rem<-0; bp<-5:15
 for (i in bp) if (Ncols%%i==0) n<-i
 if (n==0){
   rem<-10; 
   for (i in bp) if (Ncols%%i<rem) rem<-Ncols%%i
   for (i in bp) if ((Ncols-rem)%%i==0) n<-i
 }
 step<-(Ncols-rem)/n
 return(step)
}


plot.palette<-function(values=mybreaks,cols=mycol,orientation="h",mycex=2.0,
		       #mymar=c(2,6,1,2),mymgp=c(1,0.85,0),ptext=''){
		       mymar=c(3,5,0,2),mymgp=c(1,0.85,0),ptext=''){
		       #mymar=c(3,6,1,2),mymgp=c(1,1.0,0),ptext=''){
		       #mymar=c(5,4,3,4),mymgp=c(1,1.25,0),ptext=''){
# if (length(values)>10) values<-pretty(values,4) # do not use pretty as it may change min and max values!!!
 if (length(values)>5) values<-seq(min(values),max(values),length.out=5) 
 i<-array(0.1,length(values)); 
 Ncols<-length(cols)
 step<-get.bins(Ncols); n<-floor(Ncols/step)
 d<-(max(values)-min(values))/n;
 op <- par(no.readonly = TRUE) # remember old settings
 par(cex.axis=mycex,mgp=mymgp,tck=-0.2,las=1)

 if (orientation=="v"){
   par(mar=mymar)
   plot(i,values,type="n", xaxt="n", yaxt="n", xlab="",ylab="",las=2,bty="n")
   axis(4)
   usr=par("usr"); x1<-usr[1]; x2<-usr[2]; y1<-usr[3]; y2<-usr[4]
   ys<-.5*(y2+y1-d*n)
   for (k in 1:n)
     rect(x1,(d*(k-1)+ys),x2,(d*k+ys),col = cols[k*step],border=FALSE)

   rect(x1,ys,x2,(d*n+ys),border="black",col=NA)	  
 }
 if (orientation=="h"){
  par(mar=mymar,bty="n",font=2)
  plot(values,i,type="n",yaxt="n", xaxt="n", xlab="",ylab="",las=0)
  usr=par("usr"); x1<-usr[1]; x2<-usr[2]; y1<-usr[3]; y2<-usr[4]; 
  xs<-.5*(x2+x1-d*n)
  if (ptext=='')
    axis(1,at=(values+xs-values[1]),labels=values)
  if (ptext!=''){
    axis(1,at=(values[1:n]+xs-values[1]),labels=values[1:n])
    axis(1,at=values[n+1]+xs-values[1],labels=ptext)
  }
  cind<-(1:n-.5)*step
  if (length(cols)%%2) cind<-(1:n)*step
  print(cind)
  for (k in 1:n)
    rect((d*(k-1)+xs),y1,(d*k+xs),y2,col = cols[cind[k]],border=FALSE)
  rect(xs,y1,(d*n+xs),y2,border="black",col=NA)
 } 
 par(op)
#usr<-par("usr")
#message("Exiting palette: usr = ",str(usr))  
}


mypretty<-function(var,n=10){
  
  rvar<-range(var,na.rm=TRUE)
  res<-array(NA,n)
  dv<-(rvar[2]-rvar[1])/(n-1)
  res<-rvar[1]+0:(n-1)*dv
  res<-round(res,digits=2)

  return(res) 
}

palettes<-function(n=60,pname="jet",user.pal=c(1,2,3)){

 if (pname=="jet"){
  jet.colors<-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  return (jet.colors(n))
 }
 if (length(grep("brewer",pname))>0){
   #pname=brewer.<palette name>, e.g. brewer.RdBu
   pname<-substring(pname,8)
   pnmax<-brewer.pal.info[pname,1]
   if (pnmax>n) pnmax<-n
   pal<-brewer.pal(pnmax,pname)
   pal.int<-colorRampPalette(pal)
   return(pal.int(n))
 }
 if (pname=="user"){
   pal.int<-colorRampPalette(user.pal)
   return(pal.int(n))
 }

 if (pname=="lab"){
  YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
  YlOrBr.Lab <- colorRampPalette(YlOrBr, space = "Lab")
  return (YlOrBr.Lab(n))
 }

 if (pname=="terrain")
  return(terrain.colors(n))

 if (pname=="two"){
  jetlab.colors<-colorRampPalette(c("#00007F", "blue", "#007FFF","cyan","#7FFF7F","#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404"))
  return(jetlab.colors(n))
 }
 if (pname=="jetpart"){
  jetpart.colors<-colorRampPalette(c("#00007F", "blue", "#007FFF","cyan"))
  return(jetpart.colors(n))
 }
 if (pname=="labpart"){
  labpart.colors<-colorRampPalette(c("#FFFFD4","#FED98E", "#FE9929", "#D95F0E", "#993404"))
  return(labpart.colors(n))
 }
}



