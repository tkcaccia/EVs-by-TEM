# EVs-by-TEM

'''
fff
'''
link <- "/media/user/Lion/Documents/X-ray images selection/"
require("readxl")
setwd(link)

files=list.files(link)
files=files[which(!grepl( "-res", files,fixed = TRUE))]
files=files[which(!grepl( "png", files,fixed = TRUE))]
files=files[which(!grepl( "pdf", files,fixed = TRUE))]
output=gsub("\\.","-res.",files)

res=as.data.frame(read_excel("../Resolution.xlsx"))
res2=as.data.frame(read_excel("../Resolution.xlsx",sheet = 2))
res=res[res[,"Patient.ID"] %in% res2[,"Patient.ID"],]
rownames(res)=res[,1]
files=intersect(files,rownames(res))
res=res[files,]



size=list()
kk=1
for(kk in kk:length(files)){
  print(kk)
  setwd(link)
  
  c2 <- readImage(files[kk])
  
  
  
  c2=channel(c2,mode = "gray")
  c2[1:(512+128),(2048-128):2048]=median(c2)
  c2=c2-min(c2)
  c2=c2/max(c2)
  
  
  # creates an object called c2
  
  # show the image in the R graphics window
  display(c2, method = "raster")  #"raster" method means within R
  
  
  # make a brighter image by multiplying all the values by 2
  c2.b <- c2
  
  # gaussian blur
  c2.b.blur <- gblur(c2.b, sigma = 10)
  display(c2.b.blur, method = "raster")
  
  c2.b.blur[c2.b.blur<0]=0
  
  z=imageData(c2.b.blur)
  
  ma=matrix(ncol=32,nrow=32)
  for(i in 1:32){
    for(j in 1:32){
      ma[i,j]=median(z[(i-1)*64+1:64,(j-1)*64+1:64])
    }
  }
  
  display(z, method = "raster")
  display(ma, method = "raster")
  
  ll=loess(as.numeric(ma)~cbind(rep(1:32,32),rep(1:32,each=32)),span = 0.01)#,control=loess.control(surface="direct"))
  z2=matrix(predict(ll,cbind(rep(1:32,32),rep(1:32,each=32))),ncol=32)
  
  z3=matrix(ncol=2048,nrow=2048)
  for(i in 1:32){
    for(j in 1:32){
      z3[(i-1)*64+1:64,(j-1)*64+1:64]=z2[i,j]
    }
  }
  z3=gblur(z3, sigma = 100)
  
  display(z2, method = "raster")
  #z3=matrix(predict(ll,cbind(rep(1:2048,2048)/64,rep(1:2048,each=2048)/64)),ncol=2048)
  
  display(z3, method = "raster")
  
  
  
  display(4*(z-z3), method = "raster")
  
  zz=z-z3
  q=quantile(zz,prob=0.2)
  zz[zz<q]=q
  
  c2.b.blur.thres <- zz >  otsu(zz,range = range(zz)) # apply this value 
  display(c2.b.blur.thres, method = "raster")
  
  z4=gblur(c2.b.blur.thres, sigma = 5)
  
  z4 <- z4 > otsu(z4) # apply this value 
  
  display(z4, method = "raster")
  
  
  c2.b.blur.thres=fillHull(z4)
  
  
  c2.b.blur.thres.cnt <- bwlabel(c2.b.blur.thres)
  
  # show this as a coloured blobs 
  
  
  y=colorLabels(c2.b.blur.thres.cnt)
  display(y, method = "raster")
  
  rr=computeFeatures.shape(c2.b.blur.thres.cnt)
  
  w=which((rr[,"s.radius.sd"]/rr[,"s.radius.mean"]>0.4) )
  
  c2.b.blur.thres[!is.na(match(c2.b.blur.thres.cnt,w))]=0
  c2.b.blur.thres.cnt[!is.na(match(c2.b.blur.thres.cnt,w))]=0
  
  y=colorLabels(c2.b.blur.thres.cnt)
  display(y, method = "raster")
  
  
  nmask = watershed( distmap(c2.b.blur.thres), 2 )
  display(colorLabels(nmask))
  
  rr=computeFeatures.shape(nmask)
  
  w=which((rr[,"s.radius.sd"]/rr[,"s.radius.mean"]>0.15) | rr[,"s.area"]*res[kk,"nm2.pixel"]<150)
  nmask[!is.na(match(nmask,w))]=0
  display(colorLabels(nmask))
  size[[kk]]=sqrt((table(nmask)[-1]*res[kk,"nm2.pixel"])/pi)*2
  writeImage(colorLabels(nmask),output[kk],compression = "LZW")
  
  save(size,file="/home/user/Desktop/size.RData")
}







load("/home/user/Desktop/size.RData")
pat=list()
patID=unique(res[,"Patient.ID"])


number_photos=NULL
for(un in 1:length(patID)){
  w=which(patID[un]==res[,"Patient.ID"] & is.na(res[,"remove"]))
  number_photos[un]=length(w)
  pat[[un]]=as.numeric(unlist(size[w]))
  #pat[[un]][pat[[un]]>3000]=3000
}

pat=pat[-which(number_photos==0)]
patID=patID[-which(number_photos==0)]
number_photos=number_photos[-which(number_photos==0)]

for(i in 1:length(pat)){
  pdf(paste("1 - Hist-",patID[i],".pdf",sep=""),height = 4)
  hh=hist(pat[[i]],breaks = 0:60*2,xlim=c(0,120),ylim=c(0,500),axes=FALSE,xlab="Diameter (nm)",main=patID[i])
  axis(1)
  axis(2,las=2)
  counts=hh$counts
  breaks=hh$breaks[-1]-1
  ll=loess(counts~breaks,control=loess.control(surface="direct"),span = 0.15)
  x=seq(10,120,1)
  z=predict(ll,x)
  z[z<0]=0
  z[1]=0
  points(x,z,type="l",col=3,lwd=2)
  dev.off()
}




#################################################



x=seq(10,50,1)
fit=matrix(nrow=length(pat),ncol=length(x))

for(i in 1:length(pat)){
  hh=hist(pat[[i]],breaks = 0:60*2,xlim=c(0,120),ylim=c(0,500),axes=FALSE,xlab="Diameter (nm)",main=patID[i])
  
  counts=hh$counts
  breaks=hh$breaks[-1]-1
  
  ll=loess(counts~breaks,control=loess.control(surface="direct"),span = 0.15)
  
  z=predict(ll,x)
  
  z[z<0]=0
  fit[i,]=z
  
}

colnames(fit)=x
rownames(fit)=patID
##########################################################


yA = function (x) dnorm(x, mean = 16, sd =  3)
yB = function (x) dnorm(x, mean = 19, sd =  3)
yC = function (x) dnorm(x, mean = 22, sd = 3)
yD = function (x) dnorm(x, mean = 25, sd = 3)
yE = function (x) dnorm(x, mean = 28, sd = 3)
yF = function (x) dnorm(x, mean = 31, sd = 3)
yG = function (x) dnorm(x, mean = 34, sd = 3)
yH = function (x) dnorm(x, mean = 37, sd = 3)

lf=list()
lf[[1]]=yA
lf[[2]]=yB
lf[[3]]=yC
lf[[4]]=yD
lf[[5]]=yE
lf[[6]]=yF
lf[[7]]=yG
lf[[8]]=yH


number_of_functions=length(lf)

fr <- function(co) {   ## Rosenbrock Banana function
  
  sommat=rep(0,length(x))
  
  for(hh in 1:number_of_functions)
    sommat=sommat+co[hh]*lf[[hh]](x)
  
  #dnorm(x, mean = 28, sd = 2.8)
  sum((y-sommat)^2)
}



fit3=fit/number_photos

decov.fit=matrix(nrow=nrow(fit3),ncol=number_of_functions)
rownames(decov.fit)=rownames(fit3)
for(i in 1:length(pat)){
  #  png(paste("3 - Dec -",patID[i],".png",sep=""),width = 1000)#height = 4)
  pdf(paste("3 - Dec -",patID[i],".pdf",sep=""),height = 4)
  
  x=as.numeric(colnames(fit3))
  y=fit3[i,]
  oo=optim(par=rep(1,number_of_functions),
           lower = rep(0,number_of_functions),
           upper = rep(10^6,number_of_functions), fr,method="L-BFGS-B")
  
  
  plot(x,fit3[i,],type="l",axes=FALSE,main=patID[i],xlab="Diameter (nm)",ylab="Density (a.u.)")
  axis(1)
  axis(2,las=2)
  box()
  for(jj in 1:number_of_functions)
    points(x,oo$par[jj]*lf[[jj]](x),type="l",col=jj)
  
  somma=rep(0,length(x))
  
  for(jj in 1:number_of_functions)
    somma=somma+oo$par[jj]*lf[[jj]](x)
  points(x,somma,type="l",col=6,lwd=3)
  decov.fit[i,]=oo$par
  dev.off()
}
colnames(decov.fit)=c("16 nm","19 nm","22 nm","25 nm","28 nm","31 nm","34 nm","37 nm")

rownames(res2)=res2[,"Patient.ID"]
res3=res2[rownames(decov.fit),]
multi_analysis(decov.fit,as.factor(res3[,"Pathology"]))

