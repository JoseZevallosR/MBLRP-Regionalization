library(extRemes)
require(ggplot2)
library(raster)
library(devtools)
library("HyetosMinute")

#Mean
meanMBLRPM<-function(a,l,v,k,f,mx,h=1) {
  x<-(h*l*mx*v*(1+k/f))/(a-1) 
  return(x)
}
#Variance
varMBLRPM<-function(a,l,v,k,f,mx,h=1) {
  A<-(2*l*(1+k/f)*(mx^2)*(v^a))/((f^2)*((f^2)-1)*(a-1)*(a-2)*(a-3))
  B<-(2*(f^2)-2+k*f)*(f^2)*((a-3)*h*(v^(2-a))-(v^(3-a))+((v+h)^(3-a)))
  C<-k*(f*(a-3)*h*(v^(2-a))-(v^(3-a))+((v+f*h)^(3-a)))
  D<-A*(B-C)
  return(D)
}
#Covariance
covarMBLRPM<-function(a,l,v,k,f,mx,h=1,lag=1) {
  A<-(l*(1+k/f)*(mx^2)*(v^a))/((f^2)*((f^2)-1)*(a-1)*(a-2)*(a-3))
  B<-(2*(f^2)-2+k*f)*(f^2)*(((v+(lag+1)*h)^(3-a))-2*((v+lag*h)^(3-a))+((v+(lag-1)*h)^(3-a)))
  C<-k*(((v+(lag+1)*h*f)^(3-a))-(2*((v+h*lag*f)^(3-a)))+((v+(lag-1)*h*f)^(3-a))) 
  D<-A*(B-C)
  return(D)
}
#Dry probabilities
pdrMBLRPM<-function(a,l,v,k,f,h=1) {
  mt<-((1+(f*(k+f))-(0.25*f*(k+f)*(k+4*f))+((f/72)*(k+f)*(4*(k^2)+27*k*f+72*(f^2))))*v)/(f*(a-1))
  G00<-((1-k-f+1.5*k*f+(f^2)+0.5*(k^2))*v)/(f*(a-1))
  A<-(f+(k*(v/(v+(k+f)*h))^(a-1)))/(f+k)
  D<-exp(l*(-h-mt+G00*A)) 
  return(D)
}

#Bartlett-Lewis calibration
modeloBTL6p6=function(mean24,var24,cov24lag1,pdr24,var1=1,var3=1,var6=1,var12=1,var18=1,minimo=c(1.0001,0.001,0.001,0.001,0.001,0.001),maximo=c(5,0.1,5,1,1,20)){
  
  #Objective function
  fopt <- function(x) {
    a<-x[1];l<-x[2];v<-x[3];k<-x[4];f<-x[5];mx<-x[6]
    w1=1;w2=2;w3=1;w4=1;w5=1;w6=1;
    
    S1 <- w2*((varMBLRPM(a,l,v,k,f,mx,h=1)/var1)-1)^(2)
    
    S3<-w2*((varMBLRPM(a,l,v,k,f,mx,h=3)/var3)-1)^(2)
    #w3*((covarMBLRPM(a,l,v,k,f,mx,h=3,lag=1)/cov3lag1)-1)^(2)
    S6<-w2*((varMBLRPM(a,l,v,k,f,mx,h=6)/var6)-1)^(2)
    
    S12<-w2*((varMBLRPM(a,l,v,k,f,mx,h=12)/var12)-1)^(2)
    
    S18<-w2*((varMBLRPM(a,l,v,k,f,mx,h=18)/var18)-1)^(2)
    
    S24 <- w1*((meanMBLRPM(a,l,v,k,f,mx,h=24)/mean24)-1)^(2)+ w2*((varMBLRPM(a,l,v,k,f,mx,h=24)/var24)-1)^(2)+ w3*((covarMBLRPM(a,l,v,k,f,mx,h=24,lag=1)/cov24lag1)-1)^(2)+w4*((pdrMBLRPM(a,l,v,k,f,h=24)/pdr24)-1)^(2)
    
    
    S<-S24+S3+S6+S12+S18
    
    if(is.infinite(S)) {S<-10^8}
    if(is.na(S)) {S<-10^8} 
    return(S) 
  }

  # set the interior and exterior parameters bounds
  xmin <- minimo
  xmax <- maximo
  xlow <- minimo
  xup <- maximo
  
  modecal <- eas(n=6,m=30,xmin,xmax,xlow,xup,fn=fopt,maxeval=5000,ftol=1.e-10,ratio=0.99,pmut=0.95, beta=2,maxclimbs=5)
  modecal
  a<-modecal$bestpar[[1]];
  l<-modecal$bestpar[[2]];
  v<-modecal$bestpar[[3]];
  k<-modecal$bestpar[[4]];
  f<-modecal$bestpar[[5]];
  mx<-modecal$bestpar[[6]]
  # In order to use the derived parameters in the functions of HyetosR 
  # as well as in the classic version of Hyetos,please be sure that 
  # for parameters mx and sx the length units are millimeters (mm) 
  # and for parameters l, v, mx and sx the time units are days (d). 
  # For this reason, make the following unit conversions:
  
  #checar
  #l<-l*24 
  #v<-v/24 
  #mx<-mx*24
  
  # parameter set for implementation in HyetosR functions
  par <- c(a=a,l=l,v=v,k=k,f=f,mx=mx) 
  #print(modecal$bestval)
  par
}

nearpoints=function(mdist,n_neighbors=50){
  near=n_neighbors#Using 50 neigthborhs
  distance=matrix(NA,nrow = dim(mdist)[1],ncol = near)
  for (i in 1:dim(mdist)[1]){
    values=sort(mdist[i,])[2:(1+near)]
    for (j in 1:near){
      distance[i,j]=which(mdist[i,]==values[j])[1]
    }
  }
  distance
}

#Cross Validations using Inverse Distance Weigth
idwCV=function(datos,columna='a',power=2){
  x=datos
  coordinates(x) <- ~x+y
  proj4string(x)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  mdist <- distm(x)
  
  predichos=numeric()
  for (i in 1:dim(datos)[1]){
    info=datos[[columna]][-i]
    denominador=sum((1/mdist[i,-i])^power)
    predichos[i]=sum(info/mdist[i,-i]^power/denominador)
  }
  predichos
  rr=cbind(datos[c('x','y')],'var1.pred'=predichos,'observed'=datos[[columna]])
  rr$residual=rr$observed-rr[['var1.pred']]
  rr
}

#Cross Validations and recalculation of MBRLPM's parameters
repetitiveCV=function(veces=1,data,Stats,Lmin,Lmax,n_neighbors=50){
  #data contains the intial parameter estimation
  #stats is the rainfall statistics

  #Estaciones cambiantes de intervalos
  
  

  for (iter in 1:veces){
    print(paste("Number of cross validation iteration",as.character(iter)))
    
    data_help=data
    coordinates(data_help) <- ~x+y
    mdist <- distm(data_help)
    vecinos=nearpoints(mdist,n_neighbors)
    
    formulas=c(formula('a~1'),formula('l~1'),formula('v~1'),formula('k~1'),formula('f~1'),formula('mx~1'))

    for (k in 1:6){
      print(paste('Checking parameter',k))
      x <- idwCV(data,columna = c('a','l','v','k','f','mx')[k],power=2)
      
      #Checking region error
      sub=x[c('x','y','var1.pred','observed','residual')]
      sub$porcentaje=abs(sub$residual)*100/sub$observed
      sub$residual=NULL
      
      sub=cbind(sub,1:dim(sub)[1])
      names(sub)=c('x','y','var1.pred','observed','porcentaje','ID')

      #Threshold of errors
      errores=subset(sub,sub$porcentaje>=1)
      buenos=subset(sub,sub$porcentaje<1)
      
      
      n1=dim(errores)[1]
      print('Search region modification ...')
      for (j in errores$ID){
        
        cercanos=subset(buenos,  buenos$ID %in% vecinos[j,])
        if (dim(cercanos)[1]!=0){
          Lmin[k,j]=min(cercanos$observed)
          Lmax[k,j]=max(cercanos$observed)
          #print(max(cercanos$observed))
        }else{
          #checkar los valores mas cercanos
        }
        #IEsta=c(IEsta,idx)
      }
    } 
    
    n=dim(Stats)[1]
    parameters=matrix(data=NA,nrow =n,ncol = 6)
    for (i in 1:n){
      mmx=800
      momentos=Stats[i,]
      
      iter=0
      while(!(mmx>=2 & mmx<=40)){
        par=modeloBTL6p6(mean24 = momentos$mean24h,var24 = momentos$var24h,cov24lag1 = momentos$autocov24h
                         ,pdr24=momentos$dryperiod24h,var3=momentos$var3h,var1 = 1,var6=momentos$var6h,var12=momentos$var12h,var18=momentos$var18h
                         ,minimo = Lmin[,i],maximo = Lmax[,i])
        mmx=par[6]
        parameters[i,]=par
        
        iter=iter+1
        if (iter==20){break}
      }
      
    } 
    parameters=cbind(Stats[,1:2],parameters)
    names(parameters)=c('x','y','a','l','v','k','f','mx')
    data=parameters#check
    
  }
  
  list(parameters,Lmin,Lmax)
}
#################################################################
# Interpolation Functions########################################
#################################################################

RI=function(name,Estaciones,iter=1){
  #Ration interpolation
  coords=Estaciones[,c(1,2)]
  coordinates(coords)=~x+y
  proj4string(coords)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 '
  
  a=raster(name)
  names(a)='aT'
  grd.pts <- SpatialPixels(SpatialPoints((raster(name))))
  grd <- as(grd.pts, "SpatialGrid")
  
  
  a=raster::extract(a,coords)
  Estaciones$aT=a
  ratio=Estaciones[[2+iter]]/Estaciones$aT#ratio entre observado y satelite
  errores=cbind(Estaciones[,c(1,2)],err=ratio)
  coordinates(errores)=~x+y
  proj4string(errores)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 '
  proj4string(grd)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 '
  idwError=gstat::idw(formula('err~1'),errores,grd,idp=2-0)
  idwError<-idwError['var1.pred']
  gridded(idwError)=TRUE
  mapa=raster(idwError)
  mapa*raster(name)
}

RIDW=function(name,f.1=as.formula(a ~ aT),Estaciones){
  #Linear regression interpolation
  coords=Estaciones[,c(1,2)]
  coordinates(coords)=~x+y
  proj4string(coords)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 '
  
  a=raster(name)
  names(a)='aT'
  grd.pts <- SpatialPixels(SpatialPoints((raster(name))))
  grd <- as(grd.pts, "SpatialGrid")
  
  
  a=raster::extract(a,coords)
  Estaciones$aT=a
  

  # Run the regression model
  lm.1 <- lm( f.1, data=Estaciones)
  errores=cbind(Estaciones[,c(1,2)],err=lm.1$residuals)
  coordinates(errores)=~x+y
  proj4string(errores)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 '
  proj4string(grd)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 '
  idwError=gstat::idw(formula('err~1'),errores,grd,idp=2-0)
  idwError<-idwError['var1.pred']
  gridded(idwError)=TRUE
  mapa=raster(idwError)

  dat.1st <- log(raster(name))*lm.1$coefficients[2]+lm.1$coefficients[1]
  
  exp(dat.1st+mapa)
}
