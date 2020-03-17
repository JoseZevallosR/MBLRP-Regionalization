#https://mgimond.github.io/Spatial/interpolation-in-r.html
library(sp)
library(gstat)
library(raster)
library(HyetosMinute)
library(geosphere)
library(spm)
library(ggplot2)

options(warn=-1)
set.seed(123)
run=function(data,path,MaskShape,n_neighbors=20,iter=5){
  #data: contains the rainfall statistics
  #path: where to save the results
  #Maskshape: Shape form of the final results 
  
  
  n=dim(data)[1]

  #minimos y maximos de los parametros a calibrarse
  Lmin=matrix(c(1,0.001,0.001,0.001,0.0854,1),nrow = 6,ncol = n)
  Lmax=matrix(c(4,0.1,0.1,0.1,0.1,20),nrow=6,ncol=n)
  Lmin[5,]=runif(n,min = 0.01, max = 0.099)# 
  Lmax[1,]=runif(n,min = 2, max = 5)
  Lmax[2,]=runif(n,min = 0.2, max = 2)# inicializando 
  Lmax[3,]=runif(n,min = 0.2, max = 2)# inicializando 
  Lmax[4,]=runif(n,min = 0.2, max = 2)# inicializando 
  Lmax[5,]=runif(n,min = 0.2, max = 2)# inicializando 
  Lmax[6,]=runif(n,min = 5, max = 20)
  
  print('Calculating the initial parameters ...')
  #Estimacion Inicial de parametros
  parametros=matrix(data=NA,nrow =n,ncol = 6)
  
  
  for (i in 1:n){
    momentos=data[i,]
    parametros[i,]=modeloBTL6p6(mean24 = momentos$mean24h,var24 = momentos$var24h,cov24lag1 = momentos$autocov24h
                                ,pdr24=momentos$dryperiod24h,var3=momentos$var3h,var1 = 1,var6=momentos$var6h,var12=momentos$var12h,var18=momentos$var18h
                                ,minimo = Lmin[,i],maximo = Lmax[,i])
    
  }


  parameters=cbind(data[,1:2],parametros)
  names(parameters)=c('x','y','a','l','v','k','f','mx')
  
  print('Reptitive Cross Validations ...')
  info=repetitiveCV(veces = iter,parameters,data,Lmin = Lmin ,Lmax = Lmax,n_neighbors=15)
  
  
  #interpolation part
  formulas=c(formula('a~1'),formula('l~1'),formula('v~1'),formula('k~1'),formula('f~1'),formula('mx~1'))
  data=info[[1]]
  #saving the initial parameters
  write.table(data,paste0(path,'parametros1.csv'),sep = ',',row.names = F)
  coordinates(data) <- ~x+y
  proj4string(data)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  
  
  print('Saving the results ...')
  maps=stack()
  grd=readRDS('grilla.rds')
  for (k in 1:6){
    
    mapa=idw(formulas[[k]],data ,newdata=grd)
    mapa=mask(raster(mapa),MaskShape)
    maps <- addLayer(maps,mapa)
    
    writeRaster(mapa,paste0(path,'parametros','-',as.character(k),'.tif'),overwrite=TRUE)
  }
  titulo=c('a','l','v','k','f','mx')
  names(maps)=titulo
  
  maps
}

#Shape with the shape peru to make a mask
setwd("~/ILA_NORMA/Modelo/Regionalizacion/peru/")
peru=shapefile('StudyArea.shp')

#Importing the regionalization functionality
setwd("~/ILA_NORMA/Modelo/Regionalizacion")
source('RegionalizationModule.R')

#Mixed stats from gauge stations and corrected TRMM
Gauge=read.csv('ParametrosRegionMax.csv')
Gauge=Gauge[1:86,]
maps=run(Gauge,path="~/ILA_NORMA/Modelo/Regionalizacion/ParametrosEstaciones/",MaskShape=peru,iter=5,n_neighbors=10)
plot(maps)
plot(maps,col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))


#Stats of TRMM without correction
TRMM=read.csv('EstatsTRMM.csv')
coordinates(TRMM)=~x+y
proj4string(TRMM)='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
TRMM=over(peru,TRMM,returnList = TRUE)[[1]]
mapas=run(TRMM,path="~/ILA_NORMA/Modelo/Regionalizacion/TRMM_Regiones/",MaskShape=peru,n_neighbors=20,iter=3)
dev.new()
plot(maps)
plot(maps,col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))

