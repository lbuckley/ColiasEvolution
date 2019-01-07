#ESTIMATE LAMBDAS FOR COOP DATA AS FUNCTION OF ABSORPTIVITY
# APPLY TO NIWOT DATA
#USE ERBS TO PARTITION RADIATION

#CHECK WIND SPEED IN BIOPHYSICAL MODEL 

#Plant heights
heights= c(0.2, 0.5) #in m 

memory.size(max=TRUE)
memory.limit(size = 4095)

library(zoo)
library(RMAWGEN)

library( fields)
library( evd)
library( evdbayes)
library( ismev) 
library(chron) #convert dates
library(gdata)
library(maptools)
library(epicalc) 
library(RAtmosphere)
library(msm)
library(MASS)
library(SDMTools)
library(gdata) #for converting to matrix

#source biophysical model and other functions
#setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Butterflies\\Evolution\\Analysis\\")
setwd("E:\\Work\\Butterflies\\Evolution\\Analysis\\")
source("ColiasBiophysMod_27May2014_wShade.R")
source("ColiasBiophysMod_20May2014.R")

#source dtr function
source("DTRfunction.R")
#source zenith angle calculating function
source("ZenithAngleFunction.R")
#load radiation functions
source("RadiationModel_10Jun2014.R")

#microclimate model
#setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Butterflies\\Evolution\\MicroclimateModel\\")
setwd("E:\\Work\\Butterflies\\Evolution\\MicroclimateModel\\")
source("soil_temp_function_27May2014_wShade.R")
#source("soil_temp_function_noextrap_20May2014.R")
source('air_temp_at_height_z_function.R') #calculate air temp at some height z

#define geometric mean
geo_mean <- function(data) {
    log_data <- log(data)
    gm <- exp(mean(log_data[is.finite(log_data)]))
    return(gm)
}
#---------------------------------------
#LOAD CLIMATE DATA NIWOT RIDGE
#setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Butterflies\\Evolution\\ClimateData\\Niwot\\")
setwd("E:\\Work\\Butterflies\\Evolution\\ClimateData\\Niwot\\")
sites= read.csv("NiwotSites.csv")
a1.clim= read.csv("A1climate.csv")
b1.clim= read.csv("B1climate.csv")
c1.clim= read.csv("C1climate.csv")
d1.clim= read.csv("D1climate.csv")
stats= sites[,1]

#calculate as 10yr running averages
clim=c1.clim
Yr=1970
J=190
Yr.J= c(Yr, J)

clim.runmean= function(Yr.J, clim){
inds= which(clim$Julian %in% (Yr.J[2]-5):(Yr.J[2]+5) & clim$Year %in% (Yr.J[1]-5):(Yr.J[1]+5))
return( c( mean(clim[inds,"Min"], na.rm=TRUE),mean(clim[inds,"Max"], na.rm=TRUE)) )
}

Yr.J.mat= clim[,c("Year", "Julian")]
clim.out= t(apply(Yr.J.mat, MARGIN=1, FUN=clim.runmean, clim=clim))
clim$Min.ave=clim.out[,1]
clim$Max.ave=clim.out[,2]

par(mfrow=c(1,1), cex=1.3, lwd=2, mar=c(3, 3, 0.2, 0.2), mgp=c(1.5, 0.5, 0), oma=c(0,0,0,0), bty="l")
dates= as.Date(clim$Date, "%m/%d/%Y")
plot(dates,clim$Min, col="blue", type="l", ylim=range(-5,25), ylab="Temp", xlab="date")
points(dates,clim$Max, col="red", type="l")
points(dates,clim$Min.ave, col="green", type="l", lty="dashed")
points(dates,clim$Max.ave, col="orange", type="l", lty="dashed")

plot(clim$Min, col="blue", type="l", ylim=range(-5,25), ylab="Temp", xlab="date")
points(clim$Max, col="red", type="l")
points(clim$Min.ave, col="green", type="l", lty="dashed")
points(clim$Max.ave, col="orange", type="l", lty="dashed")

clim_runave=clim

#------------------------------------------
#WEATHER GENERATOR
#SEE ColiasFitnessExtremes_24Jun2014_NiwotClimate10yrAve FOR CODE
#https://github.com/ecor/RMAWGENCodeCorner

#---------------------
set.seed(1222) # Set the seed for random generation
data(trentino) # Load the dataset

#Replace data in Trentino dataset
#TEMPERATURE_MIN: Data frame containing year,month , day and daily minimum temperature in 59 stations in Trentino region
#TEMPERATURE_MAX: Data frame containing year,month , day and daily maximum temperature in 59 stations in Trentino region
#Setup: month, day, year, A1, B1, C1

#head(STATION_NAMES): Vector containing the names of the meteorological stations
#head(ELEVATION): Vector containing the elevations of the meteorological stations respectively
#head(STATION_LATLON): Matrix containing the latitude and longitude coordinates, respectively, of the meteorological stations

STATION_NAMES= sites$Site
ELEVATION= sites$Elev
STATION_LATLON= cbind(sites$Lat, sites$Lon)

date= as.Date(c1.clim$Date, "%m/%d/%Y")
c1.clim$Day= as.numeric(format(date, "%d"))
TEMPERATURE_MIN= c1.clim[,c("Year","Month","Day")]
TEMPERATURE_MAX= c1.clim[,c("Year","Month","Day")]
names(TEMPERATURE_MIN)[1:3]=c("year","month","day")
names(TEMPERATURE_MAX)[1:3]=c("year","month","day")

#Match other climate datasets
match1= match(a1.clim$Date, c1.clim$Date)
matches= na.omit(match1)
matched= which(!is.na(match1))
TEMPERATURE_MIN$A1=NA
TEMPERATURE_MAX$A1=NA
TEMPERATURE_MIN$A1[matched]<-  a1.clim$Min[matches]
TEMPERATURE_MAX$A1[matched]<-  a1.clim$Max[matches]

match1= match(b1.clim$Date, c1.clim$Date)
matches= na.omit(match1)
matched= which(!is.na(match1))
TEMPERATURE_MIN$B1=NA
TEMPERATURE_MAX$B1=NA
TEMPERATURE_MIN$B1[matched]<-  b1.clim$Min[matches]
TEMPERATURE_MAX$B1[matched]<-  b1.clim$Max[matches]

TEMPERATURE_MIN$C1= c1.clim$Min
TEMPERATURE_MAX$C1= c1.clim$Max

match1= match(d1.clim$Date, d1.clim$Date)
matches= na.omit(match1)
matched= which(!is.na(match1))
TEMPERATURE_MIN$D1=NA
TEMPERATURE_MAX$D1=NA
TEMPERATURE_MIN$D1[matched]<-  d1.clim$Min[matches]
TEMPERATURE_MAX$D1[matched]<-  d1.clim$Max[matches]

#--------------------------
#Generate temp data
year_min <- 1953
year_max <- 2011

year_min_sim <- 1953
year_max_sim <- 2011

n_GPCA_iter <- 5
n_GPCA_iteration_residuals <- 5
p <- 1
vstation <- c("C1")

generation00 <-ComprehensiveTemperatureGenerator(station=vstation,
 Tx_all=TEMPERATURE_MAX,Tn_all=TEMPERATURE_MIN,year_min=year_min,year_max=year_max,
 p=p,n_GPCA_iteration=n_GPCA_iter,n_GPCA_iteration_residuals=n_GPCA_iteration_residuals,
 sample=NULL, year_min_sim=year_min_sim,year_max_sim=year_max_sim)
#sample="monthly"

out=generation00$output
plot(unlist(out$Tx_gen), type="l")

serial_test(generation00$var)
normality_test(generation00$var)

Tx_gen= unlist(out$Tx_gen)
Tn_gen= unlist(out$Tn_gen)
#DOESN"T ACCOUNT FOR SHIFT IN MAX OVER TIME

#-----------------------
#For each year calibrate with 5 year before/after running average

for(yr in 1953:2011){
min.yrs= yr-1953
max.yrs= 2011-yr

year_min <- yr-5
year_max <- yr+5

if(min.yrs<5){
year_min <- yr-min.yrs
year_max <- yr+(10-min.yrs)
}

if(max.yrs<5){
year_min <- yr-(10-max.yrs)
year_max <- yr+max.yrs
}

year_min_sim <- yr
year_max_sim <- yr

gen <-ComprehensiveTemperatureGenerator(station=vstation,
 Tx_all=TEMPERATURE_MAX,Tn_all=TEMPERATURE_MIN,year_min=year_min,year_max=year_max,
 p=p,n_GPCA_iteration=n_GPCA_iter,n_GPCA_iteration_residuals=n_GPCA_iteration_residuals,
 sample=NULL, year_min_sim=year_min_sim,year_max_sim=year_max_sim)

out=gen$output
c1.clim[which(c1.clim$Year==yr),"Max.gen"]= unlist(out$Tx_gen)
c1.clim[which(c1.clim$Year==yr),"Min.gen"]= unlist(out$Tn_gen)
}

#TEST THAT ACCOUNTS FOR TRENDS OVER TIME
#WORKS!
c1.sum<- subset(c1.clim, c1.clim$Month<9 & c1.clim$Month>6)
plot(c1.sum$Max, c1.sum$Max.gen)
plot(c1.sum$Min, c1.sum$Min.gen)

clim.agg= aggregate(c1.sum[,5:10], by=list(c1.sum$Year), FUN=mean, na.action = na.omit)
names(clim.agg)[1]="Year"
plot(clim.agg$Year, clim.agg$Max, type="l")
points(clim.agg$Year, clim.agg$Max.gen, type="l", col="red")
plot(clim.agg$Year, clim.agg$Min, type="l")
points(clim.agg$Year, clim.agg$Min.gen, type="l", col="red")

#****************
#Simulate shift in mean

#mean_climate_Tn_sim: monthly averaged daily minimum temperatures for the simulated scenario and used by the random generator . Default is mean_climate_Tn

#Set up monthly means across years
mean_Tn=getMonthlyMean(TEMPERATURE_MIN, year_min = 1953, year_max = 2011,
  station = vstation, no_date = FALSE, origin = "1953-1-1", yearly = TRUE)
mean_Tx=getMonthlyMean(TEMPERATURE_MAX, year_min = 1953, year_max = 2011,
  station = vstation, no_date = FALSE, origin = "1953-1-1", yearly = TRUE)

#----------------------
#calculate as 10yr running averages by month

clim=c1.clim
Yr=1970
J=190
Yr.J= c(Yr, J)

clim.runmean= function(Yr.M, clim){
inds= which(clim$Julian %in% (Yr.J[2]-5):(Yr.J[2]+5) & clim$Year %in% (Yr.J[1]-5):(Yr.J[1]+5))
return( c( mean(clim[inds,"Min"], na.rm=TRUE),mean(clim[inds,"Max"], na.rm=TRUE)) )
}

Yr.J.mat= clim[,c("Year", "Julian")]
clim.out= t(apply(Yr.J.mat, MARGIN=1, FUN=clim.runmean, clim=clim))
clim$Min.ave=clim.out[,1]
clim$Max.ave=clim.out[,2]

#-----------------------

gen_shiftmean <-ComprehensiveTemperatureGenerator(station=vstation,
 Tx_all=TEMPERATURE_MAX,Tn_all=TEMPERATURE_MIN,year_min=year_min,year_max=year_max,
 p=p,n_GPCA_iteration=n_GPCA_iter,n_GPCA_iteration_residuals=n_GPCA_iteration_residuals, 
year_min_sim=year_min_sim,year_max_sim=year_max_sim, 
 mean_climate_Tn_sim=mean_Tn, mean_climate_Tx_sim=mean_Tx)

mean_Tn$'2008'

###################################3
#Put data into format for analysis

clim.genmean=c1.clim
clim.genmean$Max= clim$Max.Gen
clim.genmean$Min= clim$Min.Gen

plot(c1.clim$Max, clim.genmean$Max)
abline(a=0, b=1, col="red")

#subset to july and august
a1.clim<- subset(a1.clim, a1.clim$Month<9 & a1.clim$Month>6)
b1.clim<- subset(b1.clim, b1.clim$Month<9 & b1.clim$Month>6)
c1.clim<- subset(c1.clim, c1.clim$Month<9 & c1.clim$Month>6)
d1.clim<- subset(d1.clim, d1.clim$Month<9 & d1.clim$Month>6)
clim.genmean<- subset(clim.genmean, clim.genmean$Month<9 & clim.genmean$Month>6)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
#Set min and max temps
#Tmin=30;  # Kingsolver 1983
#Tmax=40;

#Estimate air pressure in kPa #http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
airpressure_elev<- function(h){  #H is height in meters
p= 101325* (1 - 2.25577*10^(-5)*h)^5.25588       
p= p/1000 #convert to kPa
return(p)
}

airpressures= sapply(sites$Elev, FUN=airpressure_elev) 

#Wind velocity at height z
# V_r is wind velocity at reference height

V_z<- function(V_r, z_0, z_r, z){ V_r*log((z+z_0)/z_0+1)/log((z_r+z_0)/z_0+1)}
#from Porter et al 1973

#SEE SurfaceRoughness_12May2014.R for C1 surface roughness estimate
#z=0.02 #m

#-------------------------------------------------
#PARAMETERS

#Demographic parameters
#Kingsolver 1983 Ecology 64
OviRate=0.73; # Ovipositing rates: 0.73 eggs/min (Stanton 1980) 
MaxEggs=700; # Max egg production: 700 eggs (Tabashnik 1980)
PropFlight= 0.5; # Females spend 50# of available activity time for oviposition-related
# Watt 1979 Oecologia
SurvDaily=0.6; # Daily loss rate for Colias p. eriphyle at Crested Butte, female values
#Hayes 1981, Data for Colias alexandra at RMBL
SurvMat=0.014; #1.4# survival to maturity

#pick stations and years
years= 1953:2011
days<-c(31,28,31,30,31,30,31,31,30,31,30,31) #days in months

#read/make species data
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Butterflies\\Evolution\\Data\\")
#setwd("F:\\Work\\Butterflies\\Evolution\\Data\\")
#SpecDat<-read.csv("SpecDat_23July2013.csv");
solar.abs= seq(0.4,0.7,0.01) #seq(0.4,0.7,0.001)
SpecDat= as.data.frame(solar.abs)
#SpecDat$Species="Cmeadii"
SpecDat$d=0.36
SpecDat$fur.thickness=1.46 
SpecDat= as.matrix(SpecDat)
# Solar absorptivity, proportion
# D- Thoractic diameter, cm
# Fur thickness, mm
#Fur thickness for C. eriphyle (Montrose) is 0.82mm

#Make array to store lambdas
#Lambda= survival*fecundity
#Matrix of lambdas
#dims: stats, year, Lambda 
Lambda_20cm<-array(NA, dim=c(length(years),length(stats),nrow(SpecDat),4)) #Last dimension is Lambda, FAT,Egg Viability
dimnames(Lambda_20cm)[[1]]<-years
dimnames(Lambda_20cm)[[2]]<-stats
dimnames(Lambda_20cm)[[3]]=solar.abs

Lambda_50cm<-array(NA, dim=c(length(years),length(stats),nrow(SpecDat),4)) #Last dimension is Lambda, FAT,Egg Viability
dimnames(Lambda_50cm)[[1]]<-years
dimnames(Lambda_50cm)[[2]]<-stats
dimnames(Lambda_50cm)[[3]]=solar.abs
#------------------------
#PERFORMANCE FUNCTIONS

#function for flight probability
#fl.ph<- function(x,a,b,c,d) a * exp(-0.5*(abs(x-b)/c)^d)  
#fl.ph<- function(x) 1 * exp(-0.5*(abs(x-35)/4.5)^5)  #eriphyle
#fl.ph<- function(x) 1 * exp(-0.5*(abs(x-32)/4.5)^5) #meadii
#PLACE HOLDER FUNCTION TO BROADEN AND INCLUDE FLIGHT AT LOWER TEMPERATURES
#fl.ph<- function(x) 1 * exp(-0.5*(abs(x-32)/6)^4)
fl.ph<- function(x) 1 * exp(-0.5*(abs(x-32.5)/4.7)^4)
fl.ph<- function(x) 1 * exp(-0.5*(abs(x-33.5)/5)^3.5)

#function for egg viability
#egg.viab<-function(x) ifelse(x<40, egg.viab<-1, egg.viab<- exp(-(x-40)))
#Set P_survival=0.868 at 45 (Kingsolver dissertation), log(0.868)=-5/t, t=35.32

egg.viab<-function(x) ifelse(x<40, egg.viab<-1, egg.viab<- exp(-(x-40)/35.32))
#plot(30:50, egg.viab(30:50))

#estimate flight time for Meadii using data from Ward 1977
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Butterflies\\Evolution\\Data\\")
#setwd("F:\\Work\\Butterflies\\Evolution\\Data\\")
bpop= read.csv("ButterflyPop_Ward1977.csv")
bpop= subset(bpop, bpop$Loc=="CumbPass")
bpop.tot= sum(bpop$Total)
bpop$prop= bpop$Total/bpop.tot

#Calculate weighted mean and sd for flight date
flightday.mean=wt.mean(bpop$Date, bpop$Total)
flightday.sd=wt.sd(bpop$Date, bpop$Total)

#---------------------------------------------
#Run Te calculations

for(yr.k in 1:length(years)){ #loop years
#for(yr.k in 1:length(years)){ #loop years
yr<-years[yr.k]
print(yr)

for(stat.k in 1:2){ #loop stats
#for(stat.k in 1:4){ #loop stats

if(stat.k==1)dat=c1.clim
if(stat.k==2)dat=clim.genmean
#if(stat.k==3)dat= clim.genvar
#if(stat.k==4)dat=d1.clim

lat<- sites$Lat[3]
lon<- sites$Lon[3]
elev<- sites$Elev[3]

#subset to July and August
dat2<-subset(dat, dat$Month>6 & dat$Month<9)

#subset to year
dat2<-subset(dat, dat$Year==yr)
dat2.temp= mean(dat2$Max)


if(nrow(dat2)>1 & !is.na(dat2.temp)){ #check data for dat exists

#for daylength calculations

#estimate daylength
Trise.set= suncalc(dat2$Julian, Lat = lat, Long = lon, UTC = FALSE)
set= Trise.set$sunset
rise= Trise.set$sunrise

# Calculate the diurnal temperature trend, Parton and Logan 1981
Tmat= cbind(dat2$Max, dat2$Min, rise, set )
J= dat2$Julian

#RUN EVERY TEN MINUTES
hr.dec= seq(1, 24, 1/6)

#arrays for storing temperature and radiation data
Thr.mat= matrix(NA, nrow(Tmat),length(hr.dec) )
Tsoil.mat= matrix(NA, nrow(Tmat),length(hr.dec) )
Rhr.mat= array(data=NA, dim=c(nrow(Tmat),length(hr.dec),3) ) #columns for direct, difuse, reflected radiation
Te.mat_20cm= array(data=NA, dim=c(3,nrow(Tmat),length(hr.dec), nrow(SpecDat)))
Te.mat_50cm= array(data=NA, dim=c(3,nrow(Tmat),length(hr.dec), nrow(SpecDat)))

#------------------------------
#Calculate hourly temp and radiation
for(hc in 1:length(hr.dec) ){ 
hr= hr.dec[hc]

Thr.mat[,hc]=apply(Tmat,FUN=Thour.mat, MARGIN=1, Hr=hr,  alpha=1.86, beta= -0.17, gamma=2.20)

#PICK TAU, atmospheric transmisivity, from distribution for k_t
#USES KERNAL ESTIMATE FIT TO HOURLY NREL SRRL DATA (loaded when sourcing solar radiation functions)

if(round(hr)>5 & round(hr)<21)  taus= rkde(fhat= kdes[[ round(hr)-5]], n=nrow(Tmat) )
if (!(round(hr)>5 & round(hr)<21)) taus= rep(0, nrow(Tmat))

#set negative values to zero, CHECK HOW TO CONSTRAIN
taus[taus<0]=0
#CONSTRAIN TAUs BETWEEN 7AM AND 3PM, ###CHECK
if(hr>6 & hr<16) taus[taus<0.2]=0.2

#calculate zenith in radians
psis= unlist(lapply(J, FUN=zenith.rad, Hr=hr, lat=lat, lon=lon))

#RADIATION FROM CAMPBELL & NORMAN
J.mat= cbind(J, psis, taus)
Rdat=apply(J.mat,FUN=calc.rad, MARGIN=1, lat=lat, lon=lon, elev=elev)
#returns direct, diffuse, reflected
#set negative radiation values to zero ###CHECK
Rdat[Rdat<0]=0

R_total= colSums(Rdat[1:2,])

#USE ERBS TO REPARTITION RADIATION (Olyphant 1984)
#Separate Total radiation into components
#kt is clearness index
# Models presented in Wong and Chow. 2001. Applied Energy 69(2001):1991-224
#Use Erbs et al model

#kd- diffuse fraction
kd= rep(NA, length(taus))

inds= which(taus<0.22) 
kd[inds]= 1-0.09*taus[inds]
inds= which(taus>0.22 & taus<0.8) 
kd[inds]= 0.9511 -0.1604*taus[inds] +4.388*taus[inds]^2 -16.638*taus[inds]^3 +12.336*taus[inds]^4
inds= which(taus>=0.8)
kd[inds]= 0.125 #Correction from 16.5 for Niwot from Olyphant 1984

#direct and diffuse #columns for direct, difuse, reflected radiation
Rhr.mat[,hc,1]<- R_total*(1-kd)
Rhr.mat[,hc,2]<- R_total*(kd)
Rhr.mat[,hc,3]<- Rdat[3,]

} #end hour loop

#----------------------
# RUN MICROCLIMATE MODEL

#calculate soil temp
# INPUT DATA
# air_temp: air temp (C)
# V_z: wind speed (m/s)
# TimeIn: time vector
# solrad: solar radiation 

# DEFINE PARAMETERS
# z.intervals=12 number depth intervals for soil t measurements 
# t.intervals: number time intervals
# dt: time step
# z: reference height (m) #COOP stations are 5ft=1.524m
# T_so: initial soil temp
# z_0: surface roughness length
# SSA=0.7 solar absorbtivity of soil surface
# k_so: soil conductivity
# water_content= 0.2 percetn water content
# air_pressure= 101 kPa at sea level
# rho_so= 1620 particle density of soil
# z_new= new height (m)

#data that will be used as input for the function.
#collapse matrix
temperature_vector<- unmatrix(Thr.mat,byrow=T)
time_vector<- rep(hr.dec, nrow(Thr.mat) )
wind_speed_vector<-rep(0.4, length(temperature_vector)) ####GET NIWOT WIND DATA
solrad_vector<- unmatrix(Rhr.mat[,,1],byrow=T)

#add day to begining to account for transients
temperature_vector= c(temperature_vector[1:ncol(Thr.mat)], temperature_vector)
time_vector= c(time_vector[1:ncol(Thr.mat)], time_vector)
wind_speed_vector= c(wind_speed_vector[1:ncol(Thr.mat)], wind_speed_vector)
solrad_vector= c(solrad_vector[1:ncol(Thr.mat)], solrad_vector)

#calculate soil temperatures
z_0_1=0.02 #m
z_new= 0.2 #m

Tsoil.mat= soil_temp_function_noint(z.intervals=12, z=1.524, air_temp=temperature_vector, V_z=wind_speed_vector, T_so=20, z_0=z_0_1, SSA=.7, k_so=0.293, TimeIn=time_vector, solrad=solrad_vector, water_content=.2, air_pressure=airpressures[3], rho_so=1620, z_new=z_new)
Tsoil= Tsoil.mat[,1]
Tsoil_sh= Tsoil.mat[,2]
#FAIRLY SLOW
#CURRENTLY DOES NOT ACCOUNT FOR SOIL SHADING

#TEST 
#z.intervals=12; z=1.524; air_temp=temperature_vector; V_z=wind_speed_vector; T_so=20;
#z_0=z_0_1; SSA=.7; k_so=0.293; TimeIn=time_vector; solrad=solrad_vector; water_content=.2;
#air_pressure=airpressures[stat.k]; rho_so=1620; z_new=z_new

#calculate butterfly temp at plant height
Ta_20cm= air_temp_at_height_z(z_0=z_0_1, z_r=1.524, z=0.2, T_r=temperature_vector, T_s=Tsoil)
Ta_50cm= air_temp_at_height_z(z_0=z_0_1, z_r=1.524, z=0.5, T_r=temperature_vector, T_s=Tsoil)
#Shade
Ta_20cm_sh= air_temp_at_height_z(z_0=z_0_1, z_r=1.524, z=0.2, T_r=temperature_vector, T_s=Tsoil_sh)
Ta_50cm_sh= air_temp_at_height_z(z_0=z_0_1, z_r=1.524, z=0.5, T_r=temperature_vector, T_s=Tsoil_sh)

#Turn back into matrix
Tsoil.mat= matrix(Tsoil, nrow=nrow(Thr.mat)+1, ncol=ncol(Thr.mat), byrow=T)
Ta_20cm.mat= matrix(Ta_20cm, nrow=nrow(Thr.mat)+1, ncol=ncol(Thr.mat), byrow=T)
Ta_50cm.mat= matrix(Ta_50cm, nrow=nrow(Thr.mat)+1, ncol=ncol(Thr.mat), byrow=T)
#Get rid of transient first row
Tsoil.mat= Tsoil.mat[-1,]
Ta_20cm.mat= Ta_20cm.mat[-1,]
Ta_50cm.mat= Ta_50cm.mat[-1,]
#SHADE
Tsoil.mat_sh= matrix(Tsoil_sh, nrow=nrow(Thr.mat)+1, ncol=ncol(Thr.mat), byrow=T)
Ta_20cm.mat_sh= matrix(Ta_20cm_sh, nrow=nrow(Thr.mat)+1, ncol=ncol(Thr.mat), byrow=T)
Ta_50cm.mat_sh= matrix(Ta_50cm_sh, nrow=nrow(Thr.mat)+1, ncol=ncol(Thr.mat), byrow=T)
#Get rid of transient first row
Tsoil.mat_sh= Tsoil.mat_sh[-1,]
Ta_20cm.mat_sh= Ta_20cm.mat_sh[-1,]
Ta_50cm.mat_sh= Ta_50cm.mat_sh[-1,]

#-----------------------
#loop hours
for(hc in 31:115){ 
hr= hr.dec[hc]

#Calculate zenith in degrees
psi=zenith(J, lat, lon, hr)
psi[psi>=80]=80 #set zenith position below horizon to psi=80degrees

wind_20cm= V_z(V_r= 0.4, z_0=0.02, z_r=1.524, z=0.2)
wind_20cm= rep(wind_20cm, length(psi))
wind_50cm= V_z(V_r= 0.4, z_0=0.02, z_r=1.524, z=0.5)
wind_50cm= rep(wind_50cm, length(psi))

#Vary SPECIES
for(spec.k in 1:nrow(SpecDat) ){ #loop species

#Fur thickness
delta= SpecDat[spec.k,"fur.thickness"]
#if(stat.k<3) delta= 0.85 #Value for Montrose

#Convert to Te
#20cm plant height
Tall= cbind(Ta_20cm.mat[,hc], Ta_20cm.mat_sh[,hc], Tsoil.mat[,hc], Tsoil.mat_sh[,hc], wind_20cm, Rhr.mat[,hc,1],Rhr.mat[,hc,2], psi)
Te.mat_20cm[1,,hc,spec.k]<-apply(Tall,MARGIN=1, FUN=biophys.var_sh, D=SpecDat[spec.k,"d"], delta= delta, alpha=SpecDat[spec.k,"solar.abs"])
#Set negative Te estimates to Ta ##only a few, but CHECK
inds= which(Te.mat_20cm[1,,hc,spec.k]<0)
Te.mat_20cm[1,inds,hc,spec.k]= Ta_20cm.mat[inds,hc]

#--------------
#TEST
#Temat=Tall[1,]; D=SpecDat[spec.k,"d"]; delta= SpecDat[spec.k,"fur.thickness"]; alpha=SpecDat[spec.k,"solar.abs"]

#---------------

#50cm plant height
Tall= cbind(Ta_50cm.mat[,hc], Ta_50cm.mat_sh[,hc], Tsoil.mat[,hc], Tsoil.mat_sh[,hc], wind_50cm, Rhr.mat[,hc,1],Rhr.mat[,hc,2], psi)
Te.mat_50cm[1,,hc,spec.k]<-apply(Tall,MARGIN=1, FUN=biophys.var_sh, D=SpecDat[spec.k,"d"], delta= SpecDat[spec.k,"fur.thickness"], alpha=SpecDat[spec.k,"solar.abs"])
#Set negative Te estimates to Ta ##only a few, but CHECK
inds= which(Te.mat_50cm[1,,hc,spec.k]<0)
Te.mat_50cm[1,inds,hc,spec.k]= Ta_50cm.mat[inds,hc]

##NO SHADE OPTION
#Tall= cbind(Ta_20cm.mat[,hc], Tsoil.mat[,hc], wind_20cm, Rhr.mat[,hc,1],Rhr.mat[,hc,2], psi)
#Te.mat_20cmNOSH<-apply(Tall,MARGIN=1, FUN=biophys.var, D=SpecDat[spec.k,"d"], delta= SpecDat[spec.k,"fur.thickness"], alpha=SpecDat[spec.k,"solar.abs"])

#50cm plant height
#Tall= cbind(Ta_50cm.mat[,hc], Tsoil.mat[,hc], wind_50cm, Rhr.mat[,hc,1],Rhr.mat[,hc,2], psi)
#Te.mat_50cmNOSH<-apply(Tall,MARGIN=1, FUN=biophys.var, D=SpecDat[spec.k,"d"], delta= SpecDat[spec.k,"fur.thickness"], alpha=SpecDat[spec.k,"solar.abs"])

#plot(Ta_20cm.mat[,hc], Te.mat_20cm[1,,hc,spec.k]) #, ylim=range(15,35), xlim=range(15,35) )
#points(Ta_20cm.mat[,hc], Te.mat_20cmNOSH, col="red") 
#abline(a=0, b=1)

#--------------------------------

#plot(Ta_20cm.mat[,hc], Te.mat_20cm[1,,hc,spec.k]) #, ylim=range(15,35), xlim=range(15,35) )
#points(Ta_50cm.mat[,hc], Te.mat_50cm[1,,hc,spec.k], col="red") 
#abline(a=0, b=1)
#biophys.var(Temat=Tall[1,], D=SpecDat[spec.k,"d"], delta= SpecDat[spec.k,"fur.thickness"], alpha=SpecDat[spec.k,"solar.abs"])

#DEMOGRAPHY
#Flight probability
Te.mat_20cm[2,,hc,spec.k]<-fl.ph(Te.mat_20cm[1,,hc,spec.k])
Te.mat_50cm[2,,hc,spec.k]<-fl.ph(Te.mat_50cm[1,,hc,spec.k])

#Egg viability
Te.mat_20cm[3,,hc,spec.k]<-egg.viab(Te.mat_20cm[1,,hc,spec.k])
Te.mat_50cm[3,,hc,spec.k]<-egg.viab(Te.mat_50cm[1,,hc,spec.k])

} #end species loop
} #end hour loop
#----------------------------------------
EV1=matrix(NA,2, nrow(SpecDat))

for(spec.k in 1:nrow(SpecDat) ){ #loop species

for(h in 1:2){ #loop species

if(h==1) Te.mat= Te.mat_20cm
if(h==2) Te.mat= Te.mat_20cm

#average over hours
FAT= rowSums(Te.mat[2,,,spec.k], na.rm=TRUE)/6 #Flight time, divide by 6 to return to hours
EggViab= geo_mean(Te.mat[3,,,spec.k]) #Egg viability GEOMETRIC MEAN ACROSS HOURS

#AVERAGE FAT OVER DAYS IN SEASON
FAT= mean(FAT, na.rm=TRUE)
FAT= mean(FAT, na.rm=TRUE)

#CALCULATE EGG VIABILITY OVER 5 DAY PERIOD (GEOMETRIC MEAN ACROSS HOURS)
#sample flight day from truncated normal distribution
Nind=500 #changed from 100
flightday= round(rtnorm(Nind, mean = flightday.mean, sd = flightday.sd, lower=min(J)+2, upper=max(J)-2))

Vmat= Te.mat[3,,,spec.k]
f.ind= match(flightday, J)
#calculate geometric mean of egg viability within flight period
ev.ind=sapply(1:length(f.ind), function(x)  geo_mean(c(Vmat[(f.ind[x]-2):(f.ind[x]+2),])) )
#calculate arithmetic mean of egg viability across flight period over individuals
EggViab=mean(ev.ind, na.rm=TRUE)

Eggs= FAT*60*PropFlight*OviRate*EggViab #account for Egg viability
Eggs_noViab= FAT*60*PropFlight*OviRate
EV1[1,spec.k]=FAT
EV1[2,spec.k]=EggViab

if(!is.nan(Eggs)){
MaxDay=5
Lambda1=0
for(day in 1:MaxDay){
Eggs1= min(Eggs, MaxEggs-Eggs*(day-1))  ###LIMIT MAX NUMBER EGGS
if(Eggs1<0) Eggs1=0
Lambda1= Lambda1+ SurvMat * SurvDaily^day *Eggs1;                        
}#end loop days

if(h==1){
Lambda_20cm[yr.k, stat.k,spec.k,1]= Lambda1
Lambda_20cm[yr.k, stat.k,spec.k,2]= FAT #FAT
Lambda_20cm[yr.k, stat.k,spec.k,3]= EggViab #egg viability
}
if(h==2){
Lambda_50cm[yr.k, stat.k,spec.k,1]= Lambda1
Lambda_50cm[yr.k, stat.k,spec.k,2]= FAT #FAT
Lambda_50cm[yr.k, stat.k,spec.k,3]= EggViab #egg viability
}

}#Check Eggs

#WITHOUT EGG VIAB
if(!is.nan(Eggs)){
MaxDay=5
Lambda1=0
for(day in 1:MaxDay){
Eggs1= min(Eggs_noViab, MaxEggs-Eggs_noViab*(day-1))  ###LIMIT MAX NUMBER EGGS
if(Eggs1<0) Eggs1=0
Lambda1= Lambda1+ SurvMat * SurvDaily^day *Eggs1;                        
}#end loop days

if(h==1) Lambda_20cm[yr.k, stat.k,spec.k,4]= Lambda1
if(h==2) Lambda_50cm[yr.k, stat.k,spec.k,4]= Lambda1
}#Check Eggs

} #end loop heights

} #end loop species
 
} #end check data

} #end loop stats
} #end loop years

#------------------------------------
#PUT DATA IN LONG FORMAT AND WRITE OUT

for(h in 1:2){ #loop heights

for(stat.k in 1:length(stats)){ #loop stats

if(h==1) lambda= unmatrix(Lambda_20cm[,stat.k,,1], byrow=FALSE)
if(h==2) lambda= unmatrix(Lambda_50cm[,stat.k,,1], byrow=FALSE)

dat.lamb=as.data.frame(lambda)
n.lamb=strsplit(names(lambda), ":")
n.lamb=do.call(rbind, n.lamb) 

dat.lamb$Year= n.lamb[,1]
dat.lamb$Abs= n.lamb[,2]
dat.lamb$Station= stats[stat.k]
dat.lamb$Elev= sites$Elev[stat.k]

#FAT
if(h==1) dat.lamb$FAT= unmatrix(Lambda_20cm[,stat.k,,2], byrow=FALSE)
if(h==2) dat.lamb$FAT= unmatrix(Lambda_50cm[,stat.k,,2], byrow=FALSE)

#EGG VIAB
if(h==1) dat.lamb$EggViab= unmatrix(Lambda_20cm[,stat.k,,3], byrow=FALSE)
if(h==2) dat.lamb$EggViab= unmatrix(Lambda_50cm[,stat.k,,3], byrow=FALSE)

#LAMBDA NO EGG VIAB
if(h==1) dat.lamb$lambda_noEV= unmatrix(Lambda_20cm[,stat.k,,4], byrow=FALSE)
if(h==2) dat.lamb$lambda_noEV= unmatrix(Lambda_50cm[,stat.k,,4], byrow=FALSE)

if(stat.k==1) dat.all=dat.lamb
if(stat.k>1) dat.all= rbind(dat.all, dat.lamb)
} #end loop stats

if(h==1) dat.all_20cm= dat.all
if(h==2) dat.all_50cm= dat.all

} #end loop heights

#------------------------
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Butterflies\\Evolution\\OUT\\")
#setwd("F:\\Work\\Butterflies\\Evolution\\OUT\\")
write.csv(dat.all_20cm, "LambdaLong_20cm_shade_Niwot30June_GenDataResid.csv", row.names=FALSE)
write.csv(dat.all_50cm, "LambdaLong_50cm_shade_Niwot30June_GenDataResid.csv", row.names=FALSE)

#write.csv(Lambda[,2,,1], "LambdaMat.csv")

par(mfrow=c(1,2))
c1dat=subset(dat.all_20cm, dat.all$Station== "B1")
plot(c1dat$Abs, c1dat$lambda)
plot(c1dat$Abs, c1dat$lambda_noEV, col="red")
#===================================================
#TO DO:
#WINDSPEED: Adjusted fixed value based on empirical data, 0.4m/s 

#Windspeed from dataloggers at 0.75m

setwd("F:\\Work\\Grasshoppers\\data\\PaceDatalogging\\Weather2011csv\\")
dat= read.csv("DatAll.csv")

agg= aggregate(dat$Wind, by = list(dat$site), FUN = "mean", na.rm=TRUE) 

mean(dat$Wind, na.rm=TRUE)
