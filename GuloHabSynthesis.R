
## Wolverine Habitat Predictions & Synthesis ##

# Gabe Schepens
# Yellowstone to Yukon Conservation Initiative Internship

# Manuscript 2022

#### SETUP ####

library(rgdal)
library(gdalUtils)
library(raster)
library(dplyr)

#setwd("~(...)/EnvLayers")

# area of interest
area <- shapefile("StudyArea2/StudyArea_2buffer.shp")
border <- shapefile("StudyArea2/StudyArea_20220213.shp")
# create empty raster template to match
RAST <- raster(ext = extent(area), res = 500, crs = crs(area))

# four model areas
ModAreas <- shapefile("ModelAreas/ModAreas.shp")

# environmental layers
SHRUB <- raster("Landcover2015/Shrub10kMean.tif")
CONIFER <- raster("Landcover2015/Conifer10k.tif")
MIXED <- raster("Landcover2015/Mixedwood10k.tif")

SNOW <- raster("PSSC/Snow10kMean.tif")

FSR <- raster("Roads/ResourceRdDens10km.tif")
INDUST <- raster("Roads/IndustrialDensity10k.tif")

PROTEC <- raster("ProtectedAreas/Protected10km.tif")
MARMOT <- raster("Moraines/Moraines10km.tif")

#### PREDICTIONS ####

#### Barrueto et al 2020 Model Prediction ####

vars <- stack(SNOW, SHRUB)
coefs <- c(0.3, 0.25)

# scale variables to Mean and SD of Barrueto study area
for(i in 1){
  m <- as.numeric(extract(vars[[i]], ModAreas[2,], fun = mean))
  stand <- as.numeric(extract(vars[[i]], ModAreas[2,], fun = sd ))
  vars[[i]] <- (vars[[i]] - m )/stand
}

# predict 
BarPred <- exp(sum(vars*coefs))

# CI maps
BarCIlow <- exp(sum(vars*c(0.05, -0.10)))
BarCIup <- exp(sum(vars*c(0.56, 0.51)))

#BarPred_clip <- mask(BarPred, area)
#writeRaster(BarPred, "BarruetoPred.tif")


#### Kortello et al 2019 Model Prediction ####

vars <- stack(SNOW, FSR, PROTEC, MARMOT)
coefs <- c(1.06, -1.23, 0.87, 1.28)
SE <- c(0.23, 0.26, 0.19, 0.42)

# log as in paper
vars[2:4,] <- log(vars[2:4,]+ 0.01)

# scale to kortello study area
for(i in 1:4){
  m <- as.numeric(extract(vars[[i]], ModAreas[1,], fun = mean, na.rm = T))
  stand <- as.numeric(extract(vars[[i]], ModAreas[1,], fun = sd, na.rm = T ))
  vars[[i]] <- (vars[[i]] - m )/stand
}

# predict
KortPred <- sum(vars*coefs)

# CI maps
KortCIlow <- sum(vars*(coefs - 1.98*SE))
KortCIup <- sum(vars*(coefs + 1.98*SE))
KortSE <- sum(vars*SE)


#### Heim et al 2017 Model Prediction ####

vars <- stack(SHRUB, MIXED, CONIFER, INDUST )
coefs <-  c(0.338, -0.929, 0.605, -1.243 )

# scale to heim study area
for(i in 1:4){
  m <- as.numeric(extract(vars[[i]], ModAreas[3,], fun = mean, na.rm = T))
  stand <- as.numeric(extract(vars[[i]], ModAreas[3,], fun = sd, na.rm = T ))
  vars[[i]] <- (vars[[i]] - m )/stand
}

# predict
HeimPred <- exp(sum(vars*coefs))


#### Mowat et al 2018 from open data ####

MowPred <- raster("Models/Mowat2018_Gulo_Density_Surface.tif")
MowPred <- resample(MowPred, RAST)


#### SYNTHESIS ####

# Stack four model predictions
PREDS <- stack(KortPred, BarPred, HeimPred, MowPred)
plot(PREDS)

## Habitat suitability based on equal-area value

# PERCENTILE transformation w breaks from original STUDY AREAS

PERCEN <- PREDS

for(i in 1:4){
  # extract vals from study area
  vals <- as.numeric(unlist(extract(PREDS[[i]], ModAreas[i,], na.rm = T)))
  # define 20 quantile breaks
  breaks <- quantile(vals, seq(0,1,1/20), na.rm = TRUE)
  # extrapolate inf values to each end of the "breaks" (otherwise NA on out-of-range)
  breaks[1] <- -Inf
  breaks[21] <- Inf
  # give new vals 1:20 based on quantile breaks
  PERCEN[[i]] <- cut(PREDS[[i]], breaks = breaks)
  PERCEN[[i]] <- PERCEN[[i]]/20
}

library(viridis)
plot(PERCEN, col = magma(20), main = paste(ModAreas$Name))

# visualize by original study area
par(mfrow = c(2, 2))
for(i in 1:4){
  plot(mask((PERCEN[[i]]),ModAreas[i,]) , 
       main = paste(ModAreas[i,]$Name), 
       col = magma(20))
  lines(border, col = "grey")
}
par(mfrow = c(1,1))

## Distance-weighted mean

# distance to study area rasters
DIST <- PREDS
for(i in 1:4){
  DIST[[i]] <- ModAreas[i,] %>% rasterize(RAST) %>% distance()/1000 +1
}

plot(DIST)

# weighted mean by inverse distance

ALL <- mask(weighted.mean(PERCEN, w = 1/DIST, na.rm = T), border)

plot(ALL)


#### ERROR Vis ####

#### Residuals ####
## visualize difference of synth to each prediction in its study area

for(i in 1:4){
  plot(mask((PERCEN[[i]]-ALL),ModAreas[i,]) , 
       main = paste(ModAreas[i,]$Name))
}

# hist of vals from above
par(mfrow = c(2,2))
names <- c("Kortello et al. 2019", "Barrueto et al. 2020",
           "Heim et al. 2017", "Mowat et al. 2020")
for(i in 1:4){
  hist(mask((PERCEN[[i]]-ALL),ModAreas[i,]) , 
       main = paste(names[i]),
       breaks = c(seq(-0.8, 0.8, 0.05)),
       xlim = c(-0.6, 0.6), 
       freq = F, 
       xlab = "Residual Values in Model Area")
}
par(mfrow = c(1,1))


#### Coefficient of Variation ####

library(Weighted.Desc.Stat)

# Standard deviation raster by inverse distance weight
SD <- w.sd(PERCEN, 1/DIST)

# few NAs produced from zero artifacts, change to zero manually
SD[is.na(SD)] <- 0

# coefficient of variation raster (SD/ weighted mean)
COVAR <- (SD/ALL)

hist(COVAR,
     freq = F, breaks= 20, col = viridis(20), 
     xlab = "Weighted Coefficient of Variation")
#add line at 95th percentile
abline(v = quantile(values(COVAR), 0.95, na.rm = T), lty = "dashed")

#looks pretty good. most COVAR < 1

plot(COVAR, col = viridis(20), 
     main = "Weighted Coefficient of Variation")


#### Model Agreement on "Good" Habitat ####

BEST <- PERCEN > 0.5

BESTag <- mask(sum(BEST), border)

plot(BESTag, col = magma(5), 
     main = "Model Agreeance of Areas with Habitat Suitability > 0.5")

#writeRaster(BESTag, "ModelAgreeance05.tif")

freq <- freq(BESTag)
#percent of area
freq/sum(freq[c(1:5),2])


#### Compare to Independant Data ####

##- NOT SHARED -##

## using independent denning data from Kananaskis-Ghost, Alberta Parks

DEN <- raster("KG_Denning.tif")
plot(DEN)
lines(border)

#subset habitat vals where DEN = T
denHAB <- ALL[!is.na(DEN)]

summary(denHAB)

hist(denHAB, breaks = seq(0, 1, 0.05), xlim = c(0,1), col = magma(20), freq = F, 
     xlab = "Synthesized Habitat Value",
     ylab = "Wolverine den areas\nKananaskis-Ghost Region")

## using camera data from Banff Kootenay Yoho NP, Parks Canada

bkp <- shapefile("CamData/BanffYohoKootGULOcamlocs.shp")

plot(ALL)
points(bkp)

bkpvals <- raster::extract(ALL, bkp, buffer = 100, fun = mean)

hist(bkpvals, xlim = c(0,1), breaks = seq(0, 1, 0.05),  col = magma(20), 
     xlab = "Synthesized Habitat Value", 
     ylab = "Wolverine camera detections \nBanff, Kootenay, Yoho NP")


#### POST-HOC ENVIRO AND CLIMATE ####

#### Variable correlation to the synthesized habitat raster ####

AllVARS <- stack( SHRUB, CONIFER, MIXED, SNOW, FSR, INDUST, PROTEC, MARMOT) 

#look at variable correlation to the final habitat raster

#pearson cor
cor(values(AllVARS),
    values(ALL),
    use = "na.or.complete")


### Bootstrap cor for CI ###

vals <- na.omit(data.frame(values(AllVARS), 'Hab' = values(ALL)))

# empty matching df
bootcor <- vals[1, c(1:8)]

for(i in 1:1000){
  
  # random row numbers
  row <- sample(c(1:nrow(vals)), replace = T, size = 1000)
  
  # cor on above rows, compare col 1 thru 8 to 9 (HAB)
  for(c in 1:8){
    bootcor[row,c] <- cor(vals[row, c], vals[row, 9])
  }
}

# 95th percentile values
bootCI <- data.frame("CIlow", "CIup")
for(i in 1:8){
  bootCI[i,] <- quantile(bootcor[,i], c(0.025, 0.975), na.rm = T)
}


## plot vars

plot(mask(AllVARS, area), col = viridis(20))


#### Snow Prediction ####

# Snow water equivalent model raster from https://data.pacificclimate.org/ 

SWE20 <- raster("SWE/SnowWE2020.tif")

SWE20 <- SWE20 %>% projectRaster(crs = crs(RAST)) %>% resample(RAST)

SWE20[SWE20<0] <- NA


## Class to make comparable with Modis snow cover used in modelling ##

# 10km smoothing buffer
TEN <- focalWeight(RAST, 10000,'circle')
SWEbuff <- SWE20 %>% focal(TEN)

# correlation?
cor(values(SNOW),
    values(SWEbuff),
    use = "na.or.complete")

# extract to dataframe
snow <- na.omit(data.frame("SNOW" = values(SNOW), "SWEbuff" = values(SWEbuff)))

logist <- glm(SNOW ~ SWEbuff, data = snow, family = "binomial")

summary(logist)

mod.df <- data.frame("SWEbuff" = seq(0, 1300, by = 10))

mod.df$pred <- predict(logist, newdata = mod.df, type = "response")

plot(pred ~ SWEbuff, mod.df, type = "l", 
     ylab = 'MODIS', xlab = 'Predicted Snow Water Equivalent', 
     ylim = c(0, 1))
points(snow$SNOW ~ snow$SWEbuff)
lines(pred ~ SWEbuff, mod.df, col = "red")

# Classing point = where prediction is at 0.5 (inflexion point)

mod.df[which(mod.df$pred >= 0.5),]
# SWE = 390 @ pred = 0.505

# add to plot
points(mod.df[mod.df$SWEbuff == 390,], col = "red")

## Classing snow water equivalent to binomial @ 390 ##

# check using classing matrix
library(caret)

# create binomial snow data with the cutoffs 
snowbi <- na.omit(data.frame("SNOW" = values(SNOW>0.5), "SWEbuff" = values(SWEbuff>390)))

confusionMatrix(as.factor(snowbi$SNOW), as.factor(snowbi$SWEbuff))
# classing accuracy 80.7 %



#### Appendix table: parks with 2080 snow ####

protected <- shapefile("ProtectedAreas/Protected.shp")

# use original resolution/proj SWE predictions
SWE4 <- raster("SWE/SWE_RCP4.5_May2080_EnsembleMean.tif")
SWE8 <- raster("SWE/SWE_RCP8.5_May2080_EnsembleMean.tif")

protected2 <- spTransform(protected, crs(SWE4))

plot(SWE4)
lines(protected2)

#make binary SWE rasters with 390 as binary cutoff
SWE4bi <- (SWE4 > 390)
SWE8bi <- (SWE8 > 390)

parksnow <- data.frame(park = protected$NAME)

parksnow$swe20 <- raster::extract(SWE20bi, protected, fun = mean, na.rm = T)

parksnow$swe80_4 <- raster::extract(SWE4bi, protected2, fun = mean, na.rm = T)
parksnow$swe80_8 <- raster::extract(SWE8bi, protected2, fun = mean, na.rm = T)

mean(parksnow$swe20 - parksnow$swe80, na.rm = T)

parksnow$hab <- raster::extract(ALL, protected, fun = mean, na.rm = T)

#remove parks outside of area from df

parksnow <- parksnow[!is.nan(parksnow$hab),]

#write.csv(parksnow, "ParksHabSnow2080.csv")

## Snow cover changes ##

# reproject for comparison
SWE8bi2 <- SWE8bi %>% projectRaster(crs = crs(RAST), method = 'ngb') %>% 
  resample(RAST, method = 'ngb') %>% mask(border)
SWE4bi2 <- SWE4bi %>% projectRaster(crs = crs(RAST),method = 'ngb') %>% 
  resample(RAST, method = 'ngb') %>% mask(border)
SWE20bi <- mask(SWE20bi, border)

# total snow cover values btwn 2020 and 2080

freq(SWE20bi)
## 2020 130080
freq(SWE4bi2)
## 2080 RCP4.5 118384
## 91.0%
freq(SWE8bi2)
## 2080 RCP8.5 86015
## 66.1 %



#### OTHER STATS computed for paper ####

### habitat value in parks

pro <- ALL

pro <- mask(pro, protected)

mean(values(pro), na.rm = T)
sd(values(pro), na.rm = T)

# out of parks

nopro <- ALL
nopro[!is.na(pro)]<- NA
plot(nopro)

mean(values(nopro), na.rm = T)
sd(values(nopro), na.rm = T)



