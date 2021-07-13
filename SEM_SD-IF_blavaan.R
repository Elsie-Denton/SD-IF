####
#### okay, not finding exactly the tutorial I need so I'm just going
### to start trying code and see if I can get this to work
#### goal of this analysis is to get 
#### a probit link function working with SD-IF data
### We'll start with just tillers as a predictor as that seems simplest
#### I think I'll just run tillers and a standard numeric predictor variable
#### but the text is a little light on if that is correct or not

### Elsie Denton
### 27 Dec 2019

### modified 14 Jan 2020 to use blavaan
### the count data (which is not supported by regular lavaan is causing problems with fit)


### Clear workspace

rm(list=ls())

### set working directory

setwd("C:/Users/elsie.denton/Documents/Projects/SD-IF/Data/2018/Second_Year_Survival/Final")
dir()

### load data

df<-read.csv("SD-IF_CompiledDataToShare.csv", sep=",", header=T)
head(df)


##### information on this dataset
##### Project: SD-IF, a two year studying looking at how seedlings respond to controled defoliation
##### Block: 4 blocks were used to control spatial varition in replication in this study
##### Plot: A unique identifier for each plot in this study
##### Unique: Multiple seedlings per plot were monitored in this study, this is the unique identifier for each one. 
##### Species: four letter codes for each species in study. POSE: Poa secunda, AGCR: Agropyron crisatum, PSSP: Pseudoroegneria spicata
##### Trt: treatment applied to seedlings, determined by random assignments to plots. 
#####      0 control, no manipulation
#####      30x1 thirty percent of leaf length removed once at 2 leaf stage, 
#####      30x2 thirty percent of leaf length removed once at 2 leaf stage, and again one month later
#####      70x1 seventy percent of leaf length removed once at 2 leaf stage, 
#####      70x2 seventy percent of leaf length removed once at 2 leaf stage, and again one month later
##### Fate17: did individual survive until end of first growing season (2017) 0-No, 1-Yes
##### Fate18: did individual survive until end of second growing season (2018) 0-No, 1-Yes
##### Tillers17: number of tillers produced by individual at end of the first growing season (2017)
##### Tillers18: number of tillers produced by individual at end of the second growing season (2018) 
##### Biomass18: dry mass in grams from plants harvested to ground level at the end of the second growing season (2018)
##### Inflor17: number of inflorescences produced by an individual plant at the end of the first growing season (2017)
##### Inflor18: number of inflorescences produced by an individual plant at the end of the second growing season (2018)

### dummy code treatment
#dummy code and combine with DF

library(psych)
df_dc = cbind(df, dummy.code(df$Trt))
head(df_dc)
names(df_dc)<- c("Project","Block","Plot","Unique","Species","Trt",      
"Fate17","Fate18","Tillers17","Tillers18","Biomass18","Inflor17", 
"Inflor18","t0","t30x1","t30x2","t70x1","t70x2")
### subset by species

POSE<-subset(df_dc,Species=="POSE")
AGCR<-subset(df_dc,Species=="AGCR")
PSSP<-subset(df_dc,Species=="PSSP")

#load libraries 
library(psych)
library(lavaan)
library(semPlot)

#### remind myself of format
###latent variable definition =~ is measured by
###regression ~ is regressed on
###(residual) (co)variance ~~ is correlated with
###intercept ~ 1 intercept

####
#### Actual analysis as of 14 January 2020
####

library(blavaan)

###AGCR

#lavaan model syntax
model.trt.t ='Fate18~ t30x1+t30x2+t70x1+t70x2
Fate18~Tillers17
Tillers17~t30x1+t30x2+t70x1+t70x2'

#analyze the model
#uses default and suggested prior designation
fit.trt.t.agcr=bsem(model=model.trt.t,data= AGCR,cp="srs",n.chains=3)
# show summary
summary(fit.trt.t.agcr,standardize=T,rsq=T)
fitMeasures(fit.trt.t.agcr)

#### PSSP

#lavaan model syntax
model.trt.t ='Fate18~ t30x1+t30x2+t70x1+t70x2
Fate18~Tillers17
Tillers17~t30x1+t30x2+t70x1+t70x2'
#analyze the model
fit.trt.t.pssp=bsem(model=model.trt.t,data= PSSP,cp="srs",n.chains=3)
# show summary
summary(fit.trt.t.pssp,standardize=T,rsq=T)
fitMeasures(fit.trt.t.pssp)


#### POSE


#lavaan model syntax
model.trt.t ='Fate18~ t30x1+t30x2+t70x1+t70x2
Fate18~Tillers17
Tillers17~t30x1+t30x2+t70x1+t70x2'
#analyze the model
#uses default and suggested prior designation
fit.trt.t.pose=bsem(model=model.trt.t,data= POSE,cp="srs",n.chains=3)
# show summary
summary(fit.trt.t.pose,standardize=T,rsq=T)
fitMeasures(fit.trt.t.pose)



