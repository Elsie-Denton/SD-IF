#####
##### Code to analyze Tiller, height and biomass data from SD-IF experiment
##### SD-IF stands for Seedling defoliation intensity and frequency

##### Elsie Denton
##### Date: 26 Sept 2019

### clear working directory

rm(list=ls())

### import SD-IF dataset 

setwd("C:/Users/elsie.denton/Documents/Projects/SD-IF/Data/2018/Second_Year_Survival_w_Analysis/Final")
dir()
df<-read.csv("SD-IF_CompiledDataToShare_byseedling.csv", sep=",", header=T)
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

### make sure treatment is read as a factor.

df$Trt<-as.factor(df$Trt)

### split by species
pose<-subset(df,Species=="POSE")
agcr<-subset(df,Species=="AGCR")
pssp<-subset(df,Species=="PSSP")


#### means tilllers
head(pose)
#### calculate means
library(plyr)
pose18_tillers_mean<-ddply(pose, .(Trt,Block), summarize,  mean=mean(Tillers18))
pose18_tillers_mean
ddply(pose, .(Block), summarize,  mean=mean(Tillers18),sd=sd(Tillers18),median=median(Tillers18))
ddply(pose, .(Trt), summarize,  mean=mean(Tillers18),sd=sd(Tillers18),median=median(Tillers18))
ddply(pose, .(Trt,Block), summarize,  mean=mean(Tillers18),sd=sd(Tillers18),median=median(Tillers18))

ddply(agcr, .(Block), summarize,  mean=mean(Tillers18),sd=sd(Tillers18),median=median(Tillers18))
ddply(agcr, .(Trt), summarize,  mean=mean(Tillers18),sd=sd(Tillers18),median=median(Tillers18))
ddply(agcr, .(Trt,Block), summarize,  mean=mean(Tillers18),sd=sd(Tillers18),median=median(Tillers18))

ddply(pssp, .(Block), summarize,  mean=mean(Tillers18),sd=sd(Tillers18),median=median(Tillers18))
ddply(pssp, .(Trt), summarize,  mean=mean(Tillers18),sd=sd(Tillers18),median=median(Tillers18))
ddply(pssp, .(Trt,Block), summarize,  mean=mean(Tillers18),sd=sd(Tillers18),median=median(Tillers18))

###
### initial visual exploration of data
###

### Plot the relationship between tillers and biomass in year 2

plot(pose$Tillers18,pose$Biomass18) 
plot(agcr$Tillers18,agcr$Biomass18) 
plot(pssp$Tillers18,pssp$Biomass18) 

### run some quick box plots to help visualize the data

### 2017 tillers
boxplot(Tillers17~Trt,data=pose)#no notable diff
boxplot(Tillers17~Trt,data=agcr)#increase with 2nd clipping
boxplot(Tillers17~Trt,data=pssp)#no notable diff

### 2017 inflorences
boxplot(Inflor17~Trt,data=pose)#none
boxplot(Inflor17~Trt,data=agcr)#no notable diff
boxplot(Inflor17~Trt,data=pssp)#30% only, but could easily be random

### 2018 tillers
boxplot(Tillers18~Trt,data=pose)#no notable diff
boxplot(Tillers18~Trt,data=agcr)#increase with 2nd clipping, possible increase with any clipping
boxplot(Tillers18~Trt,data=pssp)#increase 30x1, maybe

### 2018 inflorences
boxplot(Inflor18~Trt,data=pose)#none
boxplot(Inflor18~Trt,data=agcr)#possible bump in clipped twice
boxplot(Inflor18~Trt,data=pssp)#no decernable pattern

### Second Year Biomass
boxplot(Biomass18~Trt,data=pose)#no decernable pattern
boxplot(Biomass18~Trt,data=agcr)#possible bump in clipped twice
boxplot(Biomass18~Trt,data=pssp)#no decernable pattern


### test for homogeneity of variance
### reference for other tests http://www.cookbook-r.com/Statistical_analysis/Homogeneity_of_variance/

### Fligner-Killeen test used because data is highly non-normal


### 2nd year biomass <- two pass, one fail
fligner.test(Biomass18~interaction(Trt,Block),data=pose)#pass p=.06
fligner.test(Biomass18~interaction(Trt,Block),data=agcr)#fail p=0.00
fligner.test(Biomass18~interaction(Trt,Block),data=pssp)#pass p=0.09

### so, since homogenity of variance is not a given even with a lenient non-parametric test
### we'll need to use non-parametric anova, (or transform variables)

### okay, non-parametric tests don't work well for the data I have, as it doesn't handle unequal sample sizes
### we'll try a log transform
### stratch that, log transforms suck on 0s and I have a bunch of those we'll start with sqrt


### 2nd year biomass, sqrt <- two pass, one fail
fligner.test(sqrt(Biomass18)~interaction(Trt,Block),data=pose)#pass p=.06
fligner.test(sqrt(Biomass18)~interaction(Trt,Block),data=agcr)#fail p=0.01
fligner.test(sqrt(Biomass18)~interaction(Trt,Block),data=pssp)#pass p=0.09

### sqrt transform improves data in most cases

### I think a type 1 anova with weighted means makes the most sense for this data
### https://www.r-bloggers.com/r-tutorial-series-two-way-anova-with-unequal-sample-sizes/

### changed my mind and switched to a test that can handle random effects, because that is
### what block is
library(nlme)
library(lme4)


### 2nd year biomass (agcr failed variance test even after transform, but not by much, likely okay)

pose_mass_18 <- lme(sqrt(Biomass18) ~ Trt,random = ~1|Block, data=pose, method = "ML")
summary(pose_mass_18)# no difference
anova(pose_mass_18)

agcr_mass_18 <- lme(sqrt(Biomass18) ~ Trt,random = ~1|Block, data=agcr, method = "ML")
summary(agcr_mass_18)# no differences
anova(agcr_mass_18)

pssp_mass_18 <- lme(sqrt(Biomass18) ~ Trt,random = ~1|Block, data=pssp, method = "ML")
summary(pssp_mass_18)# no differences
anova(pssp_mass_18)

### get estmates for all species since there are no treatment differences

pose.mean<-ddply(pose, .(), summarize,  mean=mean(Biomass18),sd=sd(Biomass18),median=median(Biomass18))
pose.mean$l.ci<-pose.mean$mean-1.96*(pose.mean$sd/(length(pose)))
pose.mean$u.ci<-pose.mean$mean+1.96*(pose.mean$sd/(length(pose)))

pssp.mean<-ddply(pssp, .(), summarize,  mean=mean(Biomass18),sd=sd(Biomass18),median=median(Biomass18))
pssp.mean$l.ci<-pssp.mean$mean-1.96*(pssp.mean$sd/(length(pssp)))
pssp.mean$u.ci<-pssp.mean$mean+1.96*(pssp.mean$sd/(length(pssp)))

agcr.mean<-ddply(agcr, .(), summarize,  mean=mean(Biomass18),sd=sd(Biomass18),median=median(Biomass18))
agcr.mean$l.ci<-agcr.mean$mean-1.96*(agcr.mean$sd/(length(agcr)))
agcr.mean$u.ci<-agcr.mean$mean+1.96*(agcr.mean$sd/(length(agcr)))


pose.mean$.id<-"POSE"

pssp.mean$.id<-"PSSP"
agcr.mean$.id<-"AGCR"
species.biomass<-rbind(pose.mean,pssp.mean,agcr.mean)
species.biomass


### Tillers removed from the above analysis because they are count data and zero inflation
### means that they can't manage a good approximation of normality
### could try this with inflorences too, but might not help

### https://stats.idre.ucla.edu/r/faq/random-coefficient-poisson-models/
### also https://www.imachordata.com/do-not-log-transform-count-data-bitches/
agcr
library(lme4)
library(R2admb)
require(glmmADMB)

### overdispersion was a problem so switched to negative binomial model, these rutinely had better qq plots
### and residual plots than poisson models
### deviance explained from null model used as effect size, also chisq test #https://www.ssc.wisc.edu/sscc/pubs/MM/MM_TestEffects.html


### zero inflated models might be even better than neg binomial models, but Finding code for that is hard

### Tillers 2017 (PSSP a bit wanky, the others pretty great)


###
### zero inflated negative binomail model for tillers
###


library(GLMMadaptive)


z_pose_tiller_17 <- mixed_model(Tillers17 ~ Trt, random = ~ 1 | Block, data = pose,family = zi.negative.binomial(), 
zi_fixed = ~ Trt,max_coef_value=20)

summary(z_pose_tiller_17)
fixef(z_pose_tiller_17)

## qq plot
plot(ranef(z_pose_tiller_17))



### AGCR

z_agcr_tiller_17 <- mixed_model(Tillers17 ~ Trt, random = ~ 1 | Block, data = agcr,family = zi.negative.binomial(), 
zi_fixed = ~ Trt,max_coef_value=20)

summary(z_agcr_tiller_17)
fixef(z_agcr_tiller_17)

## qq plot
plot(ranef(z_agcr_tiller_17))


### PSSP

z_pssp_tiller_17 <- mixed_model(Tillers17 ~ Trt, random = ~ 1 | Block, data = pssp,family = zi.negative.binomial(), 
zi_fixed = ~ Trt,max_coef_value=20)

summary(z_pssp_tiller_17)
fixef(z_pssp_tiller_17)

## qq plot
plot(ranef(z_pssp_tiller_17))




### Tillers 2018 (zero inflated negative bionomial model)

### POSE

z_pose_tiller_18 <- mixed_model(Tillers18 ~ Trt, random = ~ 1 | Block, data = pose,family = zi.negative.binomial(), 
zi_fixed = ~ Trt,max_coef_value=20)

summary(z_pose_tiller_18)
fixef(z_pose_tiller_18)

## qq plot
plot(ranef(z_pose_tiller_18))


### AGCR


z_agcr_tiller_18 <- mixed_model(Tillers18 ~ Trt, random = ~ 1 | Block, data = agcr,family = zi.negative.binomial(), 
zi_fixed = ~ Trt,max_coef_value=20)

summary(z_agcr_tiller_18)
fixef(z_agcr_tiller_18)

### well I just found a way to test contrasts using a wald matrix which is cool
anova(z_agcr_tiller_18,L=rbind(c(0,1,0,0,0),c(0,0,1,0,0)))
anova(z_agcr_tiller_18,L=rbind(c(0,0,0,1,0),c(0,0,0,0,1)))
anova(z_agcr_tiller_18,L=rbind(c(0,1,1,0,0),c(0,0,0,1,1)))
anova(z_agcr_tiller_18,L=rbind(c(1,0,0,0,0),c(0,1,1,0,0)))
anova(z_agcr_tiller_18,L=rbind(c(1,0,0,0,0),c(0,0,0,1,1)))

## qq plot
plot(ranef(z_agcr_tiller_18))

### PSSP

z_pssp_tiller_18 <- mixed_model(Tillers18 ~ Trt, random = ~ 1 | Block, data = pssp,family = zi.negative.binomial(), 
zi_fixed = ~ Trt,max_coef_value=20)

summary(z_pssp_tiller_18)
fixef(z_agcr_tiller_18)

## qq plot
plot(ranef(z_pssp_tiller_18))



### General Conclusions ####
### No difference in biomass based on treatments
### some small differences seen in tillering based on negative binomial model with zero inflation
### not enough infloresences occured to analyze.


