#####
##### Code to analyses SD-IF data as binomial model trying to discover the probability of survival
#####
##### original author Mevin Hooten, modified by Elise Denton
##### Date: 3 Dec 2018
##### Updated: 20 July 2020

### clear working directory

rm(list=ls())

####
#### you will first need to scroll nearly to the end of this code to see where I moved the 
#### data load in. This seems strange, but I don't need to mess with it once it's loaded
#### and it was getting frustrating constantly scrolling past it.
####


####
####  Compile MCMC Algorithm
####

### set working directory


setwd("C:/Users/replace_with_appropriate_directory/")
dir()

####
#### call up data from below to use
####

AGCR.plot$palive18<-AGCR.plot$Alive.18/AGCR.plot$Count.18
PSSP.plot$palive18<-PSSP.plot$Alive.18/PSSP.plot$Count.18
POSE.plot$palive18<-POSE.plot$Alive.18/POSE.plot$Count.18


###
### 2018, second year survival
###


####
#### Set starting parameters AGCR
####

Y.AGCR=AGCR.plot$Alive.18
n.AGCR=AGCR.plot$Count.18
TRT.AGCR=AGCR.plot$Trt
S2.AGCR=sd(AGCR.plot$Alive.18/AGCR.plot$Count.18)^2
T2=2
its=20000

THETA.start.AGCR.HIGH=.9
THETA.start.AGCR.MED=.5
THETA.start.AGCR.LOW=.1
AGCR.prior=.69

####
#### Set starting parameters PSSP
####

Y.PSSP=PSSP.plot$Alive.18
n.PSSP=PSSP.plot$Count.18
TRT.PSSP=PSSP.plot$Trt
S2.PSSP=sd(PSSP.plot$Alive.18/PSSP.plot$Count.18)^2
T2=2
its=20000

THETA.start.PSSP.HIGH=.90
THETA.start.PSSP.MED=.50
THETA.start.PSSP.LOW=.10
PSSP.prior=.61

####
#### Set starting parameters POSE
####

Y.POSE=POSE.plot$Alive.18
n.POSE=POSE.plot$Count.18
TRT.POSE=POSE.plot$Trt

S2.POSE=sd(POSE.plot$Alive.18/POSE.plot$Count.18)^2
T2=2
its=20000

THETA.start.POSE.HIGH=.90
THETA.start.POSE.MED=.50
THETA.start.POSE.LOW=.10
POSE.prior=.65

####
#### load program
####


source("binom.logit.mcmc_seperatepriors_and_starting_values.R")

####
####  Fit Model to Data for each species individually, start with AGCR 
####

#### funtion (y,N,X,theta.start,theta.p,sd,theta.tune,n.mcmc)
###run for low, medium and high chains
binom.logit.out.AGCR.HIGH=binom.logit.mcmc(Y.AGCR,n.AGCR,TRT.AGCR,THETA.start.AGCR.HIGH,AGCR.prior,S2.AGCR,T2,its)
binom.logit.out.AGCR.MED=binom.logit.mcmc(Y.AGCR,n.AGCR,TRT.AGCR,THETA.start.AGCR.MED,AGCR.prior,S2.AGCR,T2,its)
binom.logit.out.AGCR.LOW=binom.logit.mcmc(Y.AGCR,n.AGCR,TRT.AGCR,THETA.start.AGCR.LOW,AGCR.prior,S2.AGCR,T2,its)


binom.logit.out.PSSP.HIGH=binom.logit.mcmc(Y.PSSP,n.PSSP,TRT.PSSP,THETA.start.PSSP.HIGH,PSSP.prior,S2.PSSP,T2,its)
binom.logit.out.PSSP.MED=binom.logit.mcmc(Y.PSSP,n.PSSP,TRT.PSSP,THETA.start.PSSP.MED,PSSP.prior,S2.PSSP,T2,its)
binom.logit.out.PSSP.LOW=binom.logit.mcmc(Y.PSSP,n.PSSP,TRT.PSSP,THETA.start.PSSP.LOW,PSSP.prior,S2.PSSP,T2,its)


binom.logit.out.POSE.HIGH=binom.logit.mcmc(Y.POSE,n.POSE,TRT.POSE,THETA.start.POSE.HIGH,POSE.prior,S2.POSE,T2,its)
binom.logit.out.POSE.MED=binom.logit.mcmc(Y.POSE,n.POSE,TRT.POSE,THETA.start.POSE.MED,POSE.prior,S2.POSE,T2,its)
binom.logit.out.POSE.LOW=binom.logit.mcmc(Y.POSE,n.POSE,TRT.POSE,THETA.start.POSE.LOW,POSE.prior,S2.POSE,T2,its)



####
####  View Trace Plot and M-H acceptance rate 
####


plot.theta<-function(dataset){
par(mfrow=c(3,2))
plot(dataset$theta1.save,type="l",ylab=expression(theta1),
  main=paste("Acceptance: ",dataset$mh.theta1/dataset$n.mcmc))
plot(dataset$theta2.save,type="l",ylab=expression(theta2),
  main=paste("Acceptance: ",dataset$mh.theta2/dataset$n.mcmc))
plot(dataset$theta3.save,type="l",ylab=expression(theta3),
  main=paste("Acceptance: ",dataset$mh.theta3/dataset$n.mcmc))
plot(dataset$theta4.save,type="l",ylab=expression(theta4),
  main=paste("Acceptance: ",dataset$mh.theta4/dataset$n.mcmc))
plot(dataset$theta5.save,type="l",ylab=expression(theta5),
  main=paste("Acceptance: ",dataset$mh.theta5/dataset$n.mcmc))
plot(dataset$s2.save,type="l",ylab=expression(s2),
  main=paste("Acceptance: ",dataset$mh.s2/dataset$n.mcmc))

}

###plot AGCR Trace plots
plot.theta(binom.logit.out.AGCR.HIGH)
plot.theta(binom.logit.out.AGCR.MED)
plot.theta(binom.logit.out.AGCR.LOW)

### plot PSSPS trace plots
plot.theta(binom.logit.out.PSSP.HIGH)
plot.theta(binom.logit.out.PSSP.MED)
plot.theta(binom.logit.out.PSSP.LOW)

### plot POSE trace plots
plot.theta(binom.logit.out.POSE.HIGH)
plot.theta(binom.logit.out.POSE.MED)
plot.theta(binom.logit.out.POSE.LOW)


####
#### Calculate Gelman and Rubin statistic
####
# m number of chains
# n number of iterations is iteration minus burn in
# x1 is values from high chain
# x2 is values from medium chain
# x3 is values from low chain
# psi is grand mean
# psi.j various are means of individual chains
# B over n winds up being the between chain variance
# W winds up being the within-chain variance
# V.hat is an estimate of the true (but unknown) variance
# R.hat, is an estimate of convergence, if near one convergence can be concluded

#set function
minus.psi <- function(x, psi){
h<-x-psi
return(h)
}

g.r.3<-function(x1,x2,x3,m,n){
psi.j1<-mean(x1)
psi.j2<-mean(x2)
psi.j3<-mean(x3)
psi.j<-c(psi.j1,psi.j2,psi.j3)
psi<-mean(c(psi.j1,psi.j2,psi.j3))
B.over.n<-(1/(m-1))*sum((minus.psi(psi.j,psi))^2)
W<-(1/(m*(n-1)))*(sum(c((sum((minus.psi(as.data.frame(x1),psi.j[1]))^2)),(sum((minus.psi(as.data.frame(x2),psi.j[2]))^2)),(sum((minus.psi
(as.data.frame(x3),psi.j[3]))^2)))))
s2.hat.plus<-((n-1)/n)*W+(B.over.n)
V.hat<-s2.hat.plus+(B.over.n*(1/3))
R.hat<-V.hat/W
return(R.hat)
}

####
#### print gelman rubin stats
####


g.r.stats<-function(ds.1,ds.2,ds.3,m,n){
gr.theta1<-g.r.3(ds.1$theta1.save,ds.2$theta1.save,ds.3$theta1.save,m,n)
gr.theta2<-g.r.3(ds.1$theta2.save,ds.2$theta2.save,ds.3$theta2.save,m,n)
gr.theta3<-g.r.3(ds.1$theta3.save,ds.2$theta3.save,ds.3$theta3.save,m,n)
gr.theta4<-g.r.3(ds.1$theta4.save,ds.2$theta4.save,ds.3$theta4.save,m,n)
gr.theta5<-g.r.3(ds.1$theta5.save,ds.2$theta5.save,ds.3$theta5.save,m,n)
gr.s2<-g.r.3(ds.1$s2.save,ds.2$s2.save,ds.3$s2.save,m,n)
gr.values<-c(gr.theta1,gr.theta2,gr.theta3,gr.theta4,gr.theta4,gr.s2)
stats.labels<-c("theta1","theta2","theta3","theta4","theta5","s2")
gr.stats<-as.data.frame(cbind(stats.labels,gr.values))
print(gr.stats)
}

### GR stats for each species
g.r.stats(binom.logit.out.AGCR.HIGH,binom.logit.out.AGCR.MED,binom.logit.out.AGCR.LOW,3,its*.2)

g.r.stats(binom.logit.out.PSSP.HIGH,binom.logit.out.PSSP.MED,binom.logit.out.PSSP.LOW,3,its*.2)

g.r.stats(binom.logit.out.POSE.HIGH,binom.logit.out.POSE.MED,binom.logit.out.POSE.LOW,3,its*.2)






####
####  violin graphs of MCMC output
####

#### note treatment must be factor I think
data_summary <- function(x) {
   m <- mean(x,na.rm=TRUE)
   ymin <- m-sd(x)*1.96
   ymax <- m+sd(x)*1.96
   return(c(y=m,ymin=ymin,ymax=ymax))
}

### write function to draw violin graphs
violin.graphs<-function(ds.1,ds.2,ds.3){
TRT.1<-rep("Control",length(ds.1$theta1.save)*3)
THETA.TRT1<-as.numeric(c(ds.1$theta1.save,ds.2$theta1.save,ds.3$theta1.save))
TRT.2<-rep("30x1",length(ds.1$theta2.save)*3)
THETA.TRT2<-as.numeric(c(ds.1$theta2.save,ds.2$theta2.save,ds.3$theta2.save))
TRT.3<-rep("30x2",length(ds.1$theta3.save)*3)
THETA.TRT3<-as.numeric(c(ds.1$theta3.save,ds.2$theta3.save,ds.3$theta3.save))
TRT.4<-rep("70x1",length(ds.1$theta4.save)*3)
THETA.TRT4<-as.numeric(c(ds.1$theta4.save,ds.2$theta4.save,ds.3$theta4.save))
TRT.5<-rep("70x2",length(ds.1$theta5.save)*3)
THETA.TRT5<-as.numeric(c(ds.1$theta5.save,ds.2$theta5.save,ds.3$theta5.save))


fit.1<-data.frame(TRT.1,THETA.TRT1)
names(fit.1)<-c("TRT","THETA")
fit.2<-data.frame(TRT.2,THETA.TRT2)
names(fit.2)<-c("TRT","THETA")
fit.3<-data.frame(TRT.3,THETA.TRT3)
names(fit.3)<-c("TRT","THETA")
fit.4<-data.frame(TRT.4,THETA.TRT4)
names(fit.4)<-c("TRT","THETA")
fit.5<-data.frame(TRT.5,THETA.TRT5)
names(fit.5)<-c("TRT","THETA")

FIT<-rbind(fit.1,fit.2,fit.3,fit.4,fit.5)

library(ggplot2)
p<-ggplot(FIT, aes(x=TRT, y=THETA, fill=TRT)) + 
  geom_violin(trim=TRUE)+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))


greycolors<-c("white","grey","black","grey","black")

q<-p+scale_fill_manual(values=greycolors) + theme_classic()+
theme(axis.ticks.length=unit(-0.25, "cm"), axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))+
labs(title="",x="Defoliation Level", y = "% Survival")

print(q)

}	


### make sure to turn on graph tracking before running this code so 
### you can see both graphs produced
### click in graphics window, Then chose history from the file menu 
### and turn on tracking with the first option
### use page up and page down to navigate between graphs

### violin graphs for each species

violin.graphs(binom.logit.out.AGCR.HIGH,binom.logit.out.AGCR.MED,binom.logit.out.AGCR.LOW)

violin.graphs(binom.logit.out.PSSP.HIGH,binom.logit.out.PSSP.MED,binom.logit.out.PSSP.LOW)

violin.graphs(binom.logit.out.POSE.HIGH,binom.logit.out.POSE.MED,binom.logit.out.POSE.LOW)



####
####  Compute Posterior Summary Statistics
####


###
### theta, aka probability of survival
###


means.out<-function(ds.1,ds.2,ds.3){

# Create the function.
#getmode <- function(v) {
#   uniqv <- unique(v)
#   uniqv[which.max(tabulate(match(v, uniqv)))]
#}

theta1<-c(ds.1$theta1.save,ds.2$theta1.save,ds.3$theta1.save)
theta2<-c(ds.1$theta2.save,ds.2$theta2.save,ds.3$theta2.save)
theta3<-c(ds.1$theta3.save,ds.2$theta3.save,ds.3$theta3.save)
theta4<-c(ds.1$theta4.save,ds.2$theta4.save,ds.3$theta4.save)
theta5<-c(ds.1$theta5.save,ds.2$theta5.save,ds.3$theta5.save)
col.1<-as.data.frame(c(mean(theta1),mean(theta2),mean(theta3),mean(theta4),mean(theta5)))
names(col.1)<-"prob.survival.mean"

col.2n3<-rbind(quantile(theta1,c(0.025,0.975)),quantile(theta2,c(0.025,0.975)),
quantile(theta3,c(0.025,0.975)),
quantile(theta4,c(0.025,0.975)),
quantile(theta5,c(0.025,0.975)))


diff<-as.data.frame(col.1-col.1[1,])
names(diff)<-"diff from control"

### probability different from control

SET<-cbind(col.1,col.2n3,diff)

return(SET)
}

#### calculate survival probabilities for each species
means.out(binom.logit.out.AGCR.HIGH,binom.logit.out.AGCR.MED,binom.logit.out.AGCR.LOW)

means.out(binom.logit.out.PSSP.HIGH,binom.logit.out.PSSP.MED,binom.logit.out.PSSP.LOW)

means.out(binom.logit.out.POSE.HIGH,binom.logit.out.POSE.MED,binom.logit.out.POSE.LOW)



####
#### sd
####


####
#### probablyity that treatments are different from each other
####

# probability that ____ is greater than distrubution 1

compare.means<-function(ds.1,ds.2,ds.3){
theta1<-c(ds.1$theta1.save,ds.2$theta1.save,ds.3$theta1.save)
theta2<-c(ds.1$theta2.save,ds.2$theta2.save,ds.3$theta2.save)
theta3<-c(ds.1$theta3.save,ds.2$theta3.save,ds.3$theta3.save)
theta4<-c(ds.1$theta4.save,ds.2$theta4.save,ds.3$theta4.save)
theta5<-c(ds.1$theta5.save,ds.2$theta5.save,ds.3$theta5.save)

#randomize samples
t1<-sample(theta1)
t2<-sample(theta2)
t3<-sample(theta3)
t4<-sample(theta4)
t5<-sample(theta5)

#probability that treatment is less than control by .2
col.1<-as.data.frame(rep(0,5))
col.1[1,]<-1-ecdf(t1-t1)(.2)
col.1[2,]<-1-ecdf(t1-t2)(0.2)
col.1[3,]<-1-ecdf(t1-t3)(0.2)
col.1[4,]<-1-ecdf(t1-t4)(0.2)
col.1[5,]<-1-ecdf(t1-t5)(0.2)
names(col.1)<-c("p.less.than.control.2")


#probability that treatment is less than control by .1
col.2<-as.data.frame(rep(0,5))
col.2[1,]<-1-ecdf(t1-t1)(.1)
col.2[2,]<-1-ecdf(t1-t2)(0.1)
col.2[3,]<-1-ecdf(t1-t3)(0.1)
col.2[4,]<-1-ecdf(t1-t4)(0.1)
col.2[5,]<-1-ecdf(t1-t5)(0.1)
names(col.2)<-c("p.less.than.control.1")

#probability that treatment is less than control by .05
col.3<-as.data.frame(rep(0,5))
col.3[1,]<-1-ecdf(t1-t1)(.05)
col.3[2,]<-1-ecdf(t1-t2)(0.05)
col.3[3,]<-1-ecdf(t1-t3)(0.05)
col.3[4,]<-1-ecdf(t1-t4)(0.05)
col.3[5,]<-1-ecdf(t1-t5)(0.05)
names(col.3)<-c("p.less.than.control.05")


#probability that treatment is greater than control by at least .05
col.4<-as.data.frame(rep(0,5))
col.4[1,]<-1-ecdf(t1-t1)(0.05)
col.4[2,]<-1-ecdf(t2-t1)(0.05)
col.4[3,]<-1-ecdf(t3-t1)(0.05)
col.4[4,]<-1-ecdf(t4-t1)(0.05)
col.4[5,]<-1-ecdf(t5-t1)(0.05)
names(col.4)<-c("p.more.at.05")

#probability that treatment is greater than control by at least .1
col.5<-as.data.frame(rep(0,5))
col.5[1,]<-1-ecdf(t1-t1)(0.1)
col.5[2,]<-1-ecdf(t2-t1)(0.1)
col.5[3,]<-1-ecdf(t3-t1)(0.1)
col.5[4,]<-1-ecdf(t4-t1)(0.1)
col.5[5,]<-1-ecdf(t5-t1)(0.1)
names(col.5)<-c("p.more.at.1")

#probability that treatment is greater than control by at least .2
col.6<-as.data.frame(rep(0,5))
col.6[1,]<-1-ecdf(t1-t1)(0.2)
col.6[2,]<-1-ecdf(t2-t1)(0.2)
col.6[3,]<-1-ecdf(t3-t1)(0.2)
col.6[4,]<-1-ecdf(t4-t1)(0.2)
col.6[5,]<-1-ecdf(t5-t1)(0.2)
names(col.6)<-c("p.more.at.2")

combine<-cbind(col.1,col.2,col.3,col.4,col.5,col.6)


print(combine)
}

#### calculate Sd for each species

compare.means(binom.logit.out.AGCR.HIGH,binom.logit.out.AGCR.MED,binom.logit.out.AGCR.LOW)
compare.means(binom.logit.out.PSSP.HIGH,binom.logit.out.PSSP.MED,binom.logit.out.PSSP.LOW)
compare.means(binom.logit.out.POSE.HIGH,binom.logit.out.POSE.MED,binom.logit.out.POSE.LOW)




####
#### contrasts of interest

#### new contrasts of interest
#### 1x vs control, 2x vs control, defoliation vs control
#### 30 vs control and 70 vs control 

#### write function to generate contrasts just listed

posterior.contrasts<-function(ds.1,ds.2,ds.3){

#shorten names for ease of entry
theta1<-c(ds.1$theta1.save,ds.2$theta1.save,ds.3$theta1.save)
theta2<-c(ds.1$theta2.save,ds.2$theta2.save,ds.3$theta2.save)
theta3<-c(ds.1$theta3.save,ds.2$theta3.save,ds.3$theta3.save)
theta4<-c(ds.1$theta4.save,ds.2$theta4.save,ds.3$theta4.save)
theta5<-c(ds.1$theta5.save,ds.2$theta5.save,ds.3$theta5.save)

#randomize samples
t1<-sample(theta1)
t2<-sample(theta2)
t3<-sample(theta3)
t4<-sample(theta4)
t5<-sample(theta5)

#proportion of (a) distribution that is less than (b) distrubution by at least .2
out.1<-as.data.frame(rep(0,7))

out.1[1,]<-1-ecdf(c(t1,t1)-c(t2,t3))(0.2)  # 30 vs control
out.1[2,]<-1-ecdf(c(t1,t1)-c(t4,t5))(0.2)  # 70 vs control

out.1[3,]<-1-ecdf(c(t4,t5)-c(t2,t3))(0.2) # 30 vs 70
out.1[4,]<-1-ecdf(c(t3,t5)-c(t2,t4))(0.2)  # 1x vs 2x
out.1[5,]<-1-ecdf(c(t1,t1)-c(t2,t4))(0.2) # 1x vs control
out.1[6,]<-1-ecdf(c(t1,t1)-c(t3,t5))(0.2) # 2x vs control
out.1[7,]<-1-ecdf(c(t1,t1,t1,t1)-c(t2,t3,t4,t5))(0.2) # defoliation vs control


#proportion of (a) distribution that is less than (b) distrubution by at least .1
out.2<-as.data.frame(rep(0,7))

out.2[1,]<-1-ecdf(c(t1,t1)-c(t2,t3))(0.1)  # 30 vs control
out.2[2,]<-1-ecdf(c(t1,t1)-c(t4,t5))(0.1)  # 70 vs control

out.2[3,]<-1-ecdf(c(t4,t5)-c(t2,t3))(0.1) # 30 vs 70
out.2[4,]<-1-ecdf(c(t3,t5)-c(t2,t4))(0.1)  # 1x vs 2x
out.2[5,]<-1-ecdf(c(t1,t1)-c(t2,t4))(0.1) # 1x vs control
out.2[6,]<-1-ecdf(c(t1,t1)-c(t3,t5))(0.1) # 2x vs control
out.2[7,]<-1-ecdf(c(t1,t1,t1,t1)-c(t2,t3,t4,t5))(0.1) # defoliation vs control

#proportion of (a) distribution that is less than (b) distrubution by at least .05
out.3<-as.data.frame(rep(0,7))

out.3[1,]<-1-ecdf(c(t1,t1)-c(t2,t3))(0.05)  # 30 vs control
out.3[2,]<-1-ecdf(c(t1,t1)-c(t4,t5))(0.05)  # 70 vs control

out.3[3,]<-1-ecdf(c(t4,t5)-c(t2,t3))(0.05) # 30 vs 70
out.3[4,]<-1-ecdf(c(t3,t5)-c(t2,t4))(0.05)  # 1x vs 2x
out.3[5,]<-1-ecdf(c(t1,t1)-c(t2,t4))(0.05) # 1x vs control
out.3[6,]<-1-ecdf(c(t1,t1)-c(t3,t5))(0.05) # 2x vs control
out.3[7,]<-1-ecdf(c(t1,t1,t1,t1)-c(t2,t3,t4,t5))(0.05) # defoliation vs control


#proportion of (a) distribution that is greater than (b) distrubution by at least .05
out.4<-as.data.frame(rep(0,7))

out.4[1,]<-1-ecdf(c(t2,t3)-c(t1,t1))(0.05)  # 30 vs control
out.4[2,]<-1-ecdf(c(t4,t5)-c(t1,t1))(0.05)  # 70 vs control

out.4[3,]<-1-ecdf(c(t2,t3)-c(t4,t5))(0.05) # 30 vs 70
out.4[4,]<-1-ecdf(c(t2,t4)-c(t3,t5))(0.05)  # 1x vs 2x
out.4[5,]<-1-ecdf(c(t2,t4)-c(t1,t1))(0.05) # 1x vs control
out.4[6,]<-1-ecdf(c(t3,t5)-c(t1,t1))(0.05) # 2x vs control
out.4[7,]<-1-ecdf(c(t2,t3,t4,t5)-c(t1,t1,t1,t1))(0.05) # defoliation vs control


#proportion of (a) distribution that is greater than (b) distrubution by at least .1
out.5<-as.data.frame(rep(0,7))

out.5[1,]<-1-ecdf(c(t2,t3)-c(t1,t1))(0.1)  # 30 vs control
out.5[2,]<-1-ecdf(c(t4,t5)-c(t1,t1))(0.1)  # 70 vs control

out.5[3,]<-1-ecdf(c(t2,t3)-c(t4,t5))(0.1) # 30 vs 70
out.5[4,]<-1-ecdf(c(t2,t4)-c(t3,t5))(0.1)  # 1x vs 2x
out.5[5,]<-1-ecdf(c(t2,t4)-c(t1,t1))(0.1) # 1x vs control
out.5[6,]<-1-ecdf(c(t3,t5)-c(t1,t1))(0.1) # 2x vs control
out.5[7,]<-1-ecdf(c(t2,t3,t4,t5)-c(t1,t1,t1,t1))(0.1) # defoliation vs control


#proportion of (a) distribution that is greater than (b) distrubution by at least .2
out.6<-as.data.frame(rep(0,7))
out.6[1,]<-1-ecdf(c(t2,t3)-c(t1,t1))(0.05)  # 30 vs control
out.6[2,]<-1-ecdf(c(t4,t5)-c(t1,t1))(0.05)  # 70 vs control

out.6[3,]<-1-ecdf(c(t2,t3)-c(t4,t5))(0.2) # 30 vs 70
out.6[4,]<-1-ecdf(c(t2,t4)-c(t3,t5))(0.2)  # 1x vs 2x
out.6[5,]<-1-ecdf(c(t2,t4)-c(t1,t1))(0.2) # 1x vs control
out.6[6,]<-1-ecdf(c(t3,t5)-c(t1,t1))(0.2) # 2x vs control
out.6[7,]<-1-ecdf(c(t2,t3,t4,t5)-c(t1,t1,t1,t1))(0.2) # defoliation vs control

out<-as.data.frame(rep(NA,7))
out<-c("30 vs con", "70 vs con",  "30 vs 70","x1 vs x2","1x vs Control","2x vs control","defoliation vs control")
names(out.1)<-"p a < b by .2"
names(out.2)<-"p a < b by .1"
names(out.3)<-"p a < b by .05"
names(out.4)<-"p a > b by .05"
names(out.5)<-"p a > b by .1"
names(out.6)<-"p a > b by .2"

compile<-cbind(out,out.1,out.2,out.3,out.4,out.5,out.6)
return(compile)

}


### calculate posterior contrasts for each species

posterior.contrasts(binom.logit.out.AGCR.HIGH,binom.logit.out.AGCR.MED,binom.logit.out.AGCR.LOW)
posterior.contrasts(binom.logit.out.PSSP.HIGH,binom.logit.out.PSSP.MED,binom.logit.out.PSSP.LOW)
posterior.contrasts(binom.logit.out.POSE.HIGH,binom.logit.out.POSE.MED,binom.logit.out.POSE.LOW)










#####
##### this is the code I was saying needs to be run first
#####



####
#### get data into number of successes per plot, number of failures format
####

### set working directory

setwd("C:/Users/replace_with_appropriate_directory/")
dir()

### load data 2018 Data

df18<-read.csv("SD-IF_2018_Secondseasonsurvival_Final.csv", sep=",", header=T)
head(df18)
dim(df18)

### reduce to only the headings I want
da18<-df18[,c(1:4,7)]
head(da18)
dim(da18)

### make a combined Plot/Color colum to be used to join datasets
da18$Plot_Color <- paste(da18$Plot,da18$Color)







### load 2017 data
### we'll use this to make sure all of the same seedlings are counted in 2017 and 2018

### set working directory

setwd("C:/Users/replace_with_appropriate_directory/")dir()

### load data

df17<-read.csv("SD-IF_2017_endofseasonsurvival_Final.csv", sep=",", header=T)
head(df17)
dim(df17)

#### reduce to only data I want for 2017

da17<-df17[,c(3:6,8,16)]
da17$Fate17<-ifelse(da17$FateF=="D",0,1)
da17<-da17[,-6]
head(da17)
dim(da17)

### make a combined Plot/Color colum
da17$Plot_Color <- paste(da17$Plot,da17$Color)

### merge data frames so we know the fate of all seedlings from 2017
library(plyr)
total<-join(da17, da18, by = "Plot_Color", type = "left", match = "all")
#total$plot<-str_replace_all(string=total$plot, pattern=" ", repl="")
#total$block<-str_replace_all(string=total$block, pattern=" ", repl="")
dim(total)
head(total)
total

### Replace all NAs for 2018 Fate with Dead
total$Fate18[is.na(total$Fate18)] <- "D"

### repalce mixed with Alive to simplify outcomes and because it was definitely alive in 2018
total$Fate_ordered<-ifelse(total$Fate18=="D","D","A")

### create nummerical representation of Alive and dead
total$Fate18_num<-ifelse(total$Fate18=="D",0,1)

### convert some comlumns to factors to avoid problems
total$Block = as.factor(total$Block)
total$Plot = as.factor(total$Plot)
total$Trt = as.factor(total$Trt)
total$Fate_ordered = as.factor(total$Fate_ordered)

head(total)
total

### remove NA for 2017 fate 1854 B1 as it was pulled up and shouldn't count toward survival

### also need to remove 1890 Y and 1880 W for same reasons, in year 2 only
total<-total[-260,]
total<-total[-168,]
total<-total[-35,]
dim(total)


#### don't remember why I wanted these ordered, but here it is
### possibly for easier coding of bayes stats later

total <- total[order(total$Species, total$Trt),]
rownames(total)<-1:nrow(total)

####
#### code treatments as numbers
####


#### create numerical representation of Alive and dead

total$Trt_num<-ifelse(total$Trt=="0",1,ifelse(total$Trt=="30x1",2,ifelse(total$Trt=="30x2",3,ifelse(total$Trt=="70x1",4,5))))
total$Trt_simp<-ifelse(total$Trt=="0",0,1)

#### make block a number again

total$Block<-as.numeric(total$Block)

dim(total)

### remove duplicate columns
total
head(total)

total<-total[,c(1:7,12:16)]

head(total)

#### rename 

####
#### seperate by species 
####


AGCR.total<-subset(total, Species=="AGCR")
head(AGCR.total)
dim(AGCR.total)
AGCR.total
AGCR.total[(50:60),]

PSSP.total<-subset(total, Species=="PSSP")
head(PSSP.total)
dim(PSSP.total)
PSSP.total
PSSP.total[(50:60),]

POSE.total<-subset(total, Species=="POSE")
head(POSE.total)
dim(POSE.total)
POSE.total
POSE.total[(50:60),]

####
#### summarize each species dataset by plot 
####

library(plyr)
AGCR.plot<-ddply(AGCR.total, .(Plot), summarize,  Alive.18=sum(Fate18_num), Count.18=length(Fate18_num),
Block=mean(Block),Trt=mean(Trt_num),Trt_simp=mean(Trt_simp),Alive.17=sum(Fate17,na.rm=TRUE), 
Count.17=length(Fate17))
AGCR.plot

PSSP.plot<-ddply(PSSP.total, .(Plot), summarize,  Alive.18=sum(Fate18_num), Count.18=length(Fate18_num),Block=mean(Block),Trt=mean(Trt_num),Trt_simp=mean(Trt_simp),Alive.17=sum(Fate17,na.rm=TRUE), Count.17=length(Fate17))
PSSP.plot

POSE.plot<-ddply(POSE.total, .(Plot), summarize,  Alive.18=sum(Fate18_num), Count.18=length(Fate18_num),Block=mean(Block),Trt=mean(Trt_num),Trt_simp=mean(Trt_simp),Alive.17=sum(Fate17,na.rm=TRUE), Count.17=length(Fate17))
POSE.plot




