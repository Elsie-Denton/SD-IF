binom.logit.mcmc <- function(y,N,X,theta.start,theta.p,sd,theta.tune,n.mcmc){
####
####  Mevin Hooten		 
####
#######################
####
####  Implements M-H sampler for binomial logit, where y~Binom(N,theta),logit(theta)~N(0,s2).
####
####  y: vector of binomial count data of length n
####  N: number of trials
####  theta: probability of success
####  s2: prior variance on logit(theta)
####  theta.tune: tuning parameter (chosen to enhance mixing) 
####
####  Example Use:
####  
####  binom.logit.mcmc(c(2,3,1,2,2),10,1,1,1000) 
####  
#######################
####
####
#######################
####  Last Updated:
####
####  20110824: Created template for later use. 
####  20110824: Coded sampler for binomial logit model.
####  20191126: added way for priors and starting values to be seperate
####
#######################

####
####  Subroutines and Libraries
####

logit <- function(x){
  log(x/(1-x))
}

logit.inv <- function(x){
  exp(x)/(1+exp(x))
}

####
####  Setup Variables 
####

n=length(y)
burn=n.mcmc*.2
mh.theta1=1
mh.theta2=1
mh.theta3=1
mh.theta4=1
mh.theta5=1
mh.s2=1

s2.save=rep(0,n.mcmc)

theta1.save=rep(0,n.mcmc)
theta2.save=rep(0,n.mcmc)
theta3.save=rep(0,n.mcmc)
theta4.save=rep(0,n.mcmc)
theta5.save=rep(0,n.mcmc)

Dbar.save=rep(0,n.mcmc)

####
####  Priors and Starting Values 
####


s2=sd

theta1=theta.start
theta2=theta.start
theta3=theta.start
theta4=theta.start
theta5=theta.start

theta.prior=theta.p

####
####  Begin Gibbs Loop 
####
  
for(k in 2:n.mcmc){
  cat(k," ");flush.console()

   ###
   ### Sample sigma
   ###
  
   #s2.star=rep(abs(rcauchy(1,0,1)),5)
#abs(rcauchy(1,0,1))
   s2.star=rep(rgamma(1,1,.01),5)
   #s2.star=rep(runif(1,0,10000),5)
   #s2.star=rep(runif(1,0,10000),5)
   s2.last=rep(s2,length(unique(X)))
#plot(density(abs(rcauchy(100,0,1))))

#plot(density(rgamma(100,1,.01)))
#rgamma(100,2,0)
  ####
  ####  Sample theta 
  ####


  ### 0% defoliation

  theta1.star=logit.inv(rnorm(1,logit(theta.prior),theta.tune)) 

#par(mfrow=c(1,1))
#plot(density(logit.inv(rnorm(1000,logit(.89),5)))) 
#hist(logit.inv(rnorm(1000,logit(.89),2)))
  
  mh1.1=sum(dbinom(y[X==1],N[X==1],theta1.star,log=TRUE))+dnorm(logit(theta1.star),0,sqrt(s2.star[1]),log=TRUE)
  mh2.1=sum(dbinom(y[X==1],N[X==1],theta1,log=TRUE))+dnorm(logit(theta1),0,sqrt(s2.last[1]),log=TRUE)
  
  mh.1=exp(mh1.1-mh2.1)

  ### 30% x 1 defoliation

  theta2.star=logit.inv(rnorm(1,logit(theta.prior),theta.tune)) 
  
  mh1.2=sum(dbinom(y[X==2],N[X==2],theta2.star,log=TRUE))+dnorm(logit(theta2.star),0,sqrt(s2.star[2]),log=TRUE)
  mh2.2=sum(dbinom(y[X==2],N[X==2],theta2,log=TRUE))+dnorm(logit(theta2),0,sqrt(s2.last[2]),log=TRUE)

  mh.2=exp(mh1.2-mh2.2)
  


  ### 30% x 2 defoliation

  theta3.star=logit.inv(rnorm(1,logit(theta.prior),theta.tune)) 
  
  mh1.3=sum(dbinom(y[X==3],N[X==3],theta3.star,log=TRUE))+dnorm(logit(theta3.star),0,sqrt(s2.star[3]),log=TRUE)
  mh2.3=sum(dbinom(y[X==3],N[X==3],theta3,log=TRUE))+dnorm(logit(theta3),0,sqrt(s2.last[3]),log=TRUE)
  mh.3=exp(mh1.3-mh2.3)
  

  ### 70% x 1 defoliation

  theta4.star=logit.inv(rnorm(1,logit(theta.prior),theta.tune)) 
  
  mh1.4=sum(dbinom(y[X==4],N[X==4],theta4.star,log=TRUE))+dnorm(logit(theta4.star),0,sqrt(s2.star[4]),log=TRUE)
  mh2.4=sum(dbinom(y[X==4],N[X==4],theta4,log=TRUE))+dnorm(logit(theta4),0,sqrt(s2.last[4]),log=TRUE)
  mh.4=exp(mh1.4-mh2.4)
  


  ### 70% x 2 defoliation

  theta5.star=logit.inv(rnorm(1,logit(theta.prior),theta.tune)) 
 
  mh1.5=sum(dbinom(y[X==5],N[X==5],theta5.star,log=TRUE))+dnorm(logit(theta5.star),0,sqrt(s2.star[5]),log=TRUE)
 mh2.5=sum(dbinom(y[X==5],N[X==5],theta5,log=TRUE))+dnorm(logit(theta5),0,sqrt(s2.last[5]),log=TRUE)
  mh.5=exp(mh1.5-mh2.5)
  
test=runif(1)
 
# if(mh.1 > runif(1)){
#    theta1=theta1.star
#    mh.theta1=mh.theta1+1
#    s2.last[1]=s2.star[1]
#  }
  
# if(mh.2 > runif(1)){
#    theta2=theta2.star
#    mh.theta2=mh.theta2+1
#    s2.last[2]=s2.star[2]
#  }

#if(mh.3 > runif(1)){
#    theta3=theta3.star
#    mh.theta3=mh.theta3+1
#    s2.last[3]=s2.star[3]
#  }

# if(mh.4 > runif(1)){
#    theta4=theta4.star
#    mh.theta4=mh.theta4+1
#    s2.last[4]=s2.star[4]
#  } 

# if(mh.5 > runif(1)){
#    theta5=theta5.star
#    mh.theta5=mh.theta5+1

#    s2.last[5]=s2.star[5]
#  } 

 
 if(mh.1 > test){
    theta1=theta1.star
    mh.theta1=mh.theta1+1
    s2.last[1]=s2.star[1]
  }
 
 if(mh.2 > test){
    theta2=theta2.star
    mh.theta2=mh.theta2+1
    s2.last[2]=s2.star[2]
  }

if(mh.3 > test){
    theta3=theta3.star
    mh.theta3=mh.theta3+1
    s2.last[3]=s2.star[3]
  }

 if(mh.4 > test){
    theta4=theta4.star
    mh.theta4=mh.theta4+1
    s2.last[4]=s2.star[4]
  } 

 if(mh.5 > test){
    theta5=theta5.star
    mh.theta5=mh.theta5+1
    s2.last[5]=s2.star[5]
  } 

 if(((mh.5+mh.4+mh.3+mh.2+mh.1)/5) > runif(1)){
    s2.last=s2.star
    mh.s2=mh.s2+1
  } 



  ###
  ### DIC Calculations 
  ###
      Dbar.tmp=rep(0,5)
      Dbar.tmp[1]=-2*sum(dbinom(y[X==1],N[X==1],theta1,log=TRUE))
      Dbar.tmp[2]=-2*sum(dbinom(y[X==2],N[X==2],theta2,log=TRUE))
      Dbar.tmp[3]=-2*sum(dbinom(y[X==3],N[X==3],theta3,log=TRUE))
      Dbar.tmp[4]=-2*sum(dbinom(y[X==4],N[X==4],theta4,log=TRUE))
      Dbar.tmp[5]=-2*sum(dbinom(y[X==5],N[X==5],theta5,log=TRUE))
    
      Dbar.save[k]=mean(Dbar.tmp)


 ####
  ####  Save Samples 
  ####

  theta1.save[k]=theta1
  theta2.save[k]=theta2
  theta3.save[k]=theta3
  theta4.save[k]=theta4
  theta5.save[k]=theta5
  s2=mean(s2.last)
  s2.save[k]=s2

}
cat("\n");flush.console()

####
#### Calculate DIC and Print to Screen
####

#### Note DIC may not be valid for this model, but I don't know what else to use

Dhat.tmp=rep(0,5)
Dhat.tmp[1]=-2*sum(dbinom(y[X==1],N[X==1],mean(theta1.save[-(1:burn)]),log=TRUE))
Dhat.tmp[2]=-2*sum(dbinom(y[X==2],N[X==2],mean(theta2.save[-(1:burn)]),log=TRUE))
Dhat.tmp[3]=-2*sum(dbinom(y[X==3],N[X==3],mean(theta3.save[-(1:burn)]),log=TRUE))
Dhat.tmp[4]=-2*sum(dbinom(y[X==4],N[X==4],mean(theta4.save[-(1:burn)]),log=TRUE))
Dhat.tmp[5]=-2*sum(dbinom(y[X==5],N[X==5],mean(theta5.save[-(1:burn)]),log=TRUE))

Dhat.tmp
Dhat=mean(Dhat.tmp)
Dbar=mean(Dbar.save[-(1:burn)])
pD=Dbar-Dhat
DIC=Dhat+2*pD
print(DIC)

####
####  Write Output 
####
 
list(theta1.save=theta1.save[-(1:burn)],mh.theta1=mh.theta1,theta2.save=theta2.save[-(1:burn)],
theta3.save=theta3.save[-(1:burn)],mh.theta3=mh.theta3,
theta4.save=theta4.save[-(1:burn)],mh.theta4=mh.theta4,
theta5.save=theta5.save[-(1:burn)],mh.theta5=mh.theta5,
s2.save=s2.save[-(1:burn)],mh.theta2=mh.theta2,mh.s2=mh.s2,
Dbar.save=Dbar.save[-(1:burn)],DIC=DIC,n.mcmc=n.mcmc)

}
