#Set up 
rm(list=ls())
library(rgenoud); library(MASS)
set.seed(123456)

source("sdmiwp.r") 
source("sdiwp.r") 

###############################Define Objects################################
g2.m.o<-NULL; g2.m.p<-NULL
psi.o<-NULL; psi.Q<-NUL; psi.p<-NULL
miwp.w<-NULL; iwp.w<-NULL
Q<-NULL; A<-NULL

n <- 500 #number of patients in each simulation
s <-1000 #number of simulations 
T=1; seq.rand <- T #if sequentially randomised 

#objects for simulation of data
psi10<-250; psi11<--1 #CD4 at time point 1
psi20<-720; psi21<--2 #CD4 at time point 2

#objects for evaluation of data 
alpha0<-250; alpha1<--1 #CD4 at time point 1
alpha2<-720; alpha3<--2 #CD4 at time point 2

###############################Simulate data################################

expit <- function(x) { exp(x)/(1+exp(x)) }

## generation of simulated data
generate<-function() {
  
  #CD4 and decision at time point one
  L1<-abs(rnorm(n,450,150)) #CD4 at time point one 
  if(seq.rand) { #if following optimal regime at time point one
    A1<-rbinom(n, 1, expit(2-0.006*L1)) #1st decision is based on L1
  } else { #if NOT following optimal regime at time point one
    A1<-rbinom(n,1,0.5) #1st decision is 50%/50%
  }
  
  #CD4 and decision at time point two
  L2<-abs(1.25*L1+rnorm(n,0,60)) #CD4 at time point two (increasing with time)
  A2<-A1 #apply first decision to the second decision i.e. if started at time point 1, then continue treatment
  if(seq.rand) { #if following optimal regime at time point two
    #2nd decision is based on L2 for those who have not started treatment 
    A2[A2==0]<-rbinom(n,1,expit(0.8-0.004*L2))[A2==0] 
  } else { #if NOT following optimal regime at time point two
    A2[A2==0]<-rbinom(n,1,0.5)[A2==0] #2nd decision is 50%/50%
  }
  
  L<-cbind(L1,L2) #combine columns of L1 and L2 ??not used later on??
  A<-cbind(A1,A2) #combine columns of A1 and A2 ??not used later on??
  
  # equation on pg 689 (1st line) - gives reduction in CD4 if A1 and A2 are not according to optimal regime
  sum.mus<-abs(psi10+psi11*L1)*(A1-as.numeric(psi10+psi11*L1>0))^2+(1-A1)*abs(psi20+psi21*L2)*(A2-as.numeric(psi20+psi21*L2>0))^2
  Y.opt<-abs(400+1.6*L1+rnorm(n,0,60)) # optimal CD4 at time point two (increasing trend)
  Y<-Y.opt-sum.mus # actual CD4 according to actual regime (sum.mus = 0 if optimal regime applied)
  
  cbind(L1,A1,L2,A2,Y) #gives a matric of the simulated values 
}

## gives s number of matrices containing simulated data
data<-lapply(1:s,function(i) generate())

################################Evaluate data###############################

m<-1000000 #number of datapoints for evaluation
x0exp<-abs(rnorm(m,450,150)) #simulated values of baseline CD4 count
x1exp<-abs(rnorm(m,1.25*x0exp,60)) #simulated values of six-month CD4 count

##evaluate final outcome under a given treatment regime
evaluate<-function(eta) {
  
  x0<-x0exp
  x1<-x1exp
  
  #assign estimates to evaluate
  eta0<-eta[1]; eta1<-eta[2]; eta2<-eta[3]; eta3<-eta[4]
  
  #optimal regimes 
  g0<-as.numeric(I(eta0+eta1*x0>0)) 
  g1<-as.numeric(g0+(1-g0)*I(eta2+eta3*x1>0))
  
  #decisions
  a0<-as.numeric(I(alpha0+alpha1*x0>0))
  a1<-as.numeric(I(alpha2+alpha3*x1>0))
  
  #final CD4 at time point two 
  y<- 400 +1.6*x0-abs(alpha0+alpha1*x0)*(a0-g0)^2-(1-a0)*abs(alpha2+alpha3*x1)*(a1-g1)^2
  
  #gives the mean CD4 at time point 2 following the given optimal regime 
  mean(y)
}
evaluate(c(250,-1,360,-1)) #mean of the CD4 counts for all individuals at time point 2 (should come out ~1120 as stated in paper)
# 250 - optimal regime is specified as time point 1 (start if CD4 <250), time point 2 (start if CD4 <360)

obqrr<-function(eta) {
  
  eta0<-eta[1]
  eta1<-eta[2]
  
  #optimal regimes 
  g0<-as.numeric(I(x0<eta0)) #start if CD4 at time point one is less than eta0
  g1<-as.numeric(g0+(1-g0)*I(x1<eta1)) #start if CD4 at time point one is less than eta1
  
  c<-as.numeric(I(a0==g0)*I(a1==g1)) #label those who are on optimal regimes
  c1<-as.numeric(I(a0!=g0)) #label those who are not on optimal regimes at time point one
  c2<-as.numeric(I(a0==g0)*I(a1!=g1)) #label those who are not on optimal regimes at time point two
  
  lamda1<-(1-g0)*ph0+g0*(1-ph0)
  lamda2<-g1*(1-(g0+(1-g0)*ph1))+(1-g1)*(g0+(1-g0)*ph1)
  
  pc<-(1-lamda1)*(1-lamda2)
  
  ym0<-g0*m.1+(1-g0)*m.0
  ym1<-g0*(m.10)+(1-g0)*(g1*m.01+(1-g1)*m.00)
  
  mean(c/pc*y+(c1-lamda1)/(1-lamda1)*ym0+(c2-lamda2*I(a0==g0))/((1-lamda1)*(1-lamda2))*ym1)
}

obqrr1<-function(eta) {
  
  eta0<-eta[1]
  eta1<-eta[2]
  
  #optimal regimes 
  g0<-as.numeric(I(x0<eta0))#start if CD4 at time point one is less than eta0
  g1<-as.numeric(g0+(1-g0)*I(x1<eta1))#start if CD4 at time point one is less than eta1

  c<-as.numeric(I(a0==g0)*I(a1==g1)) #those who are on optimal regimes
  
  lamda1<-(1-g0)*ph0+g0*(1-ph0)
  lamda2<-g1*(1-(g0+(1-g0)*ph1))+(1-g1)*(g0+(1-g0)*ph1)
  
  pc<-(1-lamda1)*(1-lamda2)
  
  mean(c/pc*y)
}

for(i in 1:s) { # for each iteration in the simulation 
  
  #pull simulated data into respective variables
  L1<-data[[i]][,1]; L2<-data[[i]][,3] 
  A1<-data[[i]][,2]; A2<-data[[i]][,4]
  Y<-data[[i]][,5]

  #rename variables
  x0<-L1 #baseline CD4
  x1<-L2 #CD4 at 6 months 
  a0<-A1 # decision at baseline 
  a1<-A2 # decision at 6months 
  y<-Y  # CD4 at 12 month 
  
  a0x0<-a0*x0 
  a00a1<-(1-a0)*a1 #a1 only matters if a0!=1 (treatment not yet initiated)
  a00x1<-(1-a0)*x1 #x1 only matters if a0!=1 (treatment not yet initiated)
  a00a1x1<-(1-a0)*a1*x1
  
  #create datasets based on treatment regimes that were followed
  dataH<-data.frame(x0,a0,x1,a1,y) #complete dataset 
  data0<-dataH[dataH[,2]==0,] #dataset including individuals who did not initiate at time 1
  data1<-dataH[dataH[,2]==1,] #dataset including individuals who initiated at time 1
  data00<-data0[data0[,4]==0,] #dataset including individuals that did not initiate at time 1 or 2
  data01<-data0[data0[,4]==1,] #dataset including individuals who did not initiate until time 2
  
  ### q-learning #####################################################################
  
  fit00<-lm(y~x0+a0+a0x0+a00x1+a00a1+a00a1x1) #fit Q-learning function (Q2 on p.689)
  
  #save coefficients from model
  beta<-summary(fit00)$coef[,1]
  
  beta0.1<-beta[1] #intercept 
  beta1.1<-beta[2] #x0
  beta2.1<-beta[3] #a0
  beta3.1<-beta[4] #a0x0
  beta4.1<-beta[5] #a00x1
  beta5.1<-beta[6] #a00a1
  beta6.1<-beta[7] #a00a1x1
  
  #treatment regime for time 2 implied by model - CD4 threshold for 2nd time point (obtain G on p 684)
  etat1.Q<-c(beta5.1/abs(beta6.1),sign(beta6.1)) #?? why divide
  
  #expected final CD4 based on the above Q-learning function Q2 (p 684 obtain Q)
  m.00<-beta0.1+beta1.1*x0+beta4.1*x1         # when a0 = 0 (beta2.1, beta3.1 not considered) AND a00 = 1 (beta4.1 and beta5.1 considered) AND a1 = 0 (beta5.1, beta6.1 are not considered) 
  m.01<-m.00+beta5.1+beta6.1*x1               # when a0 = 0(a00 = 1) AND a1 = 1
  m.10<-beta0.1+beta1.1*x0+beta2.1+beta3.1*x0 # when a0 = 1(a00 = 0) AND a1 = 0
  m.11<-m.10                                  # when a0 = 1(a00 = 0) AND a1 = 1 

  # obtain optimal outcome (p. 683 section 2 obtain Vk for the next step of Q) 
  # CD4 at time point 2 - pick the best CD4 at time point 2 depending on if they initiated treatment at time point 1 
  y0<-y # create dummy variable for the outcome at time point 1
  # for those who did not initiate treatment at time point 1
  y0[a0==0]<-ifelse(m.01[a0==0]>m.00[a0==0], # if CD4 for a0=0/a1=1 > a0=0/a1=0
                    m.01[a0==0], # then outcome at time point 1 is those who initated later 
                    m.00[a0==0]) # else outcome at time point 1 is those who did not initiate 
  # for those who intiated treatment at time point 1
  y0[a0==1]<-ifelse(m.11[a0==1]>m.10[a0==1],# if CD4 for a0=1/a1=1 > a0=1/a1=0
                    m.11[a0==1], # then outcome at time point 1 is those who initated later
                    m.10[a0==1]) # then outcome at time point 1 is those who did not initiate 
  
  #using optimised CD4 at time point 2
  fit0<-lm(y0~x0+a0+a0:x0) #fit Q-learning function (Q1 p.689)
  
  #save coefficients from model
  beta<-summary(fit0)$coef[,1]
  
  beta0.0<-beta[1] #intercept 
  beta1.0<-beta[2] #x0
  beta2.0<-beta[3] #a0
  beta3.0<-beta[4] #x0a0
  
  #estimate CD4 at time point 1 based on the above Q-learning function (Q1)
  m.0<-beta0.0+beta1.0*x0 #a0 = 0
  m.1<-beta0.0+beta1.0*x0+beta2.0+beta3.0*x0 #a0 = 1
  
  #treatment regime for time 1 implied by model
  etat0.Q<-c(beta2.0/abs(beta3.0),sign(beta3.0))
  
  eta.Q<-c(etat0.Q, etat1.Q) #combine the estimated CD4 cutoffs for time 1 and time 2 
  
  #expected value of final CD4 given our optimal regime, from our training set (Ehat in the results table in the paper)
  hatQ<-mean(ifelse(m.1>m.0, #if CD4 at time 1 for those who initiated treatment > those who did not 
                    m.1, #take the CD4 for those who initiated
                    m.0)) #else take the CD4 for those who did not
  
  #evaluate final outcome using the eta we derived from the Q learning - for validation (E in the results table in the paper)
  expY<-evaluate(eta.Q) #output is mean of the final CD4, the correct answer 
  
  summary<-c(eta.Q, hatQ, expY) 
  
  Q<-rbind(Q,t(summary)) 

  ### a-learning #####################################################################
  
  # h1 on pg. 9 
  logit1<-glm(a1~x1,family=binomial,data=data0, epsilon=1e-14)
  
  gamma<-summary(logit1)$coef[,1]
  gamma0h.1<-gamma[1] #intercept - outcome if did not initiate at time point 1
  gamma1h.1<-gamma[2] #x1 
  
  #propensity score for timepoint 2? (see bottom of pg. 9)
  ph1<-exp(gamma0h.1+gamma1h.1*x1)/(1+exp(gamma0h.1+gamma1h.1*x1)) #expit - the probabilities of starting treatmemt at time point 2  
  
  aa<-a00a1-((1-a0)*ph1) # a00a1<-(1-a0)*a1 
                         # when a0 = 1, aa = 0 
                         # when a0 = 0, a1 = 1, aa = 1 - ph1
                         # when a0 = 0, a1 = 0, aa = - ph1
                         # ?? what is aa? 
 
  x1aa<-aa*x1 

  Z2<-cbind(1,x0,a0,a0x0,a00x1,a00a1,a00a1x1) #is this the design matrix for h2 and C2 together?                                                                   
  Za1<-cbind(1,x0,a0,a0x0,a00x1,aa,x1aa)                                                      
  lZa1<-ncol(Za1)
  
  b2ahat<-as.vector(solve(t(Za1)%*%Z2)%*%t(Za1)%*%y)
  
  psi.temp<-b2ahat[6:7]
  psi.temp.1<-psi.temp
  
  eta10<-psi.temp[1]
  eta11<-psi.temp[2]
  
  eta10<-eta10/abs(eta11)
  eta11<-sign(eta11)
  etat1<-c(eta10,eta11)

  d.opt.1<-as.numeric(cbind((1-a0),a00x1) %*% psi.temp.1 >= 0)
  
  d.opt.minus.A <- d.opt.1 - a1
  Y.dash <- y + d.opt.minus.A*(cbind((1-a0),a00x1) %*% psi.temp.1)
  
  logit0<-glm(a0~x0,family=binomial, epsilon=1e-14)
  ph0<-logit0$fit
  
  aa<-a0-ph0
  a0x0<-a0*x0
  x0aa<-aa*x0
  
  Z2<-cbind(1,x0,a0,a0x0)                                                                    
  Za1<-cbind(1,x0,aa,x0aa)                                                       
  lZa1<-ncol(Za1)
  
  b2ahat<-as.vector(solve(t(Za1)%*%Z2)%*%t(Za1)%*%Y.dash)
  
  psi.temp<-b2ahat[3:4]
  psi.temp.0<-psi.temp
  
  psi<-c(psi.temp.0,psi.temp.1)
  
  eta10<-psi.temp[1]
  eta11<-psi.temp[2]
  eta10<-eta10/abs(eta11)
  eta11<-sign(eta11)
  etat0<-c(eta10,eta11)
  eta<-c(etat0,etat1)
  
  expY<-evaluate(eta)
  
  d.opt.0<-as.numeric(cbind(1,x0) %*% psi.temp.0 > 0)
  d.opt.minus.A <- d.opt.0 - a0
  hatQ<-mean(Y.dash+ d.opt.minus.A*(cbind(1,x0) %*% psi.temp.0))
  
  summary<-c(eta,hatQ,expY)
  
  A<-rbind(A,t(summary))
  
  ###############################################################################
  
  
  logit1<-glm(a1~x1,family=binomial,data=data0, epsilon=1e-14)
  
  gamma<-summary(logit1)$coef[,1]
  gamma0h.1<-gamma[1] #intercept
  gamma1h.1<-gamma[2] #x1
  
  logit0<-glm(a0~x0,family=binomial, epsilon=1e-14)
  
  gamma<-summary(logit0)$coef[,1]
  gamma0h.0<-gamma[1] #intercept
  gamma1h.0<-gamma[2] #x0
  
  ph1<-exp(gamma0h.1+gamma1h.1*x1)/(1+exp(gamma0h.1+gamma1h.1*x1)) # probability of initiating at time 2
  
  ph0<-exp(gamma0h.0+gamma1h.0*x0)/(1+exp(gamma0h.0+gamma1h.0*x0)) # probability of initiating at time 1
  
  #######################################################################
  
  etaH0<-sort(x0)[1:200]
  etaH1<-sort(x1)[1:200]
  # we fixed (η21,η11)=(−1,−1) and determined η20,η10 via a grid search; because IPWE(η), DR(η) 
  # and AIPWE(η) are step functions of η with jumps at (x1i , x2 j ) (i, j =1,...,n), 
  # we maximized in η over all (x1i,x2j)
  etaH<-expand.grid(etaH0,etaH1) 
  
  #########################################################
  idx <- which.max(apply(etaH, 1, obqrr))
  
  eta<-c(etaH[idx,][[1]],etaH[idx,][[2]])
  
  hatQ<-obqrr(eta)
  
  etat0<-c(eta[1],-1)
  etat1<-c(eta[2],-1)
  eta<-c(etat0,etat1)
  
  expY<-evaluate(eta)
  
  sd<-sdmiwp(eta) #get standard deviation 
  
  ci.l<-hatQ-1.96*sd #get lower confidence interval 
  ci.h<-hatQ+1.96*sd #get higher confidence interval 
  ci.t<-as.numeric(I(1120>ci.l)*I(1120<ci.h)) 
  
  summary<-c(eta,hatQ,expY,sd,ci.t)
  
  miwp.w<-rbind(miwp.w,t(summary))
  
  etaH<-expand.grid(etaH0,etaH1)
  ################################################
  
  idx <- which.max(apply(etaH, 1, obqrr1))
  
  eta<-c(etaH[idx,][[1]],etaH[idx,][[2]])

  hatQ<-obqrr1(eta)
  
  etat0<-c(eta[1],-1)
  etat1<-c(eta[2],-1)
  eta<-c(etat0,etat1)
  
  expY<-evaluate(eta)
  
  sd<-sdiwp(eta)
  
  ci.l<-hatQ-1.96*sd
  ci.h<-hatQ+1.96*sd
  ci.t<-as.numeric(I(1120>ci.l)*I(1120<ci.h))
  
  summary<-c(eta,hatQ,expY,sd,ci.t)
  
  iwp.w<-rbind(iwp.w,t(summary))

}

colMeans(A)
sd(A)

colMeans(Q)
sd(Q)

colMeans(miwp.w)
sd(miwp.w)

colMeans(iwp.w)
sd(iwp.w)

md<-miwp.w
round(c(mean(md[,1]),sd(md[,1]),mean(md[,2]),sd(md[,2]),mean(md[,3]),sd(md[,3]),mean(md[,4]),sd(md[,4]),mean(md[,5]),sd(md[,5]),mean(md[,6]),sd(md[,6]),mean(md[,7]),sd(md[,7]),mean(md[,8]),sd(md[,8]),mean(md[,9]),sd(md[,9]),mean(md[,11]),mean(md[,12]),mean(md[,10]),sd(md[,10])),2)

md<-Q
round(c(mean(md[,1]),sd(md[,1]),mean(md[,2]),sd(md[,2]),mean(md[,3]),sd(md[,3]),mean(md[,4]),sd(md[,4]),mean(md[,5]),sd(md[,5]),mean(md[,6]),sd(md[,6]),mean(md[,7]),sd(md[,7]),mean(md[,8]),sd(md[,8]),mean(md[,9]),sd(md[,9]),mean(md[,10]),sd(md[,10])),2)

md<-A
round(c(mean(md[,1]),sd(md[,1]),mean(md[,2]),sd(md[,2]),mean(md[,3]),sd(md[,3]),mean(md[,4]),sd(md[,4]),mean(md[,5]),sd(md[,5]),mean(md[,6]),sd(md[,6]),mean(md[,7]),sd(md[,7]),mean(md[,8]),sd(md[,8]),mean(md[,9]),sd(md[,9]),mean(md[,10]),sd(md[,10])),2)




