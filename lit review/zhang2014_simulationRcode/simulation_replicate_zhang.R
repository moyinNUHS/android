set.seed(123456)

## generation of simulated data based on optimal regime
generate<-function(n, seq.rand = 1, psi10 = 250, psi20 = 720, psi11=-1, psi21 = -1) {
  
  expit <- function(x) { exp(x)/(1+exp(x)) }
  
  #CD4 and decision at time point one
  L1<-abs(rnorm(n, 450, 150)) #CD4 at time point one 
  if(seq.rand) { #if following optimal regime at time point one
    A1<-rbinom(n, 1, expit(2-0.006*L1)) #1st decision is based on L1
  } else { #if NOT following optimal regime at time point one
    A1<-rbinom(n,1,0.5) #1st decision is 50%/50%
  }
  
  #CD4 and decision at time point two
  L2<-abs(1.25*L1 + rnorm(n,0,60)) #CD4 at time point two (increasing with time)
  A2<-A1 #apply first decision to the second decision i.e. if started at time point 1, then continue treatment
  if(seq.rand) { #if following optimal regime at time point two
    #2nd decision is based on L2 for those who have not started treatment 
    A2[A2==0]<-rbinom(n, 1 ,expit(0.8-0.004*L2))[A2==0] 
  } else { #if NOT following optimal regime at time point two
    A2[A2==0]<-rbinom(n, 1, 0.5)[A2==0] #2nd decision is 50%/50%
  }
  
  # equation on pg 689 (1st line) - gives reduction in CD4 if A1 and A2 are not according to optimal regime
  sum.mus<-abs(psi10+psi11*L1)*(A1-as.numeric(psi10+psi11*L1>0))^2+(1-A1)*abs(psi20+psi21*L2)*(A2-as.numeric(psi20+psi21*L2>0))^2
  Y.opt<-abs(400+1.6*L1+rnorm(n,0,60)) # optimal CD4 at time point two (increasing trend)
  Y<-Y.opt-sum.mus # actual CD4 according to actual regime (sum.mus = 0 if optimal regime applied)
  
  cbind(L1,A1,L2,A2,Y) #gives a matrix of the simulated values 
}


