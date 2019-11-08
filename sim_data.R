########Simulation of observational data for treatment duration of VAP#########
###############################################################################

# Based on the simulation setting used in Zhang et al. 2013 https://academic.oup.com/biomet/article-abstract/100/3/681/303040, Section 5

expit <- function(x) { exp(x) / (1 + exp(x)) } #inverse logit

set.seed(1234)

#Parameters
n = 500      #number of patients

#.......................................................#
#....................Data generation....................#
#.......................................................#

##Baseline
A0 = rep(1, n)                                      # Start treatment for everyone
L =  rnorm(n, mean = 10, sd = 5)                   # Baseline characteristics: comobidities score
C0 = rnorm(n, mean = 20, sd = 2) + 1.5 * L        # Time varying covariate: disease severity score

##Time varying, t=1 
# C1 depends on L, C0, A0
mu_C1 = C0 - 1 * A0 + 2 * L 
C1 = rnorm(n, mu_C1, 4)

# P(A1) depends on L, A0, C1
P_A1 = expit(0.3 + 0.003 * L  + 0.005 * C1)
A1 = rbinom(n , 1, P_A1)

##Time varying, t=2
# C2 depends on L, C1
mu_C2 = C1 - 1.5 * A1 + 2 * L 
C2 = rnorm(n, mu_C2, 5)

# P(A2) depends on L, C2
P_A2 = expit(0.002 * L + 0.004 * C2)
A2 = rbinom(n , 1, P_A2)
A2[which(A1 == 0)] = 0 #those who had antibiotics stopped at t=1 should have 0 in t=2

##Regime = benefits survival if at t=1, L>4 C1>6 continues antibiotics and 
#                               at t=2, L>4 C2>8 continues antibiotics

# R = rep(0, n)
# R[intersect(which(L> 0.4 & C1> 0.6 & A1 ==1), which(L> 0.4 & C2> 0.8 & A2 ==1))] = 1

##Outcome
# P_Y = expit(0.5 * L + 1 * R + 0.1 * C0 + 0.1 * C1 + 0.2 * C2) # A2 has no effect on outcome
# Y = rbinom(n, 1, P_Y)

## Note: Can the outcome be a continuous variable? Having a discrete variable complicates analysis. 
## Have to consider survival models otherwise.

## Note
#  Objective: Minimize Y and minimize time on treatment

# Assuming continuous Y, say, disease severity score at t=10
# Y is generated such that optimal regime is the following:
# If C1<=55 -> Stop treatment i.e. A1=0,A2=0,
# If C1>55 -> keep treatment i.e. A1=1
# If C2<=75 -> Stop treatment i.e. A1=1, A2=0
# If C2>75 -> Keep treatment i.e. A1=1, A2=1
# Y depends on C1, A1, C2, A2
# Y = 20 + 6*C1 + (1-A1)*|C1 - 50|*{1(C1 - 50 > 0) - A1}^2 + A1*(1-A2)*|C2 - 70|*{1(C2 - 70 > 0) - A2}^2\
#     + 5*A1 + 5*A1*A2
mu_Y = 20 + 6 * C1 + (1 - A1) * abs(C1-50) * (ifelse((C1-50)>0,1,0) - A1)**2 + A1 * (1 - A2) * abs(C2-70) * (ifelse(C2-70>0,1,0) - A2)**2 + 5*A1 + 5*A1*A2
Y = rnorm(n, mu_Y, 0.1)

# Optimal regime
# g1 = 1(C1 - 55 > 0), g2 = 1(-1 + A1 + A1*(C2-75) > 0)

##Data
d = cbind.data.frame(id = 1:n, #ID 
                     L = L,    #Baseline covariates
                     C0 = C0,  #Time varying covariates
                     C1 = C1, 
                     C2 = C2, 
                     A0 = A0,  #Treatment decisions
                     A1 = A1, 
                     A2 = A2, 
                     Y = Y)  #Outcome

d$A0=ifelse(d$A0==1, 'cont', 'stop')
d$A1=ifelse(d$A1==1, 'cont', 'stop')
d$A2=ifelse(d$A2==1, 'cont', 'stop')

if (length(which(A1==1 & A2 ==1)) + length(which(A1==1 & A2 ==0)) + length(which(A1==0 & A2 ==0)) ==n) {
  print(paste(length(which(A1==1 & A2 ==1)), 'out of' , n, 'had regime A0=1, A1=1, A2=1'))
  print(paste(length(which(A1==1 & A2 ==0)), 'out of' , n, 'had regime A0=1, A1=1, A2=0'))
  print(paste(length(which(A1==0 & A2 ==0)), 'out of' , n, 'had regime A0=1, A1=0, A2=0'))
} else {
  print ('Some patients have regimes that are not A1=0/A2=0, A1=1/A2=0, A1=1/A2=1')
}

