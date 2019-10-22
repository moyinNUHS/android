########Simulation of observational data for treatment duration of VAP#########
###############################################################################

expit <- function(x) { exp(x) / (1 + exp(x)) } #inverse logit
set.seed(1234)

#Parameters
n = 500      #number of patients

#.......................................................#
#....................Data generation....................#
#.......................................................#

##Baseline
A0 = rep(1, n)
L =  rbeta(n, shape1 = 2, shape2 = 2)                  # baseline characteristics: comobidities score
C0 = rbeta(n, shape1 = 3, shape2 = 2) + 0.15 * L       # Time varying covariate: disease severity score

##Time varying, t=1 
C1 = C0 - 0.1 + 0.2 * L - 0.15 * A0

P_A1 = expit(0.3 + 0.7 * L + 0.7 * C1)
A1 = rbinom(n , 1, P_A1)

##Time varying, t=2
C2 = C1 - 0.05 + 0.2 * L - 0 * A1 #does not depend on previous A1

P_A2 = expit(0.2 * L + 0.4 * C2 + 0.2 * A1)
A2 = rbinom(n , 1, P_A2)
A2[which(A1 == 0)] = 0 #those who had antibiotics stopped at t=1 should have 0 in t=2

##Regime = benefits survival if at t=1, L>0.4 C1>0.6 continues antibiotics and 
#                               at t=2, L>0.4 C2>0.8 continues antibiotics

R = rep(0, n)
R[intersect(which(L> 0.4 & C1> 0.6 & A1 ==1), which(L> 0.4 & C2> 0.8 & A2 ==1))] = 1

##Outcome
P_Y = expit(0.5 * L + 1 * R + 0.1 * C0 + 0.1 * C1 + 0.2 * C2) # A2 has no effect on outcome
Y = rbinom(n, 1, P_Y)

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

mean(Y[which(A1==1 & A2 ==1)]) #mortality in the regime 111
mean(Y[which(A1==1 & A2 ==0)]) #mortality in the regime 110
mean(Y[which(A1==0 & A2 ==0)]) #mortality in the regime 100
