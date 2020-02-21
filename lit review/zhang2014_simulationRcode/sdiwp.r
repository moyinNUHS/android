
sdiwp<-function(eta)

{


eta0<-eta[1]
eta1<-eta[2]
eta2<-eta[3]
eta3<-eta[4]



g0<-as.numeric(I(eta0+eta1*x0>0))

g1<-as.numeric(g0+(1-g0)*I(eta2+eta3*x1>0))



c<-as.numeric(I(a0==g0)*I(a1==g1))

c1<-as.numeric(I(a0!=g0))
c2<-as.numeric(I(a0==g0)*I(a1!=g1))


lamda1<-(1-g0)*ph0+g0*(1-ph0)

lamda2<-g1*(1-(g0+(1-g0)*ph1))+(1-g1)*(g0+(1-g0)*ph1)


pc<-(1-lamda1)*(1-lamda2)



ym<-c/pc*y


expYh3<-mean(ym)




pp0<-exp(gamma0h.0+gamma1h.0*x0)

s0<-pp0/(1+pp0)^2

pp1<-exp(gamma0h.1+gamma1h.1*x1)

s1<-pp1/(1+pp1)^2




d1<-1


d2<-mean((c/(1-lamda1)*y)/(1-lamda2)^2*(1-g0)*(2*g1-1)*s1)
d3<-mean((c/(1-lamda1)*y)/(1-lamda2)^2*(1-g0)*(2*g1-1)*s1*x1)

d4<-mean((c/(1-lamda2)*y)/(1-lamda1)^2*(2*g0-1)*s0)

d5<-mean((c/(1-lamda2)*y)/(1-lamda1)^2*(2*g0-1)*s0*x0)




r1<-t(c(d1,d2,d3,d4,d5))




d1<-0

d2<-mean(s1*I(a0==0))

d3<-mean(s1*x1*I(a0==0))

d4<-0

d5<-0

r2<-t(c(d1,d2,d3,d4,d5))




d1<-0

d2<-mean(s1*x1*I(a0==0))

d3<-mean(s1*x1*x1*I(a0==0))


d4<-0

d5<-0


r3<-t(c(d1,d2,d3,d4,d5))





d1<-0



d2<-0

d3<-0


d4<-mean(s0)

d5<-mean(s0*x0)



r4<-t(c(d1,d2,d3,d4,d5))


d1<-0

d2<-0

d3<-0

d4<-mean(s0*x0)

d5<-mean(s0*x0*x0)



r5<-t(c(d1,d2,d3,d4,d5))






h1<-rbind(r1,r2,r3,r4,r5)






p1<-ym-expYh3

p2<-(a1-ph1)*I(a0==0)

p3<-(a1-ph1)*x1*I(a0==0)


p4<-a0-ph0

p5<-(a0-ph0)*x0



p<-cbind(p1,p2,p3,p4,p5)

h2<-(t(p)%*%p)/n







invh1<-solve(h1)

sigmatrix<-invh1%*%h2%*%t(invh1)

varE<-sigmatrix[1,1]

sd3<-sqrt(varE/n)







return(sd3)




}