
sdmiwp<-function(eta)

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

ym0<-g0*m.1+(1-g0)*m.0

ym1<-g0*(g1*m.11+(1-g1)*m.10)+(1-g0)*(g1*m.01+(1-g1)*m.00)


ym<-c/pc*y+(c1-lamda1)/(1-lamda1)*ym0+(c2-lamda2*I(a0==g0))/pc*ym1


expYh3<-mean(ym)



y0<-beta0.1+beta1.1*x0+beta2.1*a0+beta3.1*a0*x0+beta4.1*(1-a0)*x1+  (1-a0)*(beta5.1+beta6.1*x1)*as.numeric(I(beta5.1+beta6.1*x1>0))

mu1<-beta0.1+beta1.1*x0+beta2.1*a0+beta3.1*a0*x0+beta4.1*(1-a0)*x1+beta5.1*(1-a0)*a1+beta6.1*(1-a0)*a1*x1
mu0<-beta0.0+beta1.0*x0+beta2.0*a0+beta3.0*a0*x0


pp0<-exp(gamma0h.0+gamma1h.0*x0)

s0<-pp0/(1+pp0)^2

pp1<-exp(gamma0h.1+gamma1h.1*x1)

s1<-pp1/(1+pp1)^2




d1<-1


d2<-mean((c/(1-lamda1)*y+(c2-I(a0==g0))/(1-lamda1)*ym1)/(1-lamda2)^2*(1-g0)*(2*g1-1)*s1)
d3<-mean((c/(1-lamda1)*y+(c2-I(a0==g0))/(1-lamda1)*ym1)/(1-lamda2)^2*(1-g0)*(2*g1-1)*s1*x1)

d4<-mean((c/(1-lamda2)*y+(c1-1)*ym0+(c2-lamda2*I(a0==g0))/(1-lamda2)*ym1)/(1-lamda1)^2*(2*g0-1)*s0)

d5<-mean((c/(1-lamda2)*y+(c1-1)*ym0+(c2-lamda2*I(a0==g0))/(1-lamda2)*ym1)/(1-lamda1)^2*(2*g0-1)*s0*x0)



d6<--mean((c2-lamda2*I(a0==g0))/pc)

d7<--mean((c2-lamda2*I(a0==g0))/pc*x0)

d8<--mean((c2-lamda2*I(a0==g0))/pc*g0)

d9<--mean((c2-lamda2*I(a0==g0))/pc*x0*g0)

d10<--mean((c2-lamda2*I(a0==g0))/pc*x1*(1-g0))

d11<--mean((c2-lamda2*I(a0==g0))/pc*g1*(1-g0))
d12<--mean((c2-lamda2*I(a0==g0))/pc*g1*x1*(1-g0))



d13<--mean((c1-lamda1)/(1-lamda1))

d14<-mean(-(c1-lamda1)/(1-lamda1)*x0)

d15<-mean(-(c1-lamda1)/(1-lamda1)*g0)

d16<-mean(-(c1-lamda1)/(1-lamda1)*g0*x0)



r1<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))




d1<-0

d2<-mean(s1*I(a0==0))

d3<-mean(s1*x1*I(a0==0))

d4<-0

d5<-0

d6<-0

d7<-0

d8<-0

d9<-0

d10<-0

d11<-0

d12<-0

d13<-0

d14<-0

d15<-0

d16<-0


r2<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))




d1<-0

d2<-mean(s1*x1*I(a0==0))

d3<-mean(s1*x1*x1*I(a0==0))


d4<-0

d5<-0

d6<-0

d7<-0

d8<-0

d9<-0

d10<-0

d11<-0

d12<-0

d13<-0

d14<-0

d15<-0

d16<-0


r3<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))





d1<-0



d2<-0

d3<-0


d4<-mean(s0)

d5<-mean(s0*x0)

d6<-0

d7<-0

d8<-0

d9<-0

d10<-0

d11<-0

d12<-0

d13<-0

d14<-0

d15<-0

d16<-0


r4<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))




d1<-0

d2<-0

d3<-0

d4<-mean(s0*x0)

d5<-mean(s0*x0*x0)

d6<-0

d7<-0

d8<-0

d9<-0

d10<-0

d11<-0

d12<-0

d13<-0

d14<-0

d15<-0

d16<-0


r5<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))





d1<-0

d2<-0

d3<-0

d4<-0

d5<-0

d6<-mean(1)

d7<-mean(1*x0)

d8<-mean(1*a0)

d9<-mean(1*a0*x0)
d10<-mean(1*(1-a0)*x1)

d11<-mean(1*a1*(1-a0))


d12<-mean(1*a1*x1*(1-a0))


d13<-0

d14<-0

d15<-0

d16<-0


r6<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))





d1<-0

d2<-0

d3<-0

d4<-0

d5<-0

d6<-mean(x0*1)

d7<-mean(x0*x0)

d8<-mean(x0*a0)

d9<-mean(x0*a0*x0)

d10<-mean(x0*(1-a0)*x1)


d11<-mean(x0*(1-a0)*a1)


d12<-mean(x0*(1-a0)*a1*x1)


d13<-0

d14<-0

d15<-0

d16<-0


r7<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))




d1<-0

d2<-0

d3<-0

d4<-0

d5<-0

d6<-mean(a0*1)

d7<-mean(a0*x0)

d8<-mean(a0*a0)

d9<-mean(a0*a0*x0)

d10<-mean(a0*(1-a0)*x1)


d11<-mean(a0*(1-a0)*a1)


d12<-mean(a0*(1-a0)*a1*x1)



d13<-0

d14<-0

d15<-0

d16<-0


r8<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))




d1<-0

d2<-0

d3<-0

d4<-0

d5<-0

d6<-mean(a0*x0*1)

d7<-mean(a0*x0*x0)

d8<-mean(a0*x0*a0)

d9<-mean(a0*x0*a0*x0)

d10<-mean(a0*x0*(1-a0)*x1)


d11<-mean(a0*x0*(1-a0)*a1)


d12<-mean(a0*x0*(1-a0)*a1*x1)



d13<-0

d14<-0

d15<-0

d16<-0


r9<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))


d1<-0

d2<-0

d3<-0

d4<-0

d5<-0

d6<-mean((1-a0)*x1*1)

d7<-mean((1-a0)*x1*x0)

d8<-mean((1-a0)*x1*a0)

d9<-mean((1-a0)*x1*a0*x0)

d10<-mean((1-a0)*x1*(1-a0)*x1)


d11<-mean((1-a0)*x1*(1-a0)*a1)


d12<-mean((1-a0)*x1*(1-a0)*a1*x1)



d13<-0

d14<-0

d15<-0

d16<-0


r10<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))



d1<-0

d2<-0

d3<-0

d4<-0

d5<-0


d6<-mean((1-a0)*a1*1)

d7<-mean((1-a0)*a1*x0)

d8<-mean((1-a0)*a1*a0)

d9<-mean((1-a0)*a1*a0*x0)

d10<-mean((1-a0)*a1*(1-a0)*x1)


d11<-mean((1-a0)*a1*(1-a0)*a1)


d12<-mean((1-a0)*a1*(1-a0)*a1*x1)
d13<-0

d14<-0

d15<-0

d16<-0


r11<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))



d1<-0

d2<-0

d3<-0

d4<-0

d5<-0


d6<-mean((1-a0)*a1*x1*1)

d7<-mean((1-a0)*a1*x1*x0)

d8<-mean((1-a0)*a1*x1*a0)

d9<-mean((1-a0)*a1*x1*a0*x0)

d10<-mean((1-a0)*a1*x1*(1-a0)*x1)


d11<-mean((1-a0)*a1*x1*(1-a0)*a1)


d12<-mean((1-a0)*a1*x1*(1-a0)*a1*x1)

d13<-0

d14<-0

d15<-0

d16<-0


r12<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))





d1<-0

d2<-0

d3<-0

d4<-0

d5<-0

d6<-mean(-1)

d7<-mean(-x0 )

d8<-mean(-a0)

d9<-mean(-a0*x0)

d10<-mean(-(1-a0)*x1)


d11<-mean(-I(beta5.1+beta6.1*x1>0)*(1-a0))

d12<-mean(-x1*I(beta5.1+beta6.1*x1>0)*(1-a0))


d13<-1

d14<-mean(1*x0)

d15<-mean(1*a0)

d16<-mean(1*a0*x0)


r13<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))




d1<-0

d2<-0

d3<-0

d4<-0

d5<-0

d6<-mean(-1*x0)

d7<-mean(-x0*x0 )

d8<-mean(-a0*x0)

d9<-mean(-a0*x0*x0)

d10<-mean(-(1-a0)*x1*x0)


d11<-mean(-I(beta5.1+beta6.1*x1>0)*(1-a0)*x0)

d12<-mean(-x1*I(beta5.1+beta6.1*x1>0)*(1-a0)*x0)



d13<-mean(x0)

d14<-mean(x0*x0)

d15<-mean(x0*a0)

d16<-mean(x0*a0*x0)


r14<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))





d1<-0

d2<-0

d3<-0

d4<-0

d5<-0



d6<-mean(-1*a0)

d7<-mean(-x0*a0 )

d8<-mean(-a0*a0)

d9<-mean(-a0*x0*a0)

d10<-mean(-(1-a0)*x1*a0)


d11<-mean(-I(beta5.1+beta6.1*x1>0)*(1-a0)*a0)

d12<-mean(-x1*I(beta5.1+beta6.1*x1>0)*(1-a0)*a0)




d13<-mean(a0)

d14<-mean(a0*x0)

d15<-mean(a0*a0)

d16<-mean(a0*a0*x0)


r15<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))



d1<-0

d2<-0

d3<-0

d4<-0

d5<-0



d6<-mean(-1*x0*a0)

d7<-mean(-x0*x0*a0)

d8<-mean(-a0*x0*a0)

d9<-mean(-a0*x0*x0*a0)

d10<-mean(-(1-a0)*x1*x0*a0)


d11<-mean(-I(beta5.1+beta6.1*x1>0)*(1-a0)*x0*a0)

d12<-mean(-x1*I(beta5.1+beta6.1*x1>0)*(1-a0)*x0*a0)


d13<-mean(a0*x0)

d14<-mean(a0*x0*x0)

d15<-mean(a0*x0*a0)

d16<-mean(a0*x0*a0*x0)


r16<-t(c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16))





h1<-rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16)






p1<-ym-expYh3

p2<-(a1-ph1)*I(a0==0)

p3<-(a1-ph1)*x1*I(a0==0)


p4<-a0-ph0

p5<-(a0-ph0)*x0

p6<-(y-mu1)

p7<-(y-mu1)*x0

p8<-(y-mu1)*a0

p9<-(y-mu1)*a0*x0

p10<-(y-mu1)*(1-a0)*x1

p11<-(y-mu1)*(1-a0)*a1

p12<-(y-mu1)*(1-a0)*a1*x1


p13<-(y0-mu0)

p14<-(y0-mu0)*x0

p15<-(y0-mu0)*a0

p16<-(y0-mu0)*a0*x0




p<-cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16)

h2<-(t(p)%*%p)/n







invh1<-solve(h1)

sigmatrix<-invh1%*%h2%*%t(invh1)

varE<-sigmatrix[1,1]

sd3<-sqrt(varE/n)







return(sd3)




}