# Solve backward Fokker-Planck for May model with Babak's pseudocode Dec 28 2019 email
# SRC 2019-12-29 + Weighted Mean 2020-01-10 + cleanup of code & notes
# see comparison of integration methods 2020-01-10; they get the same answer
# (c) S. R. Carpenter and M.S. Arani, February 2020

rm(list = ls())
graphics.off()

library(bvpSolve)

# Parameters of May model slightly modified from Carpenter et al. PNAS 2015
r =  0.75  # 0.75 is nominal
c = 1.45  # 1.35 is nominal
b = 1
K = 8.8   # 8.8 is nominal
h = 0.01 # interesting range is 0.01 to 0.09
CV = 0.1  # sigma/x when x=K
sigma = CV*K

# list for bvpSolve; if eps is specified it goes first
eps = 0.5 # documentation says 0.1 < eps < 1
plist = list(eps,r,c,b,K,h,sigma,CV)

D1 = function(x) {  # drift
  dx = r*x*(1 - (x/K)) - c*x*x/(b*b + x*x) - h*x 
  return(dx)
}

D2 = function(x) {  # diffusion
  s = sigma*x/K # note CV = sigma/K, so s = (sigma/K)*x = CV*x
  return(s)
}

# plot drift and diffusion
xvec = seq(0.01,K,length.out=100)
drift = D1(xvec)
diff = D2(xvec)

# Plot potential function
PF = rep(0,99) 
for(i in 1:99) {
  x0 = xvec[1]
  x1 = xvec[i+1]
  PF[i] = integrate(D1,lower=x0,upper=x1)$value
}
negPF = -1*PF

windows(width=7,height=10)
par(mfrow=c(3,1),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(xvec,drift,type='l',lwd=2,col='blue',xlab='x',ylab='drift')
abline(h=0,lty=2)
grid()
plot(xvec[2:100],negPF,type='l',lwd=2,col='blue',xlab='x',ylab='negative potential')
grid()
plot(xvec,diff,type='l',lwd=2,col='red',xlab='x',ylab='diffusion')

# find equilibria
xeq = rep(0,3)
print('equilibria',quote=F)
#
xroot = uniroot(D1,interval=c(0.1,1),maxiter = 100)
xeq[1] = xroot$root
print(xeq[1])
#
xroot = uniroot(D1,interval=c(1.5,3),maxiter = 100)
xeq[2] = xroot$root
print(xeq[2])
#
xroot = uniroot(D1,interval=c(4.5,7),maxiter = 100)
xeq[3] = xroot$root
print(xeq[3])

# function for solving the boundary value problem as two differential equations 
#  for T (col 1) and dT/dx (col 2)
feval2 = function(x,y,plist) {
  out1 = y[2]
  out2 = -(D1(x)*y[2]+1)/D2(x)
  return( list(c(out1,out2)) )
}

# Calculate Exit Times =======================================================

# Left basin

# set up for bvpSolve
yini = c(NA,0)
yend = c(0,NA)

# solve the left basin from x = 0 (reflecting) to x=xeq[2] (absorbing)
x = seq(0.1,xeq[2],length.out=20)  # x vector

# solve with bvpcol
trycol <- bvpcol(yini = yini, x = x, func = feval2, yend = yend, parm=plist)

windows()
par(mfrow=c(1,1),mar=c(3,4,3,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(trycol[,1],trycol[,2],type='l',lty=1,lwd=2,col='blue',xlab='x',
     ylab='Exit Time',main='left basin')
abline(v=xeq[1],lty=1,col='magenta')
abline(v=xeq[2], lty=2, col='red')

# save solution for further analysis
ETL = trycol # left exit time

# end of solution for basin 1

# Right basin 

  # set up for bvpSolve
  yini = c(0, NA)
  yend = c(NA, 0)
  
  # right basin from x=xeq[2] (absorbing) to x > xeq[3] (reflecting)
  x = seq(xeq[2],xeq[3]+0.1,length.out=20)  # x vector
  
  # solve with bvpcol
  trycol <- bvpcol(yini = yini, x = x, func = feval2, yend = yend, parm=plist)
  
  # save solution for further analysis
  ETR = trycol # right exit time
  
  windows()
  par(mfrow=c(1,1),mar=c(3,4,3,2)+0.1,cex.lab=1.5,cex.axis=1.5)
  plot(trycol[,1],trycol[,2],type='l',lty=1,lwd=2,col='blue',xlab='x',
       ylab='Exit Time',main = 'right basin')
  abline(v=xeq[2], lty=2, col='red')
  abline(v=xeq[3],lty=1,col='magenta')
  
# end of solution for basin 2
  
# CALCULATE WEIGHTS FOR AVERAGE EXIT TIME -----------------------------------------------------  
# Weight of x is p(x) from stationary distribution of Langevin equation  
# Based on appendix of Carpenter & Brock, Ecology Letters 2006 based on the book 
# 'Noise-Induced Transitions' by Horsthemke and Lefever 1984
  
# function for inner integral
finner = function (z) { 
  fx = D1(z)/(D2(z)*D2(z))
  return(fx)
}

# function for g^2 weights  
gg = function(z) {
  fx = 1/(D2(z)*D2(z))
}

# Calculate weighted average ET for left basin ===================================

x = c(0.09,ETL[,1])   # Add a point to the left so there are 20 gaps
dx = diff(x)
Tl = ETL[,2]

# Q by Stack Overflow method for double integration

innerintegral = 
  Vectorize( function(u) { gg(u)*exp((2/sigma^2)*integrate(finner, lower=x[1], upper=u)$value) } )
Qinv = integrate(innerintegral,lower=x[1],upper=x[20])$value
print('',quote=F)
print('1/Q and Q ',quote=F)
print(Qinv)
Q = 1/Qinv
print(c('Q = ',Q),quote=F)

#
wts = rep(0,20)  
for(i in 1:20) {
  ipiece = integrate(finner,lower=x[1],upper=x[i+1])
  v = (2/sigma^2)*ipiece$value
  ev = exp(v)
  opiece = ev*Q*integrate(gg,lower=x[1],upper=x[i+1])$value
  #opiece = ev*Q*gg(x[i+1])
  wts[i] = opiece
}
wts = wts/sum(wts)
print('',quote=F)
print('weights',quote=F)
print(wts)
#
meanETl = sum(Tl*wts)
print('',quote=F)
print('Mean ET for left basin',quote=F)
print(meanETl)
print('-----------------------------------------------',quote=F)

windows()
plot(x[1:20],wts,pch=19,xlab='x',ylab='weight',
     main='left basin',cex.lab=1.5,cex.axis=1.5)
grid()

# Calculate weighted average ET for right basin ===========================================================

x = c(ETR[,1],ETR[20,1]+0.2) # Add a point to the right so there are 20 gaps
dx = diff(x)
Tr = ETR[,2]

# Q by Stack Overflow method for double integration

innerintegral = 
  Vectorize( function(u) { gg(u)*exp((2/sigma^2)*integrate(finner, lower=x[1], upper=u)$value) } )
Qinv = integrate(innerintegral,lower=x[1],upper=x[20])$value
print('',quote=F)
print('1/Q and Q ',quote=F)
print(Qinv)
Q = 1/Qinv
print(c('Q = ',Q),quote=F)

#
wts = rep(0,20)  
for(i in 1:20) {
  ipiece = integrate(finner,lower=x[1],upper=x[i+1])
  v = (2/sigma^2)*ipiece$value
  ev = exp(v)
  opiece = ev*Q*integrate(gg,lower=x[1],upper=x[i+1])$value
  wts[i] = opiece
}
wts = wts/sum(wts)
print('',quote=F)
print('weights',quote=F)
print(wts)
#
meanETr = sum(Tr*wts)
print('',quote=F)
print('Mean ET for right basin',quote=F)
print(meanETr)
print('-----------------------------------------------',quote=F)

windows()
plot(x[1:20],wts,pch=19,xlab='x',ylab='weight',
     main='right basin',cex.lab=1.5,cex.axis=1.5)
grid()

