###p106 5.1###
m<-1000
x<-runif(m)
theta.hat<-mean(exp(-x))
print(theta.hat)
print(theta)
##############

###p106 5.2###
m<-10000
x<-runif(m,min=2,max=4)
theta.hat<-mean(exp(-x))
print(theta.hat)
print(exp(-2)-exp(-4))
##############

###p109 5.3###
x<-seq(.1,2.5,length=10)
m<-10000
u<-runif(m)
cdf<-numeric(length(x)) #cdf为累积分布函数，numeric表示生成向量
#cdf中储存x的10个不同值的函数值的估计值
for(i in 1:length(x)){
  g<-x[i]*exp(-(u*x[i])^2/2)
  cdf[i]<-mean(g)/sqrt(2*pi)+0.5
}
Phi<-pnorm(x)
#pnorm用于生成正态分布的密度函数
print(round(rbind(x,cdf,Phi),3))
##############

###p110 5.4###
x<-seq(.1,2.5,length=10)
m<-10000
z<-rnorm(m)
dim(x)<-length(x)
p<-apply(x,MARGIN=1,FUN=function(x,z){mean(z<x),z=z})
Phi<-pnorm(x)
print(round(rbind(x,p,Phi),3))
##############

###p112 5.5###
x<-2
m<-10000
z<-rnorm(m)
g<-(z<x)
v<-mean((g-mean(g))^2)/m
cdf<-mean(g)
c(cdf,v)
c(cdf-1.96*sqrt(v),cdf+1.96*sqrt(v))
##############

###p117 5.6###
#函数MC.Phi可以计算积分的蒙特卡洛估计量，同时还可以选择是否使用对偶抽样来计算估计值
MC.Phi<-function(x,R=10000,antithetic=TRUE){
  u<-runif(R/2)
  if(!antithetic) v<-runif(R/2) else
    v<-1-u
  cdf<-numeric(length(x))
  for(i in 1:length(x)){
    g<-x[i]*exp(-(u*x[i])^2/2)
    cdf[i]<-mean(g)/sqrt(2*pi)+0.5
  }
}
x<-seq(.1,2.5,length=10)
m<-10000
z<-rnorm(m)
dim(x)<-length(x)
p<-apply(x,MARGIN=1,FUN=function(x,z){mean(z<x)},z=z)
Phi<-pnorm(x) #理论值
print(round(rbind(x,p,Phi),3))

MC1<-MC2<-numeric(m)
x<-1.95
for(i in 1:m){
  MC1[i]<-MC.Phi(x,R=1000,anti=FALSE)
  MC2[i]<-MC.Phi(x,R=1000)
}
print(sd(MC1))
print(sd(MC2))
print((var(MC1)-var(MC2))/var(MC1)) #计算方差缩减量
##############

###p120 5.7###
m<-10000
a<--12+6*(exp(1)-1)
U<-runif(m)
T1<-exp(U)
T2<-exp(U)+a*(U-1/2)
##############

###p121 5.8###
#使用控制变量法，定义一个f(x)使得其与g(x)很接近
f<-function(u)
  exp(-.5)/(1+u^2)
g<-function(u)
  exp(-u)/(1+u^2)
set.seed(510)
u<-runif(10000)
B<-f(u)
A<-g(u)
cor(A,B) 
a<- -cov(A,B)/var(B)  #计算c.star
m<-100000
u<-runif(m)
T1<-g(u)
T2<-T1+a*(f(u)-exp(-.5)*pi/4)
##############

###p124 5.9###
set.seed(510) #随机数种子
u<-runif(10000)
f<-exp(-.5)/(1+u^2)
g<-exp(-u)/(1+u^2)
c.star<- -lm(g~f)$coeff[2]
#通过拟合回归模型进行再次估计
mu<-exp(-.5)*pi/4

u<-runif(10000)
f<-exp(-.5)/(1+u^2)
g<-exp(-u)/(1+u^2)
L<-lm(g~f)
theta.hat<-sum(L$coeff*c(1,mu))
##############

###p127 5.10###
m<-10000
theta.hat<-se<-numeric(5)
g<-function(x){
  exp(-x-log(1+x^2))*(x>0)*(x<1)
}

x<-runif(m) #使用f0
fg<-g(x)
theta.hat[1]<-mean(fg)
se[1]<-sd(fg)

x<-rexp(m,1) #使用f1
fg<-g(x)/exp(-x)
theta.hat[2]<-mean(fg)
se[2]<-sd(fg)

x<-rcauchy(m) #使用f2
i<-c(which(x>1),which(x<0))
x[i]<-2
fg<-g(x)/dcauchy(x)
theta.hat[3]<-mean(fg)
se[3]<-sd(fg)

u<-runif(m) #使用f3
x<- -log(1-u*(1-exp(-1)))
fg<-g(x)/(exp(-x)/(1-exp(-1)))
theta.hat[4]<-mean(fg)
se[4]<-sd(fg)

u<-runif(m) #使用f4
x<-tan(pi*u/4)
fg<-g(x)/(4/((1+x^2)*pi))
theta.hat[5]<-mean(fg)
se[5]<-sd(fg)

rbind(theta.hat,se)
##############

###p131 5.11###
M<-20
T2<-numeric(4)
estimates<-matrix(0,10,2)
g<-function(x){
  exp(-x-log(1+x^2))*(x>0)*(x<1)
for(i in 1:10){
  estimate[i,1]<-mean(g(runif(M)))
  T2[1]<-mean(g(runif(M/4,0,.25)))
  T2[2]<-mean(g(runif(M/4,.25,.5)))
  T2[3]<-mean(g(runif(M/4,.5,.75)))
  T2[4]<-mean(g(runif(M/4,.75,1)))
  estimate[i,2]<-mean(T2)
}
}
##############

###p133 5.12###
M<-10000
k<-10
r<-M/k
N<-50
T2<-numeric(k)
estimates<-matrix(0,N,2)
g<-function(x){
  exp(-x-log(1+x^2))*(x>0)*(x<1)
}
for(i in 1:N){
  estimates[i,1]<-mean(g(runif(M)))
  for(j in 1:k)
    T2[j]<-mean(g(runif(M/k.(j-1)/k,j/k)))
  estimates[i,2]<-mean(T2)
}
}
##############