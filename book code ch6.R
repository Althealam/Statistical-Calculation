###p136 6.1###
m<-1000
g<-numeric(m)
for (i in 1:m){
  x<-rnorm(2)
  g[i]<-abs(x[1]-x[2])
}
est<-mean(g)
se<-sqrt(sum((g-mean(g))^2))/m
##############

###p141 6.2###
#计算第一层切尾均值的抽样分布
n<-20
m<-1000
tmean<-numeric(m)
for(i in 1:m){
  x<-sort(rnorm(n))
  tmean[i]<-sum(x[2:(n-1)])/(n-2)
}
mse<-mean(tmean^2)
se<-sqrt(sum((tmean-mean(tmean))^2))/m
#计算样本中位数
n<-20
m<-1000
tmean<-numeric(m)
for(i in 1:m){
  x<-sort(rnorm(n))
  tmean[i]<-median(x)
}
mse<-mean(tmean^2)
se<-sqrt(sum((tmean-mean(tmean))^2))/m
##############

###p142 6.3###
set.seed(522)
n<-20
K<-n/2-1
mse<-matrix(0,n/2,6)
trimmed.mse<-function(n,m,k,p){
  tmean<-numeric(m)
  for(i in 1:m){
    sigma<-sample(c(1,10),size=n,replace=TRUE,prob=c(p,1-p))
    x<-sort(rnorm(n,0,sigma))
    tmean[i]<-sum(x[(k+1):(n-k)])/(n-2*k)
  }
  for(k in 0:K){
    mse[k+1,1:2]<-trimmed.mse(n=n,m=m,k=k,p=1.0)
    mse[k+1,3:4]<-trimmed.mse(n=n,m=m,k=k,p=.95)
    mse[k+1,5:6]<-trimmed.mse(n=n,m=m,k=k,p=.9)
  }
  for(k in 0:k){
    mse[k+1,1:2]<-trimmed.mse(n=n,m=m,k=k,p=1.0)
    mse[k+1,3:4]<-trimmed.mse(n=n,m=m,k=k,p=.95)
    mse[k+1,5:6]<-trimmed.mse(n=n,m=m,k=k,p=.9)
  }
}
##############

###p144 6.4###
n<-20
alpha<-.05
UCL<-replicate(1000,expr={
  x<-rnorm(n,mean=0,sd=2)
  (n-1)*var(x)/qchisq(alpha,df=n-1)
})
sum(UCL>4)
##############

###p146 6.6###
n<-20
alpha<-.05
UCL<-replicate(1000,expr={
  x<-rchisq(n,df=2)
  (n-1)*var(x)/qchisq(alpha,df=n-1)
})
sum(UCL>4)
mean(UCL>4)
##############

###p148 6.7###
#经验第一类错误率#
n<-20
alpha<-.05
mu0<-500
sigma<-100
m<-10000
p<-numeric(m)
for(j in 1:m){
  x<-rnorm(n,mu0,sigma)
  #对服从正态分布的随机样本进行检验，判断是否接受原假设
  ttest<-t.test(x,alternative="greater",mu=mu0)
  ##t.test(x,y=NULL,alternative=c("two sides","less","greater"),mu=0,paired=FALSE,var,equal=FALSE,conf.level=0.95)
  #alternative:two sides表示双侧检验,less表示右侧检验,greater表示左侧检验
  p[j]<-ttest$p.value
  #t检验的统计量p.value
}
p.hat<-mean(p<alpha)
#找出当原假设为真而拒绝原假设时的p.value，从而计算显著检验比例（观测第一类的错误率）
se.hat<-sqrt(p.hat*(1-p.hat)/m)
print(c(p.hat,se.hat))
##############

###p150 6.8###
n<-c(10,20,30,50,100,500)
cv<-qnorm(.975,0,sqrt(6/n)) #qnorm计算百分位数之下的 Z分数值

#利用sk函数计算样本偏度估计量
sk<-function(x){
  xbar<-mean(x)
  m3<-mean((x-xbar)^3)
  m2<-mean((x-xbar)^2)
  return(m3/m2^1.5)
}

p.reject<-numeric(length(n)) #计算拒绝原假设时的p值
m<-10000
for(i in 1:length(n)){
  sktests<-numeric(m)
  #检验决策以1(拒绝H0)或0(不拒绝H0)储存在向量sktests中
  for(j in 1:m){
    x<-rnorm(n[i])
    #as.integer将数据转化为整数
    sktests[j]<-as.integer(abs(sk(x))>=cv[i])
  }
  p.reject[i]<-mean(sktests)
}
##############

###p152 6.9###
n<-20
m<-1000
mu0<-500
sigma<-100
mu<-c(seq(450,650,10))
M<-length(mu)
power<-numeric(M)
for(i in 1:M){
  mu1<-mu[i]
  pvalues<-replicate(m,expr={
    x<-rnorm(n,mean=mu1,sd=sigma)
    ttest<-t.test(x,alternative="greater",mu=mu0)
    ttest$p.value
  })
  power[i]<-mean(pvalues<=.05)
  #将功效的值储存在power中
}
library(Hmisc)
plot(mu,power)
abline(v=mu0,lty=1)
abline(h=.05,lty=1)
#添加垂直误差线
se<-sqrt(power*(1-power)/m)
errbar(mu,power,yplus=power+se,yminus=power-se,xlab=bquote(theta))
lines(mu,power,lty=3)
detach(package:Hmisc)
##############

###p154 6.10###
#偏度正态性检验的功效#
alpha<-.1
n<-30
m<-2500
epsilon<-c(seq(0,.15,.01),seq(.15,1,.05))
N<-length(epsilon)
pwr<-numeric(N)
cv<-qnorm(1-alpha/2,0,sqrt(6*(n+2)/((n+1)*(n+3))))

for(j in 1:N){
  e<-epsilon[j]
  sktests<-numeric(m)
  for(i in 1:m){
    sigma<-sample(c(1,10),replace=TRUE,size=n,prob=c(1-e,e))
    x<-rnorm(n,0,sigma)
    sktests[i]<-as.integer(abs(sk[x])>=cv)
  }
  pwr[j]<-mean(sktests)
  #绘制偏度检验功效的功效曲线
  plot(epsilon,pwr,type="b",xlab=bquote(epsilon),ylim=c(0,1))
  abline(h=.1,lty=3)
  se<-sqrt(pwr*(1-pwr)/m)
  lines(epsilon,pwr+se,lty=3)
  lines(epsilon,pwr-se,lty=3)
}
##############

###p157 6.11###
#比较正态性检验的功效#
library(energy)
alpha<-.1
n<-30
m<-500
test1<-test2<-test3<-numeric(m)

cv<-qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
sim<-matrix(0,11,4)
for(i in 0:10){
  epsilon<-i*.1 #将epsilon类比p值，生成不同的p值从而生成混合正态的样本
  for(j in 1:m){
    e<-epsilon
    sigma<-sample(c(1,10),replace=TRUE,size=n,prob=c(1-e,e))
    #sample(a，n)在序列a中抽取n个元素，并将n个元素以列表的形式返回
    x<-rnorm(n,0,sigma)
    #针对不同的sigma生成n个服从正态分布的x值
    test1[j]<-as.integer(abs(sk(x)>=cv)) #偏度检验
    test2[j]<-as.integer(shapiro.test(x)$p.value<=alpha) #s-w检验
    test3[j]<-as.integer(mvnorm.etest(x,R=200)$p.value<=alpha) #Energy检验
  }
  print(c(epsilon,mean(test1),mean(test2),mean(test3)))
  sim[i+1,]<-c(epsilon,mean(test1),mean(test2),mean(test3))
}
detach(package:energy)

#绘制图形
plot(sim[,1],sim[,2],ylim=c(0,1),type="1",xlab=bquote(epsilon),ylab="power")
lines(sim[,1],sim[,3],lty=2)
lines(sim[,1],sim[,4],lty=4)
abline(h=alpha,lty=3)
legend("topright",1,c("skewness","S-W","energy"),lty=c(1,2,4),inset=.02)
##############

###p159 6.12###
x1<-rnorm(20,0,sd=1)
x2<-rnorm(20,0,sd=1.5)
y<-c(x1,x2)
group<-rep(1:2,each=length(x1))
boxplot(y~group,boxwex=.3,xlim=c(.5,2.5),main="")
points(group,y)
range(x1)
range(x2)
###p160 6.13###
maxout<-function(x,y){
  X<-x-mean(x)
  Y<-y-mean(y)
  outx<-sum(X>max(Y))+sum(X<min(Y))
  outy<-sum(Y>max(X))+sum(Y<min(X))
  return(max(c(outx,outy)))
}

n1<-n2<-20
mu1<-mu2<-0
sigma1<-sigma2<-1
m<-1000

stat<-replicate(m,expr={
  x<-rnorm(n1,mu1,sigma1)
  y<-rnorm(n2,mu2,sigma2)
  maxout(x,y)
})
print(cumsum(table(stat))/m)
print(quantile(stat,c(.8,.9,.95)))
################

###p162 6.14###
count5test<-function(x,y){
  X<-x-mean(x)
  Y<-y-mean(y)
  outx<-sum(X>max(Y))+sum(X<min(Y))
  outy<-sum(Y>max(X))+sum(Y<min(X))
  return(as.integer(max(c(outx,outy))>5))
}
n1<-n2<-20
mu1<-mu2<-0
sigma1<-sigma2<-1
m<-10000
tests<-replicate(m,expr={
  x<-rnorm(n1,mu1,sigma1)
  y<-rnorm(n2,mu2,sigma2)
  x<-x-mean(x)
  y<-y-mean(y)
  count5test(x,y)
})
alphahat<-mean(tests)
print(alphahat)
###p162 6.15###
n1<-20
n2<-30
mu1<-mu2<-0
sigma1<-sigma2<-1
m<-10000

alphahat<-mean(replicate(m,expr={
  x<-rnorm(n1,mu1,sigma1)
  y<-rnorm(n2,mu2,sigma2)
  x<-x-mean(x) #中心化
  y<-y-mean(y)
  count5test(x,y)
}))
print(alphahat)
###p162 6.16###
#使用蒙特卡洛方法估计Count Five检验的功效#
sigma1<-1
sigma<-1.5

power<-mean(replicate(m,expr={
  x<-rnorm(20,0,sigma1)
  y<-rnorm(20,0,sigma2)
  count5test(x,y)
}))
print(power)
##############