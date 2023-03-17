library(mvtnorm) #用于生成模拟数据
library(MASS)
library(pwr)
library(ICSNP)
library(Hmisc)

#################两组独立样本#################
N1<-c(15,25)
Sigma1<-matrix(c(16,-2,-2,9),byrow=TRUE,ncol=2)
y1<-round(rmvnorm(N1[1],mean=mu1,sigma=Sigma1))
y2<-round(rmvnorm(N1[2],mean=mu2,sigma=Sigma1))
mu1<-c(-4,4)
mu2<-c(3,3)
y12<-rbind(y1,y2)
y12
w1<-factor(rep(1:2,N1))

#正态性检验#
shapiro.test(y1)
shapiro.test(y2)


##########1.t检验#############
#t检验仅仅适用于检验服从正态或者近似正态分布、满足方差齐性的数据；
set.seed(2022)
n<-20
m<-1000
p.t1<-numeric(m)
alpha<-0.025
t.test(y1,y2,"two.sided")
#结论：p-value为0.01455，因此我们拒绝原假设，认为两者之间的均值不相等
#比较不同mu时的t检验的功效
for(j in 1:m){
  y10<-round(rmvnorm(N1[1],mean=mu1,sigma=Sigma1))
  y20<-round(rmvnorm(N1[2],mean=mu2,sigma=Sigma1))
  ttest<-t.test(y10,y20,alternative="two.sided")
  p.t1[j]<-ttest$p.value
}
p.hat.t1<-mean(p.t1<alpha)
p.hat.t1
#此时计算得到的Type I error为0.721
#结论：对多组样本检验时若使用t检验，则会导致type I error膨胀
#############################


########2.Hotelling T2检验##########
set.seed(200)
n<-20
m<-1000
p.ht1<-numeric(m)
alpha<-0.025
HotellingsT2(y12~w1)
#结论：计算得到的p值小于0.01，因此我们拒绝零假设，认为均值不等
#计算Hotelling T2的Type I eror
for(j in 1:m){
  y11<-round(rmvnorm(N1[1],mean=mu1,sigma=Sigma1))
  y21<-round(rmvnorm(N1[2],mean=mu2,sigma=Sigma1))
  y12<-rbind(y11,y21)
  httest<-HotellingsT2(y12~w1)
  p.ht1[j]<-httest$p.value
}
p.hat.ht1<-mean(p.ht1<alpha)
p.hat.ht1
#此时计算得到的Type I error为0.996
##################################


#########两组不独立的样本############
N2<-20
P<-2
mu<-c(95,100,85,105)
Sigma2<-15
#生成四组样本，服从均值不同的正态分布
x1t0<-rnorm(N2,mean=mu[1],sd=Sigma2)
x1t1<-rnorm(N2,mean=mu[2],sd=Sigma2)
x2t0<-rnorm(N2,mean=mu[3],sd=Sigma2)
x2t1<-rnorm(N2,mean=mu[4],sd=Sigma2)
x<-data.frame(id=factor(rep(1:N2,times=P)))
#将四组样本合并变为两组多维的样本
x1<-c(x1t0,x1t1)
x2<-c(x2t0,x2t1)
w2<-factor(rep(1:P,each=N2),labels=c('t0','t1'))
df<-aggregate(cbind(x1,x2)~id,data=x,FUN=diff)
#aggregate用于数据整形
Df<-data.matrix(df[,-1])
muH0<-c(0,0)


##正态性检验##
qqnorm(Df)
#结论：生成qq图我们发现散点大致分布在直线附近，因此正态性成立


###############1.t检验######################
set.seed(200)
t.test(Df[,1],Df[,2],"greater")
#结论：p值大于0.05，因此我们接受零假设
#t检验的Type I error#
m<-1000
p.t2<-numeric(m)
alpha<-0.05
for(j in 1:m){
  #重新进行抽样
  x1t0<-rnorm(N2,mean=mu[1],sd=Sigma2)
  x1t1<-rnorm(N2,mean=mu[2],sd=Sigma2)
  x2t0<-rnorm(N2,mean=mu[3],sd=Sigma2)
  x2t1<-rnorm(N2,mean=mu[4],sd=Sigma2)
  x10<-c(x1t0,x1t1)
  x20<-c(x2t0,x2t1)
  ttest<-t.test(x10,x20)
  p.t2[j]<-ttest$p.value
}
p.hat.t2<-mean(p.t2<alpha)
p.hat.t2
#结论：非独立样本的t检验的Type I error的值为0.16
#################################################




###############2.Hotelling t2检验##################
#test='f'基于f分布 test='chi'基于卡方分布
set.seed(2022)
HotellingsT2(Df,mu=muH0,test='f')
#结论：p值小于0.001，因此我们拒绝零假设，认为多元变量中的至少一个或者多个变量的组合在组合间表现出的均值显著不同
#Hotelling T2检验的Type I error
m<-1000
p.ht2<-numeric(m)
alpha<-0.05
for(j in 1:m){
  #重新进行抽样
  x1t0<-rnorm(N2,mean=mu[1],sd=Sigma2)
  x1t1<-rnorm(N2,mean=mu[2],sd=Sigma2)
  x2t0<-rnorm(N2,mean=mu[3],sd=Sigma2)
  x2t1<-rnorm(N2,mean=mu[4],sd=Sigma2)
  x10<-c(x1t0,x1t1)
  x20<-c(x2t0,x2t1)
  x12<-rbind(x10,x20)
  x<-data.frame(id=factor(rep(1:N2,times=P)))
  df<-aggregate(cbind(x10,x20)~id,data=x,FUN=diff)
  Df<-data.matrix(df[,-1])
  w3<-factor(rep(1:P,each=N2),labels=c('t0','t1'))
  httest<-HotellingsT2(Df,mu=muH0,test='f')
  p.ht2[j]<-httest$p.value
}
p.hat.ht2<-mean(p.ht2<alpha)
p.hat.ht2
#结论：计算得到的Type I error为0.954
####################################


##############单样本#################
n<-20
alpha<-.05
mu0<-500
sigma<-100
mu<-c(seq(450,650,10))
M<-length(mu)
power.t<-numeric(M)
power.ht<-numeric(M)
for(i in 1:M){
  mu1<-mu[i]
  pvalues.t<-replicate(m,expr={
    x<-rnorm(n,mean=mu1,sd=sigma)
    ttest<-t.test(x,alternative="greater",mu=mu0)
    ttest$p.value
  })
  power.t[i]<-mean(pvalues.t<=.05)
  pvalues.ht<-replicate(m,expr={
    y<-rnorm(n,mean=mu1,sd=sigma)
    httest<-HotellingsT2(y,mu=mu0)
    httest$p.value
  })
  power.ht[i]<-mean(pvalues.ht<=.05)
}
library(Hmisc)
plot(mu,power.t,main='Efficacy curve of t-test')
lines(mu,power.t,lty=2) #绘制t检验的功效曲线

plot(mu,power.ht,main='Efficacy curve of Hotelling T-square test')
lines(mu,power.ht,lty=2) #绘制Hotelling T2检验的功效检验
detach(package:Hmisc)
####################################
