#####chapter 6作业代码参考答案#####

####6-3####
#对n=10,20,30,40,50的样本绘制6.9中的t检验功效曲线
f.power<-function(n=10,mu=seq(450,650,10)){
  mu0<-500
  sigma<-100
  m<-1000
  M<-length(mu)
  power<-numeric(M)
  for(i in 1:M){
    mu1<-mu[i]
    pvalues<-replicate(m,expr={
      x<-rnorm(n,mean=mu1,sd=sigma)
      ttest<-t.test(x,alternative="greater",mu=mu0)
      ttest$p.value
    })
    power[i]<-mean(pvalues<=0.05)
  }
  return(power)
}

power<-matrix(NA,nrow=5,ncol=21)
plot(x=NULL,xlim=c(450,650),ylim=c(0,1),ylab="power",xlab=bquote(theta))
legend(x='bottomright',lty=1,col=c(1:5),cex=0.8,legend=c(expression(n==10),
                                                         expression(n==20),
                                                         expression(n==30),
                                                         expression(n==40)
                                                         expression(n=50)))
mu<-seq(450,650,10)
for(i in 1:5){
  power[i,]<-f.power(n=10*i)
  lines(mu,power[i,],col=i)
}
############################


####6-4####
#服从含有未知参数的对数正态分布的随机样本，构造参数mu的95%置信区间#
#使用mc方法得到置信水平的经验估计#
#分析：当X服从对数正态分布时，logX服从正态分布#
n<-20
alpha<-0.05
UCL<-replicate(1000,expr={
  x<-rlnorm(n,meanlog=0,sdlog=2)
  abs(mean(log(x)-0))<=sd(log(x))*qt(1-alpha/2,df=n-1)/sqrt(n)
  #qt函数用于生成分位数
})
sum(UCL)
mean(UCL)
#结论：经验置信水平与理论置信水平很接近
############################


####6-5####
#使用mc估计t区间对x^2(2)数据中大小为n=20的随机样本的覆盖率#
#分析：覆盖率可以理解为置信水平，与6-4类似，将总体X的分布改为卡方分布即可#
n<-20
alpha<-0.5
UCL<-replicate(1000,expr={
  x<-rchisq(n,df=2)
  abs(mean(x)-2)<=sd(x)*qt(1-alpha/2,df=n-1)/sqrt(n)
})
sum(UCL)
mean(UCL)
#结论：（1）t区间对均值的覆盖率为90.9%，低于95%
#      （2）样本服从卡方分布时，方差区间的覆盖率只有大约77%（6.6）
#因此在样本的分布改为卡方分布的情形下，方差区间的表现更加稳健
############################


####6-8####
#计算等方差F检验，对小样本、中样本、大样本比较Count Five检验和F检验的功效#
#分析：countfive检验可以利用书上的代码，F检验利用var.test函数
exercise_6_8<-function(){
  count5test<-function(x,y){
    X<-x-mean(x)
    Y<-y-mean(y)
    outx<-sum(X>max(Y))+sum(X<min(Y))
    outy<-sum(Y>max(X))+sum(Y<min(X))
    return(as.integer(max(c(outx,outy))>5))
    #as.integer用于进行数据类型的转化
  }
  n<-c(20,200,1000) #分别对应于小样本、中样本、大样本
  mu1<-mu2<-0
  sigma1<-1
  sigma2<-1.5
  m<-10000
  power1<-power2<-numeric(length(n))
  set.seed(1234)
  for(i in 1:length(n)){
    power1[i]<-mean(replicate(m,expr={
      x<-rnorm(n[i],mu1,sigma1)
      y<-rnorm(n[i],mu2,sigma2)
      x<-x-mean(x)
      y<-y-mean(y)
      count5test(x,y)
    }))
    pvalues<-replicate(m,expr={
      x<-rnorm(n[i],mu1,sigma1)
      y<-rnorm(n[i],mu2,sigma2)
      Ftest<-var.test(x,y,ratio=1,alternative="two.sided",conf.level=0.945)
      Ftest$p.value
    })
    power2[i]<-mean(pvalues<=0.055)
  }
  return(data.frame(power1,power2))
}
exercise_6_8()
############################


####6-9####
#对服从X分布的随机样本计算Gini系数，X分布取均匀分布、Bernoulli(0,1)分布#
#分析：顺序统计量可以利用sort()
gini<-function(x,mu="NULL"){
  n<-length(x)
  x.sort<-sort(x)
  coe<-2*c(1:n)-n-1
  if(mu=="NULL")
    mu<-mean(x)
  G<-sum(coe*x.sort)/(mu*n^2)
  return(G)
}
mmq.list<-function(x){
  x.mean<-mean(x,na.rm=TRUE)
  x.median<-median(x,na.rm=TRUE)
  x.deciles<-quantile(x,probs=c(1:9)/10,na.rm=TRUE)
  result<-list(x.mean,x.median,x.deciles)
  names(result)<-c("mean","median","deciles")
  return(result)
}
n<-20
m<-1000
norm.gini<-unif.gini<-binom.gini<-numeric(m)
for(i in 1:m){
  x<-rlnorm(n,meanlog=0,sdlog=1)
  y<-runif(n)
  z<-rbinom(n,size=1,prob=0.1)
  norm.gini[i]<-gini(x)
  unif.gini[i]<-gini(y)
  binom.gini[i]<-gini(z)
}
mmq.list(norm.gini)
mmq.list(unif.gini)
mmq.list(binom.gini)
draw.hist<-function(x){
  hist(x,prob=TRUE,xlim=c(max(0.75*min(x),0),min(1.25*max(x),1)),main="")
}
draw.hist(norm.gini)
draw.hist(unif.gini)
draw.hist(binom.gini[which(binom.gini>0)])
#结论：对数标准正态样本、均匀分布样本、Bernoulli(0,1)分布样本的Gini系数的平均值、中位数与十分位数分别由mmq.list(norm.gini),mmq.list(unif.gini),mmq.list(binom.gini)给出
############################


####6.A####
#当抽样总体非正态时，使用mc模拟来研究t检验的Type I error是否约等于理论显著水平alpha#
n<-20
alpha<-0.05
mu0<-1
m<-10000

#卡方分布#
p.i<-numeric(m)
for(j in 1:m){
  x<-rchisq(n,df=1) #生成服从卡方分布的随机数x，数量为n，自由度为1
  ttest.i<-t.test(x,alternative="two.sided",mu=mu0)
  p.i[j]<-ttest.i$p.value
}
p.i.hat<-mean(p.i<alpha)
se.i.hat<-sqrt(p.i.hat*(1-p.i.hat)/m)
result.i<-c(p.i.hat,se.i.hat)

#均匀分布#
p.ii<-numeric(m)
for(j in 1:m){
  y<-runif(n,min=0,max=2)
  ttest.ii<-t.test(y,alternative="two.sided",mu=mu0)
  p.ii[j]<-ttest.ii$p.value
}
p.ii.hat<-mean(p.ii<alpha)
se.ii.hat<-sqrt(p.ii.hat*(1-p.ii.hat)/m)
result.ii<-c(p.ii.hat,se.ii.hat)

#指数分布#
p.iii<-numeric(m)
for(j in 1:m){
  z<-rexp(n,rate=1) #生成参数为1的服从指数分布的n个随机数
  ttest.iii<-t.test(y,alternative="two.sided",mu=mu0)
  p.iii[j]<-ttest.iii$p.value
}
p.iii.hat<-mean(p.iii<alpha)
se.iii.hat<-sqrt(p.iii.hat*(1-p.iii.hat)/m)
result.iii<-c(p.iii.hat,se.iii.hat)

print(result.i)
print(result.ii)
print(result.iii)
#结论：三种情形下犯第一类错误的概率都接近0.05#
############################



####6.B####
#分析：程序设计思路参考课本6.9#
library(MASS)
n<-20
m<-1000
mu<-rep(0,2)
co<-seq(-1.8,1.8,0.1)
alpha<-0.05
M<-length(co)
sigma<-vector("list",M)
for(i in 1:M){
  sigma[[i]]<-matrix(c(1,co[i],co[i],4),ncol=2)
}
p<-function(method="pearson"){
  power<-numeric(M)
  for(j in 1:M){
    sigma1<-sigma[[j]]
    pvalues<-replicate(m,expr={
      X<-mvrnorm(n,mu,sigma1)
      x<-X[,1]
      y<-X[,2]
      cortest<-cor.test(x,y,method=method)
      cortest$p.value
    })
    power[j]<-mean(pvalues<=alpha)
  }
  return(power)
}
power<-matrix(rep(0,3*M),nrow=3)
power[1,]<-p()
power[2,]<-p("spearman")
power[3,]<-p("kendall")
plot(NULL,NULL,xlim=range(co),ylim=c(0,1),xlab="Cov(X,Y)",ylab="power")
for(i in 1:3){
  lines(co,power[i,],lty=i)
}
legend(x='top',lty=1:3,cex=0.8,legend=c("Pearson","Spearman","Kendall"))
p1<-function(method="pearson"){
  power<-numeric(M)
  for(j in 1:M){
    sigma1<-sigma[[j]]
    pvalues<-replicate(m,expr={
      X<-mvrnorm(n,mu,sigma1)
      x<-X[,1]
      y<-exp(X[,2])
      cortest<-cor.test(x,y,method=method)
      cortest$p.value
    })
    power[j]<-mean(pvalues<=alpha)
  }
  return(power)
}
power1<-matrix(rep(0,3*M),nrow=3)
power1[1,]<-p1()
power1[2,]<-p1("spearman")
power1[3,]<-p1("kendall")
plot(NULL,NULL,xlim=range(co),ylim=c(0,1),xlab="Cov(X,Z)",ylab="power")
for(i in 1:3){
  lines(co,power1[i,],lty=i)
}
legend(x='top',lty=1:3,cex=0.8,legend=c("Pearson","Spearman","Kendall"))

