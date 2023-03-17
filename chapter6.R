#####p164 6.3#####
#对t检验用模拟方法估计功效并绘制经验功效曲线#
n<-c(10,20,30,40,50) #不同的样本个数
m<-1000 #重复试验的次数
mu0<-500
sigma<-100
mu<-c(seq(450,650,10)) #生成不同的mu值，并将其作为x轴
M<-length(mu)
power<-numeric(M) #不同试验的功效
color<-c("red","blue","green","black","purple") #设置曲线的颜色
lty<-c(2,3,4,5,6) #设置曲线的线条形状
label<-c('','','','','')
for(k in 1:5){
  for(i in 1:M){
    mu1<-mu[i]
    pvalues<-replicate(m,expr={
      x<-rnorm(n[k],mean=mu1,sd=sigma)
      ttest<-t.test(x,alternative="greater",mu=mu0)
      ttest$p.value
    })
    power[i]<-mean(pvalues<=.05)
  }
  if(k==1){
  plot(mu,power,type='l',col="red") #当n=10时利用plot
  }
  else{
    lines(mu,power,col=color[k],lty=lty[k]) #当n=20..50时利用plot
  }
  label[k]<-paste('n=',n[k]) #图例赋值
}
legend("topleft",label,col=color,lty=1) #将图例置于左上角
#legend(x,y,lengend) 在点(x,y)添加图例
#将不同的曲线绘制在同一张图上时，往往先利用plot函数，再利用lines函数
#结论：当随着mu的增加时，不同样本数量对应的经验功效大小在增加直至稳定；
#      当保持mu不变时，随着样本数量的增加，对应的经验功效大小也增加；




#####p164 6.4#####
#使用蒙特卡洛方法得到置信水平的经验估计#
n<-20
alpha<-0.05
m<-1000
I<-numeric(n)
CL<-replicate(m,expr={
  x<-rlnorm(n) #生成对数正态分布的样本 
  #若x服从对数正态分布，则ln(x)服从正态分布
  #构造对数正态分布的关于参数mu的置信区间
  UCL<-mean(log(x))+abs(qt(alpha/2,df=n-1))*sd(log(x))/sqrt(n) #UCL表示置信上限
  LCL<-mean(log(x))-abs(qt(alpha/2,df=n-1))*sd(log(x))/sqrt(n) #LCL表示置信下限
  x>LCL&x<UCL #找到属于置信区间的x值
})
1-mean(CL)
##############
#结论：服从对数正态分布的随机样本的参数mu的经验置信水平为0.8028




#####p164 6.5#####
#t区间对随机样本的覆盖率即为置信水平#
##利用t区间进行判断##
n<-20
alpha<-.05
m<-1000
#对一个给定的概率模型进行反复抽样来确定p(LCL<x<UCL)
a<-replicate(m,expr={
  x<-rchisq(n,df=2) #生成n个服从卡方分布的随机数（自由度为df）(x由两个正态分布的随机数构成)
  #双侧检验
  UCL<-mean(x)+abs(qt(alpha/2,df=n-1))*sqrt(var(x))/sqrt(n)
  LCL<-mean(x)-abs(qt(alpha/2,df=n-1))*sqrt(var(x))/sqrt(n)
  x>LCL&x<UCL #找到属于置信区间内的x值
})
1-mean(a)  #置信水平的值

##利用方差区间进行判断##
UCL<-replicate(m,expr={
  x<-rchisq(n,df=2)
  (n-1)*var(x)/qchisq(alpha,df=n-1)
})
mean(UCL>4) #找出方差区间大于4的区间个数
#############
#结论：利用t区间得到卡方分布的随机样本的覆盖率（置信水平）为0.654
#总结：求置信水平的方法：p(LCL<x<UCL)=1-alpha，其中p的值利用重复抽样计算可得




###？？？###
#####p164 6.8#####
n<-c(10,20,50) #设置不同的样本个数
m<-1000 #重复试验的次数
power.CF<-numeric(length(n)) #设置Count Five检验对不同样本个数对应的功效
power.F<-numeric(length(n)) #设置F检验对不同样本个数对应的功效
sigma1<-1
sigma2<-1.5
alpha<-.055

#双样本等方差的Count Five检验#
count5test<-function(x,y){
  #对数据进行中心化
  X<- x-mean(x)
  Y<- y-mean(y)
  outx<-sum(X>max(Y))+sum(X<min(Y))
  outy<-sum(Y>max(X))+sum(X<min(Y))
  return(as.integer(max(c(outx,outy))>5))
}
  #拒绝H0时返回1，否则返回0

#计算countfive检验的功效
for(i in length(n)){
  n1<-n2<-n[i]
  a<-replicate(m,expr={
    x<-rnorm(n1,0,sigma1)
    y<-rnorm(n2,0,sigma2)
    CF<-count5test(x,y)
    mean(CF)
  })
  power.CF[i]<-mean(a)
}

#计算F检验的功效
for(i in length(n)){ #for循环用于计算不同的样本数量时的情况
  n1<-n2<-n[i]
  b<-replicate(m,expr={
    x<-rnorm(n1,0,sigma1)
    y<-rnorm(n2,0,sigma2)
    F<-var.test(x,y,ratio=1,"two.sided",conf.level=1-alpha)
    Ft<-F$p.value<alpha
    mean(Ft)
  })
  power.F[i]<-mean(b)
}
print(data.frame(n,power.CF,power.F))
##等方差F检验##
#F检验的作用：主要通过比较两组数据的方差，从而判断数据是否具有显著性差异(判断两总体方差是否相等，就可以用F检验)
#F检验的前提：数据满足正态分布，可以使用Shapiro-Will进行正态分布检验（本题数据为服从正态分布的随机数，因此不需要进行正态分布检验）
#F检验的函数：var.test()（alternative设置单侧/双侧检验，conf.level设置置信水平）
#############


#####p164 6.9#####
n<-20 #生成样本的个数
m<-1000 #重复试验的次数

#标准对数正态分布#
G1<-replicate(m,expr={
  g<-numeric(n) #计算各个样本的基尼系数
  for(i in 1:n){
    x<-sort(rlnorm(n)) #利用sort函数对生成的服从标准对数正态分布的随机数进行排序
    g[i]<-sum((2*i-n-1)*x[i])/(n^2*mean(x))
  }
  sum(g)
})
G1.mean<-mean(G1) #计算标准对数正态分布的均值
G1.median<-median(G1) #计算标准对数正态分布的中位数
G1.quantile<-quantile(G1,c(seq(0,1,0.1))) #计算标准对数正态分布的十分位数
hist(G1,main="Gini coefficient of standard lognormal distribution")
print(G1.quantile)

#均匀分布#
G2<-replicate(m,expr={
  g<-numeric(n) #计算各个样本的基尼系数
  for(i in 1:n){
    x<-sort(runif(n)) #利用sort函数对生成的服从均匀分布随机数进行排序
    g[i]<-sum((2*i-n-1)*x[i])/(n^2*mean(x))
  }
  sum(g)
})
G2.mean<-mean(G2) #计算均匀分布的均值
G2.median<-median(G2) #计算均匀分布的中位数
G2.quantile<-quantile(G2,c(seq(0,1,0.1))) #计算均匀分布布的十分位数
hist(G2,main="Gini coefficient of Uniformly distributed")
print(G2.quantile)

#Bernoulli(0.1)分布
G3<-replicate(m,expr={
  g<-numeric(n) #计算各个样本的基尼系数
  for(i in 1:n){
    x<-sort(rbinom(n,size=m,prob=0.1)) #利用sort函数对生成的服从标准对数正态分布的随机数进行排序
    #rbinom(n,size,prob) n为观察的数量（样本容量），size为试验次数，prob为试验成功的概率
    g[i]<-sum((2*i-n-1)*x[i])/(n^2*mean(x))
  }
  sum(g)
})
G3.mean<-mean(G3) #计算Bernoulli(0.1)分布的均值
G3.median<-median(G3) #计算Bernoulli(0.1)分布的中位数
G3.quantile<-quantile(G3,seq(0,1,0.1)) #计算Bernoulli(0.1)分布的十分位数
hist(G3,main="Gini coefficient of Bernoulli distributed")
print(G3.quantile)
#############




#####p165 6.A#####
#判断经验Type I error是否约等于理论显著水平alpha（理论上应接近）#
#Type I error的计算可以参考书本p148经验第一类错误率
n<-20
alpha<-.05
m<-10000

#卡方分布#
a<-replicate(m,expr={
  x<-rchisq(n,df=1)
  ttest1<-t.test(x,alternative="greater",mu=1)
  t1<-ttest1$p.value<alpha #<alpha为H0在显著性水平alpha下被拒绝的条件
  #此时生成逻辑值True或者False
  mean(t1)
})
errorI<-mean(a) #经验Type I error
print(errorI)

#均匀分布#
b<-replicate(m,expr={
  x<-runif(n,0,2)
  ttest2<-t.test(x,alternative="greater",mu=1)
  t2<-ttest2$p.value<alpha
  mean(t2)
})
errorI<-mean(b) #经验Type I error
print(errorI)

#指数分布#
c<-replicate(m,expr={
  x<-rexp(n,1)
  ttest3<-t.test(x,alternative="greater",mu=1)
  t3<-ttest3$p.value<alpha
  mean(t3)
})  
errorI<-mean(c) #经验Type I error
print(errorI)
#############
#总结：t.test(x,y= NULL,alternative = c("two.sided", "less", "greater"),mu = 0, paired = FALSE, var.equal = FALSE,conf.level = 0.95, ...)
#      其中x,y为需要进行检验的值，alternative为选择单侧/双侧检验，mu为进行对比的值，conf.level为显著性水平
#结论：（1）卡方分布：经验Type I error的值为0.0133，与显著性水平相差较大
#      （2）均匀分布：经验Type I error的值为0.0501，与显著性水平接近
#      （3）指数分布：经验Type I error的值为0.0202，与显著性水平相差较大
#注：可以参考书本p148例6.7的Rcode，利用for循环语句代替replicate语句




#####p165 6.B#####
library(MASS)
n<-100
alpha<-.05
mu<-c(0,1) #指定二元正态的均值
Sigmas<-matrix(c(2,1,4,5),nrow=2,ncol=2) #指定协方差矩阵
#matrix函数生成矩阵
x<-mvrnorm(n,mu,Sigmas) #抽样分布为二元正态时生成数据
cor.test(x[,1],x[,2],alternative="greater",method="pearson",conf.level=alpha) #Pearson test
cor.test(x[,1],x[,2],alternative="greater",method="kendall",conf.level=alpha) #Kendall test
cor.test(x[,1],x[,2],alternative="greater",method="spearman",conf.level=alpha) #spearman test

x<-runif(n,0,1)
y<-6*x
cor.test(x,y,alternative="greater",method="pearson",conf.level=alpha)
cor.test(x,y,alternative="greater",method="kendall",conf.level=alpha)
cor.test(x,y,alternative="greater",method="spearman",conf.level=alpha)
#总结：（1）cor.test()可以进行相关性系数的计算与检验
#           函数功能：使成对数据进行相关性检验，其中有三种方法：Pearson、Spearman、Kendall
#           函数参数：cor.test(x, y, alternative = c(“two.sided”, “less”, “greater”), method = c("pearson", "kendall", "spearman"),conf.level = 0.95)
#           参数解释：x,y是供检验的样本；alternative指定是双侧检验还是单侧检验；method为检验的方法；conf.level为检验的置信水平。
#           结果解释：rho.cor.tau为相关系数，p-value为检验p值
#      （2）基于MASS包的mvrnorm函数生成多元正态分布数据
#           mvrnorm(n=n,mu=mu, Sigma=Sigmas)  产生n个服从N（mu，Sigmas)的随机数
#结论：（1）进行pearson test时相关系数为0.3209552
#           进行Kendall test时相关系数为0.2088889
#           进行Spearman test时相关系数为0.3122352
#           显然当抽样分布为二元正态时三种检验的相关性检验（p-value）的功效都要比非参数检验的功效（rho、cor、tau）低
#      （2）当生成线性相关的二元数据时，显然非参数检验的功效高于相关性检验的功效