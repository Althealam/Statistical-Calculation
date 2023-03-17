###p193 7.2###
#使用基于自助法的水手刀法 估计 se(R)的自助法估计 的标准误差，参考"law"数据#
data(law,package="bootstrap")
B<-2000 #bootstrap重复试验的次数
n<-nrow(law) #15个学生，2个指标：LSAT与GPA
R<-numeric(B)
#生成R的bootstrap估计#
for(b in 1:B){
  i<-sample(1:n,size=n,replace=TRUE)
  LSAT<-law$LSAT[i]
  GPA<-law$GPA[i]
  R[b]<-cor(LSAT,GPA)
}
#利用Jackknife计算bootstrap估计的标准误差se.jack
se.jack<-numeric(n)
for(i in 1:n){
 R.jack<-R[-i]
 se.jack[i]<-sqrt((n-1)*mean((R.jack-mean(R.jack))^2))
}
print(sd(se.jack))
print(sqrt((n-1)*mean((se.jack-mean(se.jack))^2)))
#############
#总结：（1）sample(x,size,replace):x为抽取样本的数据来源，size为指定抽样的次数，replace=TRUE表示可以重复抽样
#      （2）data():用来加载指定的数组
#分析：利用bootstrap计算相关系数R（重复试验次数为B），再利用jackknife计算去掉其中一个R时的标准误差se.jack（重复次数为n）
#结果：基于bootstrap的jackknife估计se(R)的bootstrap估计的标准误差为0.0019左右





###p193 7.3###
#对"bootstrap"中的"law"数据的相关性统计量给出一个bootstrap法t置信区间估计
data(law,package="bootstrap")
boot.t.ci<-function(x,B=500,R=100,level=.95,statistic){
  #参数：x为样本数据，statistic为计算统计量的函数
  x<-as.matrix(x)
  n<-nrow(x)
  stat<-numeric(B) #bootstrap重复试验的统计量估计
  se<-numeric(B) #bootstrap重复试验的统计量的标准误差估计
  
  boot.se<-function(x,R,f){
    #计算标准误差估计的函数
    #参数：x为数据来源，R为重复次数，f为计算统计量的函数（statistic）
    x<-as.matrix(x)
    m<-nrow(x)
    #利用replicate函数进行bootstrap（模拟for(i in 1:B)）
    th<-replicate(R,expr={
      i<-sample(1:m,size=m,replace=TRUE) #随机抽取不同的下标i
      f(x[i,]) #计算不同下标i的统计量
    })
    return(sd(th))
  }
  
  for(b in 1:B){
    #计算bootstrap重复试验的统计量
    j<-sample(1:n,size=n,replace=TRUE) #随机抽取下标j
    y<-x[j,] 
    stat[b]<-statistic(y) #计算第b次重复试验的统计量stat[b]
    se[b]<-boot.se(y,R=R,f=statistic) #计算第b次重复试验的标准误差se[b]
  }
  stat0<-statistic(x) #观测统计量theta.hat
  t.stats<-(stat-stat0)/se #第b次重复试验的t统计量
  se0<-sd(stat)
  
  alpha<-1-level
  Qt<-quantile(t.stats,c(alpha/2,1-alpha/2),type=1) #从重复试验t[b]构成的有序样本中找出样本分位数
  names(Qt)<-rev(names(Qt))
  #names函数：给变量修改响应的元素名称
  CI<-rev(stat0-Qt*se0)
  #rev函数：用于将一个向量逆转（不是转置）
}
dat<-cbind(law$LSAT,law$GPA)
stat<-function(dat){
   cor(dat[,1],dat[,2])
}
ci<-boot.t.ci(dat,statistic=stat,B=2000,R=200)
print(ci)
#############
#分析：（参考书本p183 7.4.4自助法t区间的算法过程）
#step1:计算观测统计量theta.hat（在代码中体现为stat0）
#step2:对第b次重复试验（b=1,...B）：(bootstrap)
#      （1）从x中有放回地抽样生成第b个样本
#      （2）对第b个样本计算theta.hat[b]
#       (3)计算或者估计标准误差se.hat(theta.hat[b])(每个bootstrap样本估计不同，bootstrap估计是对目前的自助法样本x[b]重抽样，而不是x)
#      （4）计算t统计量的第b次重复试验t[b]（在代码中体现为t.stats）
#step3:从重复试验t[b]构成的有序样本中找出样本分位数t.star(alpha/2)与t.star(1-alpha/2)
#step4:计算重复试验theta.hat[b]的样本标准差se.hat[theta.hat]
#step5:计算置信限





###p193 7.4###
data(aircondit,package="boot") 
n<-nrow(aircondit)
#将数据框aircondit转化为向量mydata
mydata<-numeric(n)
for(i in 1:n){
  mydata[i]<-aircondit[i,]
}
print(mydata)
#计算lambda的极大似然估计量
lambda.hat<-1/mean(mydata)
print(lambda.hat)
#使用bootstrap对se和bias进行估计
B<-200 #重复试验次数
lambda.B<-numeric(B)
for(b in 1:B){
  i<-sample(1:n,size=n,replace=TRUE)
  hours<-aircondit$hours[i]
  lambda.B[b]<-1/mean(hours)
}
se.lambda<-sd(lambda.B) #该估计的标准误差
bias.lambda<-mean(lambda.B-lambda.hat) #该估计的偏差
print(c(se.lambda,bias.lambda))
#############
#分析：当x～E(lambda)时，lambda.hat=1/mean(x)（直接利用对数似然函数求导可得）
#总结：极大似然估计：已知某个随机样本满足某种概率分布，其中的参数未知，通过若干次试验，利用其结果推出参数的大概值
#结果：该估计的标准误差为0.00396左右，该估计的偏差为0.0010左右





###p193 7.5###
#分别使用标准正态、基本、百分位数和BCa方法计算故障间隔平均时间1/lambda的95%bootstrap置信区间#
library(boot)
data(aircondit,package="boot") #数据为1列12行，列名为hours
theta.boot<-function(dat,ind){
  #theta.boot为统计量的函数
  x<-dat[ind,1]
  mean(x)
}
dat<-cbind(aircondit$hours)
boot.obj<-boot(dat,statistic=theta.boot,R=200)
print(boot.obj)
alpha<-c(.05,.95)
#标准正态bootstrap置信区间#
print(boot.obj$t0+qnorm(alpha)*sd(boot.obj$t))
#基本bootstrap置信区间#
print(2*boot.obj$t0-quantile(boot.obj$t,rev(alpha),type=1))
#百分位数bootstrap置信区间#
print(quantile(boot.obj$t,alpha,type=6))

#BCa置信区间
boot.BCa<-function(x,th0,th,stat,conf=.95){
  #计算BCa置信区间的函数
  #参数：x=数据，stat=统计量函数，conf=分位数
  x<-as.matrix(x)
  n<-nrow(x)
  N<-1:n
  alpha<-(1+c(-conf,conf))/2
  zalpha<-qnorm(alpha)
  
  #偏差修正因子bias corrected
  z0<-qnorm(sum(th<th0)/length(th))
  
  #加速修正因子adjusted for acceleration
  th.jack<-numeric(n)
  for(i in 1:n){
    J<-N[1:(n-1)]
    th.jack[i]<-stat(x[-i,],J)
  }
  L<-mean(th.jack)-th.jack
  a<-sum(L^3)/(6*sum(L^2)^1.5)
  
  #BCa分位数的限制
  adj.alpha<-pnorm(z0+(z0+zalpha)/(1-a*(z0+zalpha)))
  limits<-quantile(th,adj.alpha,type=6)
  return(list("est"=th0,"BCa"=limits))
}
n<-nrow(aircondit)
B<-2000
x<-cbind(aircondit$hours)
theta.b<-numeric(B)
theta.hat<-mean(x)
for(b in 1:B){
  i<-sample(1:n,size=n,replace=TRUE)
  x<-aircondit$hours[i]
  theta.b[b]<-mean(x)
}
stat<-function(dat,index){
  mean(dat[index])
}
boot.BCa(x,th0=theta.hat,th=theta.b,stat=stat)
#############
#总结：（1）使用程序包"boot"中的函数"boot"和"boot.ci"来分别获取正态、基本和百分位数自助法置信区间
#      （2）quantile(x,probs=,na.rm=FALSE)：x=输入的向量或值，probs=0和1之间的值的概率，na.rm=删除NA值
#       (3) 通过生成boot.BCa函数来计算BCa置信区间





###p193 7.7###
#计算theta的样本估计并且使用bootstrap估计theta.hat的偏差bias和标准误差se#
#dat<-scale(scor, TRUE, FALSE) 
#n<-nrow(dat)
#sigma.hat<- (t(dat)%*%(dat))/n
#print(sigma.hat)
#scor数据为5列88行，列名分别为mec.vec.alg.ana.sta（数据框的第i行为第i个学生的成绩集(xi1,..,xi5)）
data(scor,package="bootstrap")
x<-cbind(scor$mec,scor$vec,scor$alg,scor$ana,scor$sta) #x为5列n行矩阵
B<-2000 
n<-nrow(scor)
theta.b<-numeric(B) #利用theta.b储存bootstrap计算得到的第一主成分解释的方差所占的比例
#计算五维成绩数据的5X5经验协方差矩阵sigma.hat（sigma的极大似然估计）
sigma.hat<-matrix(rep(0:0,times=5),5,5) #可以利用scale函数进行中心化和标准化
for(j in 1:5){
  for(k in 1:5){
    s <- 0
    for(i in 1:n){
      s<-s+((x[i,j]-mean(x[,j]))*(x[i,k]-mean(x[,k])))
    }
    sigma.hat[j,k] <- s/n
  }
}
print(sigma.hat)
#计算sigma.hat的特征值（values为特征值，vectors为特征向量）
ev<-eigen(sigma.hat)$val
#计算theta的样本估计
theta.hat<-ev[1]/sum(ev)
print(theta.hat)

#计算bootstrap估计theta的bias和se
for(b in 1:B){
  i<-sample(1:n,size=n,replace=TRUE)
  sigma.b<-matrix(rep(0:0,times=5),5,5)
  mec<-scor$mec[i]
  vec<-scor$vec[i]
  alg<-scor$alg[i]
  ana<-scor$ana[i]
  sta<-scor$sta[i]
  y<-cbind(mec,vec,alg,ana,sta) #y为5列n行的数据
  #生成sigma.hat的bootstrap估计
  for(j in 1:5){
    for(k in 1:5){
      s <- 0
      for(i in 1:n){
        s<-s+((y[i,j]-mean(y[,j]))*(y[i,k]-mean(y[,k])))
      }
      sigma.b[j,k] <- s/n
    }
  }
  ev<-eigen(sigma.b)$val
  theta.b[b]<-ev[1]/sum(ev)
}
bias.b<-mean(theta.b-theta.hat)
se.b<-sd(theta.b)
print(c(bias.b,se.b))
#############
#总结：（1）eigen()：计算特征值与特征向量
#           ev<-eigen(a) 通过ev$val访问特征值，通过ev$vec访问特征向量
#      （2）当总体X（向量）服从多元正态时，sigma的MLE（极大似然估计量）为经验协方差矩阵





###p194 7.8###
#给出数据bootstrap中的scor的theta.hat的偏差bias和标准误差se的jackknife估计
data(scor,package="bootstrap")
dat<-scale(scor, TRUE, FALSE) 
n<-nrow(dat)
sigma.hat<- (t(dat)%*%(dat))/n
print(sigma.hat)
#计算sigma.hat的特征值（values为特征值，vectors为特征向量）
ev<-eigen(sigma.hat)$val #默认为降序排列
#计算theta的样本估计theta.hat
theta.hat<-ev[1]/sum(ev)

#利用jackknife计算bias和se的估计
theta.jack<-numeric(n)
for(i in 1:n){
  #生成sigma.hat的jackknife估计
  mec<-scor$mec[-i]
  vec<-scor$vec[-i]
  alg<-scor$alg[-i]
  ana<-scor$ana[-i]
  sta<-scor$sta[-i]
  y<-cbind(mec,vec,alg,ana,sta)
  dat<-scale(y, TRUE, FALSE) 
  n<-nrow(dat)
  sigma.jack<- (t(dat)%*%(dat))/n
  ev<-eigen(sigma.jack)$val
  theta.jack[i]<-ev[1]/sum(ev)
}
bias.jack<-(n-1)*mean((theta.jack)-theta.hat)
se.jack<-sqrt((n-1)*mean((theta.jack)-mean(theta.jack))^2)
print(c(bias.jack,se.jack))
#############
#注：（1）利用jackknife计算sigma.jakc的过程与7.7计算sigma.b的过程类似
#    （2）计算经验协方差矩阵：将数据中心化后将转置矩阵与原矩阵相乘




###p194 7.10###
#将双对数模型转换为三次多项式模型重复使用“缺一法”（n折）交叉验证来选择最佳拟合模型#
library(DAAG);attach(ironslag) 
#"ironslag"("DAAG")数据中有53个含铁量的测量结果，数据为53行2列，列名分别为chemical和magnetic
a<-seq(10,40,.1)

#线性模型#
L1<-lm(magnetic~chemical)
plot(chemical,magnetic,main="Linear",pch=16)
yhat1<-L1$coef[1]+L1$coef[2]*a
lines(a,yhat1,lwd=2)

#二次多项式模型#
L2<-lm(magnetic~chemical+I(chemical^2))
plot(chemical,magnetic,main="Quadratic",pch=16)
yhat2<-L2$coef[1]+L2$coef[2]*a+L2$coef[3]*a^2
lines(a,yhat2,lwd=2)

#指数模型#
L3<-lm(log(magnetic)~chemical)
plot(chemical,magnetic,main="Exponential",pch=16)
logyhat3<-L3$coef[1]+L3$coef[2]*a
yhat3<-exp(logyhat3)
lines(a,yhat3,lwd=2)

#三次多项式模型#
L4<-lm(magnetic~chemical+I(chemical^2)+I(chemical^3))
plot(chemical,magnetic,main="Cubic polynomial",pch=16)
yhat4<-L4$coef[1]+L4$coef[2]*a+L4$coef[3]*a^2+L4$coef[4]*a^3
lines(a,yhat4,lwd=2)
  
#交叉验证法  
n<-length(magnetic)
e1<-e2<-e3<-e4<-numeric(n) #预测误差
for(k in 1:n){
    y<-magnetic[-k]
    x<-chemical[-k]
    
    #线性模型#
    J1<-lm(y~x)
    yhat1<-J1$coef[1]+J1$coef[2]*chemical[k]
    e1[k]<-magnetic[k]-yhat1
    
    #二次多项式模型#
    J2<-lm(y~x+I(x^2))
    yhat2<-J2$coef[1]+J2$coef[2]*chemical[k]+J2$coef[3]*chemical[k]^2
    e2[k]<-magnetic[k]-yhat2
    
    #指数模型#
    J3<-lm(log(y)~x)
    logyhat3<-J3$coef[1]+J3$coef[2]*chemical[k]
    yhat3<-exp(logyhat3)
    e3[k]<-magnetic[k]-yhat3
    
    #三次多项式模型#
    J4<-lm(y~x+I(x^2)+I(x^3))
    yhat4<-J4$coef[1]+J4$coef[2]*chemical[k]+J4$coef[3]*chemical[k]^2+J4$coef[4]*chemical[k]^3
    e4[k]<-magnetic[k]-yhat4 
}
#利用n折交叉验证得到预测误差估计
print(c(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2)))
#############
#总结：（1）lm(formula,data,subset,weights,method="qr)
#           参数：formula=分析类型(X~Y,X~Y+Z,X~Y*Z) data=数据源 
#           subset=拟合过程中观察的子集 weights=数据拟合的权重
#           lm函数用来拟合线性模型（回归、方差分析、协方差分析）
#      （2）拟合模型的参数可以使用lm返回的coef值
#结论：根据残差可以判断使用二次多项式模型时的残差更小，因此选择使用二次多项式模型
