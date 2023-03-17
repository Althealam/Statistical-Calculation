###p169 7.2###
#标准误差的自助法估计#
library(bootstrap) #程序包bootstrap中有法学院的数据集law
print(cor(law$LSAT,law$GPA)) #cor生成LSAT与GPA之间的相关系数
B<-200 #重复试验的次数
n<-nrow(law) #样本的容量大小
R<-numeric(B) #生成与重复试验次数大小相等的相关系数向量
#估计总体的相关系数需要使用重抽样计算得到的相关系数
for(b in 1:B){ #对第b次bootstrap重复试验
  i<-sample(1:n,size=n,replace=TRUE) #模拟抽样过程，随机抽取下标i
  LSAT<-law$LSAT[i] #找到该样本对应的LSAT
  GPA<-law$GPA[i] #找到该样本对应的GPA
  R[b]<-cor(LSAT,GPA) #计算第b个样本对应的相关系数
}
print(se.R<-sd(R)) 
hist(R,prob=TRUE) #画出关于相关系数的直方图
##############
#总结：step1:对第b次bootstrap重复试验（b=1,...,B）
#           （1）从观测样本x1,...,xn中有放回地抽样生成样本x.star(b)=(x.star[1],x.star[2],...,x.star[n])
#           （2）对第b个bootstrap样本计算第b次重复试验的统计量theta.hat[b]，其中theta.hat为(LSAT,GPA)之间的相关性R
#      step2:se(R)的bootstrap估计为重复试验的样本标准差



###p170 7.3###
#标准误差的自助法估计：boot函数
library(bootstrap)
r<-function(x,i){
  #返回第b次抽样的参数的函数
  #第一个参数为样本数据，第二个参数为指标向量
  #数据为x，指标为i，令x[i,1]取为第一个重抽样变量，x[i,2]为第二个重抽样变量
  cor(x[i,1],x[i,2])
  #x[i,1]取为第一个重抽样变量，x[i,2]取为第二个重抽样变量
}
library(boot)
obj<-boot(data=law,statistic=r,R=2000) #将结果存储在向量obj中
print(obj)
#boot(data,statistic,R,...) statistic= function(dat,ind)
##############



###p172 7.4###
#偏差的自助法估计bias(theta)#
library(bootstrap)
theta.hat<-cor(law$LSAT,law$GPA) #将每一组向量的相关系数看作参数
B<-2000 #重复试验的次数
n<-nrow(law) #2列15行的数据，列名分别为LSAT与GPA
theta.b<-numeric(B) #生成与重复试验次数大小一致的向量theta.b 
for(b in 1:B){ #对第b次bootstrao重复试验
  i<-sample(1:n,size=n,replace=TRUE) #模拟抽样过程，从观测样本中有放回的抽样生成样本
  LSAT<-law$LSAT[i]
  GPA<-law$GPA[i]
  theta.b[b]<-cor(LSAT,GPA) #将每次重复试验的相关系数保存在向量theta.b中
}
bias<-mean(theta.b-theta.hat)
#theta.hat为无偏的，即E(theta.hat)=theta.hat
##############
#总结：偏差bias的bootstrap估计使用theta.hat的bootstrap重复试验来估计theta.hat的抽样分布



###p172 7.5###
data(patch,package="bootstrap") #每一行为每一位病人的数据
#patch数据包含了一种医用贴片之后8个病人血液中某种激素的测量结果
#8行6列的数据，列名分别为subject.placebo.oldpatch.newpatch.z.y
n<-nrow(patch)
B<-2000 #重复试验次数
theta.b<-numeric(B) 
theta.hat<-mean(patch$y)/mean(patch$z) #设置我们感兴趣的参数（书本p172的公式）
#bootstrap#
for(b in 1:B){
  i<-sample(1:n,size=n,replace=TRUE) #模拟抽样过程，随机抽取下标i
  #取出我们需要的指标
  y<-patch$y[i] 
  z<-patch$z[i]
  theta.b[b]<-mean(y)/mean(z) #计算第b次重复试验的统计量
}
bias<-mean(theta.b)-theta.hat #计算偏差
se<-sd(theta.b)
print(list(est=theta.hat,bias=bias,se=se,cv=bias/se))
##############



###p174 缺一法###
#去掉向量的第i个元素的方法，常常用于水手刀法失效时（偏差的水手刀法估计）#
x<-1:5
for(i in 1:5)
  print(x[-i])
##############



###p175 7.6###
#对patch中的数据计算偏差的水手刀法估计#
data(patch,package="bootstrap")
#patch数据包含了一种医用贴片之后8个病人血液中某种激素的测量结果
#8行6列的数据，列名分别为subject.placebo.oldpatch.newpatch.z.y
n<-nrow(patch)
y<-patch$y
z<-patch$z
theta.hat<-mean(y)/mean(z)
print(theta.hat)
#利用缺一法计算偏差的水手刀法估计
theta.jack<-numeric(n)
for(i in 1:n){
  theta.jack[i]<-mean(y[-i])/mean(z[-i])}
  #常常利用算子[]来去掉第i个样本
bias<-(n-1)*(mean(theta.jack)-theta.hat)
print(bias)
##############



###p176 7.7###
#对patch中的数据计算se的jackknife法估计#
data(patch,package="bootstrap")
n<-nrow(patch)
y<-patch$y
z<-patch$z
theta.hat<-mean(y)/mean(z)
print(theta.hat)
theta.jack<-numeric(n)
for(i in 1:n){ #去掉不同的样本计算jackknife统计量(标准误差)
  theta.jack[i]<-mean(y[-i])/mean(z[-i]) #利用y[-i]或者z[-i]来表示去掉样本中第i个样本
}
se<-sqrt((n-1)*mean((theta.jack-mean(theta.jack))^2)) #书本p175公式(7.4)
print(se)
##############
#总结：[]算子给出了去掉向量第i个元素的方法



###p176 7.8###
#对1,2,...,100中的1-个整数构成的随机样本计算中位数的标准误差的水手刀法的估计#
#jackknife失效的情况
n<-10
x<-sample(1:100,size=n)
M<-numeric(n)
#缺一法#
for(i in 1:n){
  y<-x[-i]
  M[i]<-median(y)
}
Mbar<-mean(M)
print(sqrt((n-1)/n*sum((M-Mbar)^2)))
Mb<-replicate(1000,expr={
  y<-sample(x,size=n,replace=TRUE)
  median(y)
})
print(sd(Mb))
##############
#总结：在统计量不光滑的情况下可以使用弃d-jackknife法（每次重复试验时去掉d个观测值）



###p178 7.9###
#对patch数据使用基于bootstrap法的jackknife法估计标准误差
#去掉第i个观测值，估计标准误差se的算法是对每个i从剩余的n-1个观测值中重复抽样生成B次重复试验
data(patch,package="bootstrap")
n<-nrow(patch)
y<-patch$y
z<-patch$z
B<-2000
theta.b<-numeric(B) #bootstrap统计向量
indices<-matrix(0,nrow=B,ncol=n) #生成关于i的矩阵

#计算bootstrap估计的统计量
for(b in 1:B){
  i<-sample(1:n,size=n,replace=TRUE) #随机抽取下标i，每次抽取10个
  #将抽取的下标i所对应的元素值保存在y/z中
  y<-patch$y[i] 
  z<-patch$z[i]
  theta.b[b]<-mean(y)/mean(z) #计算第b次重复试验的bootstrap的统计量theta.b
  indices[b,]<-i #保存抽取的样本的下标
}

#计算se的基于bootstrap法的jackknife法的估计
se.jack<-numeric(n)
for(i in 1:n){
  keep<-(1:B)[apply(indices,MARGIN=1,FUN=function(k){!any(k==i)})] 
  #apply函数用于遍历数组中的所有维度
  #apply(X, MARGIN, FUN, ...) X：有维度的数据对象 MARGIN：=1时对行进行运算，=0时对列进行运算 FUN：待运用的函数
  #any函数：用于判断逻辑向量（至少一个成员为真时返回True，否则返回False）
  #all函数：用于判断逻辑向量（当全部成员都为真时返回True，否则返回False）
  se.jack[i]<-sd(theta.b[keep])
}

print(sd(theta.b)) #标准误差的自助法估计
print(sqrt((n-1)*mean((se.jack-mean(se.jack))^2))) #标准误差的基于自助法的水手刀法估计
##############



###p181 7.10###
#使用程序包boot中的函数boot和boot.ci来分别获取正态、基本和百分位数自助法的置信区间
library(boot)
data(patch,package="bootstrap")

theta.boot<-function(dat,ind){
  #进行重抽样的方法与进行bootstrap抽样的方法不同：
  #（自助法）重抽样：dat[int,] 
  # bootstrap抽样：boot(data,statistic,R,..)
  y<-dat[ind,1]
  z<-dat[ind,2]
  mean(y)/mean(z)
}

#运用自助法对生物等效性比计算置信区间估计
y<-patch$y
z<-patch$z
dat<-cbind(y,z)
#cbind 对矩阵进行列的合并（左右合并）
#lbind 对矩阵进行行的合并（上下合并）
boot.obj<-boot(dat,statistic=theta.boot,R=2000)
#boot(data,statistic,R,...) 从data中抽取bootstrap样本

#计算几类自助法置信区间，并与boot.ci的结果进行比较
alpha<-c(.025,.975)
#标准正态自助法置信区间#
print(boot.obj$t0+qnorm(alpha)*sd(boot.obj$t)) #qnorm为生成正态分布分位数函数
#基本自助法置信区间#
print(2*boot.obj$t0-quantile(boot.obj$t,rev(alpha),type=1)) #quantile为分位数函数 
#百分位数自助法置信区间#
print(quantile(boot.obj$t,alpha,type=6))
##############
#总结：(1)boot(data,statistic,R):data传递原数据，statistic传递统计量形式（可以为自定义函数），R传递bootstrap重复次数B
#      boot函数的输出参数为boot对象
#     （2）boot.ci():可以直接作用于这一类型的变量并且输出区间估计


###p183 7.11###
#对law数据的相关性统计量计算95%自助法置信区间估计#
library(boot)
data(law,package="bootstrap")
boot.obj<-boot(law,R=2000,statistic=function(x,i){cor(x[i,1],x[i,2])})
print(boot.ci(boot.obj,type=c("basic","norm","perc")))
##############



###p184 7.12###
#对一元或者多元样本计算自助法t置信区间的函数#
boot.t.ci<-function(x,B=500,R=100,level=.95,statistic){
  #自定义函数boot.t.ci(x,B,R,level,statistic)
  #x为样本数据 B为自助法重复试验次数 R为估计标准误差的重复试验次数 statistic为计算统计量的函数
  x<-as.matrix(x)
  n<-nrow(x)
  stat<-numeric(B) #stat=(theta.hat[1],...,theta.hat[n])
  se<-numeric(B) #se=(se.hat[1],...,se.hat[n])
  
  boot.se<-function(x,R,f){
    #x为样本数据，statistic为计算统计量的函数
    #统计f(x)的标准误差估计
    x<-as.matrix(x)
    m<-nrow(x)
    th<-replicate(R,expr={
      #bootstrap随机抽样
      i<-sample(1:m,size=m,replace=TRUE)
      f(x[i.]) #计算bootstrap统计量
    })
    return(sd(th)) #计算标准误差bootstrap估计
  }

for(b in 1:B){ #对第b次重复试验(b=1,...B)
  j<-sample(1:n,size=n,replace=TRUE)
  y<-x[j,] #从x中有放回地抽样生成第b个样本x(b)=(x1(b),...,xn(b))
  stat[b]<-statistic(y) #对第b个样本x(b)计算theta.hat[b]
  se[b]<-boot.se(y,R=R,f=statistic) #计算或者估计标准误差se(theta.hat[b])
}

stat0<-statistic(x) #stat0为观测统计量theta.hat
t.stats<-(stat-stat0)/se #计算t统计量的第b次重复试验t(b),stat为theta.hat[b]
se0<-sd(stat) 
alpha<-1-level
Qt<-quantile(t.stats,c(alpha/2,1-alpha/2),type=1) #计算置信限
names(Qt)<-rev(names(Qt))
CI<-rev(stat0-Qt*se0)
}
##############
#总结：(1)as.matrix():将一个数据框或者数据表转换为矩阵
#     （2）rev():向量逆转



###p185 7.13###
#对比率统计量计算95%的自助法t置信区间#
#应用boot.t.ci函数
dat<-cbind(patch$y,patch$z)
stat<-function(dat){
  mean(dat[,1])/mean(dat[,2])
}

boot.t.ci<-function(x,B=500,R=100,level=.95,statistic){
  #自定义函数boot.t.ci(x,B,R,level,statistic)
  #x为样本数据 B为自助法重复试验次数 R为估计标准误差的重复试验次数 statistic为计算统计量的函数
  x<-as.matrix(x)
  n<-nrow(x)
  stat<-numeric(B)
  se<-numeric(B) 
  
  boot.se<-function(x,R,f){
    #注：boot.se为局部函数
    x<-as.matrix(x)
    m<-nrow(x)
    th<-replicate(R,expr={
      i<-sample(1:m,size=m,replace=TRUE)
      f(x[i.])
    })
    return(sd(th))
  }
}

ci<-boot.t.ci(dat,statistic=stat,B=2000,R=200)
print(ci)
##############



###p186 7.14###
#BCa自助法置信区间#
boot.BCa<-function(x,th0,th,stat,conf=.95){
  #基于自助法的BCa自助法置信区间
  #th0是观察到的统计数据
  #stat是计算统计量的函数
  
  x<-as.matrix(x)
  n<-nrow(x)
  N<-1:n
  alpha<-(1+c(-conf,conf))/2
  zalpha<-qnorm(alpha) #z分位数
  
  #偏差修正因子
  z0<-qnorm(sum(th<th0)/length(th))
  #加速因子
  th.jack<-numeric(n)
  for(i in 1:n){
    J<-N[1:(n-1)]
    th.jack[i]<-stat(x[-i.],J)
  }
  L<-mean(th.jack)-th.jack
  a<-sum(L^3)/(6*sum(L^2)^1.5)
  
  #BCa置信区间限制
  adj.alpha<-qnorm(z0+(z0+zalpha)/(1-a*(z0+zalpha)))
  limits<-quantile(th,adj.alpha,type=6)
  return(list("est"=th0,"BCa"=limits))
}
##############
#总结：（1）quantile()：quantile(x, probs = , na.rm = FALSE)
#        X=输入向量或值 Probs=0和1之间的值的概率 na.rm=删除NA值



###p187 7.15###
#使用函数"boot.BCa"计算生物等效性比率统计量的BCa置信区间#
n<-nrow(patch)
B<-2000
y<-patch$y
z<-patch$z
x<-cbind(y,z)
theta.b<-numeric(B)
theta.hat<-mean(y)/mean(z)

#自助法
for(b in 1:B){
  i<-sample(1:n.size=n,replace=TRUE)
  y<-patch$y[i]
  z<-patch$z[i]
  theta.b[b]<-mean(y)/mean(z)
}
stat<-function(dat,index){
  mean(dat[index,1])/mean(dat[index,2])
}
boot.BCa(x,th0=theta.hat,th=theta.b,stat=stat)
##############



###p188 7.16###
#使用"boot.ci"的BCa自助法置信区间
boot.obj<-boot(x,statistic=stat,R=2000)
boot.ci(boot.obj,type=c("perc","bca"))
##############



###p190 7.17###
#ironslag("DAAG")数据中有53个含铁量的测量结果
#4个模型+数据的预测响应图
library(DAAG)
attach(ironslag)
a<-seq(10,40,.1)

#线性模型#
L1<-lm(magnetic~chemical)
plot(chemical,magnetic,main="Linear",pch=16)
yhat1<-L1$coef[1]+L1$coef[2]*a
lines(a,yhat1,lwd=2)

#二次模型#
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

#双对数模型#
L4<-lm(log(magnetic)~log(chemical))
plot(log(chemical),log(magnetic),main="Log-Log",pch=16)
logyhat4<-L4$coef[1]+L4$coef[2]*log(a)
lines(log(a),logyhat4,lwd=2)
##############



###p190 7.18###
n<-length(magnetic)
e1<-e2<-e3<-e4<-numeric(n)

for(k in 1:n){
  y<-magnetic[-k]
  x<-chemical[-k]
  
  J1<-lm(y~x)
  yhat1<-J1$coef[1]+J1$coef[2]*chemical[k]
  e1[k]<-magnetic[k]-yhat1
  
  J2<-lm(y~x+I(x^2))
  yhat2<-J2$coef[1]+J2$coef[2]*chemical[k]+J2$coef[3]*chemical[k]^2
  e2[k]<-magnetic[k]-yhat2
  
  J3<-lm(log(y)~x)
  logyhat3<-J3$coef[1]+J3$coef[2]*chemical[k]
  yhat3<-exp(logyhat3)
  e3[k]<-magnetic[k]-yhat3
  
  J4<-lm(log(y)~log(x))
  logyhat4<-J4$coef[1]+J4$coef[2]*log(chemical[k])
  yhat4<-exp(logyhat4)
  e4[k]<-magnetic[k]-yhat4
}

c(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2))

#模型2的残差图
par(mfrow=c(1,2))
plot(L2$fit,L2$res)
abline(0,0)
qqnorm(L2$res) 
qqline(L2$res)
par(mfrow=c(1,1))
##############
#总结：（1）QQ图（quantile-quantile plot）：用于判断某一系列的值是否符合正态分布
#      （2）qqnorm：用于判断数据是否偏离正态分布的趋势线距离
#      （3）qqline：横坐标为标准正态分布的分位数，纵坐标的输入数据的分位数，用于检验数据是否符合正态分布