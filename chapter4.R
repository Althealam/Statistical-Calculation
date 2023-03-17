######ch4 EM算法作业#####


####第一题：有限正态混合总体####
#（1）使用EM算法估计参数（参考ppt25页）#
set.seed(10000)
N<-2000
mixture.normal<-function(p1,p2,mu1,mu2,mu3,sigma){
  n1<-rbinom(1,size=N,prob=p1)
  n2<-rbinom(1,size=N,prob=p2)
  n3<-N-n1-n2
  #rbinom函数用于生成服从bernoulli分布的随机数
  X1<-rnorm(n1,mu1,sigma)
  X2<-rnorm(n2,mu2,sigma)
  X3<-rnorm(n3,mu3,sigma)
  X<-c(X1,X2,X3)
  X<-X[sample(1:N)] #X为观测数据
  plot(density(X))
  return(X) #输出随机抽样得到的X值（有限正态混合总体）
}
X=mixture.normal(0.1,0.2,1,2,3,3)

#设置初始值
p.hat10<-0.13
p.hat20<-0.23
p.hat30<-1-p.hat10-p.hat20
#利用alpha10、alpha20、alpha30分别表示p1、p2与1-p1-p2的初始值
mu10<-1.2
mu20<-2.2
mu30<-3.2
sigma10<-3.3
para<-c(p.hat10,p.hat20,p.hat30,mu10,mu20,mu30,sigma10)
tol<-1e-100 #容忍度

EM<-function(X,para){ #定义一个进行EM算法的函数EM
  para.old<-para+1 ######这个是为什么？？？？？###########
  for(j in 1:1000){
    ##Estep：
    vp10<-(para[1]*dnorm(X,para[4],para[7]))/(para[1]*dnorm(X,para[4],para[7])
                                              +para[2]*dnorm(X,para[5],para[7])
                                              +para[3]*dnorm(X,para[6],para[7]))
    vp20<-(para[2]*dnorm(X,para[5],para[7]))/(para[1]*dnorm(X,para[4],para[7])
                                              +para[2]*dnorm(X,para[5],para[7])
                                              +para[3]*dnorm(X,para[6],para[7]))
    vp30<-(para[3]*dnorm(X,para[6],para[7]))/(para[1]*dnorm(X,para[4],para[7])
                                              +para[2]*dnorm(X,para[5],para[7])
                                              +para[3]*dnorm(X,para[6],para[7]))
    phi1<-sum(vp10)
    phi2<-sum(vp20)
    phi3<-sum(vp30)
    p.hat11<-phi1/N
    p.hat21<-phi2/N
    p.hat31<-phi3/N
    #p.hat11、p.hat21、p.hat31为Esteps中的期望值的估计值
    phix1<-sum(X*vp10)
    mu11<-phix1/phi1
    phix2<-sum(X*vp20)
    mu21<-phix2/phi2
    phix3<-sum(X*vp30)
    mu31<-phix3/phi3
    phixmu1<-sum((X-mu11)^2*vp10)
    phixmu2<-sum((X-mu21)^2*vp20)
    phixmu3<-sum((X-mu31)^2*vp30)
    sigma11<-sqrt(phixmu1/phi1)
    #注：pi1中1表示进行了一次迭代后计算得到的结果
    para<-c(p.hat11,p.hat21,p.hat31,mu11,mu21,mu31,sigma11)
    ##Mstep:
    if((sqrt(sum(para-para.old)^2))/sqrt(sum(para.old^2))<tol) break 
    #判断是否小于容忍度，如果是的话则迭代结束
    para.old<-para #如果没有达到容忍度，则继续迭代
  }
  return(para)
}
EM(X,para)


#（2）使用方差估计的方法，估计（1）中/各估计的标准差/的估计
#使用bootstrao估计标准差#
B<-100
p1<-p2<-p3<-mu1<-mu2<-mu3<-sigma<-numeric(B)
para.b<-matrix()
para.b<-numeric(B)
for(b in 1:B){
  #生成bootstrap的观测值X.b
  X.b<-X[sample(1:N,size=N,replace=TRUE)]
  #有放回的从x1,...,xn中选出x1.star,...,xn.star，应用EM算法得到估计theta.hat[j]
  #使用EM算法得到估计值theta.hat[j]
  p1[b]<-EM(X.b,para)[1]
  p2[b]<-EM(X.b,para)[2]
  p3[b]<-EM(X.b,para)[3]
  mu1[b]<-EM(X.b,para)[4]
  mu2[b]<-EM(X.b,para)[5]
  mu3[b]<-EM(X.b,para)[6]
  sigma[b]<-EM(X.b,para)[7]
}
se.p1<-sd(p1)
se.p2<-sd(p2)
se.p3<-sd(p3)
se.mu1<-sd(mu1)
se.mu2<-sd(mu2)
se.mu3<-sd(mu3)
se.sigma<-sd(sigma)
print(c(se.p1,se.p2,se.p3,se.mu1,se.mu2,se.mu3,se.sigma))
##################



####第二题：灯泡寿命服从指数分布####
set.seed(2022)
exp.EM<-function(theta,m,n,t){
  #定义进行EM算法估计参数theta的函数exp.EM
  #Y=进行EM算法的样本数据，r=第二批灯泡中t时刻失效的灯泡个数个
  #n=第一批灯泡的个数，m=第二批灯泡的个数，t=失效的时刻
  
  Y0<-rexp(n,theta) #生成n个服从指数分布的随机数Y0
  Y<-Y0[sample(1:n)] #从n个服从指数分布的随机数中随机抽样生成n个灯泡的寿命Y
  Z<-rexp(m,theta)
  r<-sum(Z<=t) #找出失效时间小于等于t的灯泡

  tol<-1e-19 #设置容忍精度
  theta.old<-0.8 #设置初始值
  
  for(k in 1:1000){
    #设置进行EM算法迭代的最大次数为1000
    z1<--(r*t/(exp(theta.old*t)-1)-(m-r)*t-m/theta.old)
    theta.hat<-(m+n)/(sum(Y)+z1)
    if(sqrt((theta.hat-theta.old)^2)/sqrt(theta.old^2)<tol) break 
    #如果theta与theta.old之间的误差在容忍精度内，则跳出循环，迭代结束
    theta.old<-theta.hat #如果上述if语句不成立，则继续迭代
  }
  return(theta.hat)
}
exp.EM(1,30,20,1)

