#####ch4 EM算法作业参考答案#####

###10-1###
#（1）有限正态混合总体估计参数#
#分析：theta=(p,mu,sigma^2)为待估参数，其中p=(p1,p2,p3)，mu=(mu1,mu2,mu3)
set.seed(101)
N<-1000
n1<-rbinom(1,size=N,prob=0.1)
n2<-rbinom(1,size=N,prob=0.3)
n3<-N-n1-n2
X1<-rnorm(n1,-5,3)
X2<-rnorm(n2,3,1)
X3<-rnorm(n3,10,4)
X0<-c(X1,X2,X3)
X0<-X0[sample(1:N)] #从X中抽取N的数字
plot(density(X0),col="red") #绘制混合正态分布的图像

EM<-function(X){
  #设定迭代的初始值
  p10<-0.1
  p20<-0.3
  p30<-0.6
  mu10<-5.2
  sigma10<-3.1
  mu20<-3.3
  sigma20<-1.2
  mu30<-9.3
  sigma30<-4.4
  para<-c(p10,p20,p30,mu10,sigma10,mu20,sigma20,mu30,sigma30)
  #设定容忍精度
  tol<-1e-8 
  para.old<-para+1
  
  for(j in 1:1000){
    vp10<-(para[1]*dnorm(X,para[4],para[5]))/(para[1]*dnorm(X,para[4],para[5])+para[2]*dnorm(X,para[6],para[7])+para[3]*dnorm(X,para[8],para[9]))
    vp20<-(para[2]*dnorm(X,para[6],para[7]))/(para[1]*dnorm(X,para[4],para[5])+para[2]*dnorm(X,para[6],para[7])+para[3]*dnorm(X,para[8],para[9]))
    vp30<-(para[3]*dnorm(X,para[8],para[9]))/(para[1]*dnorm(X,para[4],para[5])+para[2]*dnorm(X,para[6],para[7])+para[3]*dnorm(X,para[8],para[9]))
  
    phi1<-sum(vp10)
    phi2<-sum(vp20)
    phi3<-sum(vp30)
    p11<-phi1/N
    p21<-phi2/N
    p31<-phi3/N
    phix1<-sum(X*vp10)
    mu11<-phix1/phi1
    phix2<-sum(X*vp20)
    mu21<-phix2/phi2
    phix3<-sum(X*vp30)
    mu31<-phix3/phi3
    phixmu1<-sum((X-mu11)^2*vp10)
    sigma11<-sqrt(phixmu1/phi1)
    phixmu2<-sum((X-mu21)^2*vp20)
    sigma21<-sqrt(phixmu2/phi2)
    phixmu3<-sum((X-mu31)^2*vp30)
    sigma31<-sqrt(phixmu3/phi3)
    
    para<-c(p11,p21,p31,mu11,sigma11,mu21,sigma21,mu31,sigma31)
    
    if((sqrt(sum((para-para.old)^2))/sqrt(sum(para.old^2)))<tol) #判断是否达到容忍度
    {break} #如果达到容忍度，则迭代结束
    para.old<-para #如果没有达到容忍度，则继续迭代
  }
  return(list(estimate=para,iter=j,tol=tol))
}
print(list(estimate=para,iter=j,tol=tol))

#（2）利用方差估计的某个方法估计各个估计的标准差的估计#
B<-50 #设定bootstrap重复实验次数
n<-N
theta.b<-matrix(0,nrow=B,ncol=9)
theta.b[1,]<-EM(X0)$estimate
for(b in 2:B){
  Xnew<-X0[i]
  theta.b[b,]<-EM(Xnew)$estimate
}
print(theta.b)
sqrt(diag(cov(theta.b)))
#######################



####10-2####
#估计灯泡寿命的参数theta#
set.seed(101)
N<-1000
n<-600
m<-400
t<-5
r<-160
#设定样本服从theta=10的指数分布
x1<-rexp(n,1/10)
x2<-rexp(m,1/10)
x0<-c(x1,x2)
x0<-x0[sample(1:N)]
plot(density(x0))

theta0<-9
para<-theta0
tol<-1e-8
para.old<-para+1

for(j in 1:1000){
  ht<-exp(-t*para)/(1-exp(-t*para))
  theta1<-(m+n)/(sum(x1)+(m-r)*(t+1/para)+r*(1/para-t*ht))
  para<-theta1
  
  if(sqrt(sum((para-para.old)^2))/sqrt(sum(para.old^2))<tol)
  {break}
  para.old<-para
}
list(estimate=para,iter=j,tol=tol)
#######################