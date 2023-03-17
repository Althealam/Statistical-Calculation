#####chapter3 参考答案#####



#####3-3#####
#利用逆变换法模拟一个服从Pareto(2,2)分布的样本#
#分析过程：利用逆变换法由累积分布函数求出逆变换的函数
#          由此模拟出一个服从Pareto(2,2)分布的样本
n<-1000
u<-runif(n) #生成服从均匀分布的n个随机数u
x<-2/sqrt(1-u) #逆变换法生成服从分布函数为F(x)的随机数x
hist(x,prob=TRUE,xlim=c(0,40),main=bquote(f(x)==1-(2/x)^2),col="pink")
#绘制样本密度直方图（使用hist函数）
#xlim与ylim由于设置x轴与y轴的范围，main用于设置直方图标题，col用于设置直方图颜色
y<-seq(2,30,0.01) 
lines(y,8/(y^3),col="steelblue") #绘制理论密度曲线f(x)使用lines函数
############
#结论：经验分布与理论分布基本吻合。



#####3-4#####
#分析过程：
#方法一：使用接受拒绝法生成服从Rayleigh(sigma)分布的随机变量#
#算法：Step1:生成一个随机变量x，使其服从正态分布N(sigma,sigma^2)；
#      Step2:生成一个随机变量u，使其服从U(0,1)
#      Step3:如果u<ex*exp(1/2-x/sigma)/sigma，则令y=x，否则拒绝x，返回步骤2
#方法二：参考hz#
#选取一系列sigma
sigma<-c(1,2,4,8,16,32)
for(i in 1:length(sigma)){
  #设置种子以保证伪随机数的一致性
  set.seed(i)
  #直方图的标题
  title<-c("Histogram of Rayleigh","parameter is",sigma[i])
  #生成两个正态分布
  x<-rnorm(1000,0,sigma[i])
  y<-rnorm(1000,0,sigma[i])
  #生成Rayleigh分布的随机数
  z<-sqrt(x^2+y^2)
  #绘制直方图
  hist(z,prob=TRUE,breaks=seq(0,6*sigma[i],length.out=20),main=title,col="skyblue")
  #绘制Rayleigh密度函数
  x1<-seq(0,6*sigma[i],length.out=100000)
  y1<-(x1/sigma[i]^2)*exp(-(x1^2)/(2*sigma[i]^2))
  lints(x1,y1,col="steelblue")
}
############
#结论：样本密度直方图与密度曲线大致吻合，生成效果较好



#####3-9#####
n<-1000 #随机数个数
y<-rbeta(n,2,2) #令Y=(X+1)/2，Y~Be(2,2)，由Y产生随机数
x<-2*y-1 #将Y产生的随机数结果回代
hist(x,prob=TRUE,main=expression(f(x)==(3/4)(1-x^2)),col="skyblue")
#由此得到直方图
z<-seq(-1,1,0.01) #将结果进行拟合
lines(z,(3/4)*(1-z^2),col="steelblue")

n2<-1000 #通过题目给出的公式定义得到相同的结果图像并进行拟合
u<-vector(mode="numeric",length=1000)
for(i in 1:n2){
  u1<-runif(n2,-1,1)
  u2<-runif(n2,-1,1)
  u3<-runif(n2,-1,1)
  ifelse(abs(u3[i])>=abs(u2[i]) && abs(u3[i])>=abs(u1[i]),u[i]<-u2[i],u[i]<-u3[i])
}
hist(u,prob=TRUE,main=expression(f(x)==(3/4)(1-x^2)))
z<-seq(-1,1,0.01)
lines(z,(3/4)*(1-z^2))
############
#结论：样本密度直方图与密度曲线大致吻合，而且与理论直方图图像较为接近，说明生成效果好。



#####3-11#####
#观察混合变量的经验分布是否是双峰分布#
#分析过程：模拟混合变量的步骤：
#Step1:生成一个随机数u~U(0,1)
#Step2:如果u<p1，生成一个来自N(0,1)分布的x，否则生成一个来自N(3,1)分布的x
f<-function(n,mu1,mu2,sigma1,sigma2,p1){
  x1<-rnorm(n,mu1,sigma1)
  x2<-rnorm(n,mu2,sigma2)
  u<-runif(n)
  k<-as.integer(u<p1)
  x<-k*x1+(1-k)*x2
  return(x)
}
n<-1000
mu1<-0
mu2<-3
sigma1<-1
sigma2<-1
#density()：
#混合变量的概率为0.5时
x0.5=f(n,mu1,mu2,sigma1,sigma2,0.5) 
hist(x0.5,prob=TRUE,xlim=c(-4,6),ylim=c(0,0.23),main="Histogram of x p1=0.5")
lines(density(x0.5))
par(mfrow=c(2,2)) #par函数用于设置或获取图形参数

#混合变量的概率为0.4时
x0.4=f(n,mu1,mu2,sigma1,sigma2,0.4)
hist(x0.4,prob=TRUE,xlim=c(-4,6),ylim=c(0,0.28),main="Histogram of x p1=0.4")
lines(density(x0.4))

#混合变量的概率为0.6时
x0.6=f(n,mu1,mu2,sigma1,sigma2,0.6)
hist(x0.6,prob=TRUE,xlim=c(-4,6),ylim=c(0,0.28),main="Histogram of x p1=0.6")
lines(density(x0.6))

#混合变量的概率为0.25时
x0.25=f(n,mu1,mu2,sigma1,sigma2,0.25)
hist(x0.25,prob=TRUE,xlim=c(-4,6),ylim=c(0,0.28),main="Histogram of x p1=0.25")
lines(density(x0.25))

#混合变量的概率为0.75时
x0.75=f(n,mu1,mu2,sigma1,sigma2,0.75)
hist(x0.75,prob=TRUE,xlim=c(-4,6),ylim=c(0,0.28),main="Histogram of x p1=0.75")
lines(density(x0.75))

par(mfrow=c(1,1))
############
#总结：（1）par函数可以用于实现一页多图
#      （2）mfrow=c(3,5)表示15张图，3行5列
#结论：（1）经验直方图与理论密度曲线拟合程度较好；
#      （2）对于混合变量的经验分布，可以从直方图中很明显的看出来两个混合变量



#####3-13#####
#分析过程：（1）先生成服从Gamma(r,beta)分布的观测值lambda
#          （2）再由y=lambda*exp(-lambda*y)可生成指数-伽马混合变量的观测值
n<-1000
r<-4
bata<-2
lambda<-rgamma(n,r,bata) #rgamma用于生成服从伽马分布的随机数
y<-rexp(n,lambda) #rexp用于生成服从指数分布的随机数，其中参数服从gamma分布
hist(y,prob=TRUE,breaks=20,xlim=c(0,10))
lines(density(y,from=0))
x>0
x<-seq(0,10,0.05)
lines(x,r*beta^r/(beta+x)^(r+1),col="red") #pareto密度曲线
############



#####3-14#####
#分析过程：仿照课本例3.18中的rmvn.choleski函数
rmvn.chol<-function(n,mu,Sigma){
  d<-length(mu)
  Q<-chol(Sigma) #chol函数用于进行cholesky分解
  Z<-matrix(rnorm(n*d),nrow=n,ncol=d)
  X<-Z%*%Q+matrix(mu,n,d,byrow=TRUE) #"%*%"用于矩阵相乘
  return(X)
}
 n<-200
 mu<-0:2
 Sigma<-matrix(c(1,-0.5,0.5,-0.5,1,-0.5,0.5,-0.5,1),nrow=3)
 X<-rmvn.chol(n,mu,Sigma)
 pairs(X)
 #总结：（1）pairs函数：绘制矩阵散点图
 #      （2）矩阵散点图结果分析：可以显示两个数据之间的相关性。当线性相关时可以看出来散点排列在一条直线上下
 #结论：每一对边缘分布的联合分布都大致显示出多元正态分布的椭圆对称性
 #      并且第一个分量与第三个分量呈现正相关，其余两对呈现负相关，这与协方差矩阵吻合
 
 
 #####3-16#####
 #函数norm.std用于将学生的考试成绩样本(X1,X2)(闭卷)和(X3,X4,X5)(开卷)标准化
 norm.std<-function(X)
 {
   n<-nrow(X)
   xbar<-apply(X,2,mean)
   J<-rep(1,n)
   S<-cov(X)
   Q<-chol(S)
   Z<-(X-J%*%t(Xbar))%*%solve(Q)
   return(Z)
 }
 #函数rm.minor用于提出非对角元中微小的量（一般小于10^-10量级）
 rm.minor<-function(X,min=1e-10){
   X[which(abs(X)<min)]<-0
   return(X)
 }
 mydata<-read.csv("scor.csv")
 X<-as.matrix(mydata)
 X12<-X[,1:2]
 Z12<-norm.std(X12) #标准化(X1,X2)
 X345<-X[,3:5]
 Z345<-norm.std(X345) #标准化(X3,X4,X5)
 Z<-cbind(Z12,Z345)
 rm.minor(apply(Z,2,mean))
 rm.minor(cov(Z))

 