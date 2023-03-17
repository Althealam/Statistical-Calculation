#####chapter5 参考答案#####


###5-3###
#利用u(0,5)抽样计算蒙特卡洛估计量theta.hat并估计方差，并利用指数分布抽样计算新的蒙特卡洛估计量theta.star#
#1.unirform(0,0.5)
set.seed(10000)
m<-10000
x<-runif(m,0,0.5)
y1<-exp(-x) #计算x属于(0,0.5)时的函数值
theta.hat<-0.5*mean(y1) #利用(a,b)区间上的积分计算法，0.5为b-a的值
var.theta.hat.est<-0.5^var(y1)/m #计算theta.hat的方差估计值（est表示估计值）

#2.exp
n<-10000
y2<-rexp(n,2) #生成服从指数分布的n个随机数
theta.star<-length(which(y2>=0 & y2<0.5))/n #利用示性函数，如果满足y>=0&y<0.5时示性函数的值为1
#相当于hist-miss，如果y2落入区间(0,0.5)内则计为1
var.theta.star.est<-theta.star*(1-theta.star)/n #计算theta.star的方差的估计值

c(theta.hat,theta.star)
c(var.theta.hat.est,var.theta.star.est)
###########
#结论：theta.hat与theta.star接近theta的理论值，当m=n=10000时，var(theta.hat)<<var(theta.star)



###5-6、5-7###
#利用对偶变量法计算cov(exp(U),exp(1-U))和var(exp(U)+exp(1-U))，其中U～u(0,1)#
#分别使用对偶变量法和简单的蒙特卡洛方法来估计theta#
#计算使用对偶变量法得到的方差缩减百分比的经验估计，并和理论值进行比较#
#分析过程：cov(exp(U),exp(1-U))=E(e)-E(exp(U))E(exp(1-U))=e-(e-1)^2
#          var(exp(U))=E(exp(2U))-(E(exp(U)))^2
#1.简单的蒙特卡洛变量法#
m<-10000
x<-runif(m)
y1<-exp(x)
theta.hat<-mean(y1)
var.theta.hat.est<-var(y1)/m

#2.对偶变量法#
n<-10000
R<-1000
theta.star<-numeric(R)
for(i in 1:R){
  u<-runif(n/2)
  v<-1-u
  u<-c(u,v)
  y2<-exp(u)
  theta.star[i]<-mean(y2)
}
var.theta.star.est<-var(theta.star) #计算对偶变量法计算得到的theta.star的方差的估计量
theta.star<-mean(theta.star) #计算利用对偶变量法计算得到的theta.star
c(theta.hat,theta.star)
c(var.theta.hat.est,var.theta.star.est)

#3.计算方差缩减百分比
(var.theta.hat.est-var.theta.star.est)/var.theta.hat.est
###########
#结论：模拟的方差缩减百分比与理论推导的十分接近



###5-13###
#生成函数图像并计算重要抽样法的估计值和方差#
#分析过程：分别取函数正态分布N(2^1/2,1)的密度函数与exp(-x)
#选择这两个函数的原因：支撑在固定的定义域上并且接近于g(x)
g<-function(x){
  (x^2*exp(-x^2/2)/sqrt(2*pi))*(x>1) #设置当x>1时有这个函数成立
}
x<-seq(1.05,5,0.05) #初始值为1.05，结束值为5，间隔为0.05
y<-g(x) #计算当x属于上述区间时的函数值
par(mfrow=c(1,2)) #在同一个页面上生成两张图，1行2列
par(pin=c(2,1))

#绘制经验密度直方图
plot(x,y,type='l',ylim=c(0,0.5),ylab="") #type设置曲线类型(散点、虚线、...)，ylim设置y轴范围，ylab设置y轴的标签
lines(x,dnorm(x,mean=sqrt(2),1),lty=2) #绘制正态分布N(2^1/2,1)的理论密度曲线
legend(x='topright',legend=c(expression(g),expression(f_1),expression(f_2)),lty=c(1,2,3),bty='n',cex=0.8)
plot(x,y/dnorm(x,mean=sqrt(2),1),type='l',lty=2,ylim=c(0,2),ylab="")
lines(x,g(x)/exp(-x),lty=3)
legend(x='topright',legend=c(expression(g/f_1),expression(g/f_2)),lty=c(2,3),bty='n',cex=0.8)

m<-10000
theta.hat<-se<-numeric(2)

x<-rnorm(m.mean=sqrt(2),1) #f_1
i<-c(which(x>1),which(x<0)) #满足条件x>1或者x<0的i值
x[i]<-2
g.f<-g(x)/dnorm(x,mean=sqrt(2),1) #g/f的值，其中dnorm(x,mean=sqrt(2),1)
theta.hat[1]<-mean(g.f)
se[1]<-sd(g.f)

x<-rexp(m,1) #f_2
g.f<-g(x)/exp(-x) #g/f
theta.hat[2]<-mean(g.f)
se[2]<-sd(g.f)
###########
#总结：R中绘制图形函数线条属性的参数：
#     （1）lty=指定线条类型；（2）lwd=指定线条宽度；（3）pch=指定线条的点形状；（4）cex=指定线条的连接点形状大小；



#5-15#
#得到5-13中的分层重要抽样估计并和5-10中的结果进行比较#
g<-function(x){ #计算g(x)
  exp(-x)/(1+x^2)
}

n<-5
f<-function(x,n,j){
  exp(-x)/(exp((1-j)/n)-exp(-j/n))
}

m<-10000
theta.hat<-se<-numeric(n)
stra.est<-function(m,n,j){
  u<-runif(m)
  x<-j/n-log(u+exp(1/n)*(1-u))
  g.f<-g(x)/f(x,n,j)
  theta.hat<-mean(g.f)
  se<-sd(g.f)
  return(c(theta.hat,se))
}

stra.est1<-function(j){
  stra.est(m,n,j)
}

j<-seq(1,n)
stra.result<-as.data.frame(lapply(j,stra.est1),row.names=c("theta.i_est","sd_i"),col.names=seq(1,5))
apply(stra.result,1,sum)
###########
#结论：运行结果部分的theta_i.est与例5.10中的结果很接近，并且标准差sd_i要更小。
