###chapter 4书本代码###

####例1:椒花蛾####
moth<-function(p,n.obs){
  #p=(pc,pt) 等位基因的概率
  #n.obs=(nc,ni,nt) 可观测到的数据
  n<-sum(n.obs)
  nc<-n.obs[1] #黑色可观测数据
  ni<-n.obs[2] #中间可观测数据
  nt<-n.obs[3] #浅色可观测数据
  ntt<-nt #ntt是不可观测数据，nt为可观测数据，且ntt=nt
  
  cat(p,"\n") #cat函数用于打印字符串
  pct<-pit<-ptt<-rep(0,20) #不可观测等位基因的概率，先将其设置为0向量
  #进行EM迭代#
  pct[1]<-p[1]
  pit[1]<-p[2]
  ptt[1]<-1-p[1]-p[2]
  for(i in 2:20){#迭代次数
    
    pc.old<-pct[i-1]
    pi.old<-pit[i-1]
    pt.old<-ptt[i-1]
    #注：pc.pld,pi.old,pt.old为最初设定的值，同时为上一次迭代时的最大值点
    
    den<-pc.old^2+2*pc.old*pi.old+2*pc.old*pt.old
    
    #利用条件期望（参考ppt的p7）
    ##E步：计算期望值##
    ncc<-nc*pc.old^2/den #ncc.i-1
    nci<-2*nc*pc.old*pi.old/den #nci.i-1
    nii<-ni*pi.old^2/(pi.old^2+2*pi.old*pt.old) #nii.i-1
    nit<-2*ni*pi.old*pt.old/(pi.old^2+2*pi.old*pt.old) #nit.i-1
    nct<-2*nc*pc.old*pt.old/den #nct.i-1
    
    ##M步：关于theta的最大化##
    #利用迭代公式（参考ppt的p8）
    #Q(theta|theta.t)为在观测到X=x，以及在参数theta=theta.t的条件下完全数据的联合对数似然函数的期望
    #对该期望求导则可以得到以下迭代公式
    pct[i]<-(2*ncc+nci+nct)/(2*n)
    pit[i]<-(2*nii+nit+nci)/(2*n)
    ptt[i]<-(2*ntt+nct+nit)/(2*n)
    #注：pct[i],pit[i],ptt[i]为第i次迭代时的最大值点
  }
  return(list(pct=pct,pit=pit,ptt=ptt))
}

#导入数据
n.obs<-c(85,196,341)
p<-c(1/3,1/3) #假定p的初始值（EM算法的初始步骤都需要自己假定初值）
a<-moth(p,n.obs)
pct<-a$pct
pit<-a$pit
ptt<-a$ptt

rcc=sqrt((diff(pct)^2+diff(pit)^2)/(pct[-20]^2+pit[-20]^2))
rcc=c(0,rcc)

d1=(pct[-1]-pct[20])/(pct[-20]-pct[20])
d1=c(d1,0)
d2=(pit[-1]-pit[20])/(pit[-20]-pit[20])
d2=c(d2,0)

print(cbind(pct,pit,rcc,d1,d2)[1:9,],digits=5)
#############
#总结EM算法：
#（1）E步：计算Q(theta|theta.t)=E[L(theta|Y)|x,theta.t]
#（2）M步：关于theta的最大化Q(theta|theta.t)，并记theta.t+1表示此时的最大值点
#（3）返回到E步，直至收敛准则达到



####例2:Bayes后验众数####
theta=0.3
n<-50
x<-drop(rmultinom(1,n,c((2+theta)/4,(1-theta)/2,theta/4)))
a<-.01
b<-.01

th<-rep(0,20)
th[1]<-0.2
for(t in 2:20){
  num<-x[1]*th[t-1]/(2+th[t-1])+x[3]+a-1
  den<-x[1]*th[t-1]/(2+th[t-1])+x[2]+x[3]+a+b-2
  th[t]<-num/den
}
rcc<-sqrt((diff(th)^2+diff(th)^2)/(th[-20]^2+th[-20]^2))
print(cbind(th,c(0,rcc)),digits=5)
#############



####例3:混合分布中参数的估计问题####
set.seed(100)
N<-2000
n1<-rbinom(1,size=N,prob=0.3) #生成概率为0.3，数量为1，试验次数为N的服从二项分布的随机数
n2<-N-n1 #利用二项分布的性质可以知道n1+n2=N
X1<-rnorm(n1,1,3) #生成数量为n1，均值为1，方差为3的服从生态分布的随机数
X2<-rnorm(n2,10,1) #生成数量为n2，均值为10，方差为1的服从正态分布的随机数
X<-c(X1,X2)
X<-X[sample(1:N)] #利用[]算子，从X中获得自己想要的数据（从X中抽样N次）
plot(density(X))

#设定初始值
alpha10<-0.4
alpha20<-0.6
mu10<-1.2
sigma10<-3.1
mu20<-9.3
sigma20<-1.2
para<-c(alpha10,alpha20,mu10,sigma10,mu20,sigma20) #待估参数的初始值
tol<-1e-20 #用于判断终止条件的容忍限度
para.old<-para+1

for(j in 1:1000){
  #可以使用count计算EM算法总共进行了几次
  #count<-count+1
  vp10<-(para[1]*dnorm(X,para[3],para[4]))/(para[1]*dnorm(X,para[3],para[4])+para[2]*dnorm(X,para[5],para[6]))
  vp20<-(para[2]*dnrom(X,para[5],para[6]))/(para[1]*dnorm(X,para[3],para[4])+para[2]*dnorm(X,para[5],para[6]))
  #dnorm返回标准正态分布密度函数在某点处的函数值（即代入具体的自变量值，计算因变量的值）
  phi1<-sum(vp10)
  phi2<-sum(vp20)
  alpha11<-phi1/N
  alpha21<-phi2/N
  phix1<-sum(X*vp10)
  mu11<-phix1/phi1
  phix2<-sum(X*vp20)
  mu21<-sum((X-mu11)^2*vp10)
  phixmu1<-sum((X-mu11)^2*vp10)
  sigma11<-sqrt(phixmu1/phi1)
  phixmu2<-sum((X-mu21)^2*vp20)
  sigma21<-sqrt(phixmu2/phi2)
  
  para<-c(alpha11,alpha21,mu11,sigma11,mu21,sigma21)
  
  if(sqrt(sum((para-para.old)^2))/sqrt(sum(para.old^2))<tol) break 
  #终止条件，如果新参数与旧参数的相对值小于容忍限度时则结束
  #sqrt(sum((para-para.old)^2))与sqrt(sum(para.old^2))分别表示欧式范数
  para.old<-para #如果不满足终止条件，则继续EM算法，将新的参数值赋于para.old
}
print(list(estimate=para,iter=j,tol=tol))
#总结：（1）二项分布：在n次独立重复的Bernoulli试验中，设每次试验中事件A发生的概率为p，用X表示n重Bernoulli试验中事件A发生的次数
#      （2）rbinom(n,size,prob) n=观察的次数，size=试验次数，prob=每次试验成功的概率
#           从给定样本生成给定概率的所需数量的随机数
#      （3）dnorm 输入x轴上的数值，输出的是该点正态分布概率密度的函数值
#           eg:dnorm(z)：标准正态分布密度函数f(x)在x=z处的函数值

