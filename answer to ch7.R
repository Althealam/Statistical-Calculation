####ch7 作业代码参考答案#####

####7-2####
#参考law的数据，使用基于bootstrap的jackknife方法估计se(R)的bootstrap估计的标准误差#
#算法：（1）通过bootstrap方法得到样本大小为B的样本x1.star,...,xB.star
#      （2）令J(i)表示bootstrap样本中的不包含xi的样本指标，B(i)表示不含xi的样本个数
#      （3）丢掉B-B(i)个含有xi的v样本后其余样本计算一个bootstrap重复
#      （4）得到标准差se估计量的jackknife估计
library(bootstrap)
#初始化
n<-nrow(law)
B<-2000 #bootstrap重复次数
R<-numeric(B)
indices<-matrix(0,nrow=B,ncol=n) #行数目为B，列数目为n
#利用bootstrap计算se(R)
for(b in 1:B){
  i<-sample(1:n,size=n,replace=TRUE) #随机抽取指标
  LSAT<-law$LSAT[i]
  GPA<-law$GPA[i]
  R[b]<-cor(LSAT,GPA)
  indices[b,]<-i
}

#基于bootstrap的jackknife方法
se.jack<-numeric(n)
for(i in 1:n){
  keep<-(1:B)[apply(indices,MARGIN=1,
                    FUN=function(k){!any(k==i)})]
  se.jack[i]<-sd(R[keep])
}
print(sd(R)) #标准差
print(sqrt((n-1)*mean((se.jack-mean(se.jack))^2))) #计算标准误差
#################


####7-3####
#对7-2中的相关性统计量给出一个自助法t置信区间的估计#
#分析：100(1-alpha)%的置信区间可以利用theta.hat与se(theta.hat)表示
boot.t.ci<-function(x,B=500,R=100,level=.95,statistic){
  x<-as.matrix(x)
  n<-nrow(x)
  stat<-numeric(B)
  se<-numeric(B)
  
  boot.se<-function(x,R,fun){
    #利用bootstrap估计f(x)的标准误差
    x<-as.matrix(x)
    m<-nrow(x)
    th<-replicate(R,expr={
      i<-sample(1:m,size=m,replace=TRUE) #利用sample随机抽取指标
      fun(x[i,]) #对抽取的x使用函数fun
    })
    return(sd(th))
  }
  
  for(b in 1:B){
    j<-sample(1:n,size=n,replace=TRUE)
    y<-x[j,]
    stat[b]<-statistic(y)
    se[b]<-boot.se(y,R=R,fun=statistic)
  }
  
  #估计
  stat0<-statistix(x) #原始数据
  t.stats<-(stat-stat0)/se
  se0<-sd(stat)
  alpha<-1-level
  
  Qt<-quantile(t.stats,c(alpha/2,1-alpha/2),type=1)
  names(Qt)<-rev(names(Qt))
  CI<-rev(stat0-Qt*se0)
  return(CI)
}

library(bootstrap)
stat<-function(x){
  cor(x[,1],x[,2])
}
ci<-boot.t.ci(law,statistic=stat,B=2000,R=200)
print(ci)

