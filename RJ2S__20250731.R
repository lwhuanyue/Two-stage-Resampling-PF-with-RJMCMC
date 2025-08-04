# Code for 2 stage resampling method with RJMCMC
# Version 2024 Oct
# 
library(stats)
require(ReacTran)
require(fields)
library(nleqslv)
# 01.Pre functions --------------------------------------------------------
taper <- function(dist) Wendland(dist,scale*2,1,1)
#For plotting the bands.
plot.PI=function(ens,color,probs=c(.1,.9),xgrid=1:nrow(ens)) {
  quants=apply(ens,1,quantile,probs)
  polygon(c(xgrid,rev(xgrid)),c(quants[1,],rev(quants[2,])),col=adjustcolor(color,alpha.f=0.3),border=NA)
}
#Func: 1-d advection function
model2 <- function (t, C, parms,...) {
  # Parameter:
  #   K
  #   c.x
  #   c.v
  #   DIS
  K   <- parms[1]
  c.x <- parms[1+(1:(K+2))] #
  c.v <- parms[K+3+(1:(K+1))]
  p.dis <- parms[2*K+5]
  v1 <- rep(0,n)
  for (k in 1:(K+1)){
    v1[(c.x[k]+1):c.x[k+1]] <- rep(c.v[k],length((c.x[k]+1):c.x[k+1]))
  }
  smw <- 30 # smooth window
  if(K != 0 ){
    for (k in  1:K){
      xa <- max(1,c.x[k+1]-smw)
      xb <- min(n,c.x[k+1]+smw)
      va <- c.v[k]
      vb <- c.v[k+1]
      v1[xa:xb] <- va + (vb-va)*(xa:xb-xa)/(xb-xa)
    }
  }
  dC <- advection.1D(C=C, v = v1, dx = dx, C.up = C[n],
                     C.down = C[1], ...)$dC  +
    tran.1D(C = C,  D = p.dis, C.up = C[n],
            dx = dx)$dC
  return(list(dC))
}

# ~~~~~~~~~~~~~~~~~  The settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


set.seed(20250731)
lambda_p <- 2 # prior mean
kmax <- 3
N <- 80
times <- 0:600
sig.eps0 <- 0.2 #2

L <- 400 #
n <- L+1
#Gamma 
g_alpha <- 0.4
g_beta  <- 0.95
kalmangain <- 1; plot_t = 0
dx = 0.2
G = 1
# 02. The prior function --------------------------------------------------

#2.1 --- The prior of the model index. ---
# 'the true model following the k-th model' means there are k points, 
#      dividing the space into k+1 parts.


#用sample做截断泊松
dpois(c(0:kmax), lambda = lambda_p) #
sample(c(0:kmax),N,replace=TRUE,dpois(c(0:kmax), lambda = lambda_p)) 

#location [0,L]~(2k+1) then k
sample.locs <- function(k,L){
  rs.unif <- runif(2*k+1,0,1)
  vec <- sort(ceiling(rs.unif*L))
  dup_positions <- which(duplicated(vec) | duplicated(vec, fromLast = TRUE))
  while (length(dup_positions)!=0) {
    vec[dup_positions[1]] <-  vec[dup_positions[1]]-1
    dup_positions <- which(duplicated(vec) | duplicated(vec, fromLast = TRUE))
    #vec;dup_positions
  }
  return(vec[seq(2, length(vec), 2)])
}
##set.seed(1)
sample.locs(5,100)

#k+1 velocities
#rgamma(k+1,shape = 2,rate = 1)

 
EnPara <- matrix(NA, ncol = N, nrow = 2*kmax+2)
EnPara[1,] <- sample(c(0:kmax),N,replace=TRUE,dpois(c(0:kmax), lambda = lambda_p)) 
##set.seed(111)
for(i in 1:N){
  kk <- EnPara[1,i]
  if(kk>0){
    c.x <- sample.locs(kk,L)
    EnPara[1+(1:(kk)),i] <-   c.x}
  #重新定义 gamma分布的参数
  c.v <- rgamma(kk+1,shape=g_alpha,rate=g_beta)
  EnPara[kk+1+(1:(kk+1)),i] <- c.v  
}

#对于状态向量的先验


#~~~~~~~~~~~~~~~~~~~~~~~~~Add for Part 2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Other ensemble of parameters.


# 03.prior -----------------------------------------------------------------


#true
#~~~~ 3.1 ~~~~ The same setting with the previous part.



#Pre work
### true forecast distribution
locs=1:n
distmat=rdist(locs,locs)  #
scale=20
sig.eps=.2
kappa=3


lik.norm=function(sigma,y,x){exp(sum(dnorm(y,x,sigma,log=FALSE)))}  # log=true/false??
distmat=rdist(1:n,1:n); scale=100; sig.eps=.5 ;kappa=3; df.t=2

Obs.site <- sort(sample(c(1:L),ceiling(n/10))) #obs location
m.obs <- length(Obs.site)                     # m
m <- m.obs
### true forecast distribution
#Will be rehanced in this way.
forecast.mean.true=rep(0,n)+.2
forecast.cov.true=exp(-(distmat/scale)^0.8)
forecast.chol.true=t(chol(forecast.cov.true))

K.true=forecast.cov.true[,Obs.site]%*%
  solve(forecast.cov.true[Obs.site,Obs.site]+sig.eps^2*diag(m.obs))



t.obs <- (0:50)*10 # Leave 500-600 for forecast.

EnPara.l <- EnPara                              # Parameter for lik guiding
EnPv1d <- runif(N, min = 0.2, max = 0.7)        # Parameter for Normal method

#set.seed(1234)
##!!!!!初始向量的矩阵

#Enstat <- forecast.mean.true+forecast.chol.true%*%matrix(rnorm(n*N),nrow=n)+3
#Enstat <- matrix(NA, nrow = n, ncol = N)
#for(i in 1:N){Enstat[,i] <- Cini+rep(rnorm(1,0,0.2),n)}

#plot(Enstat[,1])
#for(i in 1:N){lines(Enstat[,i])}


#目前有N个state和parameter的粒子



# 04. DA process ----------------------------------------------------------

#首先定义提议改变分布的的函数
# 给定当前k和固定的参数c 保证减少间断点数的变化 c取小一点
Prop.Index<- function(k,c=0.3,lambda = lambda_p,kmax=kmax){
  # k: is the index of current model
  # c: the adjusting parameter for probability
  # lambda: the p of poisson. 
  if(k<kmax){bk <- c*min(1,dpois(k+1, lambda = lambda_p)/dpois(k, lambda = lambda_p))}else{bk <- 0}
  if(k>0){dk <- c*min(1,dpois(k-1, lambda = lambda_p)/dpois(k, lambda = lambda_p))}else{dk <- 0}
  etak <- (1-bk-dk)/2
  phik <- (1-bk-dk)/2
  if(k==0){phik <- 0}else{phik <- (1-bk-dk)/2}
  pp <- c(bk,dk,etak,phik)/sum(c(bk,dk,etak,phik))
  return(list('p'=pp,'index'=sample(c(1:4),1,FALSE,pp)))
}
Prop.Index(0,c=0.3,lambda = lambda_p,kmax=kmax)
vec <- rep(0,200)
data <- unlist(sapply(vec,Prop.Index,lambda = lambda_p,kmax=kmax)[2,])

#测试k和c组合下4种提案类型的概率
hist(data,
     col = 'pink',
     breaks = seq(0,5,1),#
     density = FALSE,#添加阴影线，此处不是密度线，默认斜率为45°，没有此项作出的图为实心
     freq = FALSE#
)

lines(density(data),
      col='blue',
      lty=2,
      lwd=2
)


#For all the proposal

pro.type <- Prop.Index(k=1,lambda = lambda_p,kmax=kmax)$index




# 05. 真实的参数以及真实参数下的时空场 ----------------------------------------------------
#Para.true <- c(2,0.5,0.7,0.3,101,301)
#x <- seq(0,100, by=0.25) #Loc
#n <- n;print(n)
#m <- 200 # number of observations

#dx <- x[2]-x[1]

t.obs <- seq(0,times[length(times)], by = 10)





#set.seed(1)
#Cini=as.numeric(forecast.mean.true+forecast.chol.true%*%rnorm(n))+3
sx <- seq(0,1,by=0.0025)
Cini <- 80*sin(60*sx/pi)*sx*(2/3-sx)*(1-sx)*exp(-2*sx);#Cini[330:401]<- rep(0,72)
plot(Cini)
Cini[401]
#Real evolution for observations. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
K = 2 # If 2 change points, then 3 stage.
c.x <- c(0,101,251,n) #注意:左边加0,右边加n
#c.v <- c(0.1,0.6,0.1)
c.v <- c(0.7,0.2,0.4)
c.x.true <- c.x; c.v.true <- c.v;
#c.v <- c(2,3,2)
para <- c(K, c.x, c.v, 0.0)
out <- ode.1D(Cini, times, model2, para, method = "adams", 
              names = c("C"), dimens = n)
image(out,xlab = 'times',ylab = 'positions')
dim(out)
plot(out[601,2:402],type='l');lines(Cini,col='blue')

# out2 <- matrix(NA,nrow = length(times),ncol = n)
# out2[1,] <- Cini
# for (i in 1: (length(times)-1)){ #i=1
#   out2[i+1,] <-  ode.1D(out2[i,], 
#                         (i-1):i, #0:1,# 
#                         model2, para, method = "adams",
#                         names = c("C"), dimens = n)[2,(1:n)+1]
#   print(i)
# }
# image(out2,xlab = 'times',ylab = 'positions')
# dim(out);dim(out2)
#既然这样，直接用多次的作为真实数据 (不好 差太多)
#out[,1+(1:n)] <- out2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~比较的初始设定~~~~~~~~~~~~~~~~~~~#
Enstat <- matrix(NA, nrow = n, ncol = N)                           #
for(i in 1:N){Enstat[,i] <- Cini*(1+rep(rnorm(1,0,0.1),n))}          #
Enstat_1S  <- Enstat # Classical MCMC   1个断点                            #
Enstat_2S  <- Enstat # Classical MCMC 同真实数据一样2个断点， 3个steps     #
Enstat_3S  <- Enstat # Classical MCMC   3个断点， 4个steps                 # 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


plot(Enstat[,1])
for(i in 1:N){lines(Enstat[,i])}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~比较的参数集合~~~~~~~~~~~~~~~~~~~#
#

# Classical MCMC   1个断点                                                 #
EnPara_1S <- matrix(NA, ncol = N, nrow = 2*kmax+2)
EnPara_1S[1,] <- rep(1,N) 
#set.seed(111)
for(i in 1:N){
  kk <- EnPara_1S[1,i]
  if(kk>0){
    c.x <- sample.locs(kk,L)
    EnPara_1S[1+(1:(kk)),i] <-   c.x}
  c.v <- rgamma(kk+1,shape = 2,rate = 1)
  EnPara_1S[kk+1+(1:(kk+1)),i] <- c.v  
}

# Classical MCMC 同真实数据一样2个断点， 3个steps                          #
#
EnPara_2S <- matrix(NA, ncol = N, nrow = 2*kmax+2)
EnPara_2S[1,] <- rep(2,N) 

for(i in 1:N){
  kk <- EnPara_2S[1,i]
  if(kk>0){
    c.x <- sample.locs(kk,L)
    EnPara_2S[1+(1:(kk)),i] <-   c.x}
  c.v <- rgamma(kk+1,shape = 2,rate = 1)
  EnPara_2S[kk+1+(1:(kk+1)),i] <- c.v  
}


# Classical MCMC   3个断点， 4个steps                                      # 
#
# Classical MCMC   1个断点                                                 #
EnPara_3S <- matrix(NA, ncol = N, nrow = 2*kmax+2)
EnPara_3S[1,] <- rep(3,N) 

for(i in 1:N){
  kk <- EnPara_3S[1,i]
  if(kk>0){
    c.x <- sample.locs(kk,L)
    EnPara_3S[1+(1:(kk)),i] <-   c.x}
  c.v <- rgamma(kk+1,shape = 2,rate = 1)
  EnPara_3S[kk+1+(1:(kk+1)),i] <- c.v  
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



# 06.Resampling 重采样步骤 ----------------------------------------------------------------

# 6.0 先定义提议函数.


#PPF1: Propose Parameter Function --- Birth

Birth_velocity <- function(v0,c_o1,c_o2,c_n,mu_b){
  # 输入原来的速度v0
  # 输入原来的两个端点 c_o1 c_o2 #注意是字母old
  # 输入新提议的点 c_n
  # 均匀分布的随机数 mu_b
  
  A <- (c_n - c_o1)#/L
  B <- (c_o2 - c_n)#/L
  C <- (c_o2 - c_o1)*log(v0)#/L
  D <- (1-mu_b)/mu_b
  
  equations_type_A <- function(vars) {
    x1 <- vars[1]
    x2 <- vars[2]
    eq1 <- A * log(x1) + B * log(x2) - C  # 第一方程
    eq2 <- x1 / x2 - D                    # 第二方程
    c(eq1, eq2)
  }
  
  # 设置初始值（估计值）
  start_vals <- c(0.5, 0.5)
  # 使用nleqslv求解方程组
  solution <- nleqslv(start_vals, equations_type_A)
  return(c(max(min(solution$x[1],3),0.05),max(min(solution$x[2],3),0.05)) )
}

PPF1 <- function(para_0){
  # drawing a point uniformly from [0, L]
  #ipos <- sample(0:k_0,1) # i-th position. 0 is the first step
  
  #
  #para_0 <- c(2,100,200,0.1,0.3,0.5)
  
  k_0 <- para_0[1]
  
  
  if(k_0 != 0){
    c.x <- para_0[1+(1:(k_0))]
    c.v <- para_0[1+k_0+(1:(k_0+1))]
    
    #与之前的抽样方法不同，之前的是先选一个区间，现在是直接从[0,L]上选
    ipos <- max(ceiling(runif(1,0,1)*L),L-1)
    #ipos <- 316
    # 然后判断新提议的点是否重复
    
    while (any(abs(ipos - c.x) <= 2)) {
      ipos <- sample(0:L, 1)  # 重新从区间[0, L]中抽取d
    }
    
    insert_pos <- which(c.x > ipos)[1]  # 找到第一个大于d的位置
    if (is.na(insert_pos)) {
      # 如果ipos大于所有元素，插入到最后
      c.x <- c(c.x, ipos)
      insert_pos <- length(c.x)
    } else {
      # 否则将ipos插入到相应位置
      c.x <- append(c.x, ipos, after = insert_pos - 1)
    }
    
    # 输出结果
    #cat("插入后的向量:", c.x, "\n")
    #cat("d的位置:", insert_pos, "\n")
    
    
    #u <-runif(1,0,1)
    u <-runif(1,0.4,0.6)
    v0=c.v[insert_pos]
    if(k_0==0){v0=c.v}
    if(insert_pos>1){c_o1=c.x[insert_pos-1]}else{c_o1=0}
    if(insert_pos<k_0+1){c_o2=c.x[insert_pos+1]}else{c_o2=n}
    prp.v <- Birth_velocity(v0,c_o1,c_o2,ipos,u)
    #prp.v <- c(2.56580150,0.05366406)
    c.v[insert_pos] <- prp.v[1]
    c.v <- append(c.v,prp.v[2],after=insert_pos)
    #cat(paste(c('The proposed paramter:',k_0+1,c.x,c.v)))
  }else{ #k==0
    c.v <- para_0[1+k_0+(1:(k_0+1))]
    ipos <- ceiling(runif(1,0,1)*L)
    #u <-runif(1,0,1)
    u <-runif(1,0.4,0.6)
    v0=c.v
    c_o1=0
    c_o2=n
    prp.v <- Birth_velocity(v0,c_o1,c_o2,ipos,u)
    c.x <- ipos
    c.v <- prp.v
  }
  para_new <- c(k_0+1,c.x,c.v)
  
  
  ratio_f <- ((2*k_0+2)*(2*k_0+3)/L**2)*(ipos-c_o1)*(c_o2-ipos)/(c_o2-c_o1)
  ratio_v <- (g_beta)**g_alpha/gamma(g_alpha)*(prp.v[1]*prp.v[2]/v0)**(g_alpha-1)*exp(g_beta*(prp.v[1]+prp.v[2]-v0))
  ratio_p <- Prop.Index(k_0+1,c=0.3,lambda = lambda_p,kmax=kmax)$p[2]*L/
    (Prop.Index(k_0,c=0.3,lambda = lambda_p,kmax=kmax)$p[1]*(k_0+1))
  ratio_r <- dpois(k_0+1, 2, log = FALSE)/dpois(k_0, 2, log = FALSE)
  jaco <- (prp.v[1]+prp.v[2])**2/v0
  ratio <- ratio_f*ratio_v*ratio_p*ratio_r*jaco
  #ratio <- 1
  return(list('para'=para_new, 'ratio'=ratio))
}

para_0 <- c(2,100,200,0.2,0.3,0.5)
PPF1(c(2,100,200,0.2,0.4,0.6))
PPF1(c(2,100,200,0.2,0.3,0.5))
PPF1(c(1,200,0.2,0.6))
PPF1(c(0,0.8))

#PPF2: Propose Parameter Function --- Death
PPF2_2 <- function(para_0){
  # drawing a point uniformly from [0, L]
  #ipos <- sample(0:k_0,1) # i-th position. 0 is the first step
  
  #
  #para_0 <- c(2,100,200,0.1,0.3,0.5)
  
  k_0 <- para_0[1]
  c.x <- para_0[1+(1:(k_0))]
  c.v <- para_0[1+k_0+(1:(k_0+1))]
  
  ipos <- sample(c(1:k_0),1)
  c.x.new <- c.x[-ipos]
  
  v_o1 <- c.v[ipos]
  v_o2 <- c.v[ipos+1]
  if(ipos > 1){c_o1 <- c.x[ipos-1]}else{c_o1 <-0}    # 去除点的前一个点
  if(ipos < k_0){c_o2 <- c.x[ipos+1]}else{c_o2 <- n} # 去除点的后一个点
  c_s <- c.x[ipos] #要去除的点
  
  v_new <- exp(((c_s-c_o1)*log(v_o1)+(c_o2-c_s)*log(v_o2))/(c_o2-c_o1)) 
  
  c_0 <- c.x[ipos]
  c.x <- c.x[-ipos]
  c.v <- c.v[-ipos]
  c.v[ipos] <- v_new
  para_new <-c(k_0-1,c.x,c.v)
  
  prp.v <- c(v_o1,v_o2)
  v0 <- v_new
  c_s <- ipos
  
  
  ratio_f <- ((2*k_0)*(2*k_0+1)/L**2)*(c_0-c_o1)*(c_o2-c_0)/(c_o2-c_o1)
  ratio_v <- (g_beta)**g_alpha/gamma(g_alpha)*(prp.v[1]*prp.v[2]/v0)**(g_alpha-1)*exp(g_beta*(prp.v[1]+prp.v[2]-v0))
  ratio_p <- Prop.Index(k_0,c=0.3,lambda = lambda_p,kmax=kmax)$p[2]*L/
    (Prop.Index(k_0-1,c=0.3,lambda = lambda_p,kmax=kmax)$p[1]*(k_0))
  ratio_r <- dpois(k_0, 2, log = FALSE)/dpois(k_0-1, 2, log = FALSE)
  jaco <- (prp.v[1]+prp.v[2])**2/v0
  ratio <- 1/(ratio_f*ratio_v*ratio_p*ratio_r*jaco)
  
  #ratio <- 1
  return(list('para'=para_new,'ratio'=ratio))
}


#PPF2的rate用PPF1的rate计算就好了
PPF1.fix <- function(para_0,ipos1,v_1){
  # drawing a point uniformly from [0, L]
  #ipos <- sample(0:k_0,1) # i-th position. 0 is the first step
  
  #
  #para_0 <- c(2,100,200,0.1,0.3,0.5)
  
  k_0 <- para_0[1]
  
  
  if(k_0 != 0){
    c.x <- para_0[1+(1:(k_0))]
    c.v <- para_0[1+k_0+(1:(k_0+1))]
    
    #与之前的抽样方法不同，之前的是先选一个区间，现在是直接从[0,L]上选
    #ipos <- ceiling(runif(1,0,1)*L)
    ipos <- ipos1
    # 然后判断新提议的点是否重复
    
    while (any(abs(ipos - c.x) <= 2)) {
      ipos <- sample(0:L, 1)  # 重新从区间[0, L]中抽取d
    }
    
    insert_pos <- which(c.x > ipos)[1]  # 找到第一个大于d的位置
    if (is.na(insert_pos)) {
      # 如果ipos大于所有元素，插入到最后
      c.x <- c(c.x, ipos)
      insert_pos <- length(c.x)
    } else {
      # 否则将ipos插入到相应位置
      c.x <- append(c.x, ipos, after = insert_pos - 1)
    }
    
    # 输出结果
    #cat("插入后的向量:", c.x, "\n")
    #cat("d的位置:", insert_pos, "\n")
    
    
    #u <-runif(1,0,1)
    u <-runif(1,0.4,0.6)
    v0=c.v[insert_pos]
    if(k_0==0){v0=c.v}
    if(insert_pos>1){c_o1=c.x[insert_pos-1]}else{c_o1=0}
    if(insert_pos<k_0+1){c_o2=c.x[insert_pos+1]}else{c_o2=n}
    prp.v <- v_1
    #prp.v <- Birth_velocity(v0,c_o1,c_o2,ipos,u)
    #prp.v <- c(2.56580150,0.05366406)
    c.v[insert_pos] <- prp.v[1]
    c.v <- append(c.v,prp.v[2],after=insert_pos)
    #cat(paste(c('The proposed paramter:',k_0+1,c.x,c.v)))
  }else{ #k==0
    c.v <- para_0[1+k_0+(1:(k_0+1))]
    #ipos <- ceiling(runif(1,0,1)*L)
    #u <-runif(1,0,1)
    ipos <- ipos1
    u <-runif(1,0.4,0.6)
    v0=c.v
    c_o1=0
    c_o2=n
    #prp.v <- Birth_velocity(v0,c_o1,c_o2,ipos,u)
    prp.v <- v_1
    c.x <- ipos
    c.v <- prp.v
  }
  para_new <- c(k_0+1,c.x,c.v)
  
  
  ratio_f <- ((2*k_0+2)*(2*k_0+3)/L**2)*(ipos-c_o1)*(c_o2-ipos)/(c_o2-c_o1)
  ratio_v <- (g_beta)**g_alpha/gamma(g_alpha)*(prp.v[1]*prp.v[2]/v0)**(g_alpha-1)*exp(g_beta*(prp.v[1]+prp.v[2]-v0))
  ratio_p <- Prop.Index(k_0+1,c=0.3,lambda = lambda_p,kmax=kmax)$p[2]*L/
    (Prop.Index(k_0,c=0.3,lambda = lambda_p,kmax=kmax)$p[1]*(k_0+1))
  ratio_r <- dpois(k_0+1, 2, log = FALSE)/dpois(k_0, 2, log = FALSE)
  jaco <- (prp.v[1]+prp.v[2])**2/v0
  ratio <- ratio_f*ratio_v*ratio_p*ratio_r*jaco
  #ratio <- 1
  return(list('para'=para_new, 'ratio'=ratio))
}
PPF2 <- function(para_0){
  # drawing a point uniformly from [0, L]
  #ipos <- sample(0:k_0,1) # i-th position. 0 is the first step
  
  #
  #para_0 <- c(2,100,200,0.1,0.3,0.5)
  
  k_0 <- para_0[1]
  c.x <- para_0[1+(1:(k_0))]
  c.v <- para_0[1+k_0+(1:(k_0+1))]
  
  ipos <- sample(c(1:k_0),1)
  c.x.new <- c.x[-ipos]
  
  v_o1 <- c.v[ipos]
  v_o2 <- c.v[ipos+1]
  if(ipos > 1){c_o1 <- c.x[ipos-1]}else{c_o1 <-0}    # 去除点的前一个点
  if(ipos < k_0){c_o2 <- c.x[ipos+1]}else{c_o2 <- n} # 去除点的后一个点
  c_s <- c.x[ipos] #要去除的点
  
  v_new <- exp(((c_s-c_o1)*log(v_o1)+(c_o2-c_s)*log(v_o2))/(c_o2-c_o1)) 
  
  c_0 <- c.x[ipos]
  c.x <- c.x[-ipos]
  c.v <- c.v[-ipos]
  c.v[ipos] <- v_new
  para_new <-c(k_0-1,c.x,c.v)
  
  prp.v <- c(v_o1,v_o2)
  v0 <- v_new
  c_s <- c_0
  ratio <- 1/PPF1.fix(para_new,c_0,prp.v)$ratio
  #new version
  return(list('para'=para_new,'ratio'=ratio))
}

para_0 <- c(3,100,150,200,0.2,0.4,0.6,0.8)
PPF2(c(3,100,150,200,0.2,0.4,0.6,0.8))$ratio
PPF2_2(c(3,100,150,200,0.2,0.4,0.6,0.8))$ratio
PPF2(c(2,100,150,0.4,0.6,0.8))
#PPF3: Propose Parameter Function --- Propose a new velocity.
PPF3 <- function(para_0){
  k_0 <- para_0[1]
  c.x <- para_0[1+(1:(k_0))]
  c.v <- para_0[1+k_0+(1:(k_0+1))]
  
  ipos <- sample(c(1:(k_0+1)),1)
  
  #mu_v <- runif(1,-1/2,1/2)
  mu_v <- runif(1,-1/5,1/5)
  v_new <- exp(mu_v)*c.v[ipos]
  v_old <- c.v[ipos]
  c.v[ipos] <- v_new
  if(k_0>0){para_new <-c(k_0,c.x,c.v)}else{para_new <-c(k_0,c.v)}
  
  ratio <- ((v_new/v_old)**g_alpha)*exp(-g_beta*(v_new-v_old))
  #ratio <- 1
  return(list('para'=para_new,'ratio'=ratio))
}
PPF3(c(0,1.8))

vec_v <- rep(NA,100)
for(i in 1:100){vec_v[i] <- PPF3(c(0,2.5))$para[2]}
plot(c(1:100),vec_v,type='l')
# Propose Parameter Function --- Propose a new location.
PPF4 <- function(para_0){
  k_0 <- para_0[1]
  c.x <- para_0[1+(1:(k_0))]
  c.v <- para_0[1+k_0+(1:(k_0+1))]
  
  ipos <- sample(c(1:k_0),1)
  if(ipos>1){c_left <- c.x[ipos-1]}else{c_left <- 1}
  if(ipos<k_0){c_right <- c.x[ipos+1]}else{c_right <- n}
  new_x <- ceiling(runif(1,c_left,c_right-1.5))
  old_x <- c.x[ipos]
  c.x[ipos] <- new_x
  ratio <- (c_right-old_x)*(old_x-c_left)/(c_right-old_x)*(old_x-c_left)
  para_new <-c(k_0,c.x,c.v)
  #ratio <- 1
  return(list('para'=para_new,'ratio'=ratio))
}
para_0 <- c(2, 264, 400,   0.1997230,   0.4621581,   0.3372146)
PPF4(para_0)
#PPF0: Propose Parameter Function --- Type
Prp.Para <- function(para_0,prpd.index){
  k <- para_0[1]
  if(prpd.index ==  1){prpd.para<- PPF1(para_0)}
  if(prpd.index ==  2){prpd.para<- PPF2(para_0)}
  if(prpd.index ==  3){prpd.para<- PPF3(para_0)}
  if(prpd.index ==  4){prpd.para<- PPF4(para_0)}
  return(prpd.para)
}

#prpd.index <- Prop.Index(k,c=0.3,lambda=lambda_p,kmax=kmax)$index
#Prp.Para(para_0,prpd.index)

#现在先有三个矩阵，上一时刻的状态向量，这一时刻的状态向量，以及参数的矩阵
#set.seed(1)

DA_Process <- function(Enstat,EnPara,Mthd){
  En.s.t_1 <- Enstat # 表示前一时刻的状态向量 x_(t-1)
  En.s.t <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
  En.lik <- rep(NA,N)
  
  #Def for Bootstrap Resampling.
  Mid.En.s.t_1 <- En.s.t_1
  Mid.En.s.t   <- En.s.t
  Mid.En.lik   <- En.lik
  Mid.EnPara   <- EnPara
  for (i.t in 1:(length(t.obs)-1)){   # Loop of time. i.t=1 i.t=2
    t1 <- t.obs[i.t]; t2 <- t.obs[i.t+1]
    sta.obs <- out[t2+1,-1] *(1+ sig.eps*rnorm(n))
    #sta.obs <- out[t2+1,-1] + sig.eps*rnorm(n) #Perturbed observations.
    if(0){
      #i.t < length(t.obs)/2
      #后面把*9改成*18
      sig.eps <- (10-(i.t-1)/length(t.obs)*9)*sig.eps0}else{ 
        sig.eps <- sig.eps0}
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #Loop_1:       Propagate          #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    
    
    # 在此步骤中，把状态向量(En.s.t_1) propagate到新时刻(En.s.t)
    for (i.n in 1:N){
      k_0 <- EnPara[1,i.n]
      para_0 <- EnPara[1:(2*k_0+2),i.n]; print(para_0)
      if(k_0 == 0){
        p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
      }else{
        p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
      }
      # 这是原来的 p4evo0 <- c(k_0,0,para_0[2*(1:k_0)+1], n, para_0[2*(0:k_0)+2],0.0)
      sta.ini <- En.s.t_1[,i.n]
      sta.fore0 <- ode.1D(sta.ini, t1:t2, model2, p4evo0, method = "adams", 
                          names = c("C"), dimens = n)[t2-t1+1,1:n+1]
      lik0 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore0[Obs.site])
      En.lik[i.n] <- lik0
      En.s.t[,i.n] <- sta.fore0
    }
    
    if(plot_t){
      plot(out[t2+1,-1],type='l',lwd = 3,ylab = paste('t=',t2))
      points(Obs.site,sta.obs[Obs.site])
      for(i.n in 1:N){lines(En.s.t[,i.n],col = 'red')}
    } 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #Loop_2:       Bootstrap          #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    if(max(En.lik) != 0){
      Fathers <- sample(c(1:N),N,En.lik,replace = TRUE)
    }else{Fathers <- c(1:N)}
    #Fathers <- c(1:N)
    # 重采样 代替 注意要用中间变量！！！！ 
    # 只考虑一半
    II2 <-  which(En.lik < median(En.lik)) 
    for (i.n in II2) {
      Mid.En.s.t_1[,i.n] <- En.s.t_1[,Fathers[i.n]]
      Mid.En.s.t[,i.n]   <- En.s.t[,Fathers[i.n]]
      Mid.En.lik[i.n]   <- En.lik[Fathers[i.n]]
      Mid.EnPara[,i.n]   <- EnPara[,Fathers[i.n]]
    }
    Mid.En.s.t_1[,setdiff(c(1:N),II2)] <- En.s.t_1[,setdiff(c(1:N),II2)]
    Mid.En.s.t[,setdiff(c(1:N),II2)]   <- En.s.t[,setdiff(c(1:N),II2)]
    Mid.En.lik[setdiff(c(1:N),II2)]   <- En.lik[setdiff(c(1:N),II2)]
    Mid.EnPara[,setdiff(c(1:N),II2)]   <- EnPara[,setdiff(c(1:N),II2)] 
    #替换完成后，再把矩阵都代回去
    En.s.t_1 <- Mid.En.s.t_1
    En.s.t   <- Mid.En.s.t
    En.lik   <- Mid.En.lik
    EnPara   <- Mid.EnPara
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #Loop_3:  RJMCMC Re-sampling      #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    for (i.n in 1:N){   # Loop of ensemble member. i.n = 2
      #因为第一步的重采样，重复第一步的输入参数步骤
      for (gg in 1:G) {
        k_0 <- EnPara[1,i.n]
        para_0 <- EnPara[1:(2*k_0+2),i.n]; print(para_0)
        if(k_0 == 0){
          p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
        }else{
          p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
        }
        #接下来就是提议新的参数
        if (Mthd == 1){
          prpd.index <- Prop.Index(k_0,c=0.3,lambda=lambda_p,kmax=kmax)$index
        }else{
          prpd.index <- sample(c(3,4),1)  
        }
        
        propose1 <- Prp.Para(para_0,prpd.index) 
        para_1 <- propose1$para; print(para_1) # 这时候要定义函数了，提一个新的函数
        k_1 <- para_1[1]
        if(k_1 == 0){
          p4evo1 <- c(k_1,0,n,para_1[1+k_1+(1:(k_1+1))],0.0)
        }else{
          p4evo1 <- c(k_1,0,para_1[1+(1:(k_1))], n, para_1[1+k_1+(1:(k_1+1))],0.0)
        }
        sta.ini <- En.s.t_1[,i.n]
        sta.fore1 <- ode.1D(sta.ini, t1:t2, model2, p4evo1, method = "adams", 
                            names = c("C"), dimens = n)[t2-t1+1,1:n+1]
        
        
        
        lik1 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore1[Obs.site])
        if(plot_t){lines(sta.fore1,col='blue')}
        # 比较 En.lik[i.n] 和 lik0
        
        #注意两点，一个是雅各比，一个是初始的向量别错了
        
        #然后再加if
        
        #~~~~~~~~~~~~~~~~End loop # ~~~~~~~~~~
        
        #原来未更新的部分
        lik0 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore0[Obs.site])
        lik1 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore1[Obs.site])
        if((lik1*propose1$ratio)>lik0*runif(1,0,1)){
          EnPara[,i.n] <- rep(NA,2*kmax+2)
          EnPara[1:(2*k_1+2),i.n] <- para_1
          sta.fore0 <- sta.fore1
        }
      }
      En.s.t[,i.n] <- sta.fore0
      # K*(obs- perturbed forecast.)
      if(kalmangain){En.s.t[,i.n] <- sta.fore0 + 
        K.true %*% as.matrix(sta.obs[Obs.site] - sta.fore0[Obs.site] - sig.eps*rnorm(m))} 
    }
    # for (i.n in 1:N){while(max(Enstat[,i.n])> 10){
    #     i.new <- sample(c(1:N),1)
    #     Enstat[,i.n] <- Enstat[,i.new]
    #     EnPara[,i.n] <- EnPara[,i.new]
    #   }}
    En.s.t_1 <- En.s.t
  }
  return(list('EnS'=En.s.t_1,'EnP'=EnPara))
}

#Res1_Rj <- DA_Process(Enstat,EnPara,1)
#Res1_2S <- DA_Process(Enstat_2S,EnPara_2S,0)
Enstat_2S



En.s.t_1 <- Enstat # 表示前一时刻的状态向量 x_(t-1)
En.s.t <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
En.lik <- rep(NA,N)

#Def for Bootstrap Resampling.
Mid.En.s.t_1 <- En.s.t_1
Mid.En.s.t   <- En.s.t
Mid.En.lik   <- En.lik
Mid.EnPara   <- EnPara

# 1.0 The main loop. -------------------------------------------------------
#2025-07-31
Indice.m <- matrix(NA, ncol =length(t.obs)-1,nrow = kmax+1)
Loca.m  <- matrix(NA, ncol =length(t.obs)-1,nrow = 2)
index_k <- as.vector(table(factor(EnPara[1, ], levels = 0:3)))
#Indice.m[,1] <- index_k
#EnPara_2S <- EnPara
# For easier index, the dim of observations is n, not m.
for (i.t in 1:(length(t.obs)-1)){   # Loop of time. i.t=1 i.t=2

  t1 <- t.obs[i.t]; t2 <- t.obs[i.t+1]
  sta.obs <- out[t2+1,-1] *(1+ sig.eps*rnorm(n))
  #sta.obs <- out[t2+1,-1] + sig.eps*rnorm(n) #Perturbed observations.
  if(0){
    #i.t < length(t.obs)/2
    #后面把*9改成*18
    sig.eps <- (10-(i.t-1)/length(t.obs)*9)*sig.eps0}else{ 
      sig.eps <- sig.eps0}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Loop_1:       Propagate          #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  
  
  # 在此步骤中，把状态向量(En.s.t_1) propagate到新时刻(En.s.t)
  for (i.n in 1:N){
    k_0 <- EnPara[1,i.n]
    para_0 <- EnPara[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    }else{
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    # 这是原来的 p4evo0 <- c(k_0,0,para_0[2*(1:k_0)+1], n, para_0[2*(0:k_0)+2],0.0)
    sta.ini <- En.s.t_1[,i.n]
    sta.fore0 <- ode.1D(sta.ini, t1:t2, model2, p4evo0, method = "adams", 
                        names = c("C"), dimens = n)[t2-t1+1,1:n+1]
    lik0 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore0[Obs.site])
    En.lik[i.n] <- lik0
    En.s.t[,i.n] <- sta.fore0
  }
  
  if(plot_t){
    plot(out[t2+1,-1],type='l',lwd = 3,ylab = paste('t=',t2))
    points(Obs.site,sta.obs[Obs.site])
    for(i.n in 1:N){lines(En.s.t[,i.n],col = 'red')}
  } 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Loop_2:       Bootstrap          #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  if(max(En.lik) != 0){
    Fathers <- sample(c(1:N),N,En.lik,replace = TRUE)
  }else{Fathers <- c(1:N)}
  #Fathers <- c(1:N)
  # 重采样 代替 注意要用中间变量！！！！ 
  # 只考虑一半
  II2 <-  which(En.lik < median(En.lik)) 
  for (i.n in II2) {
    Mid.En.s.t_1[,i.n] <- En.s.t_1[,Fathers[i.n]]
    Mid.En.s.t[,i.n]   <- En.s.t[,Fathers[i.n]]
    Mid.En.lik[i.n]   <- En.lik[Fathers[i.n]]
    Mid.EnPara[,i.n]   <- EnPara[,Fathers[i.n]]
  }
  Mid.En.s.t_1[,setdiff(c(1:N),II2)] <- En.s.t_1[,setdiff(c(1:N),II2)]
  Mid.En.s.t[,setdiff(c(1:N),II2)]   <- En.s.t[,setdiff(c(1:N),II2)]
  Mid.En.lik[setdiff(c(1:N),II2)]   <- En.lik[setdiff(c(1:N),II2)]
  Mid.EnPara[,setdiff(c(1:N),II2)]   <- EnPara[,setdiff(c(1:N),II2)] 
  #替换完成后，再把矩阵都代回去
  En.s.t_1 <- Mid.En.s.t_1
  En.s.t   <- Mid.En.s.t
  En.lik   <- Mid.En.lik
  EnPara   <- Mid.EnPara
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Loop_3:  RJMCMC Re-sampling      #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  for (i.n in 1:N){   # Loop of ensemble member. i.n = 2
    #因为第一步的重采样，重复第一步的输入参数步骤
    for (gg in 1:G) {
      k_0 <- EnPara[1,i.n]
      para_0 <- EnPara[1:(2*k_0+2),i.n]; print(para_0)
      if(k_0 == 0){
        p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
      }else{
        p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
      }
      #接下来就是提议新的参数
      prpd.index <- Prop.Index(k_0,c=0.3,lambda=lambda_p,kmax=kmax)$index
      propose1 <- Prp.Para(para_0,prpd.index) 
      para_1 <- propose1$para; print(para_1) # 这时候要定义函数了，提一个新的函数
      k_1 <- para_1[1]
      if(k_1 == 0){
        p4evo1 <- c(k_1,0,n,para_1[1+k_1+(1:(k_1+1))],0.0)
      }else{
        p4evo1 <- c(k_1,0,para_1[1+(1:(k_1))], n, para_1[1+k_1+(1:(k_1+1))],0.0)
      }
      sta.ini <- En.s.t_1[,i.n]
      sta.fore1 <- ode.1D(sta.ini, t1:t2, model2, p4evo1, method = "adams", 
                          names = c("C"), dimens = n)[t2-t1+1,1:n+1]
      
      
      
      lik1 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore1[Obs.site])
      if(plot_t){lines(sta.fore1,col='blue')}
      # 比较 En.lik[i.n] 和 lik0
      
      #注意两点，一个是雅各比，一个是初始的向量别错了
      
      #然后再加if
      
      #~~~~~~~~~~~~~~~~End loop # ~~~~~~~~~~
      
      #原来未更新的部分
      lik0 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore0[Obs.site])
      lik1 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore1[Obs.site])
      if((lik1*propose1$ratio)>lik0*runif(1,0,1)){
        EnPara[,i.n] <- rep(NA,2*kmax+2)
        EnPara[1:(2*k_1+2),i.n] <- para_1
        sta.fore0 <- sta.fore1
      }
    }
    En.s.t[,i.n] <- sta.fore0
    # K*(obs- perturbed forecast.)
    if(kalmangain){En.s.t[,i.n] <- sta.fore0 + 
      K.true %*% as.matrix(sta.obs[Obs.site] - sta.fore0[Obs.site] - sig.eps*rnorm(m))} 
  }
  # for (i.n in 1:N){while(max(Enstat[,i.n])> 10){
  #     i.new <- sample(c(1:N),1)
  #     Enstat[,i.n] <- Enstat[,i.new]
  #     EnPara[,i.n] <- EnPara[,i.new]
  #   }}
  En.s.t_1 <- En.s.t
  #indice
  index_k <- as.vector(table(factor(EnPara[1, ], levels = 0:3)))
  Indice.m[,i.t] <- index_k
  # 统计 [80,120] 区间内的元素个数
  count_75_125 <- sum(EnPara >= 75 & EnPara <= 125, na.rm = TRUE)
  
  # 统计 [230,270] 区间内的元素个数
  count_225_275 <- sum(EnPara >= 225 & EnPara <= 275, na.rm = TRUE)
  Loca.m[,i.t] <- c(count_75_125,count_225_275)
  
}

#Show the results~~~~~~~~~
lines(EnPara[1,],type = 'b',col='red')
print(sum(EnPara[1,]== 1))

tt <- 1+t.obs[i.t+1]
m.true  <- out[tt,-1]
mean.rj <- rowMeans(En.s.t_1)
par(mfrow=c(1,1),mar=c(2,2,0,0.5)+2.5)
plot(m.true,type ='l', ylim = c(min(En.s.t_1),7),
     cex.axis = 2, cex.lab =2, cex.main =2,
     main='The pointwise posterior 80% credible intervals',xlab = 'Position',ylab = 'Value')
for (i.n in 1:N){lines(En.s.t_1[,i.n],col = 'lightcoral',lwd=2)}
plot.PI(En.s.t_1,'blue')
i.t
#for (i.n in which(EnPara[1,]==3)){lines(Enstat[,i.n],col = 'lightcoral',lwd=2)}
#for (i.n in which(EnPara[1,]==2)){lines(Enstat[,i.n],col = 'darkorange',lwd=2)}
#for (i.n in which(EnPara[1,]==0)){lines(Enstat[,i.n],col = 'blue',lwd=2)}
lines(m.true, lwd = 2)

sum(EnPara[1,]==3)
sum(EnPara[1,]==2)
sum(EnPara[1,]==1)
sum(EnPara[1,]==0)
en2 <- which(EnPara[1,]==2)
apply(EnPara[,en2],1,mean)

par(mfrow=c(1,2),mar=c(2,2,0,0.5)+2.5)
for(i in 1:kmax){
  if(i==1){plot(density(EnPara[3,which(EnPara[1,]==2)]),ylim=c(0,0.008),
                main = 'The denmsities of the locations.',
                xlab='Locations')}
  for(j in 1:i){
    lines(density(EnPara[1+j,which(EnPara[1,]==i)]),col=i,pch=i,lwd=3)
  }
  print(apply(EnPara[,which(EnPara[1,]==i)],1,mean))}
for(i in 2:3){abline(v=c.x.true[i], col = "red", lty="dashed")}

#画图

# library
library(ggplot2)
library(dplyr)
library(hrbrthemes)

kmax <- 3
N.en <- c()        #表示第几个模型对应有多少个粒子
for (i in 1:kmax){N.en <- c(N.en,sum(EnPara[1,]==i))}
d.m.index <- c()     #表示第几个模型
d.para.order <- c()  #表示是第几个断点
d.location <- c()    #该粒子对应的location
for(i in 1:kmax){
  for(j in 1:i){
    d.m.index <- c(d.m.index,rep(paste("Model ", i),N.en[i]))
    d.para.order <- c(d.para.order, rep(j,N.en[i]))
    d.location <- c(d.location,EnPara[1+j,which(EnPara[1,]==i)])
  }
}



#Try
par(mfrow=c(2,1))
par(mar=c(0,5,3,3))
hist(d.location[d.m.index=='Model  2'&d.para.order==1] , breaks=4, 
     xaxt="n", las=1,
     cex.axis = 2, cex.lab = 2,ylim=c(0,40),
     xlim=c(0,400), col=rgb(1,0,0,0.5), xlab="", 
     ylab="Number of particles", main="" )
hist(d.location[d.m.index=='Model  2'&d.para.order==2], 
     cex.axis = 2, cex.lab = 2,
     breaks=10, xlim=c(0,400), col="tomato3", add=T)
par(mar=c(5,5,0,3))
hist(d.location[d.m.index=='Model  1'&d.para.order==1] , main="" , xlim=c(0,400), ylab="Number of particles", 
     cex.axis = 2, cex.lab = 2,
     xlab="Locations", ylim=c(40,0), las=1 , col=rgb(0,0,1,0.5), breaks=10)


#Try
par(mfrow=c(2,1))
par(mar=c(0,5,3,3))
hist(d.location[d.m.index=='Model  2'&d.para.order==1] , breaks=4, 
     xaxt="n", las=1,
     cex.axis = 2, cex.lab = 2,
     xlim=c(0,400), ylim=c(0,40), col=rgb(1,0,0,0.5), xlab="", 
     ylab="Number of particles", main="" )
hist(d.location[d.m.index=='Model  2'&d.para.order==2], 
     cex.axis = 2, cex.lab = 2,
     breaks=10, xlim=c(0,400), col="tomato3", add=T)
par(mar=c(5,5,0,3))
hist(d.location[d.m.index=='Model  3'&d.para.order==1] , main="" , xlim=c(0,400), ylab="Number of particles", 
     xlab="Locations", ylim=c(40,0), 
     cex.axis = 2, cex.lab = 2,
     las=1 , col='yellow', breaks=10)
hist(d.location[d.m.index=='Model  3'&d.para.order==2] , main="" , xlim=c(0,400), ylab="Number of particles", 
     xlab="Locations", ylim=c(40,0),
     cex.axis = 2, cex.lab = 2,
     las=1 , col='green', breaks=10,add =T)
hist(d.location[d.m.index=='Model  3'&d.para.order==3] , main="" , xlim=c(0,400), ylab="Number of particles", 
     xlab="Locations", ylim=c(40,0),
     cex.axis = 2, cex.lab = 2,
     las=1 , col='darkgreen', breaks=10,add =T)




# 3.0 比较之前的程序  断点数=1 -------------------------------------------------------------
En.s.t_1_1S <- Enstat_1S # 表示前一时刻的状态向量 x_(t-1)
En.s.t_1S <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
En.lik_1S <- rep(NA,N)

#Def for Bootstrap Resampling.
Mid.En.s.t_1_1S <- En.s.t_1_1S
Mid.En.s.t_1S   <- En.s.t_1S
Mid.En.lik_1S   <- En.lik_1S
Mid.EnPara_1S   <- EnPara_1S
for (i.t in 1:(length(t.obs)-1)){   # Loop of time. i.t=1 i.t=2
  t1 <- t.obs[i.t]; t2 <- t.obs[i.t+1]
  sta.obs <- out[t2+1,-1] *(1+ sig.eps*rnorm(n))
  #sta.obs <- out[t2+1,-1] + sig.eps*rnorm(n) #Perturbed observations.
  if(0){
    #i.t < length(t.obs)/2
    #后面把*9改成*18
    sig.eps <- (10-(i.t-1)/length(t.obs)*9)*sig.eps0}else{ 
      sig.eps <- sig.eps0}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Loop_1:       Propagate          #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # 在此步骤中，把状态向量(En.s.t_1) propagate到新时刻(En.s.t)
  for (i.n in 1:N){
    k_0 <- EnPara_1S[1,i.n]
    para_0 <- EnPara_1S[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    }else{
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    # 这是原来的 p4evo0 <- c(k_0,0,para_0[2*(1:k_0)+1], n, para_0[2*(0:k_0)+2],0.0)
    sta.ini <- En.s.t_1_1S[,i.n]
    sta.fore0 <- ode.1D(sta.ini, t1:t2, model2, p4evo0, method = "adams", 
                        names = c("C"), dimens = n)[t2-t1+1,1:n+1]
    lik0 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore0[Obs.site])
    En.lik_1S[i.n] <- lik0
    En.s.t_1S[,i.n] <- sta.fore0
  }
  if(plot_t){
    plot(out[t2+1,-1],type='l',lwd = 3,ylab = paste('t=',t2))
    points(Obs.site,sta.obs[Obs.site])
    for(i.n in 1:N){lines(En.s.t_1S[,i.n],col = 'red')}
  } 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Loop_2:       Bootstrap          #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  if(max(En.lik_1S) != 0){
    Fathers <- sample(c(1:N),N,En.lik_1S,replace = TRUE)
  }else{Fathers <- c(1:N)}
  #Fathers <- c(1:N)
  # 重采样 代替 注意要用中间变量！！！！ 
  # 只考虑一半
  II2 <-  which(En.lik_1S < median(En.lik_1S)) 
  for (i.n in II2) {
    Mid.En.s.t_1_1S[,i.n] <- En.s.t_1_1S[,Fathers[i.n]]
    Mid.En.s.t_1S[,i.n]   <- En.s.t_1S[,Fathers[i.n]]
    Mid.En.lik_1S[i.n]   <- En.lik_1S[Fathers[i.n]]
    Mid.EnPara_1S[,i.n]   <- EnPara_1S[,Fathers[i.n]]
  }
  Mid.En.s.t_1_1S[,setdiff(c(1:N),II2)] <- En.s.t_1_1S[,setdiff(c(1:N),II2)]
  Mid.En.s.t_1S[,setdiff(c(1:N),II2)]   <- En.s.t_1S[,setdiff(c(1:N),II2)]
  Mid.En.lik_1S[setdiff(c(1:N),II2)]   <- En.lik_1S[setdiff(c(1:N),II2)]
  Mid.EnPara_1S[,setdiff(c(1:N),II2)]   <- EnPara_1S[,setdiff(c(1:N),II2)] 
  #替换完成后，再把矩阵都代回去
  En.s.t_1_1S <- Mid.En.s.t_1_1S
  En.s.t_1S   <- Mid.En.s.t_1S
  En.lik_1S   <- Mid.En.lik_1S
  EnPara_1S   <- Mid.EnPara_1S
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Loop_3:  RJMCMC Re-sampling      #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  for (i.n in 1:N){   # Loop of ensemble member. i.n = 2
    #因为第一步的重采样，重复第一步的输入参数步骤
    k_0 <- EnPara_1S[1,i.n]
    para_0 <- EnPara_1S[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    }else{
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    #接下来就是提议新的参数
    prpd.index <- sample(c(3,4),1)
    propose1 <- Prp.Para(para_0,prpd.index) 
    para_1 <- propose1$para; print(para_1) # 这时候要定义函数了，提一个新的函数
    k_1 <- para_1[1]
    if(k_1 == 0){
      p4evo1 <- c(k_1,0,n,para_1[1+k_1+(1:(k_1+1))],0.0)
    }else{
      p4evo1 <- c(k_1,0,para_1[1+(1:(k_1))], n, para_1[1+k_1+(1:(k_1+1))],0.0)
    }
    sta.ini <- En.s.t_1_1S[,i.n]
    sta.fore1 <- ode.1D(sta.ini, t1:t2, model2, p4evo1, method = "adams", 
                        names = c("C"), dimens = n)[t2-t1+1,1:n+1]
    
    
    lik1 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore1[Obs.site])
    
    # 比较 En.lik[i.n] 和 lik0
    
    #注意两点，一个是雅各比，一个是初始的向量别错了
    
    #然后再加if
    
    #~~~~~~~~~~~~~~~~End loop # ~~~~~~~~~~
    
    #原来未更新的部分
    lik0 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore0[Obs.site])
    lik1 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore1[Obs.site])
    if((lik1*propose1$ratio)>lik0*runif(1,0,1)){
      EnPara_1S[,i.n] <- rep(NA,2*kmax+2)
      EnPara_1S[1:(2*k_1+2),i.n] <- para_1
      sta.fore0 <- sta.fore1
    }
    # K*(obs- perturbed forecast.)
    En.s.t_1S[,i.n] <- sta.fore0
    if(kalmangain){En.s.t_1S[,i.n] <- sta.fore0 + 
      K.true %*% as.matrix(sta.obs[Obs.site] - sta.fore0[Obs.site] - sig.eps*rnorm(m))
    }
  }
  # for (i.n in 1:N){while(max(Enstat[,i.n])> 10){
  #     i.new <- sample(c(1:N),1)
  #     Enstat[,i.n] <- Enstat[,i.new]
  #     EnPara[,i.n] <- EnPara[,i.new]
  #   }}
  En.s.t_1_1S <- En.s.t_1S
}


# 4.0 比较之前的程序  断点数=2 -------------------------------------------------------------
En.s.t_1_2S <- Enstat_2S # 表示前一时刻的状态向量 x_(t-1)
En.s.t_2S <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
En.lik_2S <- rep(NA,N)

#Def for Bootstrap Resampling.
Mid.En.s.t_1_2S <- En.s.t_1_2S
Mid.En.s.t_2S   <- En.s.t_2S
Mid.En.lik_2S   <- En.lik_2S
Mid.EnPara_2S   <- EnPara_2S
for (i.t in 1:(length(t.obs)-1)){   # Loop of time. i.t=1 i.t=2
  t1 <- t.obs[i.t]; t2 <- t.obs[i.t+1]
  sta.obs <- out[t2+1,-1] *(1+ sig.eps*rnorm(n))
  #sta.obs <- out[t2+1,-1] + sig.eps*rnorm(n) #Perturbed observations.
  if(0){
    #i.t < length(t.obs)/2
    #后面把*9改成*18
    sig.eps <- (10-(i.t-1)/length(t.obs)*9)*sig.eps0}else{ 
      sig.eps <- sig.eps0}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Loop_1:       Propagate          #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # 在此步骤中，把状态向量(En.s.t_1) propagate到新时刻(En.s.t)
  for (i.n in 1:N){
    k_0 <- EnPara_2S[1,i.n]
    para_0 <- EnPara_2S[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    }else{
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    # 这是原来的 p4evo0 <- c(k_0,0,para_0[2*(1:k_0)+1], n, para_0[2*(0:k_0)+2],0.0)
    sta.ini <- En.s.t_1_2S[,i.n]
    sta.fore0 <- ode.1D(sta.ini, t1:t2, model2, p4evo0, method = "adams", 
                        names = c("C"), dimens = n)[t2-t1+1,1:n+1]
    lik0 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore0[Obs.site])
    En.lik_2S[i.n] <- lik0
    En.s.t_2S[,i.n] <- sta.fore0
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Loop_2:       Bootstrap          #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  if(max(En.lik_2S) != 0){
    Fathers <- sample(c(1:N),N,En.lik_2S,replace = TRUE)
  }else{Fathers <- c(1:N)}
  #Fathers <- c(1:N)
  # 重采样 代替 注意要用中间变量！！！！ 
  # 只考虑一半
  II2 <-  which(En.lik_2S < median(En.lik_2S)) 
  for (i.n in II2) {
    Mid.En.s.t_1_2S[,i.n] <- En.s.t_1_2S[,Fathers[i.n]]
    Mid.En.s.t_2S[,i.n]   <- En.s.t_2S[,Fathers[i.n]]
    Mid.En.lik_2S[i.n]   <- En.lik_2S[Fathers[i.n]]
    Mid.EnPara_2S[,i.n]   <- EnPara_2S[,Fathers[i.n]]
  }
  Mid.En.s.t_1_2S[,setdiff(c(1:N),II2)] <- En.s.t_1_2S[,setdiff(c(1:N),II2)]
  Mid.En.s.t_2S[,setdiff(c(1:N),II2)]   <- En.s.t_2S[,setdiff(c(1:N),II2)]
  Mid.En.lik_2S[setdiff(c(1:N),II2)]   <- En.lik_2S[setdiff(c(1:N),II2)]
  Mid.EnPara_2S[,setdiff(c(1:N),II2)]   <- EnPara_2S[,setdiff(c(1:N),II2)] 
  #替换完成后，再把矩阵都代回去
  En.s.t_1_2S <- Mid.En.s.t_1_2S
  En.s.t_2S   <- Mid.En.s.t_2S
  En.lik_2S   <- Mid.En.lik_2S
  EnPara_2S   <- Mid.EnPara_2S
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Loop_3:  RJMCMC Re-sampling      #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  for (i.n in 1:N){   # Loop of ensemble member. i.n = 2
    #因为第一步的重采样，重复第一步的输入参数步骤
    k_0 <- EnPara_2S[1,i.n]
    para_0 <- EnPara_2S[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    }else{
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    #接下来就是提议新的参数
    prpd.index <- sample(c(3,4),1)
    #prpd.index <- Prop.Index(k_0,c=0.3,lambda=lambda_p,kmax=kmax)$index
    propose1 <- Prp.Para(para_0,prpd.index) 
    para_1 <- propose1$para; print(para_1) # 这时候要定义函数了，提一个新的函数
    k_1 <- para_1[1]
    if(k_1 == 0){
      p4evo1 <- c(k_1,0,n,para_1[1+k_1+(1:(k_1+1))],0.0)
    }else{
      p4evo1 <- c(k_1,0,para_1[1+(1:(k_1))], n, para_1[1+k_1+(1:(k_1+1))],0.0)
    }
    sta.ini <- En.s.t_1_2S[,i.n]
    sta.fore1 <- ode.1D(sta.ini, t1:t2, model2, p4evo1, method = "adams", 
                        names = c("C"), dimens = n)[t2-t1+1,1:n+1]
    
    
    
    lik1 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore1[Obs.site])
    
    # 比较 En.lik[i.n] 和 lik0
    
    #注意两点，一个是雅各比，一个是初始的向量别错了
    
    #然后再加if
    
    #~~~~~~~~~~~~~~~~End loop # ~~~~~~~~~~
    
    #原来未更新的部分
    lik0 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore0[Obs.site])
    lik1 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore1[Obs.site])
    if((lik1*propose1$ratio)>lik0*runif(1,0,1)){
      EnPara_2S[,i.n] <- rep(NA,2*kmax+2)
      EnPara_2S[1:(2*k_1+2),i.n] <- para_1
      sta.fore0 <- sta.fore1
    }
    # K*(obs- perturbed forecast.)
    En.s.t_2S[,i.n] <- sta.fore0
    if(kalmangain){En.s.t_2S[,i.n] <- sta.fore0 + 
      K.true %*% as.matrix(sta.obs[Obs.site] - sta.fore0[Obs.site] - sig.eps*rnorm(m))
    }
  }
  # for (i.n in 1:N){while(max(Enstat[,i.n])> 10){
  #     i.new <- sample(c(1:N),1)
  #     Enstat[,i.n] <- Enstat[,i.new]
  #     EnPara[,i.n] <- EnPara[,i.new]
  #   }}
  En.s.t_1_2S <- En.s.t_2S
}


# 5.0 断点数==3 --------------------------------------------------------------
En.s.t_1_3S <- Enstat_3S # 表示前一时刻的状态向量 x_(t-1)
En.s.t_3S <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
En.lik_3S <- rep(NA,N)

#Def for Bootstrap Resampling.
Mid.En.s.t_1_3S <- En.s.t_1_3S
Mid.En.s.t_3S   <- En.s.t_3S
Mid.En.lik_3S   <- En.lik_3S
Mid.EnPara_3S   <- EnPara_3S
for (i.t in 1:(length(t.obs)-1)){   # Loop of time. i.t=1 i.t=2
  t1 <- t.obs[i.t]; t2 <- t.obs[i.t+1]
  sta.obs <- out[t2+1,-1] *(1+ sig.eps*rnorm(n))
  #sta.obs <- out[t2+1,-1] + sig.eps*rnorm(n) #Perturbed observations.
  if(0){
    #i.t < length(t.obs)/2
    #后面把*9改成*18
    sig.eps <- (10-(i.t-1)/length(t.obs)*9)*sig.eps0}else{ 
      sig.eps <- sig.eps0}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Loop_1:       Propagate          #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # 在此步骤中，把状态向量(En.s.t_1) propagate到新时刻(En.s.t)
  for (i.n in 1:N){
    k_0 <- EnPara_3S[1,i.n]
    para_0 <- EnPara_3S[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    }else{
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    # 这是原来的 p4evo0 <- c(k_0,0,para_0[2*(1:k_0)+1], n, para_0[2*(0:k_0)+2],0.0)
    sta.ini <- En.s.t_1_3S[,i.n]
    sta.fore0 <- ode.1D(sta.ini, t1:t2, model2, p4evo0, method = "adams", 
                        names = c("C"), dimens = n)[t2-t1+1,1:n+1]
    lik0 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore0[Obs.site])
    En.lik_3S[i.n] <- lik0
    En.s.t_3S[,i.n] <- sta.fore0
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Loop_2:       Bootstrap          #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  if(max(En.lik_3S) != 0){
    Fathers <- sample(c(1:N),N,En.lik_3S,replace = TRUE)
  }else{Fathers <- c(1:N)}
  #Fathers <- c(1:N)
  # 重采样 代替 注意要用中间变量！！！！ 
  # 只考虑一半
  II2 <-  which(En.lik_3S < median(En.lik_3S)) 
  for (i.n in II2) {
    Mid.En.s.t_1_3S[,i.n] <- En.s.t_1_3S[,Fathers[i.n]]
    Mid.En.s.t_3S[,i.n]   <- En.s.t_3S[,Fathers[i.n]]
    Mid.En.lik_3S[i.n]   <- En.lik_3S[Fathers[i.n]]
    Mid.EnPara_3S[,i.n]   <- EnPara_3S[,Fathers[i.n]]
  }
  Mid.En.s.t_1_3S[,setdiff(c(1:N),II2)] <- En.s.t_1_3S[,setdiff(c(1:N),II2)]
  Mid.En.s.t_3S[,setdiff(c(1:N),II2)]   <- En.s.t_3S[,setdiff(c(1:N),II2)]
  Mid.En.lik_3S[setdiff(c(1:N),II2)]   <- En.lik_3S[setdiff(c(1:N),II2)]
  Mid.EnPara_3S[,setdiff(c(1:N),II2)]   <- EnPara_3S[,setdiff(c(1:N),II2)] 
  #替换完成后，再把矩阵都代回去
  En.s.t_1_3S <- Mid.En.s.t_1_3S
  En.s.t_3S   <- Mid.En.s.t_3S
  En.lik_3S   <- Mid.En.lik_3S
  EnPara_3S   <- Mid.EnPara_3S
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Loop_3:  RJMCMC Re-sampling      #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  for (i.n in 1:N){   # Loop of ensemble member. i.n = 2
    #因为第一步的重采样，重复第一步的输入参数步骤
    k_0 <- EnPara_3S[1,i.n]
    para_0 <- EnPara_3S[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    }else{
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    #接下来就是提议新的参数
    prpd.index <- sample(c(3,4),1)
    propose1 <- Prp.Para(para_0,prpd.index) 
    para_1 <- propose1$para; print(para_1) # 这时候要定义函数了，提一个新的函数
    k_1 <- para_1[1]
    if(k_1 == 0){
      p4evo1 <- c(k_1,0,n,para_1[1+k_1+(1:(k_1+1))],0.0)
    }else{
      p4evo1 <- c(k_1,0,para_1[1+(1:(k_1))], n, para_1[1+k_1+(1:(k_1+1))],0.0)
    }
    sta.ini <- En.s.t_1_3S[,i.n]
    sta.fore1 <- ode.1D(sta.ini, t1:t2, model2, p4evo1, method = "adams", 
                        names = c("C"), dimens = n)[t2-t1+1,1:n+1]
    
    
    
    lik1 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore1[Obs.site])
    
    # 比较 En.lik[i.n] 和 lik0
    
    #注意两点，一个是雅各比，一个是初始的向量别错了
    
    #然后再加if
    
    #~~~~~~~~~~~~~~~~End loop # ~~~~~~~~~~
    
    #原来未更新的部分
    lik0 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore0[Obs.site])
    lik1 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore1[Obs.site])
    if((lik1*propose1$ratio)>lik0*runif(1,0,1)){
      EnPara_3S[,i.n] <- rep(NA,2*kmax+2)
      EnPara_3S[1:(2*k_1+2),i.n] <- para_1
      sta.fore0 <- sta.fore1
    }
    # K*(obs- perturbed forecast.)
    En.s.t_3S[,i.n] <- sta.fore0
    if(kalmangain){En.s.t_3S[,i.n] <- sta.fore0 + 
      K.true %*% as.matrix(sta.obs[Obs.site] - sta.fore0[Obs.site] - sig.eps*rnorm(m))
    }
  }
  # for (i.n in 1:N){while(max(Enstat[,i.n])> 10){
  #     i.new <- sample(c(1:N),1)
  #     Enstat[,i.n] <- Enstat[,i.new]
  #     EnPara[,i.n] <- EnPara[,i.new]
  #   }}
  En.s.t_1_3S <- En.s.t_3S
}
# 画图 ----------------------------------------------------------------------
En.s.t_1 <- Res1_Rj$EnS
En.s.t_1_2S <- Res1_2S$EnS
par(mfrow=c(1,1),mar=c(2,2,0,0.5)+2.5)
plot(m.true,type ='l',
     #ylim = c(min(En.s.t_1)-0.1,6),
     ylim = c(-7,6),
     cex.axis = 2, cex.lab =2, cex.main =2,
     main='The pointwise posterior 80% credible intervals',xlab = 'Locations',ylab = 'State')

plot.PI(En.s.t_1,'blue')
plot.PI(En.s.t_1_1S,'green')
plot.PI(En.s.t_1_2S,'red')
plot.PI(En.s.t_1_3S,'purple')
i.t
mean_rj  <- rowMeans(En.s.t_1);lines(mean_rj,col='green',lwd=3)
mean_1S  <- rowMeans(En.s.t_1_1S);lines(mean_1S,col='yellow',lwd=3)
mean_2S  <- rowMeans(En.s.t_1_2S);lines(mean_2S,col='red',lwd=3)
mean_3S  <- rowMeans(En.s.t_1_3S);lines(mean_3S,col='purple',lwd=3)
lines(m.true, lwd = 3)
legend("topleft",legend=c("The new method","The true velocities"),
       col=c("blue","black"), bty = "n", 
       cex=2,    y.intersp = 0.5, x.intersp = 0.1, xjust = 0,           
       lty=1,lwd=3)
legend("bottomleft",legend=c("pMCMC-1","pMCMC-2","pMCMC-3"),
       col=c("green","red","purple"), bty = "n", 
       cex=2,    y.intersp = 0.5, x.intersp = 0.1, xjust = -0.2,           
       lty=1,lwd=3)
#计算MSE
mean((mean_rj-m.true)**2)
mean((mean_1S-m.true)**2)
mean((mean_2S-m.true)**2)
mean((mean_3S-m.true)**2) #可能是与粒子数有关？？

# 6.0 增加预测的比较 -------------------------------------------------------------

#首先看逻辑
# 预测的话先有 最后时刻的ensemble of states parameters
# 和真实的 states, parameters
En.s <- En.s.t_1
En.p <- EnPara

En.s <- En.s.t_1_1S
En.p <- EnPara_1S
s.true <- out[length(times),1+(1:n)]
p.true <- para
times.p <- times[length(times)]+c(0:50)
out.pre <- ode.1D(s.true, times.p, model2, p.true, method = "adams", 
                  names = c("C"), dimens = n)[length(times.p),1+(1:n)]
En.predict <- matrix(NA,ncol = N, nrow <- n)
for(i.n in  1:N){
  k_0 <- En.p[1,i.n]
  para_0 <- En.p[1:(2*k_0+2),i.n]; print(para_0)
  if(k_0 == 0){
    p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
  }else{
    p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
  }
  En.predict[,i.n] <- ode.1D(En.s[,i.n], times.p, model2, p4evo0, method = "adams", 
                             names = c("C"), dimens = n)[length(times.p),1+(1:n)]
}
plot(out.pre, type='l')
for(i.n in 1:N){lines(En.predict[,i.n],col='red')}
plot.PI(En.predict,'blue')
lines(out.pre,lwd=3)

#
En.s <- En.s.t_1[,which(EnPara[1,]==2)]
En.p <- EnPara[,which(EnPara[1,]==2)]
Pre_En <- function(En.s,En.p,s.true,p.true,times.p){
  out.pre <- ode.1D(s.true, times.p, model2, p.true, method = "adams", 
                    names = c("C"), dimens = n)[length(times.p),1+(1:n)]
  N <- dim(En.s)[2]
  En.predict <- matrix(NA,ncol = N, nrow <- n)
  for(i.n in  1:N){
    k_0 <- En.p[1,i.n]
    para_0 <- En.p[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    }else{
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    En.predict[,i.n] <- ode.1D(En.s[,i.n], times.p, model2, p4evo0, method = "adams", 
                               names = c("C"), dimens = n)[length(times.p),1+(1:n)]
  }
  plot(out.pre, type='l',ylab = 'The state',xlab = 'Locations',cex.axis = 1.5, cex.lab = 1.5)
  for(i.n in 1:N){lines(En.predict[,i.n],col='lightblue')}
  plot.PI(En.predict,'blue')
  lines(out.pre,lwd=3,col='red')
  mean_pre  <- rowMeans(En.predict);
  lines(mean_pre,col='darkblue',lwd=2)
  MSPE <- mean((mean_pre-out.pre)**2)
  return(list('En'=En.predict,'MSPE'=MSPE,'Pre'=mean_pre))
}


Res_r <- Pre_En(En.s.t_1,EnPara,out[length(times),1+(1:n)],para,times.p)
Res_1 <- Pre_En(En.s.t_1_1S,EnPara_1S, out[length(times),1+(1:n)],para,times.p)
Res_2 <- Pre_En(En.s.t_1_2S,EnPara_2S,out[length(times),1+(1:n)],para,times.p)
Res_3 <-  Pre_En(En.s.t_1_3S,EnPara_3S,out[length(times),1+(1:n)],para,times.p)
Res_rB <- Pre_En(En.s.t_1[,which(EnPara[1,]==2)],EnPara[,which(EnPara[1,]==2)],
                 out[length(times),1+(1:n)],para,times.p)
Res_rB <- Pre_En(En.s.t_1[,which(EnPara[1,]==2)],EnPara[,which(EnPara[1,]==2)],
                 out[length(times),1+(1:n)],para,times.p)
Res_r$MSPE
Res_2$MSPE
Res_rB$MSPE
Res_1$MSPE
Res_3$MSPE


#以上程序是能跑通的，结果也很好 2024.11.13

#以上的程序都整理在了11.14的那个func.R中
#11.14的两个程序跑出的结果有问题 还需再看 或者试试多一点观测


# 函数，输入一个参数返回每一点的速度值 ------------------------------------------------------


# para_0
# k_0 <- para_0[1]
# c.x <- c(1,para_0[1+(1:(k_0))],n)
# c.v <- para_0[1+k_0+(1:(k_0+1))]
# v_vec <- rep(NA, n)
# for (i.k in 1:(k_0+1)) {
#   v_vec[c.x[i.k]:c.x[i.k+1]] <- rep(c.v[i.k],c.x[i.k+1]-c.x[i.k]+1)}

#~~~~~~~~~~~~~~~~~~~画速度参数的预测
V_Points <- function(para_0,n){
  k_0 <- para_0[1]
  c.x <- c(1,para_0[1+(1:(k_0))],n)
  c.v <- para_0[1+k_0+(1:(k_0+1))]
  v_vec <- rep(NA, n)
  for (i.k in 1:(k_0+1)) {
    v_vec[c.x[i.k]:c.x[i.k+1]] <- rep(c.v[i.k],c.x[i.k+1]-c.x[i.k]+1)}
  
  return(v_vec)
}
V_Plot <- function(En.P,col1,col2){
  V_matrix <-matrix(NA,ncol = N,nrow = n)
  for (i.n in 1:N) {V_matrix[,i.n] <- V_Points(En.P[,i.n],n)}
  plot.PI(V_matrix,col1)
  lines(rowMeans(V_matrix),lwd=3,col=col2)
}

plot(V_Points(c(K,101,251,c.v.true),n),col='black',ylim=c(0.1,2),
     cex.axis = 1.5, cex.lab = 1.5,
     xlab ='Locations',ylab ='The pointwise velocities')
V_Plot(EnPara,'blue','darkblue')
V_Plot(EnPara_2S,'orange','red')
V_Plot(EnPara_1S,'green','darkgreen')
V_Plot(EnPara_3S,'pink','purple')
legend("topright",legend=c("The new method","pMCMC-1","pMCMC-2","pMCMC-3","The true velocities"),
       col=c("blue","green","red","purple","black"),                
       lty=1,lwd=3)
points(V_Points(c(K,101,251,c.v.true),n))
# V_matrix <-matrix(NA,ncol = N,nrow = n)
# for (i.n in 1:N) {V_matrix[,i.n] <- V_Points(EnPara[,i.n],n)}
# plot(rowMeans(V_matrix),type='l',lwd=2,col='darkblue',ylim=c(0.1,0.8))
# plot.PI(V_matrix,'blue')
# lines(V_Points(c(K,101,251,c.v.true),n),col='red',lwd=2)


#这个程序就不要动了 放在论文中
# 把观测数据从1/10改成1/5就是文章中的第80个观测点的结论。



# 查看前 5 列
print(Indice.m[, 1:5])


library(tidyr)  # 用于 pivot_longer
library(dplyr)  # 用于数据处理

# 转换为数据框并添加类别标签
df <- as.data.frame(t(Indice.m))  # 转置矩阵（60x4）
colnames(df) <- c("Count_0", "Count_1", "Count_2", "Count_3")  # 列名
df$Time <- 1:60  # 添加时间列

# 转换为长格式
df_long <- df %>%
  pivot_longer(
    cols = -Time,
    names_to = "Category",
    values_to = "Count"
  )

# 查看整理后的数据
head(df_long)


library(ggplot2)

ggplot(df_long, aes(x = Time, y = Count, color = Category)) +
  geom_line(lwd = 1) +  # 绘制曲线
  geom_point(size = 1.5, alpha = 0.6) +  # 添加点（可选）
  scale_color_manual(
    values = c("Count_0" = "blue", "Count_1" = "red", "Count_2" = "green", "Count_3" = "purple"),
    labels = c("0", "1", "2", "3")  # 图例标签
  ) +
  coord_cartesian(ylim = c(0, 80))+
  labs(
    x = "Time step t",
    y = "Number of particles",
    color = "Model"
  ) +
  theme_minimal(base_size = 14) +  # 基础字体大小设为14
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # 标题字体加大
    axis.title = element_text(size = 16),  # 坐标轴标题大小
    axis.text = element_text(size = 14),   # 坐标轴文字大小
    legend.title = element_text(size = 16), # 图例标题大小
    legend.text = element_text(size = 14),  # 图例文字大小
    legend.position = "right"  # 图例在右侧
  )






# 图2 ----------------------------------------------------------------------
# 统计 [80,120] 区间内的元素个数
count_75_125 <- sum(EnPara >= 75 & EnPara <= 125, na.rm = TRUE)

# 统计 [230,270] 区间内的元素个数
count_225_275 <- sum(EnPara >= 225 & EnPara <= 275, na.rm = TRUE)

# 输出结果
cat("Number of elements in [80,120]:", count_80_120, "\n")
cat("Number of elements in [230,270]:", count_230_270, "\n")



# 查看前 5 列
print(Loca.m[, 1:5])


library(tidyr)  # 用于 pivot_longer
library(dplyr)  # 用于数据处理

# 转换为数据框并添加类别标签
df <- as.data.frame(t(Loca.m))  # 转置矩阵（60x4）
colnames(df) <- c( "Loca_1", "Loca_2")  # 列名
df$Time <- 1:60  # 添加时间列

# 转换为长格式
df_long <- df %>%
  pivot_longer(
    cols = -Time,
    names_to = "Category",
    values_to = "Count"
  )

# 查看整理后的数据
head(df_long)


library(ggplot2)

ggplot(df_long, aes(x = Time, y = Count, color = Category)) +
  geom_line(lwd = 1) +  # 绘制曲线
  geom_point(size = 1.5, alpha = 0.6) +  # 添加点（可选）
  scale_color_manual(
    values = c("Loca_1"= "lightblue", "Loca_2"= "orange"),
    labels = c("Loca_1", "Loca_2")  # 图例标签
  ) +
  coord_cartesian(ylim = c(0, 80))+
  labs(
    x = "Time step t",
    y = "Number of particles",
    color = "spatial point"
  ) +
  theme_minimal(base_size = 14) +  # 基础字体大小设为14
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # 标题字体加大
    axis.title = element_text(size = 16),  # 坐标轴标题大小
    axis.text = element_text(size = 14),   # 坐标轴文字大小
    legend.title = element_text(size = 16), # 图例标题大小
    legend.text = element_text(size = 14),  # 图例文字大小
    legend.position = "right"  # 图例在右侧
  )

