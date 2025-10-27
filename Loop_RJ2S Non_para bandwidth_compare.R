R_path <- 'E:/Study-2024/09.RJ_2stage_resample_Oct24/4. Second_review CODE/'

#source(paste0(R_path, '.R'))
source(paste0(R_path, 'RJ2S_basic function and setting.R'))
source(paste0(R_path, 'Func_RJ2S  add_perturbation.R'))
source(paste0(R_path, 'Func_RJ2S  Para KDE.R'))
source(paste0(R_path, 'Func_RJ2S  test acceptance and ess.R'))
source(paste0(R_path, 'Func_RJ2S Para estimation proposal.R'))
source(paste0(R_path, 'Func_RJ2S Non_para estimation proposal.R'))

# The loop for different bandwiths.

# 先选取N=50

# the three  processes are 1 time, 5 times, and 10 times.

N_Loop <- 20
Res_ALL <- matrix(NA,ncol = 8,nrow = N_Loop)
Mtx.ESS_1 <-matrix(NA,ncol= N_Loop,nrow = length(t.obs)-1)
Mtx.ESS_2 <-matrix(NA,ncol= N_Loop,nrow = length(t.obs)-1)
Mtx.ESS_3 <-matrix(NA,ncol= N_Loop,nrow = length(t.obs)-1)
Mtx.ACP_1 <-matrix(NA,ncol= N_Loop,nrow = length(t.obs)-1)
Mtx.ACP_2 <-matrix(NA,ncol= N_Loop,nrow = length(t.obs)-1)
Mtx.ACP_3 <-matrix(NA,ncol= N_Loop,nrow = length(t.obs)-1)

#单跑一次的，固定样本数，固定观测点数的N_Loop次重复
#i.all <- 10
for(i.all in 1:N_Loop){
  
  set.seed(i.all+1006)
  lambda_p <- 2 #泊松先验的均值
  kmax <- 3
  N <- 60
  times <- 0:600
  sig.eps0 <- 0.2 #2
  
  L <- 400 #
  n <- L+1
  #Gamma分布的参数
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
  
  #对间断点的位置选择，从[0,L]中抽取(2k+1)个点，选择中间的k个点
  
  set.seed(i.all+1)
  sample.locs(5,100)
  
  #k+1段速度的先验分布
  #rgamma(k+1,shape = 2,rate = 1)
  
  # 则可以对任意一个k，抽取相应的参数
  
  #从先验分布中抽取样本 
  EnPara <- matrix(NA, ncol = N, nrow = 2*kmax+2)
  EnPara[1,] <- sample(c(0:kmax),N,replace=TRUE,dpois(c(0:kmax), lambda = lambda_p)) 
  set.seed(i.all+111)
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
  
  
  # 03.初始条件 -----------------------------------------------------------------
  
  
  #对于初始条件（真实条件）
  #~~~~ 3.1 ~~~~ The same setting with the previous part.
  
  
  
  #Pre work
  ### true forecast distribution
  locs=1:n
  distmat=rdist(locs,locs)  #
  scale=20
  sig.eps=.2
  kappa=3
  
  
  
  distmat=rdist(1:n,1:n); scale=100; sig.eps=.5 ;kappa=3; df.t=2
  
  Obs.site <- sort(sample(c(1:L),ceiling(n/10))) #观测点的位置 改观测数
  #Obs.site <- sort(sample(c(1:L),20)) 
  m.obs <- length(Obs.site)                     #观测数
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
  
  set.seed(i.all+1234)
  ##!!!!!初始向量的矩阵
  
  #Enstat <- forecast.mean.true+forecast.chol.true%*%matrix(rnorm(n*N),nrow=n)+3
  #Enstat <- matrix(NA, nrow = n, ncol = N)
  #for(i in 1:N){Enstat[,i] <- Cini+rep(rnorm(1,0,0.2),n)}
  
  #plot(Enstat[,1])
  #for(i in 1:N){lines(Enstat[,i])}
  
  
  #目前有N个state和parameter的粒子
  
  
  
  # 04. DA process ----------------------------------------------------------
  
  
  
  
  
  
  
  
  # 05. 真实的参数以及真实参数下的时空场 ----------------------------------------------------
  #Para.true <- c(2,0.5,0.7,0.3,101,301)
  #x <- seq(0,100, by=0.25) #Loc
  #n <- n;print(n)
  #m <- 200 # number of observations
  
  #dx <- x[2]-x[1]
  
  t.obs <- seq(0,times[length(times)], by = 10)
  
  
  
  
  
  set.seed(i.all+1)
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
  Enstat_NP  <- Enstat # Classical MCMC   1个断点                            #
  Enstat_PE  <- Enstat # Classical MCMC 同真实数据一样2个断点， 3个steps     #
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  
  plot(Enstat[,1])
  for(i in 1:N){lines(Enstat[,i])}
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~比较的参数集合~~~~~~~~~~~~~~~~~~~#
  #
  
  # Classical MCMC   1个断点                                                 #
  EnPara_NP <- EnPara
  
  # Classical MCMC 同真实数据一样2个断点， 3个steps                          #
  #
  EnPara_PE <- EnPara
  
  
  
  
  
  
  # 06.Resampling 重采样步骤 ----------------------------------------------------------------
  
  # 6.0 先定义提议函数.
  
  
  
  
  #prpd.index <- Prop.Index(k,c=0.3,lambda=lambda_p,kmax=kmax)$index
  #Prp.Para(para_0,prpd.index)
  
  #现在先有三个矩阵，上一时刻的状态向量，这一时刻的状态向量，以及参数的矩阵
  set.seed(i.all+1)
  
  
  
  
  
  En.s.t_1 <- Enstat # 表示前一时刻的状态向量 x_(t-1)
  En.s.t <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
  En.lik <- rep(NA,N)
  
  #Def for Bootstrap Resampling.
  Mid.En.s.t_1 <- En.s.t_1
  Mid.En.s.t   <- En.s.t
  Mid.En.lik   <- En.lik
  Mid.EnPara   <- EnPara
  
  # 1.0 The main loop. -------------------------------------------------------
  
  
  
  
  
  
  # 3.0 non parameter estimation -------------------------------------------------------------
  En.s.t_1_NP <- Enstat_NP # 表示前一时刻的状态向量 x_(t-1)
  En.s.t_NP <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
  En.lik_NP <- rep(NA,N)
  
  #Def for Bootstrap Resampling.
  Mid.En.s.t_1_NP <- En.s.t_1_NP
  Mid.En.s.t_NP   <- En.s.t_NP
  Mid.En.lik_NP   <- En.lik_NP
  Mid.EnPara_NP   <- EnPara_NP
  set.seed(1014+i.all)
  Res1_NP <- DA_Process_NonP(Enstat_NP,EnPara_NP,1,0.1) #2nd
  En.s.t_1_NP <- Res1_NP$EnS
  EnPara_NP <- Res1_NP$EnP
  
  En.s.t_1_PE <- Enstat_PE # 表示前一时刻的状态向量 x_(t-1)
  En.s.t_PE <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
  En.lik_PE <- rep(NA,N)
  
  Mid.En.s.t_1_PE <- En.s.t_1_PE
  Mid.En.s.t_PE   <- En.s.t_PE
  Mid.En.lik_PE   <- En.lik_PE
  Mid.EnPara_PE   <- EnPara_PE
  set.seed(1014+i.all)
  Res1_PE <- DA_Process_NonP(Enstat_PE,EnPara_PE,1,10) #3rd
  En.s.t_1_PE <- Res1_PE$EnS
  EnPara_PE <- Res1_PE$EnP  
  
  set.seed(1014+i.all)
  Res1_Rj <- DA_Process_NonP(Enstat,EnPara,1,1)
  
  EnPara    <- Res1_Rj$EnP
  En.s.t_1  <- Res1_Rj$EnS
  
  mean_rj  <- rowMeans(En.s.t_1);#lines(mean_rj,col='green',lwd=3)
  mean_NP  <- rowMeans(En.s.t_1_NP);#lines(mean_NP,col='yellow',lwd=3)
  mean_PE  <- rowMeans(En.s.t_1_PE);#lines(mean_PE,col='red',lwd=3)
  
  
  #计算MSE
  m.true <- out[length(times),1+(1:n)]
  Res_ALL[i.all,1] <- mean((mean_rj-m.true)**2)
  Res_ALL[i.all,2] <- mean((mean_NP-m.true)**2)
  Res_ALL[i.all,3] <- mean((mean_PE-m.true)**2)
  
  
  # 6.0 增加预测的比较 -------------------------------------------------------------
  
  #首先看逻辑
  # 预测的话先有 最后时刻的ensemble of states parameters
  # 和真实的 states, parameters
  En.s <- En.s.t_1
  En.p <- EnPara
  
  En.s <- En.s.t_1_NP
  En.p <- EnPara_NP
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
  
  
  
  Res_r <- Pre_En(En.s.t_1,EnPara,out[length(times),1+(1:n)],para,times.p)
  Res_1 <- Pre_En(En.s.t_1_NP,EnPara_NP, out[length(times),1+(1:n)],para,times.p)
  Res_2 <- Pre_En(En.s.t_1_PE,EnPara_PE,out[length(times),1+(1:n)],para,times.p)
  
  
  Res_ALL[i.all,5] <- Res_r$MSPE
  Res_ALL[i.all,7] <- Res_2$MSPE
  Res_ALL[i.all,6] <- Res_1$MSPE
  
  Mtx.ESS_1[,i.all] <- Res1_Rj$Acp[,10]
  Mtx.ACP_1[,i.all] <- Res1_Rj$Acp[,9]/N
  
  Mtx.ESS_2[,i.all] <- Res1_NP$Acp[,10]
  Mtx.ACP_2[,i.all] <- Res1_NP$Acp[,9]/N
  
  Mtx.ESS_3[,i.all] <- Res1_PE$Acp[,10]
  Mtx.ACP_3[,i.all] <- Res1_PE$Acp[,9]/N
  
}

colMeans(Res_ALL)
ESS_1RJ <- rowMeans(Mtx.ESS_1);ESS_2NP <- rowMeans(Mtx.ESS_2);ESS_3PE <- rowMeans(Mtx.ESS_3)
ACP_1RJ <- rowMeans(Mtx.ACP_1);ACP_2NP <- rowMeans(Mtx.ACP_2);ACP_3PE <- rowMeans(Mtx.ACP_3)
library(ggplot2)
library(tidyr)

# Create data frame
time <- 1:60
data <- data.frame(
  Time = time,
  #Preset_Distribution = c(27.48723, 27.52074, 27.11971, 26.05647, 25.41210, 26.18622, 25.97759, 25.99319, 25.98209, 28.09217, 25.82977, 25.41876, 25.89369, 26.37012, 27.14707, 27.65339, 26.27867, 26.46956, 24.95643, 26.71892, 26.30998, 25.20318, 26.26254, 24.98950, 25.87268, 25.86762, 24.87031, 25.04287, 25.67605, 25.51372, 26.51410, 27.46069, 26.52473, 25.72909, 25.58654, 25.29118, 25.10681, 26.76137, 27.90858, 25.98649, 27.95671, 26.40052, 25.69673, 25.98161, 25.97946, 26.98091, 27.38491, 26.67978, 25.01190, 26.03781, 27.74172, 26.83884, 26.43189, 27.17804, 27.58353, 26.64462, 26.45022, 24.78959, 27.03945, 26.27829),
  #Nonparametric_Estimation = c(29.05354, 27.52514, 28.13478, 26.87407, 26.81944, 25.67685, 26.04109, 26.04457, 26.64802, 24.47181, 25.32974, 28.14573, 24.99087, 26.40190, 27.82323, 24.46112, 26.79432, 24.81219, 25.36175, 25.16935, 26.23749, 28.51134, 26.56864, 26.18657, 27.35499, 26.67095, 25.70665, 26.77505, 28.00448, 24.88653, 25.15180, 25.03451, 24.94383, 28.70464, 26.61396, 25.94946, 25.10737, 26.20445, 25.92530, 26.83071, 26.27338, 27.81369, 27.62229, 28.26863, 26.56641, 26.29265, 27.22119, 26.61098, 27.76441, 26.67274, 26.57731, 27.40297, 27.41839, 26.39433, 27.48980, 26.72722, 26.21521, 26.32102, 25.80490, 26.82570),
  #Parametric_Estimation = c(26.99103, 27.06380, 26.51209, 26.31654, 27.23992, 26.07361, 26.69149, 25.87445, 26.20207, 25.56512, 25.68182, 27.53572, 26.16147, 25.75235, 26.67771, 26.51088, 26.06696, 26.80573, 26.36600, 25.54872, 26.44845, 25.78267, 26.28029, 27.26288, 25.87598, 25.12601, 26.27198, 26.41723, 26.66368, 25.52146, 25.90565, 24.88061, 25.95565, 25.98976, 25.86181, 25.00147, 25.69163, 26.14205, 26.68700, 26.41599, 26.65495, 25.07200, 27.37560, 26.17269, 25.97759, 25.97783, 26.05663, 26.75010, 27.26916, 26.04473, 25.60774, 25.92472, 27.21170, 25.32477, 24.63850, 25.75034, 26.63036, 26.42073, 26.23881, 24.46844)
  Preset_Distribution =      rowMeans(Mtx.ESS_1[,c(1:18)]),
  Nonparametric_Estimation = rowMeans(Mtx.ESS_2[,c(1:18)]),
  Parametric_Estimation =    rowMeans(Mtx.ESS_3[,c(1:18)])
)

# Convert to long format
data_long <- pivot_longer(data, cols = -Time, names_to = "Method", values_to = "ESS")

# Create plot
ggplot(data_long, aes(x = Time, y = ESS, color = Method)) +
  geom_line(lwd = 1.2) +  # Thicker lines
  geom_point(size = 1.5) +
  scale_color_manual(
    values = c("Preset_Distribution" = "#E41A1C", 
               "Nonparametric_Estimation" = "#377EB8", 
               "Parametric_Estimation" = "#4DAF4A"),
    labels = c("Original bandwidths", "0.1 times bandwidths", "10 times bandwidths")
  ) +
  labs(
    title = "ESS over Time for Three Estimation Methods",
    x = "Time",
    y = "ESS",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10)
  )



library(ggplot2)
library(tidyr)

# Create data frame
time <- 1:60
data <- data.frame(
  Time = time,
  #Preset_Distribution = c(0.3983333, 0.4016667, 0.4550000, 0.4941667, 0.4875000, 0.5050000, 0.5125000, 0.4866667, 0.5150000, 0.4866667, 0.5041667, 0.4925000, 0.4675000, 0.5258333, 0.4975000, 0.4866667, 0.4916667, 0.4516667, 0.4675000, 0.4725000, 0.4750000, 0.4525000, 0.4775000, 0.4791667, 0.5058333, 0.4933333, 0.4883333, 0.5241667, 0.5358333, 0.5250000, 0.5025000, 0.5191667, 0.5041667, 0.5283333, 0.5133333, 0.5033333, 0.4883333, 0.4983333, 0.4808333, 0.4825000, 0.4683333, 0.4975000, 0.5158333, 0.4641667, 0.4708333, 0.4900000, 0.4766667, 0.5025000, 0.4608333, 0.5066667, 0.4850000, 0.5433333, 0.5291667, 0.5233333, 0.5550000, 0.5216667, 0.5591667, 0.5716667, 0.5175000, 0.5450000),
  #Nonparametric_Estimation = c(0.10250000, 0.07250000, 0.04333333, 0.04750000, 0.04750000, 0.06583333, 0.05500000, 0.06000000, 0.06333333, 0.03833333, 0.04416667, 0.05833333, 0.05333333, 0.04750000, 0.06666667, 0.04916667, 0.05916667, 0.05083333, 0.04083333, 0.05750000, 0.06000000, 0.04583333, 0.06250000, 0.04583333, 0.05083333, 0.04500000, 0.04583333, 0.05166667, 0.03750000, 0.03916667, 0.04083333, 0.05250000, 0.03500000, 0.03833333, 0.05333333, 0.03666667, 0.05083333, 0.05083333, 0.04250000, 0.05250000, 0.05500000, 0.05416667, 0.05000000, 0.05166667, 0.04333333, 0.05416667, 0.05083333, 0.05333333, 0.04500000, 0.03583333, 0.06000000, 0.04916667, 0.05000000, 0.04083333, 0.05666667, 0.04583333, 0.04250000, 0.05083333, 0.04583333, 0.04666667),
  #Parametric_Estimation = c(0.4316667, 0.3925000, 0.4675000, 0.4866667, 0.4733333, 0.5016667, 0.5300000, 0.4916667, 0.5225000, 0.5008333, 0.4508333, 0.5075000, 0.4966667, 0.5000000, 0.4566667, 0.4408333, 0.4750000, 0.4700000, 0.4916667, 0.4658333, 0.4691667, 0.4683333, 0.4533333, 0.4425000, 0.5125000, 0.4933333, 0.5241667, 0.5358333, 0.5058333, 0.5408333, 0.5475000, 0.5175000, 0.5125000, 0.5325000, 0.5141667, 0.5058333, 0.5391667, 0.5200000, 0.5000000, 0.4791667, 0.4750000, 0.4475000, 0.4783333, 0.4800000, 0.4883333, 0.4508333, 0.4391667, 0.5100000, 0.4375000, 0.4766667, 0.4883333, 0.4775000, 0.5250000, 0.5258333, 0.5216667, 0.5075000, 0.5183333, 0.5066667, 0.4875000, 0.5191667),
  Preset_Distribution =      rowMeans(Mtx.ACP_1),#[,c(1:18)]
  Nonparametric_Estimation = rowMeans(Mtx.ACP_2),
  Parametric_Estimation =    rowMeans(Mtx.ACP_3)
)

# Convert to long format
data_long <- pivot_longer(data, cols = -Time, names_to = "Method", values_to = "Acceptance_Rate")

# Create plot
ggplot(data_long, aes(x = Time, y = Acceptance_Rate, color = Method)) +
  geom_line(lwd = 1.2) +
  geom_point(size = 1.5) +
  scale_color_manual(
    values = c("Preset_Distribution" = "#E41A1C", 
               "Nonparametric_Estimation" = "#377EB8", 
               "Parametric_Estimation" = "#4DAF4A"),
    #labels = c("Original bandwidths", "5 times bandwidths", "10 times bandwidths")
    labels = c("Original bandwidths", "0.1 times bandwidths", "10 times bandwidths")
  ) +
  labs(
    title = "Average Acceptance Rate over Time for Three Estimation Methods",
    x = "Time",
    y = "Average Acceptance Rate",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))  # Format y-axis as percentage


# notes <- 'N=60, Loop=20, bandwidth=1, 0.1, 10.'
# # save
# save(Mtx.ESS_1, Mtx.ESS_2, Mtx.ESS_3,
#      Mtx.ACP_1, Mtx.ACP_2, Mtx.ACP_3,
#      Res_ALL,   notes,
#      file = "BW_N=60.RData")