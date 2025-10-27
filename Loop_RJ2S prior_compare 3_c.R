
library(ggplot2)
library(dplyr)
library(tidyr)

R_path <- 'E:/Study-2024/09.RJ_2stage_resample_Oct24/4. Second_review CODE/'

#source(paste0(R_path, '.R'))
source(paste0(R_path, 'RJ2S_basic function and setting.R'))
source(paste0(R_path, 'Func_RJ2S  add_perturbation.R'))
source(paste0(R_path, 'Func_RJ2S  Para KDE.R'))
source(paste0(R_path, 'Func_RJ2S  test acceptance and ess.R'))
source(paste0(R_path, 'Func_RJ2S Para estimation proposal.R'))
source(paste0(R_path, 'Func_RJ2S Non_para estimation proposal.R'))

# The loop for different paramter

N_Loop <-20 # The repeat time
#N.Vec <- c(40,60,80)
c_vec <- c(0.1,0.3,0.5)
Res_List <- list()

iN = 1; i.all = 1
for (iN in 1:length(c_vec)){
  Res_ALL <- matrix(NA,ncol = 15,nrow = N_Loop) #记录四个方法的两个误差指标
  #单跑一次的，固定样本数，固定观测点数的N_Loop次重复
  Prop.Index<- function(k,c=c_vec[iN],lambda = lambda_p,kmax=kmax){
    # k: is the index of current model
    # c: the adjusting parameter for probability
    # lambda: the p of poisson. 
    if(k<kmax){bk <- c_vec[iN]*min(1,dpois(k+1, lambda = lambda_p)/dpois(k, lambda = lambda_p))}else{bk <- 0}
    if(k>0){dk <- c_vec[iN]*min(1,dpois(k-1, lambda = lambda_p)/dpois(k, lambda = lambda_p))}else{dk <- 0}
    etak <- (1-bk-dk)/2
    phik <- (1-bk-dk)/2
    if(k==0){phik <- 0}else{phik <- (1-bk-dk)/2}
    pp <- c(bk,dk,etak,phik)/sum(c(bk,dk,etak,phik))
    return(list('p'=pp,'index'=sample(c(1:4),1,FALSE,pp)))
  }
  for(i.all in 1:N_Loop){
    
    # 这里是一次模拟的开始
    set.seed(i.all+1006)
    lambda_p <- 2 #泊松先验的均值
    kmax <- 3
    #N <- 120
    N <- 80
    times <- 0:600
    sig.eps0 <- 0.2 #2
    
    L <- 400 #
    n <- L+1
    #Gamma分布的参数
    g_alpha <- 0.40
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
    
    set.seed(100+i.all+1)
    sample.locs(5,100)
    
    #k+1段速度的先验分布
    #rgamma(k+1,shape = 2,rate = 1)
    
    # 则可以对任意一个k，抽取相应的参数
    
    #从先验分布中抽取样本 
    EnPara <- matrix(NA, ncol = N, nrow = 2*kmax+2)
    EnPara[1,] <- sample(c(0:kmax),N,replace=TRUE,dpois(c(0:kmax), lambda = lambda_p)) 
    set.seed(100+i.all+111)
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
    
    set.seed(100+i.all+1234)
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
    
    
    
    
    
    set.seed(100+i.all+1)
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
    # Enstat_1S  <- Enstat # Classical MCMC   1个断点                            #
    # Enstat_2S  <- Enstat # Classical MCMC 同真实数据一样2个断点， 3个steps     #
    # Enstat_3S  <- Enstat # Classical MCMC   3个断点， 4个steps                 # 
    # Enstat_SM  <- Enstat # SMC 方法 同真实数据一样2个断点， 3个steps     #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    
    plot(Enstat[,1])
    for(i in 1:N){lines(Enstat[,i])}
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~比较的参数集合~~~~~~~~~~~~~~~~~~~#
    #
    
    # # Classical MCMC   1个断点                                                 #
    # EnPara_1S <- matrix(NA, ncol = N, nrow = 2*kmax+2)
    # EnPara_1S[1,] <- rep(1,N) 
    # set.seed(100+i.all+111)
    # for(i in 1:N){
    #   kk <- EnPara_1S[1,i]
    #   if(kk>0){
    #     c.x <- sample.locs(kk,L)
    #     EnPara_1S[1+(1:(kk)),i] <-   c.x}
    #   c.v <- rgamma(kk+1,shape = 2,rate = 1)
    #   EnPara_1S[kk+1+(1:(kk+1)),i] <- c.v  
    # }
    # 
    # # Classical MCMC 同真实数据一样2个断点， 3个steps                          #
    # #
    # EnPara_2S <- matrix(NA, ncol = N, nrow = 2*kmax+2)
    # EnPara_2S[1,] <- rep(2,N) 
    # 
    # for(i in 1:N){
    #   kk <- EnPara_2S[1,i]
    #   if(kk>0){
    #     c.x <- sample.locs(kk,L)
    #     EnPara_2S[1+(1:(kk)),i] <-   c.x}
    #   c.v <- rgamma(kk+1,shape = 2,rate = 1)
    #   EnPara_2S[kk+1+(1:(kk+1)),i] <- c.v  
    # }
    # 
    # 
    # # Classical MCMC   3个断点， 4个steps                                      # 
    # #
    # # Classical MCMC   1个断点                                                 #
    # EnPara_3S <- matrix(NA, ncol = N, nrow = 2*kmax+2)
    # EnPara_3S[1,] <- rep(3,N) 
    # 
    # for(i in 1:N){
    #   kk <- EnPara_3S[1,i]
    #   if(kk>0){
    #     c.x <- sample.locs(kk,L)
    #     EnPara_3S[1+(1:(kk)),i] <-   c.x}
    #   c.v <- rgamma(kk+1,shape = 2,rate = 1)
    #   EnPara_3S[kk+1+(1:(kk+1)),i] <- c.v  
    # }
    # 
    # 
    # EnPara_SM <- EnPara_2S
    # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 
    # 
    # 
    # # 06.Resampling 重采样步骤 ----------------------------------------------------------------
    # 
    # # 6.0 先定义提议函数.
    # 
    # 
    # 
    # 
    # #prpd.index <- Prop.Index(k,c=0.3,lambda=lambda_p,kmax=kmax)$index
    # #Prp.Para(para_0,prpd.index)
    # 
    # #现在先有三个矩阵，上一时刻的状态向量，这一时刻的状态向量，以及参数的矩阵
    # set.seed(100+i.all+1)
    # 
    # 
    # 
    # 
    # 
    # En.s.t_1 <- Enstat # 表示前一时刻的状态向量 x_(t-1)
    # En.s.t <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
    # En.lik <- rep(NA,N)
    # 
    # #Def for Bootstrap Resampling.
    # Mid.En.s.t_1 <- En.s.t_1
    # Mid.En.s.t   <- En.s.t
    # Mid.En.lik   <- En.lik
    # Mid.EnPara   <- EnPara
    
    # 1.0 The main loop. -------------------------------------------------------
    
    #EnPara_2S <- EnPara
    # For easier index, the dim of observations is n, not m.
    
    runtime_RJ <- system.time({
      Res1_Rj <- DA_Process(Enstat,EnPara,1)
    })
    EnPara    <- Res1_Rj$EnP
    En.s.t_1  <- Res1_Rj$EnS
    
    
    
    
    
    # # 3.0 比较之前的程序  断点数=1 -------------------------------------------------------------
    # En.s.t_1_1S <- Enstat_1S # 表示前一时刻的状态向量 x_(t-1)
    # En.s.t_1S <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
    # En.lik_1S <- rep(NA,N)
    # 
    # 
    # 
    # 
    # runtime_1S <- system.time({
    #   Res1_1S <- DA_Process(Enstat_1S,EnPara_1S,0)
    # })
    # 
    # En.s.t_1_1S <- Res1_1S$EnS
    # EnPara_1S <- Res1_1S$EnP
    # 
    # 
    # # 4.0 比较之前的程序  断点数=2 -------------------------------------------------------------
    # En.s.t_1_2S <- Enstat_2S # 表示前一时刻的状态向量 x_(t-1)
    # En.s.t_2S <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
    # En.lik_2S <- rep(NA,N)
    # 
    # 
    # 
    # 
    # 
    # runtime_2S <- system.time({
    #   Res1_2S <- DA_Process(Enstat_2S,EnPara_2S,0)
    # })
    # 
    # En.s.t_1_2S <- Res1_2S$EnS
    # EnPara_2S <- Res1_2S$EnP
    # 
    # # 5.0 断点数==3 --------------------------------------------------------------
    # En.s.t_1_3S <- Enstat_3S # 表示前一时刻的状态向量 x_(t-1)
    # En.s.t_3S <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
    # En.lik_3S <- rep(NA,N)
    # 
    # 
    # 
    # 
    # 
    # runtime_3S <- system.time({
    #   Res1_3S <- DA_Process(Enstat_3S,EnPara_3S,0)
    # })
    # 
    # En.s.t_1_3S <- Res1_3S$EnS
    # EnPara_3S <- Res1_3S$EnP
    # 
    # 
    # # 比较SMC方法~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 
    # En.s.t_1_SM <- Enstat_SM # 表示前一时刻的状态向量 x_(t-1)
    # En.s.t_SM <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
    # En.lik_SM <- rep(NA,N)
    # 
    # 
    # 
    # 
    # 
    # runtime_SM <- system.time({
    #   Res1_SM <- DA_Process_SMC(Enstat_SM,EnPara_SM,0)
    # })
    # 
    # En.s.t_1_SM <- Res1_SM$EnS
    # EnPara_SM <- Res1_SM$EnP
    
    #~~~~~~~~~~~~~~~~~~结果~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    mean_rj  <- rowMeans(En.s.t_1);#lines(mean_rj,col='green',lwd=3)
    # mean_1S  <- rowMeans(En.s.t_1_1S);#lines(mean_1S,col='yellow',lwd=3)
    # mean_2S  <- rowMeans(En.s.t_1_2S);#lines(mean_2S,col='red',lwd=3)
    # mean_3S  <- rowMeans(En.s.t_1_3S);#lines(mean_3S,col='purple',lwd=3)
    
    #计算MSE
    m.true <- out[length(times),1+(1:n)]
    Res_ALL[i.all,1] <- mean((mean_rj-m.true)**2)
    # Res_ALL[i.all,2] <- mean((mean_1S-m.true)**2)
    # Res_ALL[i.all,3] <- mean((mean_2S-m.true)**2)
    # Res_ALL[i.all,4] <- mean((mean_3S-m.true)**2) #可能是与粒子数有关？？
    
    
    # mean_SM  <- rowMeans(En.s.t_1_SM)
    # Res_ALL[i.all,5] <- mean((mean_SM-m.true)**2)
    
    
    Res_ALL[i.all,11] <- runtime_RJ[3]
    # Res_ALL[i.all,12] <- runtime_1S[3]
    # Res_ALL[i.all,13] <- runtime_2S[3]
    # Res_ALL[i.all,14] <- runtime_3S[3]
    # Res_ALL[i.all,15] <- runtime_SM[3]
    
    # 6.0 增加预测的比较 -------------------------------------------------------------
    
    #首先看逻辑
    # 预测的话先有 最后时刻的ensemble of states parameters
    # 和真实的 states, parameters
    # En.s <- En.s.t_1
    # En.p <- EnPara
    # 
    # En.s <- En.s.t_1_1S
    # En.p <- EnPara_1S
    # s.true <- out[length(times),1+(1:n)]
    # p.true <- para
    # times.p <- times[length(times)]+c(0:50)
    # out.pre <- ode.1D(s.true, times.p, model2, p.true, method = "adams", 
    #                   names = c("C"), dimens = n)[length(times.p),1+(1:n)]
    # En.predict <- matrix(NA,ncol = N, nrow <- n)
    # for(i.n in  1:N){
    #   k_0 <- En.p[1,i.n]
    #   para_0 <- En.p[1:(2*k_0+2),i.n]; print(para_0)
    #   if(k_0 == 0){
    #     p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    #   }else{
    #     p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    #   }
    #   En.predict[,i.n] <- ode.1D(En.s[,i.n], times.p, model2, p4evo0, method = "adams", 
    #                              names = c("C"), dimens = n)[length(times.p),1+(1:n)]
    # }
    # plot(out.pre, type='l')
    # for(i.n in 1:N){lines(En.predict[,i.n],col='red')}
    # plot.PI(En.predict,'blue')
    # lines(out.pre,lwd=3)
    # 
    # #
    # En.s <- En.s.t_1[,which(EnPara[1,]==2)]
    # En.p <- EnPara[,which(EnPara[1,]==2)]
    
    
    
    Res_r <- Pre_En(En.s.t_1,EnPara,out[length(times),1+(1:n)],para,times.p)
    # Res_1 <- Pre_En(En.s.t_1_1S,EnPara_1S, out[length(times),1+(1:n)],para,times.p)
    # Res_2 <- Pre_En(En.s.t_1_2S,EnPara_2S,out[length(times),1+(1:n)],para,times.p)
    # Res_3 <-  Pre_En(En.s.t_1_3S,EnPara_3S,out[length(times),1+(1:n)],para,times.p)
    
    Res_ALL[i.all,6] <- Res_r$MSPE
    # Res_ALL[i.all,8] <- Res_2$MSPE
    # Res_ALL[i.all,7] <- Res_1$MSPE
    # Res_ALL[i.all,9] <- Res_3$MSPE
    
    # Res_S <-  Pre_En(En.s.t_1_SM,EnPara_SM,out[length(times),1+(1:n)],para,times.p)
    # Res_ALL[i.all,10] <- Res_S$MSPE
    
  }
  
  # 将结果存储到列表中
  Res_List[[iN]] <- Res_ALL
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
  # 动态创建变量名，并赋值
  assign(paste0("Res_c", c_vec[iN]), Res_ALL)
  save(Res_ALL,  file = paste0("E:/Study-2024/09.RJ_2stage_resample_Oct24/NewRes_PARA_c",
                               c_vec[iN],".RData"))
}
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

Mse_Ave_c <- matrix(NA,nrow = 2,ncol =length(c_vec))
for (iN in 1:length(c_vec)){
  Mse_Ave_c[1,iN] <- mean(Res_List[[iN]][,1])
  Mse_Ave_c[2,iN] <- mean(Res_List[[iN]][,6])
}
results <- data.frame(
  par.c = c_vec,
  MSE =  Mse_Ave_c[1,],
  MSPE = Mse_Ave_c[2,]
)

# 转换成长格式
results_long <- results %>%
  pivot_longer(cols = c(MSE, MSPE), 
               names_to = "Metric", 
               values_to = "Value")

# 绘制折线图
ggplot(results_long, aes(x = par.c, y = Value, color = Metric, group = Metric)) +
  geom_line(lwd = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("MSE" = "#E69F00", "MSPE" = "#56B4E9")) +
  labs(title = "MSE and MSPE across Different c Values",
       x = "c Value",
       y = "Mean Value",
       color = "Metric") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(0.1, 0.5, by = 0.2)) +
  scale_y_continuous(limits = c(0, max(results_long$Value) * 1.1))
