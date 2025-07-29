# 2024.11.14
# The functions of RJMC-2 Step Re-sampling PF method. 

library(stats)
require(ReacTran)
require(fields)
library(nleqslv)
#画图
library(ggplot2)
library(dplyr)
library(hrbrthemes)

# 01. Pre Functions -------------------------------------------------------
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
  c.x <- parms[1+(1:(K+2))] #注意区分: 演化方程内部的是因为para前后加了0和n
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
lik.norm=function(sigma,y,x){exp(sum(dnorm(y,x,sigma,log=FALSE)))}

# 02. Prior Functions -----------------------------------------------------

#对间断点的位置选择，从[0,L]中抽取(2k+1)个点，选择中间的k个点
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

# 03. Proposal Functions --------------------------------------------------
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
} #另一种思路 互相验证结果一致性

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

#PPF0: Propose Parameter Function --- Type
Prp.Para <- function(para_0,prpd.index){
  k <- para_0[1]
  if(prpd.index ==  1){prpd.para<- PPF1(para_0)}
  if(prpd.index ==  2){prpd.para<- PPF2(para_0)}
  if(prpd.index ==  3){prpd.para<- PPF3(para_0)}
  if(prpd.index ==  4){prpd.para<- PPF4(para_0)}
  return(prpd.para)
}



# DA_Process --------------------------------------------------------------

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

# Prediction --------------------------------------------------------------
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
  plot(out.pre, type='l')
  for(i.n in 1:N){lines(En.predict[,i.n],col='lightblue')}
  plot.PI(En.predict,'blue')
  lines(out.pre,lwd=3,col='red')
  mean_pre  <- rowMeans(En.predict);
  lines(mean_pre,col='darkblue',lwd=2)
  MSPE <- mean((mean_pre-out.pre)**2)
  return(list('En'=En.predict,'MSPE'=MSPE,'Pre'=mean_pre))
}


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
