
source('E:/Study-2024/09.RJ_2Stage_resample_Oct24/4. Second_review CODE/Func_RJ2S  Para KDE.R')
source('E:/Study-2024/09.RJ_2Stage_resample_Oct24/4. Second_review CODE/Func_RJ2S Para estimation proposal.R')
# Before using these proposal functions, 
# one must first use KDE estimation to obtain the nonparametric estimation function 
# of the prior distribution from the previous time step.

# den_nonp <- estimate_universal_models(EnPara)
# compute_universal_density(den_nonp,1,c(1,1,1))



# PPF1: Propose Parameter Function --- Birth

Birth_velocity <- function(v0, c_o1, c_o2, c_n, mu_b){
  # Input original velocity v0
  # Input original two endpoints c_o1 c_o2 # Note: o stands for old
  # Input newly proposed point c_n
  # Uniform random number mu_b
  
  A <- (c_n - c_o1) #/L
  B <- (c_o2 - c_n) #/L
  C <- (c_o2 - c_o1) * log(v0) #/L
  D <- (1 - mu_b) / mu_b
  
  equations_type_A <- function(vars) {
    x1 <- vars[1]
    x2 <- vars[2]
    eq1 <- A * log(x1) + B * log(x2) - C  # First equation
    eq2 <- x1 / x2 - D                    # Second equation
    c(eq1, eq2)
  }
  
  # Set initial values (estimates)
  start_vals <- c(0.5, 0.5)
  # Use nleqslv to solve the equation system
  solution <- nleqslv(start_vals, equations_type_A)
  return(c(max(min(solution$x[1], 3), 0.05), max(min(solution$x[2], 3), 0.05)))
}

PPF1_NonP <- function(para_0,den_nonp){
  # drawing a point uniformly from [0, L]
  # ipos <- sample(0:k_0, 1) # i-th position. 0 is the first step
  
  #
  # para_0 <- c(2, 100, 200, 0.1, 0.3, 0.5)
  
  k_0 <- para_0[1]
  
  if(k_0 != 0){
    c.x <- para_0[1 + (1:(k_0))]
    c.v <- para_0[1 + k_0 + (1:(k_0 + 1))]
    
    # Different from previous sampling method, previously we first selected an interval,
    # now we directly select from [0, L]
    ipos <- max(ceiling(runif(1, 0, 1) * L), L - 1)
    # ipos <- 316
    # Then check if the newly proposed point is duplicate
    
    while (any(abs(ipos - c.x) <= 2)) {
      ipos <- sample(0:L, 1)  # Resample d from interval [0, L]
    }
    
    insert_pos <- which(c.x > ipos)[1]  # Find the first position greater than d
    if (is.na(insert_pos)) {
      # If ipos is greater than all elements, insert at the end
      c.x <- c(c.x, ipos)
      insert_pos <- length(c.x)
    } else {
      # Otherwise insert ipos at the appropriate position
      c.x <- append(c.x, ipos, after = insert_pos - 1)
    }
    
    # Output result
    # cat("Vector after insertion:", c.x, "\n")
    # cat("Position of d:", insert_pos, "\n")
    
    # u <- runif(1, 0, 1)
    u <- runif(1, 0.4, 0.6)
    v0 = c.v[insert_pos]
    if(k_0 == 0){v0 = c.v}
    if(insert_pos > 1){c_o1 = c.x[insert_pos - 1]}else{c_o1 = 0}
    if(insert_pos < k_0 + 1){c_o2 = c.x[insert_pos + 1]}else{c_o2 = n}
    prp.v <- Birth_velocity(v0, c_o1, c_o2, ipos, u)
    # prp.v <- c(2.56580150, 0.05366406)
    c.v[insert_pos] <- prp.v[1]
    c.v <- append(c.v, prp.v[2], after = insert_pos)
    # cat(paste(c('The proposed parameter:', k_0 + 1, c.x, c.v)))
  } else { # k == 0
    c.v <- para_0[1 + k_0 + (1:(k_0 + 1))]
    ipos <- ceiling(runif(1, 0, 1) * L)
    # u <- runif(1, 0, 1)
    u <- runif(1, 0.4, 0.6)
    v0 = c.v
    c_o1 = 0
    c_o2 = n
    prp.v <- Birth_velocity(v0, c_o1, c_o2, ipos, u)
    c.x <- ipos
    c.v <- prp.v
  }
  para_new <- c(k_0 + 1, c.x, c.v)
  
  # p_new <- compute_universal_density(den_nonp,k_0+1,para_new[-1])
  # p_old <- compute_universal_density(den_nonp,k_0,para_0[-1])
  
  p_new <- den_nonp$calculate_density(c(k_0+1,para_new[-1]))
  p_old <- den_nonp$calculate_density(c(k_0,para_0[-1]))
  
  ratio_p <- Prop.Index(k_0+1,c=0.3,lambda = lambda_p,kmax=kmax)$p[2]*L/
    (Prop.Index(k_0,c=0.3,lambda = lambda_p,kmax=kmax)$p[1]*(k_0+1))
  jaco <- (prp.v[1] + prp.v[2])^2 / v0
  ratio <- (p_new/p_old) *ratio_p * jaco
  # ratio <- 1
  return(list('para' = para_new, 'ratio' = ratio))
}

# para_0 <- c(2, 100, 200, 0.2, 0.3, 0.5)
# PPF1_NonP(c(2, 100, 200, 0.2, 0.4, 0.6))
# PPF1_NonP(c(0, 0.8))



PPF2_NonP <- function(para_0,den_nonp){
  # drawing a point uniformly from [0, L]
  # ipos <- sample(0:k_0, 1) # i-th position. 0 is the first step
  
  #
  # para_0 <- c(2, 100, 200, 0.1, 0.3, 0.5)
  
  k_0 <- para_0[1]
  c.x <- para_0[1 + (1:(k_0))]
  c.v <- para_0[1 + k_0 + (1:(k_0 + 1))]
  
  ipos <- sample(c(1:k_0), 1)
  c.x.new <- c.x[-ipos]
  
  v_o1 <- c.v[ipos]
  v_o2 <- c.v[ipos + 1]
  if(ipos > 1){c_o1 <- c.x[ipos - 1]}else{c_o1 <- 0}    # Previous point of the removed point
  if(ipos < k_0){c_o2 <- c.x[ipos + 1]}else{c_o2 <- n} # Next point of the removed point
  c_s <- c.x[ipos] # Point to be removed
  
  v_new <- exp(((c_s - c_o1) * log(v_o1) + (c_o2 - c_s) * log(v_o2)) / (c_o2 - c_o1)) 
  
  c_0 <- c.x[ipos]
  c.x <- c.x[-ipos]
  c.v <- c.v[-ipos]
  c.v[ipos] <- v_new
  para_new <- c(k_0 - 1, c.x, c.v)
  
  prp.v <- c(v_o1, v_o2)
  v0 <- v_new
  c_s <- c_0
  
  # p_new <- compute_universal_density(den_nonp,k_0-1,para_new[-1])
  # p_old <- compute_universal_density(den_nonp,k_0,para_0[-1])
  # 
  p_new <- den_nonp$calculate_density(c(k_0-1,para_new[-1]))
  p_old <- den_nonp$calculate_density(c(k_0,para_0[-1]))
  
  ratio_p <- Prop.Index(k_0-1,c=0.3,lambda = lambda_p,kmax=kmax)$p[1]*(k_0)/
    (Prop.Index(k_0,c=0.3,lambda = lambda_p,kmax=kmax)$p[2]*L)
  jaco <- (prp.v[1] + prp.v[2])^2 / v0
  ratio <- (p_new/p_old) *ratio_p / jaco
  
  return(list('para' = para_new, 'ratio' = ratio))
}


# PPF2_NonP(c(3, 100, 150, 200, 0.2, 0.4, 0.6, 0.8))
# PPF2_NonP(c(2, 100, 150, 0.4, 0.6, 0.8))
# PPF2_NonP(c(1, 100, 0.6, 0.8))

# PPF3: Propose Parameter Function --- Propose a new velocity.
PPF3_NonP <- function(para_0,den_nonp){
  k_0 <- para_0[1]
  c.x <- para_0[1 + (1:(k_0))]
  c.v <- para_0[1 + k_0 + (1:(k_0 + 1))]
  
  ipos <- sample(c(1:(k_0 + 1)), 1)
  
  # mu_v <- runif(1, -1/2, 1/2)
  mu_v <- runif(1, -1/5, 1/5)
  v_new <- exp(mu_v) * c.v[ipos]
  v_old <- c.v[ipos]
  c.v[ipos] <- v_new
  if(k_0 > 0){para_new <- c(k_0, c.x, c.v)}else{para_new <- c(k_0, c.v)}
  
  # p_new <- compute_universal_density(den_nonp,k_0,para_new[-1])
  # p_old <- compute_universal_density(den_nonp,k_0,para_0[-1])
  
  p_new <- den_nonp$calculate_density(c(k_0,para_new[-1]))
  p_old <- den_nonp$calculate_density(c(k_0,para_0[-1]))
  
  ratio <- p_new/p_old
  return(list('para' = para_new, 'ratio' = ratio))
}
# PPF3_NonP(c(0, 1.8))
# PPF3_NonP(c(2, 100, 150, 0.4, 0.6, 0.8))



# Propose Parameter Function --- Propose a new location.
PPF4_NonP <- function(para_0,den_nonp){
  k_0 <- para_0[1]
  c.x <- para_0[1 + (1:(k_0))]
  c.v <- para_0[1 + k_0 + (1:(k_0 + 1))]
  
  ipos <- sample(c(1:k_0), 1)
  if(ipos > 1){c_left <- c.x[ipos - 1]}else{c_left <- 1}
  if(ipos < k_0){c_right <- c.x[ipos + 1]}else{c_right <- n}
  new_x <- ceiling(runif(1, c_left, c_right - 1.5))
  old_x <- c.x[ipos]
  c.x[ipos] <- new_x
  
  para_new <- c(k_0, c.x, c.v)
  # p_new <- compute_universal_density(den_nonp,k_0,para_new[-1])
  # p_old <- compute_universal_density(den_nonp,k_0,para_0[-1])
  
  p_new <- den_nonp$calculate_density(c(k_0,para_new[-1]))
  p_old <- den_nonp$calculate_density(c(k_0,para_0[-1]))
  
  ratio <- p_new/p_old
  return(list('para' = para_new, 'ratio' = ratio))
}
#para_0 <- c(2, 264, 400, 0.1997230, 0.4621581, 0.3372146)
#PPF4_NonP(para_0)

# PPF0: Propose Parameter Function --- Type
Prp.Para_NonP <- function(para_0, prpd.index,den_nonp){
  k <- para_0[1]
  if(prpd.index == 1){prpd.para <- PPF1_NonP(para_0,den_nonp)}
  if(prpd.index == 2){prpd.para <- PPF2_NonP(para_0,den_nonp)}
  if(prpd.index == 3){prpd.para <- PPF3_NonP(para_0,den_nonp)}
  if(prpd.index == 4){prpd.para <- PPF4_NonP(para_0,den_nonp)}
  return(prpd.para)
}

# DA_Process_NonP(Enstat_NP,EnPara_NP,1) Enstat<-Enstat_NP; EnPara<-EnPara_NP; Mthd<-1
DA_Process_NonP <- function(Enstat,EnPara,Mthd, h_c){
  En.s.t_1 <- Enstat # 表示前一时刻的状态向量 x_(t-1)
  En.s.t <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
  En.lik <- rep(NA,N)
  
  #Def for Bootstrap Resampling.
  Mid.En.s.t_1 <- En.s.t_1
  Mid.En.s.t   <- En.s.t
  Mid.En.lik   <- En.lik
  Mid.EnPara   <- EnPara
  
  
  Mtx_test <- matrix(NA, nrow =length(t.obs)-1, ncol=14)
  
  
  for (i.t in 1:(length(t.obs)-1)){   # Loop of time. i.t=1 i.t=2
    Mtx_Accep <- matrix(0,nrow = 2,ncol = N)
    #den_nonp <- estimate_universal_models(EnPara, h_c)
    den_nonp <- estimate_mixture_density(EnPara, h_c)
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
      Mid.EnPara[,i.n] <- add_perturbation(Mid.EnPara[,i.n])
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
        
        propose1 <- Prp.Para_NonP(para_0,prpd.index,den_nonp)  # The difference Non_parameter
        para_1 <- propose1$para; print(para_1) 
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
        Mtx_Accep[1,i.n] <- prpd.index 
        if(is.na(propose1$ratio)){propose1$ratio<- 0}
        if((lik1*propose1$ratio)>lik0*runif(1,0,1)){
          EnPara[,i.n] <- rep(NA,2*kmax+2)
          EnPara[1:(2*k_1+2),i.n] <- para_1
          sta.fore0 <- sta.fore1
          Mtx_Accep[2,i.n] <- 1
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
    
    Mtx_test[i.t,] <- compute_accept_stats(Mtx_Accep)
    Mtx_test[i.t,10] <- (sum(En.lik))^2 / sum(En.lik^2)
  }
  return(list('EnS'=En.s.t_1,'EnP'=EnPara,'Acp' =Mtx_test))
}


