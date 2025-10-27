# Code for 2 stage resampling method with RJMCMC
# Version 2025 Oct
# 
library(stats)
require(ReacTran)
require(fields)
library(nleqslv)

source('E:/Study-2024/09.RJ_2stage_resample_Oct24/4. Second_review CODE/Func_RJ2S  Para KDE.R')

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
N <- 40
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


lambda_estimate <- mean(EnPara[1,])


# 从EnPara矩阵中提取各模型对应的数据
extract_model_data <- function(EnPara) {
  # 获取模型索引（第一行）
  model_indices <- as.integer(EnPara[1, ])
  
  # 初始化存储各模型数据的列表
  model_data <- list(
    model_0 = numeric(0),
    model_1 = numeric(0),
    model_2 = numeric(0),
    model_3 = numeric(0)
  )
  
  # 遍历每一列，根据模型索引提取对应的行数据
  for (i in 1:ncol(EnPara)) {
    model_index <- model_indices[i]
    
    if (model_index == 0) {
      # 模型0：取第二行
      if (!is.na(EnPara[2, i])) {
        model_data$model_0 <- c(model_data$model_0, EnPara[2, i])
      }
    } 
    else if (model_index == 1) {
      # 模型1：取第3-4行
      rows_1 <- 3:4
      valid_values <- na.omit(as.numeric(EnPara[rows_1, i]))
      model_data$model_1 <- c(model_data$model_1, valid_values)
    }
    else if (model_index == 2) {
      # 模型2：取第4-6行
      rows_2 <- 4:6
      valid_values <- na.omit(as.numeric(EnPara[rows_2, i]))
      model_data$model_2 <- c(model_data$model_2, valid_values)
    }
    else if (model_index == 3) {
      # 模型3：取第5-8行
      rows_3 <- 5:8
      valid_values <- na.omit(as.numeric(EnPara[rows_3, i]))
      model_data$model_3 <- c(model_data$model_3, valid_values)
    }
  }
  
  return(model_data)
}

# 显示各模型数据的统计信息
summarize_model_data <- function(model_data) {
  cat("Model Data Summary\n")
  cat("==================\n")
  
  for (model_name in names(model_data)) {
    samples <- model_data[[model_name]]
    cat(paste("\n", model_name, ":\n"))
    cat("  Sample size:", length(samples), "\n")
    
    if (length(samples) > 0) {
      cat("  Range: [", round(min(samples), 4), ", ", round(max(samples), 4), "]\n", sep = "")
      cat("  Mean:", round(mean(samples), 4), "\n")
      cat("  Standard deviation:", round(sd(samples), 4), "\n")
      cat("  First 10 values:", paste(round(head(samples, 10), 4), collapse = ", "), "\n")
    } else {
      cat("  No data available\n")
    }
  }
}

# 提取数据
model_data <- extract_model_data(EnPara)

# 显示汇总信息
summarize_model_data(model_data)

# 访问特定模型的数据
cat("\n=== Accessing Individual Model Data ===\n")
cat("Model 0 data (first 10 values):\n")
print(head(model_data$model_0, 10))

cat("\nModel 1 data (first 10 values):\n")
print(head(model_data$model_1, 10))

cat("\nModel 2 data (first 10 values):\n")
print(head(model_data$model_2, 10))

cat("\nModel 3 data (first 10 values):\n")
print(head(model_data$model_3, 10))

# 可视化各模型数据的分布
plot_model_distributions <- function(model_data) {
  par(mfrow = c(2, 2))
  
  for (i in 1:4) {
    model_name <- paste0("model_", i-1)
    samples <- model_data[[model_name]]
    
    if (length(samples) > 0) {
      hist(samples, breaks = 20, main = paste("Distribution of", model_name),
           xlab = "Values", col = "lightblue", border = "white")
      box()
    } else {
      plot(0, 0, type = "n", xlab = "", ylab = "", main = paste(model_name, "- No Data"))
      text(0, 0, "No data available")
    }
  }
  par(mfrow = c(1, 1))
}

# 绘制分布图
plot_model_distributions(model_data)

# 生成详细的数据报告
generate_data_report <- function(model_data) {
  cat("\n=== Detailed Data Report ===\n")
  
  report_df <- data.frame(
    Model = character(),
    Sample_Size = numeric(),
    Min = numeric(),
    Max = numeric(),
    Mean = numeric(),
    Median = numeric(),
    SD = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (model_name in names(model_data)) {
    samples <- model_data[[model_name]]
    
    if (length(samples) > 0) {
      new_row <- data.frame(
        Model = model_name,
        Sample_Size = length(samples),
        Min = round(min(samples), 4),
        Max = round(max(samples), 4),
        Mean = round(mean(samples), 4),
        Median = round(median(samples), 4),
        SD = round(sd(samples), 4)
      )
      report_df <- rbind(report_df, new_row)
    } else {
      new_row <- data.frame(
        Model = model_name,
        Sample_Size = 0,
        Min = NA,
        Max = NA,
        Mean = NA,
        Median = NA,
        SD = NA
      )
      report_df <- rbind(report_df, new_row)
    }
  }
  
  print(report_df)
  return(report_df)
}

# 生成报告
data_report <- generate_data_report(model_data)

# 保存提取的数据到变量（便于后续分析）
model_0_samples <- model_data$model_0
model_1_samples <- model_data$model_1
model_2_samples <- model_data$model_2
model_3_samples <- model_data$model_3

all_model_samples <- c(model_0_samples, model_1_samples, model_2_samples, model_3_samples)



# 估计gamma分布参数的完整函数
estimate_gamma_parameters <- function(samples) {
  # 移除NA值和确保为正数（gamma分布要求x>0）
  samples_clean <- samples[!is.na(samples) & samples > 0]
  
  if (length(samples_clean) == 0) {
    stop("No valid positive samples available for gamma distribution estimation")
  }
  
  cat("Gamma Distribution Parameter Estimation\n")
  cat("=======================================\n")
  cat("Sample size:", length(samples_clean), "\n")
  cat("Sample range: [", min(samples_clean), ", ", max(samples_clean), "]\n", sep = "")
  cat("Sample mean:", round(mean(samples_clean), 4), "\n")
  cat("Sample variance:", round(var(samples_clean), 4), "\n\n")
  
  # 方法1: 矩估计法
  # Gamma(shape = α, rate = β) 
  # 均值 μ = α/β, 方差 σ² = α/β²
  # 所以: α = μ²/σ², β = μ/σ²
  mean_val <- mean(samples_clean)
  var_val <- var(samples_clean)
  
  shape_moment <- mean_val^2 / var_val
  rate_moment <- mean_val / var_val
  scale_moment <- 1 / rate_moment  # scale参数 = 1/rate
  
  # 方法2: 最大似然估计 (使用MASS包)
  if (require(MASS, quietly = TRUE)) {
    fit_mle <- fitdistr(samples_clean, "gamma")
    shape_mle <- fit_mle$estimate["shape"]
    rate_mle <- fit_mle$estimate["rate"]
    scale_mle <- 1 / rate_mle
    mle_se <- fit_mle$sd
  } else {
    shape_mle <- NA
    rate_mle <- NA
    scale_mle <- NA
    mle_se <- c(shape = NA, rate = NA)
    warning("MASS package not available, MLE estimation skipped")
  }
  
  # 方法3: 使用对数矩估计（对偏态分布更稳健）
  log_mean <- mean(log(samples_clean))
  log_var <- var(log(samples_clean))
  
  # 对数矩估计公式
  shape_log <- (1 + sqrt(1 + 4 * log_var / 3)) / (4 * log_var)
  rate_log <- shape_log / exp(log_mean + digamma(shape_log) - log(shape_log))
  
  # 显示结果
  cat("Parameter Estimates:\n")
  cat("-------------------\n")
  cat("Moment Estimation:\n")
  cat("  Shape (α):", round(shape_moment, 4), "\n")
  cat("  Rate (β):", round(rate_moment, 4), "\n")
  cat("  Scale (θ = 1/β):", round(scale_moment, 4), "\n\n")
  
  if (!is.na(shape_mle)) {
    cat("Maximum Likelihood Estimation:\n")
    cat("  Shape (α):", round(shape_mle, 4), "±", round(mle_se["shape"], 4), "\n")
    cat("  Rate (β):", round(rate_mle, 4), "±", round(mle_se["rate"], 4), "\n")
    cat("  Scale (θ = 1/β):", round(scale_mle, 4), "\n\n")
  }
  
  cat("Log-Moment Estimation:\n")
  cat("  Shape (α):", round(shape_log, 4), "\n")
  cat("  Rate (β):", round(rate_log, 4), "\n")
  cat("  Scale (θ = 1/β):", round(1/rate_log, 4), "\n\n")
  
  # 拟合优度检验
  perform_goodness_of_fit <- function(samples, shape, rate) {
    # Kolmogorov-Smirnov检验
    ks_test <- ks.test(samples, "pgamma", shape = shape, rate = rate)
    
    return(list(
      ks_statistic = ks_test$statistic,
      ks_pvalue = ks_test$p.value,
      fit_quality = ifelse(ks_test$p.value > 0.05, "Good fit", "Poor fit")
    ))
  }
  
  # 使用矩估计结果进行拟合优度检验
  gof <- perform_goodness_of_fit(samples_clean, shape_moment, rate_moment)
  
  cat("Goodness of Fit Test:\n")
  cat("---------------------\n")
  cat("Kolmogorov-Smirnov test:\n")
  cat("  Statistic:", round(gof$ks_statistic, 4), "\n")
  cat("  p-value:", round(gof$ks_pvalue, 4), "\n")
  cat("  Conclusion:", gof$fit_quality, "\n")
  
  return(list(
    samples = samples_clean,
    sample_size = length(samples_clean),
    sample_statistics = list(
      mean = mean_val,
      variance = var_val,
      min = min(samples_clean),
      max = max(samples_clean)
    ),
    moment_estimation = list(
      shape = shape_moment,
      rate = rate_moment,
      scale = scale_moment
    ),
    mle_estimation = list(
      shape = shape_mle,
      rate = rate_mle,
      scale = scale_mle,
      standard_errors = mle_se
    ),
    log_moment_estimation = list(
      shape = shape_log,
      rate = rate_log,
      scale = 1/rate_log
    ),
    goodness_of_fit = gof
  ))
}

# 可视化gamma分布拟合
plot_gamma_fit <- function(samples, shape_est, rate_est, title = "Gamma Distribution Fit") {
  par(mfrow = c(2, 2))
  
  # 1. 直方图与拟合密度曲线
  hist(samples, probability = TRUE, breaks = 30,
       main = paste(title, "\nHistogram with Fitted Gamma"),
       xlab = "Value", col = "lightblue", border = "white")
  
  x_range <- seq(0, max(samples) * 1.2, length.out = 100)
  lines(x_range, dgamma(x_range, shape = shape_est, rate = rate_est), 
        col = "red", lwd = 2)
  
  legend("topright", legend = c("Data", "Fitted Gamma"), 
         fill = c("lightblue", "red"))
  
  # 2. Q-Q图
  theoretical_quantiles <- qgamma(ppoints(length(samples)), 
                                  shape = shape_est, rate = rate_est)
  plot(theoretical_quantiles, sort(samples),
       xlab = "Theoretical Gamma Quantiles", 
       ylab = "Sample Quantiles",
       main = paste(title, "\nQ-Q Plot"))
  abline(0, 1, col = "red")
  
  # 3. 累积分布函数比较
  ecdf_empirical <- ecdf(samples)
  plot(ecdf_empirical, main = paste(title, "\nCDF Comparison"),
       xlab = "Value", ylab = "Cumulative Probability")
  lines(x_range, pgamma(x_range, shape = shape_est, rate = rate_est), 
        col = "red", lwd = 2)
  legend("bottomright", legend = c("Empirical", "Theoretical"), 
         col = c("black", "red"), lwd = 2)
  
  # 4. 概率-概率图 (P-P图)
  empirical_probs <- ecdf_empirical(samples)
  theoretical_probs <- pgamma(samples, shape = shape_est, rate = rate_est)
  plot(theoretical_probs, empirical_probs,
       xlab = "Theoretical Probabilities", 
       ylab = "Empirical Probabilities",
       main = paste(title, "\nP-P Plot"),
       xlim = c(0, 1), ylim = c(0, 1))
  abline(0, 1, col = "red")
  
  par(mfrow = c(1, 1))
}

# 主分析函数
analyze_gamma_distribution <- function(samples) {
  # 参数估计
  results <- estimate_gamma_parameters(samples)
  
  # 使用矩估计结果进行可视化（因为MLE可能不可用）
  plot_gamma_fit(results$samples, 
                 results$moment_estimation$shape,
                 results$moment_estimation$rate,
                 "Gamma Distribution Analysis")
  
  return(results)
}

# 运行分析
gamma_results <- analyze_gamma_distribution(all_model_samples)

# 生成预测区间
generate_prediction_intervals <- function(shape, rate, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  intervals <- qgamma(probs, shape = shape, rate = rate)
  names(intervals) <- paste0(probs * 100, "%")
  return(intervals)
}

# 生成预测区间（使用矩估计参数）
prediction_intervals <- generate_prediction_intervals(
  gamma_results$moment_estimation$shape,
  gamma_results$moment_estimation$rate
)

cat("\nPrediction Intervals (using moment estimation):\n")
print(round(prediction_intervals, 4))

# 简洁版本：如果只需要基本的矩估计
estimate_gamma_simple <- function(samples) {
  samples_clean <- samples[!is.na(samples) & samples > 0]
  mean_val <- mean(samples_clean)
  var_val <- var(samples_clean)
  
  shape <- mean_val^2 / var_val
  rate <- mean_val / var_val
  scale <- 1 / rate
  
  return(list(shape = shape, rate = rate, scale = scale))
}

# 快速估计
simple_estimate <- estimate_gamma_simple(all_model_samples)

Est_para <- function(EnPara){
  
  lambda_estimate <- mean(EnPara[1,])
  
  model_data <- extract_model_data(EnPara)
  model_0_samples <- model_data$model_0
  model_1_samples <- model_data$model_1
  model_2_samples <- model_data$model_2
  model_3_samples <- model_data$model_3
  
  all_model_samples <- c(model_0_samples, model_1_samples, model_2_samples, model_3_samples)
  simple_estimate <- estimate_gamma_simple(all_model_samples)
  return(c(lambda_estimate,simple_estimate$shape,simple_estimate$rate))
}
#Est_para(EnPara)
#Enstat<-Enstat_2S;EnPara<-EnPara_2S;Mthd<-1
DA_Process_EstP <- function(Enstat,EnPara,Mthd){
  En.s.t_1 <- Enstat # 表示前一时刻的状态向量 x_(t-1)
  En.s.t <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
  En.lik <- rep(NA,N)
  
  #Def for Bootstrap Resampling.
  Mid.En.s.t_1 <- En.s.t_1
  Mid.En.s.t   <- En.s.t
  Mid.En.lik   <- En.lik
  Mid.EnPara   <- EnPara
  
  #Add the test 
  
  Mtx_test <- matrix(NA, nrow =length(t.obs)-1, ncol=14)
  
  for (i.t in 1:(length(t.obs)-1)){   # Loop of time. i.t=1 i.t=2
    Mtx_Accep <- matrix(0,nrow = 2,ncol = N)
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
        

        
        
        
        #原来未更新的部分
        lik0 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore0[Obs.site])
        lik1 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore1[Obs.site])
        Mtx_Accep[1,i.n] <- prpd.index 
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
    est_para <- Est_para(EnPara) 
    lambda_p <- est_para[1] # prior mean

    #Gamma 
    g_alpha <- est_para[2]
    g_beta  <- est_para[3]
    
    
    Mtx_test[i.t,] <- compute_accept_stats(Mtx_Accep)
    Mtx_test[i.t,10] <- (sum(En.lik))^2 / sum(En.lik^2)
  }
  return(list('EnS'=En.s.t_1,'EnP'=EnPara,'Acp' =Mtx_test))
}
DA_Process <- function(Enstat,EnPara,Mthd){
  En.s.t_1 <- Enstat # 表示前一时刻的状态向量 x_(t-1)
  En.s.t <- matrix(NA,nrow = n, ncol = N)   # 表示新时刻的状态向量 x_t
  En.lik <- rep(NA,N)
  
  #Def for Bootstrap Resampling.
  Mid.En.s.t_1 <- En.s.t_1
  Mid.En.s.t   <- En.s.t
  Mid.En.lik   <- En.lik
  Mid.EnPara   <- EnPara
  
  
  #Add the test 
  
  Mtx_test <- matrix(NA, nrow =length(t.obs)-1, ncol=14)
  
  
  for (i.t in 1:(length(t.obs)-1)){   # Loop of time. i.t=1 i.t=2
    Mtx_Accep <- matrix(0,nrow = 2,ncol = N)
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
        Mtx_Accep[1,i.n] <- prpd.index
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