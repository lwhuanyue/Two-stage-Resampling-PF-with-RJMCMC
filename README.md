##  Particle Filter Estimation with Two-step RJMCMC Resampling for Sequential Hierarchical Bayesian Model 

### Overview
This repository provides the implementation code for our novel  with **Two-Stage Resampling Particle Filter** for estimating **Sequential Hierarchical Bayesian Model (SHBM)**, addressing key challenges in Data Assimilation (DA) involving uncertain dynamic model structures and varying parameter dimensions.

### Key Features
- **Sequential Hierarchical Bayesian Model (SHBM)**: A flexible framework integrating state space modeling with hierarchical parameter structures to handle model uncertainty
- **Two-Stage Resampling Scheme**:
  - **Stage 1**: Bootstrap particle filter for efficient preliminary resampling
  - **Stage 2**: RJMCMC-based resampling for particle diversification in trans-dimensional spaces
- **Simultaneous Capabilities**:
  - Automatic model structure identification
  - Parameter estimation in variable-dimensional spaces
  - State assimilation and prediction

### Problem Solved
Traditional DA methods often suffer from:
- **Filter divergence** due to incorrect model assumptions
- **Particle impoverishment** in high-dimensional parameter spaces
- **Inability to handle** varying parameter dimensions across different model structures

Our method effectively overcomes these limitations through the two-stage resampling approach, maintaining particle diversity while preventing filter divergence.

### Applications
- Advection equation models
- Lorenz 96 systems
- Spatiotemporal process studies
- Systems with structural uncertainties and parameter dimension variations

### Reference
For theoretical foundations, methodological details, and experimental results, please refer to our accompanying publication.

## Installation

You can install the development version of SIStree from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lwhuanyue/Two-stage-Resampling-PF-with-RJMCMC")
```




## 📚 File Description

### 🎯 Core Algorithm Files
- **Illustration_RJ2S.R** - A minimal, reproducible, end-to-end workflow that demonstrates the proposed method, explains how to use the codebase, and generates all relevant results.
- **RJ2S_basic_function_and_setting.R** - Provides basic functions for PF with two-stage RJMCMC resampling method, including proposal functions for four transition types and the main DA process function `DA_Process(Enstat, EnPara, Mthd)`
- **R2Resample_Lorenz96_Main.R** - The main/master script for performing data assimilation with the Lorenz96 model.
  
### 📊 Analysis Scripts
- **Loop_R2S_prior_compare2_alpha_and_c.R** - Sensitivity analysis for parameter alpha, comparing MSEs and MSPEs of the proposed method under different alpha values through independent repeated simulations
- **Loop_RU2S_prior_compare_3_c.R** - Sensitivity analysis for parameter c, comparing MSEs and MSPEs of the proposed method under different c values through independent repeated simulations
- **Loop_R2S_Non_para_method_compare.R** - Comparison of three parameter distribution estimation methods, evaluating Acceptance rate, ESS, MSEs, and MSPEs through independent repeated simulations
- **Loop_R2S_Non_para_bandwidth_compare.R** - Bandwidth analysis, comparing Acceptance rate, ESS, MSEs, and MSPEs for different bandwidths in parameter estimation method through independent repeated simulations

### 📈 Plotting Scripts
- **Plot_RJ2s_prior_compare_3_c.R** - Generates plots for parameter c sensitivity analysis based on results from "Loop RU2S prior compare 3_c.R", corresponding to the right subplot of Figure 2 in the associated paper
- **Plot_R2s_prior_compare_2_alpha.R** - Generates plots for parameter alpha sensitivity analysis based on results from "Loop R2S prior_compare2 alpha and c.R", corresponding to the left subplot of Figure 2 in the associated paper

### 🔧 Utility Functions
- **Func_R2S_test_acceptance_and_ess.R** - Contains functions to calculate monitoring metrics such as acceptance rate and ESS for each of the four transition types at each time step
- **Func_R2S_Para_KDE.R** - Kernel density estimation functions used for non-parametric estimation
- **Func_R2S_add_perturbation.R** - Contains functions to add perturbations to duplicated particles after the first resampling step
- **Func_RU2S_Para_estimation_proposal.R** - Parameter distribution estimation methods and corresponding DA main function `DA_Process_EstP`
- **Func_RU2S_Non_para_estimation_proposal.R** - Non-parametric distribution function estimation methods and corresponding DA main function `DA_Process_NonP`
- **Resample_Func_L96.R** - Core functions and utilities for the Lorenz96 model data assimilation.

## Usage

We show how to use the code with examples.

## 🚀 Quick Start

### 1. Load the necessary scripts
```r
R_path <- ''

source(paste0(R_path, 'RJ2S_basic function and setting.R'))
source(paste0(R_path, 'Func_RJ2S  add_perturbation.R'))
source(paste0(R_path, 'Func_RJ2S  Para KDE.R'))
source(paste0(R_path, 'Func_RJ2S  test acceptance and ess.R'))
source(paste0(R_path, 'Func_RJ2S Para estimation proposal.R'))
source(paste0(R_path, 'Func_RJ2S Non_para estimation proposal.R'))

```

### 2. Run the DA process

Generating Spatiotemporal Field Data - 1D Advection Equation
```
out <- ode.1D(Cini, times, model2, para, method = "adams", names = c("C"), dimens = n) #The true states over time.
image(out,xlab = 'times',ylab = 'positions')   # The figure of true states over time.
```

Generate the initial ensemble of state vectors and parameters. For the definition of the parameter space, please refer to Sections 2 and 4 of the associated paper.
```
set.seed(1028) # Fix a random seed
N <- 80 # The size of ensemble
# Generate the initial ensemble of particles of states and parameter --------
Enstat <- matrix(NA, nrow = n, ncol = N)                   # The initial ensemble of particles of states
for(i in 1:N){Enstat[,i] <- Cini*(1+rep(rnorm(1,0,0.1),n))}          

EnPara <- matrix(NA, ncol = N, nrow = 2*kmax+2)            # The initial ensemble of particles of parameters
EnPara[1,] <- sample(c(0:kmax),N,replace=TRUE,dpois(c(0:kmax), lambda = lambda_p)) # The index of model.
for(i in 1:N){
  kk <- EnPara[1,i]
  if(kk>0){
    c.x <- sample.locs(kk,L)
    EnPara[1+(1:(kk)),i] <-   c.x}
  c.v <- rgamma(kk+1,shape=g_alpha,rate=g_beta)
  EnPara[kk+1+(1:(kk+1)),i] <- c.v  
}

```

The ensemble of state vectors and parameters is updated sequentially with the acquisition of observations via the proposed method. 
The observations are postulated to be the authentic spatiotemporal field data corrupted by additive Gaussian noise. 
This DA process is encapsulated within the function `DA_process` in our codebase.

```
En.s.t_1 <- Enstat # Represents state vector from previous time step x_(t-1)
En.s.t <- matrix(NA,nrow = n, ncol = N)   # Represents state vector at new time step x_t
En.lik <- rep(NA,N)

#Def for Bootstrap Resampling.
Mid.En.s.t_1 <- En.s.t_1
Mid.En.s.t   <- En.s.t
Mid.En.lik   <- En.lik
Mid.EnPara   <- EnPara
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
    #Later change *9 to *18
    sig.eps <- (10-(i.t-1)/length(t.obs)*9)*sig.eps0}else{ 
      sig.eps <- sig.eps0}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Loop_1:       Propagate          #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  
  
  # In this step, propagate the state vector (En.s.t_1) to the new time (En.s.t)
  for (i.n in 1:N){
    k_0 <- EnPara[1,i.n]
    para_0 <- EnPara[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    }else{
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    # This is the original: p4evo0 <- c(k_0,0,para_0[2*(1:k_0)+1], n, para_0[2*(0:k_0)+2],0.0)
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
  # Resampling replacement. Note: must use intermediate variables!!!!
  # Only consider half
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
  #After replacement is complete, assign the matrices back
  En.s.t_1 <- Mid.En.s.t_1
  En.s.t   <- Mid.En.s.t
  En.lik   <- Mid.En.lik
  EnPara   <- Mid.EnPara
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #Loop_3:  RJMCMC Re-sampling      #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  for (i.n in 1:N){   # Loop of ensemble member. i.n = 2
    #Due to the first resampling step, repeat the input parameter step from the first step
    for (gg in 1:G) {
      k_0 <- EnPara[1,i.n]
      para_0 <- EnPara[1:(2*k_0+2),i.n]; print(para_0)
      if(k_0 == 0){
        p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
      }else{
        p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
      }
      #Next, propose new parameters
      prpd.index <- Prop.Index(k_0,c=0.3,lambda=lambda_p,kmax=kmax)$index
      propose1 <- Prp.Para(para_0,prpd.index) 
      para_1 <- propose1$para; print(para_1) # Now need to define the function, propose a new function
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
      # Compare En.lik[i.n] and lik0
      
      #Note two points: one is the Jacobian, the other is not to mistake the initial vector
      
      #Then add if
      
      #~~~~~~~~~~~~~~~~End loop # ~~~~~~~~~~
      
      #Original unupdated part
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
  # Count the number of elements in the interval [80,120]
  count_75_125 <- sum(EnPara >= 75 & EnPara <= 125, na.rm = TRUE)
  
  # Count the number of elements in the interval [230,270]
  count_225_275 <- sum(EnPara >= 225 & EnPara <= 275, na.rm = TRUE)
  Loca.m[,i.t] <- c(min(count_75_125,N),min(count_225_275,N))
  
}

tt <- 1+t.obs[i.t+1]
m.true  <- out[tt,-1]
mean.rj <- rowMeans(En.s.t_1)
```

Plot the result figures that are used in the associated article. The code for this can be found in `Illustration_RJ2S .R`.

```
par(mfrow=c(1,1),mar=c(2,2,0,0.5)+2.5)
plot(m.true,type ='l', ylim = c(min(En.s.t_1,m.true),7),
     cex.axis = 2, cex.lab =2, cex.main =2,
     main='The pointwise posterior 80% credible intervals',xlab = 'Position',ylab = 'Value')
for (i.n in 1:N){lines(En.s.t_1S[,i.n],col = 'lightcoral',lwd=2)}
plot.PI(En.s.t_1S,'blue')
# for (i.n in 1:N){lines(En.s.t_1[,i.n],col = 'lightcoral',lwd=2)}
# plot.PI(En.s.t_1,'blue')
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

# Plotting

# Load libraries
kmax <- 3
N.en <- c()        # Represents how many particles correspond to each model
for (i in 1:kmax){N.en <- c(N.en,sum(EnPara[1,]==i))}
d.m.index <- c()     # Represents the model index
d.para.order <- c()  # Represents which breakpoint order
d.location <- c()    # The location corresponding to that particle
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
     xlab="Locations", ylim=c(40,0), las=1 , col=rgb(0,0,1,0.5), breaks=2)


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
     las=1 , col='yellow', breaks=5)
hist(d.location[d.m.index=='Model  3'&d.para.order==2] , main="" , xlim=c(0,400), ylab="Number of particles", 
     xlab="Locations", ylim=c(40,0),
     cex.axis = 2, cex.lab = 2,
     las=1 , col='green', breaks=5,add =T)
hist(d.location[d.m.index=='Model  3'&d.para.order==3] , main="" , xlim=c(0,400), ylab="Number of particles", 
     xlab="Locations", ylim=c(40,0),
     cex.axis = 2, cex.lab = 2,
     las=1 , col='darkgreen', breaks=5,add =T)

# Check first 5 columns
print(Indice.m[, 1:5])


library(tidyr)  # For pivot_longer
library(dplyr)  # For data processing

# Convert to data frame and add category labels
df <- as.data.frame(t(Indice.m))  # Transpose matrix (60x4)
colnames(df) <- c("Count_0", "Count_1", "Count_2", "Count_3")  # Column names
df$Time <- 1:60  # Add time column

# Convert to long format
df_long <- df %>%
  pivot_longer(
    cols = -Time,
    names_to = "Category",
    values_to = "Count"
  )

# View the tidied data
head(df_long)


library(ggplot2)

ggplot(df_long, aes(x = Time, y = Count, color = Category)) +
  geom_line(lwd = 1) +  # Draw curves
  geom_point(size = 1.5, alpha = 0.6) +  # Add points (optional)
  scale_color_manual(
    values = c("Count_0" = "blue", "Count_1" = "red", "Count_2" = "green", "Count_3" = "purple"),
    labels = c("0", "1", "2", "3")  # Legend labels
  ) +
  coord_cartesian(ylim = c(0, 80))+
  labs(
    x = "Time step t",
    y = "Number of particles",
    color = "Model"
  ) +
  theme_minimal(base_size = 14) +  # Set base font size to 14
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Increase title font
    axis.title = element_text(size = 16),  # Axis title size
    axis.text = element_text(size = 14),   # Axis text size
    legend.title = element_text(size = 16), # Legend title size
    legend.text = element_text(size = 14),  # Legend text size
    legend.position = "right"  # Legend on the right
  )






# 3.1 Number of particles ----------------------------------------------------------------------
# Count the number of elements in the interval [80,120]
count_75_125 <- sum(EnPara >= 75 & EnPara <= 125, na.rm = TRUE)

# Count the number of elements in the interval [230,270]
count_225_275 <- sum(EnPara >= 225 & EnPara <= 275, na.rm = TRUE)





# Check first 5 columns
print(Loca.m[, 1:5])



# Convert to data frame and add category labels
df <- as.data.frame(t(Loca.m))  # Transpose matrix (60x4)
colnames(df) <- c( "Loca_1", "Loca_2")  # Column names
df$Time <- 1:60  # Add time column

# Convert to long format
df_long <- df %>%
  pivot_longer(
    cols = -Time,
    names_to = "Category",
    values_to = "Count"
  )

# View the tidied data
head(df_long)



ggplot(df_long, aes(x = Time, y = Count, color = Category)) +
  geom_line(lwd = 1) +  # Draw curves
  geom_point(size = 1.5, alpha = 0.6) +  # Add points (optional)
  scale_color_manual(
    values = c("Loca_1"= "lightblue", "Loca_2"= "orange"),
    labels = c("Loca_1", "Loca_2")  # Legend labels
  ) +
  coord_cartesian(ylim = c(0, 80))+
  labs(
    x = "Time step t",
    y = "Number of particles",
    color = "spatial point"
  ) +
  theme_minimal(base_size = 14) +  # Set base font size to 14
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Increase title font
    axis.title = element_text(size = 16),  # Axis title size
    axis.text = element_text(size = 14),   # Axis text size
    legend.title = element_text(size = 16), # Legend title size
    legend.text = element_text(size = 14),  # Legend text size
    legend.position = "right"  # Legend on the right
  )
```

<p align="center">
  <img src="https://github.com/lwhuanyue/Two-stage-Resampling-PF-with-RJMCMC/blob/Images/Hist_21.png" width="45%" />
  <img src="https://github.com/lwhuanyue/Two-stage-Resampling-PF-with-RJMCMC/blob/Images/Hist_23.png" width="45%" /> 
</p>
Figure : The results of our method for model identification and location parameter estimation with $N=80$. The upper half of both subplots shows the bar plots of the position parameters corresponding to the particles with \( k = 2 \). The light red and red bars represent the positions of the first and second discrete spatial points, respectively. In the left subplot, the blue bar plot represents the position parameters corresponding to the particles with \( k = 1 \). The bottom-right subplot shows the bar plot of the position parameters for particles with \( k = 3 \), where yellow, light green, and dark green bars correspond to the positions of the three discrete spatial points.
<p align="center">
  <img src="https://github.com/lwhuanyue/Two-stage-Resampling-PF-with-RJMCMC/blob/Images/k_time.png" width="45%" />
  <img src="https://github.com/lwhuanyue/Two-stage-Resampling-PF-with-RJMCMC/blob/Images/l_time.png" width="45%" /> 
</p>
Figure : Particle identification performance for model indices and discrete spatial point detection with ensemble size $N=80$. The left subfigure displays the temporal evolution of particle counts identifying model indices, showing that the majority of particles recognize models with either $k=2$ or $k=3$.
The right subfigure presents the time-varying counts of particles detecting two discrete spatial points, demonstrating that nearly all particles successfully identify both discrete spatial points as time progresses. To quantify this, the number of particles falling within the intervals [75,125] and [225,275] are respectively used to represent the number of particles identifying the two breakpoints.
<p align="center">
  <img src="https://github.com/lwhuanyue/Two-stage-Resampling-PF-with-RJMCMC/blob/Images/Pre_velocities.png" width="45%" />
</p>
Figure : Estimated velocities. The true pointwise velocity is shown by the black line, alongside the posterior means (solid lines) and pointwise posterior 80$\%$ credible intervals (shaded bands). These intervals represent the regions between the upper $0.1$ quantile and lower $0.1$ quantile of the $N$ particles, as obtained using various methods. Notably, fixing $k = 1$ (green line) or $k = 3$ (purple line) causes significant deviations in the velocity estimates. In contrast, the proposed approach yields the most accurate results, closely matching the true velocity.

The figures and tables in the related article can be reproduced using the `Illustration_RJ2S.R script` and other code.
