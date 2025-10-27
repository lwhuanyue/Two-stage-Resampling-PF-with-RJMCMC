R_path <- 'E:/Study-2024/09.RJ_2stage_resample_Oct24/4. Second_review CODE/'

#source(paste0(R_path, '.R'))
source(paste0(R_path, 'RJ2S_basic function and setting.R'))
source(paste0(R_path, 'Func_RJ2S  add_perturbation.R'))
source(paste0(R_path, 'Func_RJ2S  Para KDE.R'))
source(paste0(R_path, 'Func_RJ2S  test acceptance and ess.R'))
source(paste0(R_path, 'Func_RJ2S Para estimation proposal.R'))
source(paste0(R_path, 'Func_RJ2S Non_para estimation proposal.R'))


# The initial settings --------------------------------------------------------------------

# In this part, we demonstrate an application of our proposed two-stage resampling method
# in the data assimilation process for the one-dimensional shallow water equations.
# We also present: four types of RJMCMC moves (birth/death/split/merge), 
#                  ESS traces, 
#                  acceptance rates per move type, 
#   and finally estimate the posterior distribution of parameters using 
#   both non-parametric and parametric methods.

# The evolution function
#Func: 1-d advection function as the evolution function, which is M(·) in the corresponding paper.
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

# There are some settings

set.seed(1028) # Fix a random seed

#The hyperparameter for the prior distribution of the parameter.
#Poisson distribution
lambda_p <- 2 # lamda, also the prior mean
kmax <- 3     # truncated  Poisson
#Gamma distribution
g_alpha <- 0.4  #alpha
g_beta  <- 0.95 #beta

N <- 80 # The size of ensemble       

#On the specification of spatio-temporal fields
times <- 0:600  # Time steps
sig.eps0 <- 0.2 # Evolution model error sigma^2 
L <- 400        # The length of the spatial domain, then: 0,1,...,399,400
n <- L+1        # n: dimension of the state x
dx = 0.2        # The step size of x 
G = 1           # Loop of the RJMCMC resampling

t.obs <- seq(0,times[length(times)], by = 10) # Observations are obtained at every ten time steps.


# The true parameter and spatial-temporal field ---------------------------

#Real evolution for observations. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sx <- seq(0,1,by=0.0025)
Cini <- 80*sin(60*sx/pi)*sx*(2/3-sx)*(1-sx)*exp(-2*sx); #Initial state, x_0
K = 2 # Two change-points in the field. If 2 change points, then 3 stage.
c.x <- c(0,101,251,n) # The location of change-points
c.v <- c(0.7,0.2,0.4) # The piecewise velocities. The velocity value in every segment.
c.x.true <- c.x; c.v.true <- c.v;
para <- c(K, c.x, c.v, 0.0) # The true parameter.
out <- ode.1D(Cini, times, model2, para, method = "adams", names = c("C"), dimens = n) #The true states over time.
image(out,xlab = 'times',ylab = 'positions')   # The figure of true states over time.
dim(out)                                       # The dimension of the output.
plot(out[601,2:402],type='l');lines(Cini,col='blue') #The initial state x_0 (blue) and state at 601 time step (black).


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
#The DA_process function
DA_Process <- function(Enstat, EnPara, Mthd){
  # Initialize state vectors and likelihood
  En.s.t_1 <- Enstat  # Previous state vector x_(t-1)
  En.s.t <- matrix(NA, nrow = n, ncol = N)   # New state vector x_t
  En.lik <- rep(NA, N)
  
  # Define variables for Bootstrap Resampling
  Mid.En.s.t_1 <- En.s.t_1
  Mid.En.s.t <- En.s.t
  Mid.En.lik <- En.lik
  Mid.EnPara <- EnPara
  
  # Initialize test matrix
  Mtx_test <- matrix(NA, nrow = length(t.obs)-1, ncol = 14)
  
  # Time loop
  for (i.t in 1:(length(t.obs)-1)){
    print(paste('Time step t =',i.t))
    Mtx_Accep <- matrix(0, nrow = 2, ncol = N)
    t1 <- t.obs[i.t]
    t2 <- t.obs[i.t+1]
    
    # Generate perturbed observations
    sta.obs <- out[t2+1, -1] * (1 + sig.eps * rnorm(n))
    
    # Adjust observation error if needed
    if(0){
      # i.t < length(t.obs)/2
      # Later change *9 to *18
      sig.eps <- (10 - (i.t-1)/length(t.obs)*9) * sig.eps0
    } else { 
      sig.eps <- sig.eps0
    }
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Loop_1: State Propagation      #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    # Propagate state vector from En.s.t_1 to En.s.t
    for (i.n in 1:N){
      k_0 <- EnPara[1, i.n]
      para_0 <- EnPara[1:(2*k_0+2), i.n]
      
      # Prepare parameters for evolution
      if(k_0 == 0){
        p4evo0 <- c(k_0, 0, n, para_0[1+k_0+(1:(k_0+1))], 0.0)
      } else {
        p4evo0 <- c(k_0, 0, para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))], 0.0)
      }
      
      # Initial state and forward propagation
      sta.ini <- En.s.t_1[, i.n]
      sta.fore0 <- ode.1D(sta.ini, t1:t2, model2, p4evo0, method = "adams", 
                          names = c("C"), dimens = n)[t2-t1+1, 1:n+1]
      
      # Calculate likelihood
      lik0 <- lik.norm(sig.eps, sta.obs[Obs.site], sta.fore0[Obs.site])
      En.lik[i.n] <- lik0
      En.s.t[, i.n] <- sta.fore0
    }
    
    # Optional plotting
    if(plot_t){
      plot(out[t2+1, -1], type = 'l', lwd = 3, ylab = paste('t =', t2))
      points(Obs.site, sta.obs[Obs.site])
      for(i.n in 1:N){
        lines(En.s.t[, i.n], col = 'red')
      }
    }
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Loop_2: Bootstrap Resampling   #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    # Perform resampling based on likelihood
    if(max(En.lik) != 0){
      Fathers <- sample(1:N, N, En.lik, replace = TRUE)
    } else {
      Fathers <- 1:N
    }
    
    # Update ensemble members with lower likelihood
    II2 <- which(En.lik < median(En.lik))
    for (i.n in II2) {
      Mid.En.s.t_1[, i.n] <- En.s.t_1[, Fathers[i.n]]
      Mid.En.s.t[, i.n] <- En.s.t[, Fathers[i.n]]
      Mid.En.lik[i.n] <- En.lik[Fathers[i.n]]
      Mid.EnPara[, i.n] <- EnPara[, Fathers[i.n]]
    }
    
    # Keep ensemble members with higher likelihood unchanged
    Mid.En.s.t_1[, setdiff(1:N, II2)] <- En.s.t_1[, setdiff(1:N, II2)]
    Mid.En.s.t[, setdiff(1:N, II2)] <- En.s.t[, setdiff(1:N, II2)]
    Mid.En.lik[setdiff(1:N, II2)] <- En.lik[setdiff(1:N, II2)]
    Mid.EnPara[, setdiff(1:N, II2)] <- EnPara[, setdiff(1:N, II2)]
    
    # Update main variables
    En.s.t_1 <- Mid.En.s.t_1
    En.s.t <- Mid.En.s.t
    En.lik <- Mid.En.lik
    EnPara <- Mid.EnPara
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Loop_3: RJMCMC Resampling      #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    for (i.n in 1:N){   # Ensemble member loop
      # Reinitialize parameters after resampling
      for (gg in 1:G) {
        k_0 <- EnPara[1, i.n]
        para_0 <- EnPara[1:(2*k_0+2), i.n]
        
        if(k_0 == 0){
          p4evo0 <- c(k_0, 0, n, para_0[1+k_0+(1:(k_0+1))], 0.0)
        } else {
          p4evo0 <- c(k_0, 0, para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))], 0.0)
        }
        
        # Propose new parameters
        if (Mthd == 1){
          prpd.index <- Prop.Index(k_0, c = 0.3, lambda = lambda_p, kmax = kmax)$index
        } else {
          prpd.index <- sample(c(3, 4), 1)
        }
        
        propose1 <- Prp.Para(para_0, prpd.index)
        para_1 <- propose1$para
        k_1 <- para_1[1]
        
        if(k_1 == 0){
          p4evo1 <- c(k_1, 0, n, para_1[1+k_1+(1:(k_1+1))], 0.0)
        } else {
          p4evo1 <- c(k_1, 0, para_1[1+(1:(k_1))], n, para_1[1+k_1+(1:(k_1+1))], 0.0)
        }
        
        # Propagate with proposed parameters
        sta.ini <- En.s.t_1[, i.n]
        sta.fore1 <- ode.1D(sta.ini, t1:t2, model2, p4evo1, method = "adams", 
                            names = c("C"), dimens = n)[t2-t1+1, 1:n+1]
        
        lik1 <- lik.norm(sig.eps, sta.obs[Obs.site], sta.fore1[Obs.site])
        
        if(plot_t){
          lines(sta.fore1, col = 'blue')
        }
        
        # Calculate likelihoods for comparison
        lik0 <- lik.norm(sig.eps, sta.obs[Obs.site], sta.fore0[Obs.site])
        lik1 <- lik.norm(sig.eps, sta.obs[Obs.site], sta.fore1[Obs.site])
        
        Mtx_Accep[1, i.n] <- prpd.index
        
        # Metropolis-Hastings acceptance step
        if((lik1 * propose1$ratio) > lik0 * runif(1, 0, 1)){
          EnPara[, i.n] <- rep(NA, 2*kmax+2)
          EnPara[1:(2*k_1+2), i.n] <- para_1
          sta.fore0 <- sta.fore1
          Mtx_Accep[2, i.n] <- 1
        }
      }
      
      En.s.t[, i.n] <- sta.fore0
      
      # Apply Kalman gain if enabled
      if(kalmangain){
        En.s.t[, i.n] <- sta.fore0 + 
          K.true %*% as.matrix(sta.obs[Obs.site] - sta.fore0[Obs.site] - sig.eps * rnorm(m))
      }
    }
    
    # Update state for next time step
    En.s.t_1 <- En.s.t
    
    # Compute acceptance statistics
    Mtx_test[i.t, ] <- compute_accept_stats(Mtx_Accep)
    Mtx_test[i.t, 14] <- (sum(En.lik))^2 / sum(En.lik^2)
  }
  
  return(list('EnS' = En.s.t_1, 'EnP' = EnPara, 'Acp' = Mtx_test))
}


# The example. ------------------------------------------------------------

Res <- DA_Process(Enstat, EnPara,1) 
#Enstat: The initial ensemble of particles of states
#Enara:  The initial ensemble of particles of parameters
#Mthd = 1: the proposed 2-stage RJMCMC resampling method.
# The observations are inside the DA_Process function, by adding a Gaussian noise to the true:
#      sta.obs <- out[t2+1,-1] *(1+ sig.eps*rnorm(n))
# The results for the first five time points, 
# associated with the four types of proposal acceptance rates.
head(Res$Acp,5)
#The four types are: (a) adding a discontinuity point,
#                    (b) removing a discontinuity point, 
#                    (c) changing one of the velocity, 
#                and (d) changing the position of one discontinuity point.
# Each row in Res$Acp represents the result at each time step.
#'   \item [1-4] Count of proposal types (a)-(d) 
#'   \item [5-8] Count of accepted proposal types (a)-(d) (%)
#'   \item [9] Sum of accepted proposal
#'   \item [10-13] Acceptance rates for types (a)-(d)
#'   \item [14] Effective sample size (ESS)



# Plot the RJMCMC move proposals and  acceptance rates per move type --------
library(ggplot2)
library(tidyr)
library(dplyr)
 
#~~~~ 1.

# Extract first four columns and add time step
plot_data <- data.frame(
  Time = 1:nrow(Res$Acp),
  Type_a = Res$Acp[, 1],
  Type_b = Res$Acp[, 2], 
  Type_c = Res$Acp[, 3],
  Type_d = Res$Acp[, 4]
)

plot_data_long <- pivot_longer(plot_data, 
                               cols = -Time, 
                               names_to = "Transition_Type", 
                               values_to = "Proposal_Count")

ggplot(plot_data_long, aes(x = Time, y = Proposal_Count, color = Transition_Type)) +
  geom_line(lwd = 1.2) +
  geom_point(size = 1.2) +
  labs(title = "Proposal Counts of Four Transition Types",
       x = "Time step",
       y = "Proposal Count",
       color = "Transition") +
  theme_minimal() +
  scale_color_manual(values = c("Type_a" = "#E41A1C", 
                                "Type_b" = "#377EB8", 
                                "Type_c" = "#4DAF4A", 
                                "Type_d" = "#984EA3"))+
  theme(text = element_text(size = 16)) 


#~~~~ 2.
plot_data_accept <- data.frame(
  Time = 1:nrow(Res$Acp),
  Type_a = Res$Acp[, 5],  
  Type_b = Res$Acp[, 6],  
  Type_c = Res$Acp[, 7],  
  Type_d = Res$Acp[, 8]   
)

plot_data_accept_long <- pivot_longer(plot_data_accept, 
                                      cols = -Time, 
                                      names_to = "Transition_Type", 
                                      values_to = "Acceptance_Count")

ggplot(plot_data_accept_long, aes(x = Time, y = Acceptance_Count, color = Transition_Type)) +
  geom_line(lwd = 1.2) +
  geom_point(size = 1.2) +
  labs(title = "Acceptance Counts of Four Transition Types",
       x = "Time step",
       y = "Acceptance Count",
       color = "Transition") +
  theme_minimal() +
  scale_color_manual(values = c("Type_a" = "#E41A1C", 
                                "Type_b" = "#377EB8", 
                                "Type_c" = "#4DAF4A", 
                                "Type_d" = "#984EA3")) +
  theme(text = element_text(size = 16))


#~~~~ 2. 
plot_data_acceptance_rate <- data.frame(
  Time = 1:nrow(Res$Acp),
  Type_a = Res$Acp[, 10],  
  Type_b = Res$Acp[, 11],  
  Type_c = Res$Acp[, 12],  
  Type_d = Res$Acp[, 13]   
)


plot_data_rate_long <- pivot_longer(plot_data_acceptance_rate, 
                                    cols = -Time, 
                                    names_to = "Transition_Type", 
                                    values_to = "Acceptance_Rate")


ggplot(plot_data_rate_long, aes(x = Time, y = Acceptance_Rate, color = Transition_Type)) +
  geom_line(lwd = 1.2) +
  geom_point(size = 1.2) +
  labs(title = "Acceptance Rates of Four Transition Types",
       x = "Time step",
       y = "Acceptance Rate",
       color = "Transition") +
  theme_minimal() +
  scale_color_manual(values = c("Type_a" = "#E41A1C", 
                                "Type_b" = "#377EB8", 
                                "Type_c" = "#4DAF4A", 
                                "Type_d" = "#984EA3")) +
  theme(text = element_text(size = 16)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))  




# posterior summaries of parameter -----------------------------------------

# The states.
Enstat_post <- Res$EnS
m.true  <- out[601,-1]
mean <- rowMeans(Enstat_post)
par(mfrow=c(1,1),mar=c(2,2,0,0.5)+2.5)
plot(m.true,type ='l', ylim = c(min(En.s.t_1)-1,6),
     cex.axis = 2, cex.lab =2, cex.main =2,
     main='The pointwise posterior 80% credible intervals',xlab = 'Position',ylab = 'Value')
for (i.n in 1:N){lines(Enstat_post[,i.n],col = 'lightcoral',lwd=2)}
plot.PI(Enstat_post,'blue');lines(m.true,col='black',lwd=2);lines(mean,col='blue',lwd=1.2)


# #Part 2 Comparison ------------------------------------------------------


#  comparison with pMCMC methods ------------------------------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Initial Comparison Setup ~~~~~~~~~~~~~~~~~~~#

#' Here, we compare our method with the pMCMC method. 
#' Enstat, Enstat_1S, Enstat_2S, and Enstat_3S represent the initial state ensembles 
#' for the proposed method and the pMCMC methods with k=1, 2, and 3 breakpoints respectively, 
#' which are set to the same ensemble here.
Enstat <- matrix(NA, nrow = n, ncol = N)                                           #
for(i in 1:N){Enstat[,i] <- Cini*(1+rep(rnorm(1,0,0.1),n))}                        #
Enstat_1S  <- Enstat # Classical MCMC with 1 breakpoint                            #
Enstat_2S  <- Enstat # Classical MCMC with 2 breakpoints (same as real data), 3 steps#
Enstat_3S  <- Enstat # Classical MCMC with 3 breakpoints, 4 steps                  # 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Comparison Parameter Ensembles ~~~~~~~~~~~~~~~~~~~#
#' The ensembls EnPara_1S, EnPara_2S, and EnPara_3S represent the initial parameter 
#'  corresponding to k=1, 2, and 3 breakpoints, respectively.

# Classical MCMC with 1 breakpoint                                                 #
EnPara_1S <- matrix(NA, ncol = N, nrow = 2*kmax+2)
EnPara_1S[1,] <- rep(1,N) 

for(i in 1:N){
  kk <- EnPara_1S[1,i]
  if(kk>0){
    c.x <- sample.locs(kk,L)
    EnPara_1S[1+(1:(kk)),i] <-   c.x}
  c.v <- rgamma(kk+1,shape = 2,rate = 1)
  EnPara_1S[kk+1+(1:(kk+1)),i] <- c.v  
}

# Classical MCMC with 2 breakpoints (same as real data), 3 steps                          #
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


# Classical MCMC with 3 breakpoints, 4 steps                                      # 
#
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

# Comare 1.1 Compared pMCMC Method: Number of Breakpoints = 1 -------------------------------------------------------------
En.s.t_1_1S <- Enstat_1S # Represents the state vector at previous time step x_(t-1)
En.s.t_1S <- matrix(NA,nrow = n, ncol = N)   # Represents the state vector at new time step x_t
En.lik_1S <- rep(NA,N)

# Definitions for Bootstrap Resampling.
Mid.En.s.t_1_1S <- En.s.t_1_1S
Mid.En.s.t_1S   <- En.s.t_1S
Mid.En.lik_1S   <- En.lik_1S
Mid.EnPara_1S   <- EnPara_1S
for (i.t in 1:(length(t.obs)-1)){   # Time Loop. i.t=1 i.t=2
  t1 <- t.obs[i.t]; t2 <- t.obs[i.t+1]
  sta.obs <- out[t2+1,-1] *(1+ sig.eps*rnorm(n))
  #sta.obs <- out[t2+1,-1] + sig.eps*rnorm(n) # Perturbed observations.
  if(0){
    #i.t < length(t.obs)/2
    # Later change *9 to *18
    sig.eps <- (10-(i.t-1)/length(t.obs)*9)*sig.eps0
  } else { 
    sig.eps <- sig.eps0
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Loop_1:       Propagate         #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # In this step, propagate the state vector (En.s.t_1) to the new time (En.s.t)
  for (i.n in 1:N){
    k_0 <- EnPara_1S[1,i.n]
    para_0 <- EnPara_1S[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    } else {
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    # This was the original: p4evo0 <- c(k_0,0,para_0[2*(1:k_0)+1], n, para_0[2*(0:k_0)+2],0.0)
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
  # Loop_2:       Bootstrap         #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  if(max(En.lik_1S) != 0){
    Fathers <- sample(c(1:N),N,En.lik_1S,replace = TRUE)
  } else {Fathers <- c(1:N)}
  # Fathers <- c(1:N)
  # Resampling replacement. Caution: Must use intermediate variables!!!!
  # Only consider the lower half
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
  # After replacement is complete, assign the matrices back
  En.s.t_1_1S <- Mid.En.s.t_1_1S
  En.s.t_1S   <- Mid.En.s.t_1S
  En.lik_1S   <- Mid.En.lik_1S
  EnPara_1S   <- Mid.EnPara_1S
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Loop_3:  RJMCMC Re-sampling     #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  for (i.n in 1:N){   # Loop over ensemble members. i.n = 2
    # Due to the resampling in the first step, repeat the input parameter step from the first step
    k_0 <- EnPara_1S[1,i.n]
    para_0 <- EnPara_1S[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    } else {
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    # Next, propose new parameters
    prpd.index <- sample(c(3,4),1)
    propose1 <- Prp.Para(para_0,prpd.index) 
    para_1 <- propose1$para; print(para_1) # Now we need to define the function, propose a new function
    k_1 <- para_1[1]
    if(k_1 == 0){
      p4evo1 <- c(k_1,0,n,para_1[1+k_1+(1:(k_1+1))],0.0)
    } else {
      p4evo1 <- c(k_1,0,para_1[1+(1:(k_1))], n, para_1[1+k_1+(1:(k_1+1))],0.0)
    }
    sta.ini <- En.s.t_1_1S[,i.n]
    sta.fore1 <- ode.1D(sta.ini, t1:t2, model2, p4evo1, method = "adams", 
                        names = c("C"), dimens = n)[t2-t1+1,1:n+1]
    
    lik1 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore1[Obs.site])
    
    # Compare En.lik[i.n] and lik0
    
    # Pay attention to two points: one is the Jacobian, the other is not to mistake the initial vector
    
    # Then add if
    
    #~~~~~~~~~~~~~~~~End loop # ~~~~~~~~~~
    
    # Original unupdated part
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


# Comare 1.2 Compared pMCMC Method: Number of Breakpoints = 2 -------------------------------------------------------------
En.s.t_1_2S <- Enstat_2S # Represents the state vector at previous time step x_(t-1)
En.s.t_2S <- matrix(NA,nrow = n, ncol = N)   # Represents the state vector at new time step x_t
En.lik_2S <- rep(NA,N)

# Definitions for Bootstrap Resampling.
Mid.En.s.t_1_2S <- En.s.t_1_2S
Mid.En.s.t_2S   <- En.s.t_2S
Mid.En.lik_2S   <- En.lik_2S
Mid.EnPara_2S   <- EnPara_2S
for (i.t in 1:(length(t.obs)-1)){   # Time Loop. i.t=1 i.t=2
  t1 <- t.obs[i.t]; t2 <- t.obs[i.t+1]
  sta.obs <- out[t2+1,-1] *(1+ sig.eps*rnorm(n))
  #sta.obs <- out[t2+1,-1] + sig.eps*rnorm(n) # Perturbed observations.
  if(0){
    #i.t < length(t.obs)/2
    # Later change *9 to *18
    sig.eps <- (10-(i.t-1)/length(t.obs)*9)*sig.eps0
  } else { 
    sig.eps <- sig.eps0
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Loop_1:       Propagate         #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # In this step, propagate the state vector (En.s.t_1) to the new time (En.s.t)
  for (i.n in 1:N){
    k_0 <- EnPara_2S[1,i.n]
    para_0 <- EnPara_2S[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    } else {
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    # This was the original: p4evo0 <- c(k_0,0,para_0[2*(1:k_0)+1], n, para_0[2*(0:k_0)+2],0.0)
    sta.ini <- En.s.t_1_2S[,i.n]
    sta.fore0 <- ode.1D(sta.ini, t1:t2, model2, p4evo0, method = "adams", 
                        names = c("C"), dimens = n)[t2-t1+1,1:n+1]
    lik0 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore0[Obs.site])
    En.lik_2S[i.n] <- lik0
    En.s.t_2S[,i.n] <- sta.fore0
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Loop_2:       Bootstrap         #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  if(max(En.lik_2S) != 0){
    Fathers <- sample(c(1:N),N,En.lik_2S,replace = TRUE)
  } else {Fathers <- c(1:N)}
  # Fathers <- c(1:N)
  # Resampling replacement. Caution: Must use intermediate variables!!!!
  # Only consider the lower half
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
  # After replacement is complete, assign the matrices back
  En.s.t_1_2S <- Mid.En.s.t_1_2S
  En.s.t_2S   <- Mid.En.s.t_2S
  En.lik_2S   <- Mid.En.lik_2S
  EnPara_2S   <- Mid.EnPara_2S
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Loop_3:  RJMCMC Re-sampling     #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  for (i.n in 1:N){   # Loop over ensemble members. i.n = 2
    # Due to the resampling in the first step, repeat the input parameter step from the first step
    k_0 <- EnPara_2S[1,i.n]
    para_0 <- EnPara_2S[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    } else {
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    # Next, propose new parameters
    prpd.index <- sample(c(3,4),1)
    # prpd.index <- Prop.Index(k_0,c=0.3,lambda=lambda_p,kmax=kmax)$index
    propose1 <- Prp.Para(para_0,prpd.index) 
    para_1 <- propose1$para; print(para_1) # Now we need to define the function, propose a new function
    k_1 <- para_1[1]
    if(k_1 == 0){
      p4evo1 <- c(k_1,0,n,para_1[1+k_1+(1:(k_1+1))],0.0)
    } else {
      p4evo1 <- c(k_1,0,para_1[1+(1:(k_1))], n, para_1[1+k_1+(1:(k_1+1))],0.0)
    }
    sta.ini <- En.s.t_1_2S[,i.n]
    sta.fore1 <- ode.1D(sta.ini, t1:t2, model2, p4evo1, method = "adams", 
                        names = c("C"), dimens = n)[t2-t1+1,1:n+1]
    
    lik1 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore1[Obs.site])
    
    # Compare En.lik[i.n] and lik0
    
    # Pay attention to two points: one is the Jacobian, the other is not to mistake the initial vector
    
    # Then add if
    
    #~~~~~~~~~~~~~~~~End loop # ~~~~~~~~~~
    
    # Original unupdated part
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


# Comare 1.3 Compared pMCMC Method: Number of Breakpoints = 3 --------------------------------------------------------------
En.s.t_1_3S <- Enstat_3S # Represents the state vector at previous time step x_(t-1)
En.s.t_3S <- matrix(NA,nrow = n, ncol = N)   # Represents the state vector at new time step x_t
En.lik_3S <- rep(NA,N)

# Definitions for Bootstrap Resampling.
Mid.En.s.t_1_3S <- En.s.t_1_3S
Mid.En.s.t_3S   <- En.s.t_3S
Mid.En.lik_3S   <- En.lik_3S
Mid.EnPara_3S   <- EnPara_3S
for (i.t in 1:(length(t.obs)-1)){   # Time Loop. i.t=1 i.t=2
  t1 <- t.obs[i.t]; t2 <- t.obs[i.t+1]
  sta.obs <- out[t2+1,-1] *(1+ sig.eps*rnorm(n))
  #sta.obs <- out[t2+1,-1] + sig.eps*rnorm(n) # Perturbed observations.
  if(0){
    #i.t < length(t.obs)/2
    # Later change *9 to *18
    sig.eps <- (10-(i.t-1)/length(t.obs)*9)*sig.eps0
  } else { 
    sig.eps <- sig.eps0
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Loop_1:       Propagate         #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # In this step, propagate the state vector (En.s.t_1) to the new time (En.s.t)
  for (i.n in 1:N){
    k_0 <- EnPara_3S[1,i.n]
    para_0 <- EnPara_3S[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    } else {
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    # This was the original: p4evo0 <- c(k_0,0,para_0[2*(1:k_0)+1], n, para_0[2*(0:k_0)+2],0.0)
    sta.ini <- En.s.t_1_3S[,i.n]
    sta.fore0 <- ode.1D(sta.ini, t1:t2, model2, p4evo0, method = "adams", 
                        names = c("C"), dimens = n)[t2-t1+1,1:n+1]
    lik0 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore0[Obs.site])
    En.lik_3S[i.n] <- lik0
    En.s.t_3S[,i.n] <- sta.fore0
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Loop_2:       Bootstrap         #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  if(max(En.lik_3S) != 0){
    Fathers <- sample(c(1:N),N,En.lik_3S,replace = TRUE)
  } else {Fathers <- c(1:N)}
  # Fathers <- c(1:N)
  # Resampling replacement. Caution: Must use intermediate variables!!!!
  # Only consider the lower half
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
  # After replacement is complete, assign the matrices back
  En.s.t_1_3S <- Mid.En.s.t_1_3S
  En.s.t_3S   <- Mid.En.s.t_3S
  En.lik_3S   <- Mid.En.lik_3S
  EnPara_3S   <- Mid.EnPara_3S
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Loop_3:  RJMCMC Re-sampling     #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  for (i.n in 1:N){   # Loop over ensemble members. i.n = 2
    # Due to the resampling in the first step, repeat the input parameter step from the first step
    k_0 <- EnPara_3S[1,i.n]
    para_0 <- EnPara_3S[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    } else {
      p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    }
    # Next, propose new parameters
    prpd.index <- sample(c(3,4),1)
    propose1 <- Prp.Para(para_0,prpd.index) 
    para_1 <- propose1$para; print(para_1) # Now we need to define the function, propose a new function
    k_1 <- para_1[1]
    if(k_1 == 0){
      p4evo1 <- c(k_1,0,n,para_1[1+k_1+(1:(k_1+1))],0.0)
    } else {
      p4evo1 <- c(k_1,0,para_1[1+(1:(k_1))], n, para_1[1+k_1+(1:(k_1+1))],0.0)
    }
    sta.ini <- En.s.t_1_3S[,i.n]
    sta.fore1 <- ode.1D(sta.ini, t1:t2, model2, p4evo1, method = "adams", 
                        names = c("C"), dimens = n)[t2-t1+1,1:n+1]
    
    lik1 <- lik.norm(sig.eps,sta.obs[Obs.site],sta.fore1[Obs.site])
    
    # Compare En.lik[i.n] and lik0
    
    # Pay attention to two points: one is the Jacobian, the other is not to mistake the initial vector
    
    # Then add if
    
    #~~~~~~~~~~~~~~~~End loop # ~~~~~~~~~~
    
    # Original unupdated part
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
# Plotting ----------------------------------------------------------------------
En.s.t_1 <- Res$EnS

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
# Calculate MSE
mean((mean_rj-m.true)**2)
mean((mean_1S-m.true)**2)
mean((mean_2S-m.true)**2)
mean((mean_3S-m.true)**2)

# Predictions:
# 2.1 Adding Prediction Comparison -------------------------------------------------------------

# First examine the logic
# For prediction, we start with the final time step's ensemble of states and parameters
# and the true states and parameters
En.s <- En.s.t_1
En.p <- EnPara

En.s <- En.s.t_1_1S
En.p <- EnPara_1S
s.true <- out[length(times), 1+(1:n)]
p.true <- para
times.p <- times[length(times)] + c(0:50)
out.pre <- ode.1D(s.true, times.p, model2, p.true, method = "adams", 
                  names = c("C"), dimens = n)[length(times.p), 1+(1:n)]
En.predict <- matrix(NA, ncol = N, nrow = n)
for(i.n in 1:N){
  k_0 <- En.p[1,i.n]
  para_0 <- En.p[1:(2*k_0+2),i.n]; print(para_0)
  if(k_0 == 0){
    p4evo0 <- c(k_0, 0, n, para_0[1+k_0+(1:(k_0+1))], 0.0)
  } else {
    p4evo0 <- c(k_0, 0, para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))], 0.0)
  }
  En.predict[,i.n] <- ode.1D(En.s[,i.n], times.p, model2, p4evo0, method = "adams", 
                             names = c("C"), dimens = n)[length(times.p), 1+(1:n)]
}
plot(out.pre, type = 'l')
for(i.n in 1:N){lines(En.predict[,i.n], col = 'red')}
plot.PI(En.predict, 'blue')
lines(out.pre, lwd = 3)

#

Pre_En <- function(En.s, En.p, s.true, p.true, times.p){
  arg_name <- as.character(match.call()$En.s)
  out.pre <- ode.1D(s.true, times.p, model2, p.true, method = "adams", 
                    names = c("C"), dimens = n)[length(times.p), 1+(1:n)]
  N <- dim(En.s)[2]
  En.predict <- matrix(NA, ncol = N, nrow = n)
  for(i.n in 1:N){
    k_0 <- En.p[1,i.n]
    para_0 <- En.p[1:(2*k_0+2),i.n]; print(para_0)
    if(k_0 == 0){
      p4evo0 <- c(k_0, 0, n, para_0[1+k_0+(1:(k_0+1))], 0.0)
    } else {
      p4evo0 <- c(k_0, 0, para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))], 0.0)
    }
    En.predict[,i.n] <- ode.1D(En.s[,i.n], times.p, model2, p4evo0, method = "adams", 
                               names = c("C"), dimens = n)[length(times.p), 1+(1:n)]
  }
  plot(out.pre, type = 'l', ylab = 'The state', xlab = 'Locations', cex.axis = 1.5, cex.lab = 1.5,
       main = paste("redictions by:", arg_name))
  for(i.n in 1:N){lines(En.predict[,i.n], col = 'lightblue')}
  plot.PI(En.predict, 'blue')
  lines(out.pre, lwd = 3, col = 'red')
  mean_pre <- rowMeans(En.predict);
  lines(mean_pre, col = 'darkblue', lwd = 2)
  MSPE <- mean((mean_pre - out.pre)^2)
  return(list('En' = En.predict, 'MSPE' = MSPE, 'Pre' = mean_pre))
}

Res_r <- Pre_En(En.s.t_1, Res$EnP, out[length(times), 1+(1:n)], para, times.p)
Res_1 <- Pre_En(En.s.t_1_1S, EnPara_1S, out[length(times), 1+(1:n)], para, times.p)
Res_2 <- Pre_En(En.s.t_1_2S, EnPara_2S, out[length(times), 1+(1:n)], para, times.p)
Res_3 <- Pre_En(En.s.t_1_3S, EnPara_3S, out[length(times), 1+(1:n)], para, times.p)

Res_r$MSPE
Res_2$MSPE
Res_1$MSPE
Res_3$MSPE


# 2.2 Velocity ------------------------------------------------------------

V_Points <- function(para_0,n){
  k_0 <- para_0[1]
  if(k_0 >0){c.x <- c(1,para_0[1+(1:(k_0))],n)
  c.v <- para_0[1+k_0+(1:(k_0+1))]
  v_vec <- rep(NA, n)
  for (i.k in 1:(k_0+1)) {
    v_vec[c.x[i.k]:c.x[i.k+1]] <- rep(c.v[i.k],c.x[i.k+1]-c.x[i.k]+1)}
  }else{
      v_vec <- rep(para_0[2],n)
    }
  return(v_vec)
}
V_Plot <- function(En.P,col1,col2){
  V_matrix <-matrix(NA,ncol = N,nrow = n)
  for (i.n in 1:N) {V_matrix[,i.n] <- V_Points(En.P[,i.n],n)}
  plot.PI(V_matrix,col1)
  lines(rowMeans(V_matrix),lwd=3,col=col2)
}
par(mfrow=c(1,1),mar=c(2,2,0,0.5)+2.5)
plot(V_Points(c(K,101,251,c.v.true),n),col='black',ylim=c(0.1,2),
     cex.axis = 1.5, cex.lab = 1.5,
     xlab ='Locations',ylab ='The pointwise velocities')
V_Plot(Res$EnP,'blue','darkblue')
V_Plot(EnPara_2S,'orange','red')
V_Plot(EnPara_1S,'green','darkgreen')
V_Plot(EnPara_3S,'pink','purple')
legend("topright",legend=c("The new method","pMCMC-1","pMCMC-2","pMCMC-3","The true velocities"),
       col=c("blue","green","red","purple","black"),                
       lty=1,lwd=3)
points(V_Points(c(K,101,251,c.v.true),n))



# 3. Plot the results in the paper ----------------------------------------
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

