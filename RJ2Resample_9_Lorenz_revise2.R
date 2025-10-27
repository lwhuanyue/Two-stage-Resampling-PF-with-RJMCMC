# 2024-11-25 
source('E:/Study-2024/09.RJ_2stage_resample_Oct24/4. Second_review CODE/RJ2Resample_Func_L96_revise2.R')


# These are our constants
Ns <- 40  # Number of variables
FF <- 2  # Forcing


L96 <- function(time, state, parms,...) {
  # Setting up vector
  N <- length(state)
  d <- rep(0, N)
  FF <- parms
  # Loops over indices (with operations and R modulo handling edge cases)
  for (i in 1:N) {
    index.l <- c(i-2,i-1,i,i+1) %% N
    index.l[index.l == 0] <- N
    d[i] <- (state[index.l[4]] - state[index.l[1]]) * state[index.l[2]] - state[index.l[3]] + FF
  }
  return(list(d))
}

x0 <- rep(FF, Ns)  # Initial state (equilibrium)
x0[1] <- x0[1] + 0.01  # Add small perturbation to the first variable
#x0[20] <- x0[1] + 0.01
#x0 <- 1+0.2*sin(2*pi*c(1:Ns)/40);plot(x0,type ='l')
delta_t <- 0.03
t <- seq(0, 30, 0.03)

x <- ode(y = x0, times = t, func = L96, parms=FF,method = "ode45")

# Plot the first three variables
# scatter3D(x[,2], x[,3], x[,4], 
#           phi = 2, type = "l",
#           col = ramp.col(col=c("cyan","red"),
#                          n=length(x[, 3])),ticktype = "detailed", 
#           lwd = 2, bty = "g")
# image(x[,2:30])
# plot(x[,5],type= 'l')
# 
# # 1.1 The true model and observations. ------------------------------------
# #Firstly, define the model3
# L96_M2 <- function(time, state, parms,...) {
#   # Setting up vector
#   N <- length(state)
#   d <- rep(0, N)
#   FF <- parms
#   # Loops over indices (with operations and R modulo handling edge cases)
#   for (i in 1:N) {
#     index.l <- c(i-2,i-1,i,i+1,i+2) %% N
#     index.l[index.l == 0] <- N
#     d[i] <- (1)*( (state[index.l[2]] - state[index.l[5]]) * state[index.l[4]] - state[index.l[3]]  + FF)
#   }
#   return(list(d))
# }
# 
# Nonlmodel <- function(x,parms,t){ #t now is the vector of time points.
#   FF1  <- parms[1]
#   IND <- parms[2]
#   x0 <- x
#   if(IND==1){res <- ode(y = x0, times = t, func = L96, parms=FF1,method = "ode45")} else{
#     res <- ode(y = x0, times = t, func = L96_M2, parms=FF1,method = "ode45")
#   }
#   return(res)
# }
# 
# 
# x0 <- rep(1, Ns)  # Initial state (equilibrium)
# x0[1] <- x0[1] + 0.01
# mu.true1 <- c(2, 1); mu.true2 <- c(1, 2); delta=0.05 # Case 1
# #mu.true1 <- c(6, 1); mu.true2 <- c(4, 2); delta=0.01 # Case 2
# 
# mu.prior <- c(6, 1)
# mu.update <- mu.prior
# t.mid <- 500
# true.out1 <- Nonlmodel(x0, mu.true1, t[1:t.mid])
# Cini2 <- as.vector(true.out1[t.mid,])
# true.out2 <- Nonlmodel(Cini2[(1:Ns)+1], mu.true2, t[(t.mid):length(t)])
# true.out <- matrix(NA, ncol = Ns,nrow = length(t))
# true.out[1:t.mid,] <- true.out1[,(1:Ns)+1]
# true.out[(t.mid+1):length(t),] <- true.out2[2:(length(t)-t.mid+1),(1:Ns)+1]
# image(true.out)
# plot(true.out[,33],type= 'l')
# plot(true.out[1,],type= 'l')
# res <- ode(y = x0, times = t, func = L96, parms=FF,method = "ode45")
# plot(true.out[,4]);lines(res[,5],col ='yellow',lwd=3); lines(x[,5],col ='red');
# lines(true.out[,40],col = 'green')
# plot(true.out[,1]);lines(res[,2],col ='yellow',lwd=3); lines(x[,2],col ='red');
# scatter3D(true.out[,2], true.out[,3], true.out[,1], 
#           phi = 2, type = "l",
#           col = ramp.col(col=c("cyan","red"),
#                          n=length(x[, 3])),ticktype = "detailed", 
#           lwd = 2, bty = "g")



# Piecewise Lorenz96 --------------------------------------------------------------

# Change FF to FF[i]
# Here parms is (k,cx,cv)
L96_RJ <- function(time, state, parms,...) {
  # Setting up vector
  n <- length(state)
  d <- rep(0, n)
  #FF <- parms
  if(parms[1]!=0){FF <- V_Points(parms,n)}else{FF <- rep(parms[2],n)}
  # Loops over indices (with operations and R modulo handling edge cases)
  for (i in 1:n) {
    index.l <- c(i-2,i-1,i,i+1) %% n
    index.l[index.l == 0] <- n
    d[i] <- (state[index.l[4]] - state[index.l[1]]) * state[index.l[2]] - state[index.l[3]] + FF[i]
  }
  return(list(d))
}

x2 <- ode(y = x0, times = t, func = L96_RJ, parms= c(0,2),method = "ode45")
x2-x



# 2. Main program below --------------------------------------------------------------

set.seed(1028)
set.seed(10)
lambda_p <- 1 # Mean of Poisson prior
kmax <- 2
N <- 100


t.obs <- c((0:100)*10);  Obs.site <- c(1:40)


sig.eps0 <- 0.1 #2

L <- 39 #
n <- L+1
# Parameters for Gamma distribution
g_alpha <- 0.4
g_beta  <- 0.95
kalmangain <- 1; plot_t = 0; #Mthd = 1 
dx = 0.2
G = 1
# 02. The prior function --------------------------------------------------

#2.1 --- The prior of the model index. ---
# 'the true model following the k-th model' means there are k points, 
#      dividing the space into k+1 parts.


# Use sample for truncated Poisson
dpois(c(0:kmax), lambda = lambda_p) #
sample(c(0:kmax),N,replace=TRUE,dpois(c(0:kmax), lambda = lambda_p)) 

# For breakpoint location selection, sample (2k+1) points from [0,L], choose the middle k points

set.seed(1)
sample.locs(5,100)

# Prior distribution for k+1 segment speeds
#rgamma(k+1,shape = 2,rate = 1)

# Then for any given k, sample corresponding parameters

# Sample from prior distribution
EnPara <- matrix(NA, ncol = N, nrow = 2*kmax+2)
EnPara[1,] <- sample(c(0:kmax),N,replace=TRUE,dpois(c(0:kmax), lambda = lambda_p)) 
set.seed(111)
for(i in 1:N){
  kk <- EnPara[1,i]
  if(kk>0){
    c.x <- sample.locs(kk,L)
    EnPara[1+(1:(kk)),i] <-   c.x}
  # Redefine gamma distribution parameters
  c.v <- rgamma(kk+1,shape=g_alpha,rate=g_beta)
  EnPara[kk+1+(1:(kk+1)),i] <- c.v  
}



#~~~~~~~~~~~~~~~~~~~~~~~~~Add for Part 2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Other ensemble of parameters.


# 03.Initial conditions -----------------------------------------------------------------


#For initial conditions (true conditions)
#~~~~ 3.1 ~~~~ The same setting with the previous part.



#Pre work
### true forecast distribution
locs=1:n
distmat=rdist(locs,locs)  #
scale=20
sig.eps=.1
kappa=3



distmat=rdist(1:n,1:n); scale=100; sig.eps=.1 ;kappa=3; df.t=2

#Obs.site <- sort(sample(c(1:L),ceiling(n/2))) #Observation point locations Change observation count
m.obs <- length(Obs.site)                     #Number of observations
m <- m.obs
### true forecast distribution
#Will be rehanced in this way.
forecast.mean.true=rep(0,n)+.2
forecast.cov.true=exp(-(distmat/scale)^0.8)
forecast.chol.true=t(chol(forecast.cov.true))

K.true=forecast.cov.true[,Obs.site]%*%
  solve(forecast.cov.true[Obs.site,Obs.site]+sig.eps^2*diag(m.obs))



#t.obs <- (0:60)*10 # Leave 500-600 for forecast.

EnPara.l <- EnPara                              # Parameter for lik guiding
EnPv1d <- runif(N, min = 0.2, max = 0.7)        # Parameter for Normal method

set.seed(1234)
##!!!!!Initial vector matrix

#Enstat <- forecast.mean.true+forecast.chol.true%*%matrix(rnorm(n*N),nrow=n)+3
#Enstat <- matrix(NA, nrow = n, ncol = N)
#for(i in 1:N){Enstat[,i] <- Cini+rep(rnorm(1,0,0.2),n)}

#plot(Enstat[,1])
#for(i in 1:N){lines(Enstat[,i])}


#Currently have N particles of state and parameter



# 04. DA process ----------------------------------------------------------








# 05. True parameters and spatiotemporal field under true parameters ----------------------------------------------------
#Para.true <- c(2,0.5,0.7,0.3,101,301)
#x <- seq(0,100, by=0.25) #Loc
#n <- n;print(n)
#m <- 200 # number of observations

#dx <- x[2]-x[1]
times <- 0:1000
t.obs <- seq(0,times[length(times)], by = 10)





set.seed(1)
#Cini=as.numeric(forecast.mean.true+forecast.chol.true%*%rnorm(n))+3
sx <- seq(0,1,by=0.0025)
Cini <- x0
plot(Cini)

#Real evolution for observations. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

para <- c(1, 20, 1, 3)
#para <- c(0,2)
out <- ode.1D(Cini, times*delta_t, L96_RJ, para, method = "ode45", 
              names = c("C"), dimens = n)
image(out,xlab = 'times',ylab = 'positions')
dim(out)
plot(out[601,2:41],type='l');lines(Cini,col='blue')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Modified until here on 11/25 evening

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
#Since it's like this, directly use multiple runs as real data (not good, too different)
#out[,1+(1:n)] <- out2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Initial setup for comparison~~~~~~~~~~~~~~~~~~~#
Enstat <- matrix(NA, nrow = n, ncol = N)                           #
for(i in 1:N){Enstat[,i] <- Cini*(1+rep(rnorm(1,0,0.1),n))}          #
Enstat_1S  <- Enstat # Classical MCMC   1 breakpoint                            #
Enstat_2S  <- Enstat # Classical MCMC Same as real data 2 breakpoints, 3 steps     #
Enstat_3S  <- Enstat # Classical MCMC   3 breakpoints, 4 steps                 # 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


plot(Enstat[,1])
for(i in 1:N){lines(Enstat[,i])}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Parameter sets for comparison~~~~~~~~~~~~~~~~~~~#
#

# Classical MCMC   1 breakpoint                                                 #
EnPara_1S <- matrix(NA, ncol = N, nrow = 2*kmax+2)
EnPara_1S[1,] <- rep(1,N) 
set.seed(111)
for(i in 1:N){
  kk <- EnPara_1S[1,i]
  if(kk>0){
    c.x <- sample.locs(kk,L)
    EnPara_1S[1+(1:(kk)),i] <-   c.x}
  c.v <- rgamma(kk+1,shape = 2,rate = 1)
  EnPara_1S[kk+1+(1:(kk+1)),i] <- c.v  
}

# Classical MCMC Same as real data 2 breakpoints, 3 steps                          #
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


# Classical MCMC   0 breakpoints, 4 steps                                      # 
#
# Classical MCMC   1 breakpoint                                                 #
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



# 06.Resampling step ----------------------------------------------------------------

# 6.0 First define proposal functions.




#prpd.index <- Prop.Index(k,c=0.3,lambda=lambda_p,kmax=kmax)$index
#Prp.Para(para_0,prpd.index)

#Now have three matrices: state vector from previous time step, state vector from current time step, and parameter matrix
set.seed(1)





En.s.t_1 <- Enstat # Represents state vector from previous time step x_(t-1)
En.s.t <- matrix(NA,nrow = n, ncol = N)   # Represents state vector at new time step x_t
En.lik <- rep(NA,N)

#Def for Bootstrap Resampling.
Mid.En.s.t_1 <- En.s.t_1
Mid.En.s.t   <- En.s.t
Mid.En.lik   <- En.lik
Mid.EnPara   <- EnPara

# 1.0 The main loop. -------------------------------------------------------

#EnPara_2S <- EnPara
# For easier index, the dim of observations is n, not m.
#Enstat=Enstat_0S;EnPara=EnPara_0S;Mthd=0
DA_Process_Lorenz96 <- function(Enstat,EnPara,Mthd){
  En.s.t_1 <- Enstat # Represents state vector from previous time step x_(t-1)
  En.s.t <- matrix(NA,nrow = n, ncol = N)   # Represents state vector at new time step x_t
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
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #Loop_1:       Propagate          #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    
    
    # In this step, propagate the state vector (En.s.t_1) to the new time (En.s.t)
    for (i.n in 1:N){
      k_0 <- EnPara[1,i.n]
      para_0 <- EnPara[1:(2*k_0+2),i.n]; print(para_0)
      # if(k_0 == 0){
      #   p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
      # }else{
      #   p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
      # }
      p4evo0 <- para_0 
      # This is the original: p4evo0 <- c(k_0,0,para_0[2*(1:k_0)+1], n, para_0[2*(0:k_0)+2],0.0)
      sta.ini <- En.s.t_1[,i.n]
      sta.fore0 <- ode.1D(sta.ini, (t1:t2)*delta_t, L96_RJ, p4evo0, method = "ode45", 
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
        # if(k_0 == 0){
        #   p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
        # }else{
        #   p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
        # }
        #Next, propose new parameters
        if (Mthd == 1){
          prpd.index <- Prop.Index(k_0,c=0.3,lambda=lambda_p,kmax=kmax)$index
        }else{
          prpd.index <- sample(c(3,4),1)  
        }
        if(Mthd == 0|k_0==0){prpd.index <- 3}
        propose1 <- Prp.Para(para_0,prpd.index) 
        para_1 <- propose1$para; print(para_1) # Now need to define the function, propose a new function
        
        k_1 <- para_1[1]
        if(para_1[1]!=0&para_1[2]==0){para_1[1+(1:(k_1))]<-para_1[1+(1:(k_1))]+1}
        # if(k_1 == 0){
        #   p4evo1 <- c(k_1,0,n,para_1[1+k_1+(1:(k_1+1))],0.0)
        # }else{
        #   p4evo1 <- c(k_1,0,para_1[1+(1:(k_1))], n, para_1[1+k_1+(1:(k_1+1))],0.0)
        # }
        p4evo1 <- para_1
        sta.ini <- En.s.t_1[,i.n]
        sta.fore1 <- ode.1D(sta.ini, (t1:t2)*delta_t, L96_RJ, p4evo1, method = "ode45", 
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
  }
  if(plot_t){lines(out[t2+1,-1],lwd = 3)
    points(Obs.site,sta.obs[Obs.site])} 
  return(list('EnS'=En.s.t_1,'EnP'=EnPara))
}

Res1_Rj <- DA_Process_Lorenz96(Enstat,EnPara,1)

EnPara    <- Res1_Rj$EnP
En.s.t_1  <- Res1_Rj$EnS



# 2.0 The compared method -------------------------------------------------

En.s.t_1_2S <- Enstat_2S # Represents state vector from previous time step x_(t-1)
En.s.t_2S <- matrix(NA,nrow = n, ncol = N)   # Represents state vector at new time step x_t
En.lik_2S <- rep(NA,N)

#Def for Bootstrap Resampling.
Mid.En.s.t_1_2S <- En.s.t_1_2S
Mid.En.s.t_2S   <- En.s.t_2S
Mid.En.lik_2S   <- En.lik_2S
Mid.EnPara_2S   <- EnPara_2S

Res1_2S <- DA_Process_Lorenz96(Enstat_2S,EnPara_2S,0)
En.s.t_1_2S <- Res1_2S$EnS
EnPara_2S <- Res1_2S$EnP

mean_rj  <- rowMeans(En.s.t_1)
mean_2S  <- rowMeans(En.s.t_1_2S)
m.true <- out[length(times),1+(1:n)]
#Res_ALL[i.all,1] <- mean((mean_rj-m.true)**2)
#Res_ALL[i.all,3] <- mean((mean_2S-m.true)**2)



# 3. The Loop -------------------------------------------------------------
En.s <- En.s.t_1
En.p <- EnPara
s.true <- out[length(times),1+(1:n)]
p.true <- para
times.p <- times[length(times)]+c(0:100)

Pre_En_Lorenz96 <- function(En.s,En.p,s.true,p.true,times.p){
  out.pre <- ode.1D(s.true, times.p*delta_t, L96_RJ, p.true, method = "ode45", 
                    names = c("C"), dimens = n)[length(times.p),1+(1:n)]
  N <- dim(En.s)[2]
  En.predict <- matrix(NA,ncol = N, nrow <- n)
  for(i.n in  1:N){
    k_0 <- En.p[1,i.n]
    para_0 <- En.p[1:(2*k_0+2),i.n]; print(para_0)
    # if(k_0 == 0){
    #   p4evo0 <- c(k_0,0,n,para_0[1+k_0+(1:(k_0+1))],0.0)
    # }else{
    #   p4evo0 <- c(k_0,0,para_0[1+(1:(k_0))], n, para_0[1+k_0+(1:(k_0+1))],0.0)
    # }
    p4evo0 <- para_0
    En.predict[,i.n] <- ode.1D(En.s[,i.n], times.p*delta_t, L96_RJ, 
                               p4evo0, method = "ode45", 
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
#For initial conditions (true conditions)
#~~~~ 3.1 ~~~~ The same setting with the previous part.



#Pre work
### true forecast distribution
locs=1:n
distmat=rdist(locs,locs)  #
scale=20
sig.eps=.1
kappa=3



distmat=rdist(1:n,1:n); scale=100; sig.eps=.1 ;kappa=3; df.t=2

#Obs.site <- sort(sample(c(1:L),ceiling(n/10))) #Observation point locations Change observation count
m.obs <- length(Obs.site)                     #Number of observations
m <- m.obs
### true forecast distribution
#Will be rehanced in this way.
forecast.mean.true=rep(0,n)+.2
forecast.cov.true=exp(-(distmat/scale)^0.8)
forecast.chol.true=t(chol(forecast.cov.true))

K.true=forecast.cov.true[,Obs.site]%*%
  solve(forecast.cov.true[Obs.site,Obs.site]+sig.eps^2*diag(m.obs))



#t.obs <- (0:60)*10 # Leave 500-600 for forecast.

EnPara.l <- EnPara                              # Parameter for lik guiding
EnPv1d <- runif(N, min = 0.2, max = 0.7)        # Parameter for Normal method


##!!!!!Initial vector matrix

#Enstat <- forecast.mean.true+forecast.chol.true%*%matrix(rnorm(n*N),nrow=n)+3
#Enstat <- matrix(NA, nrow = n, ncol = N)
#for(i in 1:N){Enstat[,i] <- Cini+rep(rnorm(1,0,0.2),n)}

#plot(Enstat[,1])
#for(i in 1:N){lines(Enstat[,i])}


#Currently have N particles of state and parameter



# 04. DA process ----------------------------------------------------------








# 05. True parameters and spatiotemporal field under true parameters ----------------------------------------------------
#Para.true <- c(2,0.5,0.7,0.3,101,301)
#x <- seq(0,100, by=0.25) #Loc
#n <- n;print(n)
#m <- 200 # number of observations

#dx <- x[2]-x[1]
times <- 0:1000
t.obs <- seq(0,times[length(times)], by = 10)




#Cini=as.numeric(forecast.mean.true+forecast.chol.true%*%rnorm(n))+3
sx <- seq(0,1,by=0.0025)
Cini <- x0
plot(Cini)

#Real evolution for observations. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

para <- c(1, 20, 1, 3)
out <- ode.1D(Cini, times*delta_t, L96_RJ, para, method = "ode45", 
              names = c("C"), dimens = n)
#image(out,xlab = 'times',ylab = 'positions')
dim(out)
plot(out[601,2:41],type='l');lines(Cini,col='blue')

# First select N=40, 60(Server), 80 100
N_Loop <-20
Res_ALL <- matrix(NA,ncol = 6, nrow = N_Loop)
Res_mean_Para <- matrix(NA, ncol = 2+4+2+4, nrow = N_Loop)
Res_Para <- list()
# Single run, fixed sample size, fixed observation points, N_Loop repetitions
i.all <- 5
for(i.all in 1:N_Loop){
  # Sample from prior distribution
  set.seed(i.all+20241126)
  EnPara <- matrix(NA, ncol = N, nrow = 2*kmax+2)
  EnPara[1,] <- sample(c(0:kmax),N,replace=TRUE,dpois(c(0:kmax), lambda = lambda_p)) 
  
  for(i in 1:N){
    kk <- EnPara[1,i]
    if(kk>0){
      c.x <- sample.locs(kk,L)
      EnPara[1+(1:(kk)),i] <-   c.x}
    # Redefine gamma distribution parameters
    c.v <- rgamma(kk+1,shape=g_alpha,rate=g_beta)
    EnPara[kk+1+(1:(kk+1)),i] <- c.v  
  }
  
  # Prior for state vector
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~Add for Part 2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #          Other ensemble of parameters.
  
  
  # 03.Initial conditions -----------------------------------------------------------------
  
  
  
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Modified until here on 11/25 evening
  
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
  # Since it's like this, directly use multiple runs as real data (not good, too different)
  #out[,1+(1:n)] <- out2
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Initial setup for comparison~~~~~~~~~~~~~~~~~~~#
  Enstat <- matrix(NA, nrow = n, ncol = N)                           #
  for(i in 1:N){Enstat[,i] <- Cini*(1+rep(rnorm(1,0,0.1),n))}          #
  Enstat_1S  <- Enstat # Classical MCMC   1 breakpoint                            #
  Enstat_0S  <- Enstat # Classical MCMC Same as real data 0 breakpoints,      #
  #Enstat_3S  <- Enstat # Classical MCMC   3 breakpoints, 4 steps                 # 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  
  plot(Enstat[,1])
  for(i in 1:N){lines(Enstat[,i])}
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Parameter sets for comparison~~~~~~~~~~~~~~~~~~~#
  #
  
  
  
  # Classical MCMC Same as real data 2 breakpoints, 3 steps                          #
  #
  EnPara_0S <- matrix(NA, ncol = N, nrow = 2*kmax+2)
  EnPara_0S[1,] <- rep(0,N) 
  EnPara_1S <- matrix(NA, ncol = N, nrow = 2*kmax+2)
  EnPara_1S[1,] <- rep(1,N)
  
  for(i in 1:N){
    kk <- EnPara_0S[1,i]
    if(kk>0){
      c.x <- sample.locs(kk,L)
      EnPara_0S[1+(1:(kk)),i] <-   c.x}
    c.v <- rgamma(kk+1,shape = 2,rate = 1)
    EnPara_0S[kk+1+(1:(kk+1)),i] <- c.v  
  }
  
  for(i in 1:N){
    kk <- EnPara_1S[1,i]
    if(kk>0){
      c.x <- sample.locs(kk,L)
      EnPara_1S[1+(1:(kk)),i] <-   c.x}
    c.v <- rgamma(kk+1,shape = 2,rate = 1)
    EnPara_1S[kk+1+(1:(kk+1)),i] <- c.v  
  }
  
  
  
  
  
  
  
  
  
  
  En.s.t_1 <- Enstat # Represents state vector from previous time step x_(t-1)
  En.s.t <- matrix(NA,nrow = n, ncol = N)   # Represents state vector at new time step x_t
  En.lik <- rep(NA,N)
  
  #Def for Bootstrap Resampling.
  Mid.En.s.t_1 <- En.s.t_1
  Mid.En.s.t   <- En.s.t
  Mid.En.lik   <- En.lik
  Mid.EnPara   <- EnPara
  
  # 1.0 The main loop. -------------------------------------------------------
  
  #EnPara_2S <- EnPara
  # For easier index, the dim of observations is n, not m.
  
  Res1_Rj <- DA_Process_Lorenz96(Enstat,EnPara,1)
  
  EnPara    <- Res1_Rj$EnP
  En.s.t_1  <- Res1_Rj$EnS
  
  
  
  # 2.0 The compared method -------------------------------------------------
  
  #K=0
  En.s.t_1_0S <- Enstat_0S # Represents state vector from previous time step x_(t-1)
  En.s.t_0S <- matrix(NA,nrow = n, ncol = N)   # Represents state vector at new time step x_t
  En.lik_0S <- rep(NA,N)
  
  #Def for Bootstrap Resampling.
  Mid.En.s.t_1_0S <- En.s.t_1_0S
  Mid.En.s.t_0S   <- En.s.t_0S
  Mid.En.lik_0S   <- En.lik_0S
  Mid.EnPara_0S   <- EnPara_0S
  
  Res1_0S <- DA_Process_Lorenz96(Enstat_0S,EnPara_0S,0)
  En.s.t_1_0S <- Res1_0S$EnS
  EnPara_0S <- Res1_0S$EnP
  
  
  #K=1
  En.s.t_1_1S <- Enstat_1S # Represents state vector from previous time step x_(t-1)
  En.s.t_1S <- matrix(NA,nrow = n, ncol = N)   # Represents state vector at new time step x_t
  En.lik_1S <- rep(NA,N)
  
  #Def for Bootstrap Resampling.
  Mid.En.s.t_1_1S <- En.s.t_1_1S
  Mid.En.s.t_1S   <- En.s.t_1S
  Mid.En.lik_1S   <- En.lik_1S
  Mid.EnPara_1S   <- EnPara_1S
  
  Res1_1S <- DA_Process_Lorenz96(Enstat_1S,EnPara_1S,0)
  En.s.t_1_1S <- Res1_1S$EnS
  EnPara_1S <- Res1_1S$EnP
  
  mean_rj  <- rowMeans(En.s.t_1)
  mean_0S  <- rowMeans(En.s.t_1_0S)
  mean_1S  <- rowMeans(En.s.t_1_1S)
  m.true <- out[length(times),1+(1:n)]
  Res_ALL[i.all,1] <- mean((mean_rj-m.true)**2)
  Res_ALL[i.all,2] <- mean((mean_0S-m.true)**2)
  Res_ALL[i.all,3] <- mean((mean_1S-m.true)**2)
  
  #Predictions
  Res_r <- Pre_En_Lorenz96(En.s.t_1,EnPara,out[length(times),1+(1:n)],para,times.p)
  Res_0 <- Pre_En_Lorenz96(En.s.t_1_0S,EnPara_0S,out[length(times),1+(1:n)],para,times.p)
  Res_1 <- Pre_En_Lorenz96(En.s.t_1_1S,EnPara_1S,out[length(times),1+(1:n)],para,times.p)
  
  Res_ALL[i.all,4] <- Res_r$MSPE
  Res_ALL[i.all,5] <- Res_0$MSPE
  Res_ALL[i.all,6] <- Res_1$MSPE
  
  Res_Para[[i.all]] <- rbind(EnPara,EnPara_0S,EnPara_1S)
  # Generate a long matrix result
  # Respectively for RJ k=0 and k=1
  if(sum(EnPara[1,]==0)==0){Res_mean_Para[i.all, 1:2] <- rep(NA,2)
  }else{
    Res_mean_Para[i.all, 1:2] <- rowMeans(EnPara[1:2,which(EnPara[1,]==0)])
  }
  if(sum(EnPara[1,]==1)==0){Res_mean_Para[i.all, 2+1:4] <- rep(NA,4)
  }else{
    Res_mean_Para[i.all, 2+1:4] <- rowMeans(EnPara[1:4,which(EnPara[1,]==1)])
  }
  Res_mean_Para[i.all, 6+1:2] <- rowMeans(EnPara_0S[1:2,])
  Res_mean_Para[i.all, 8+1:4] <- rowMeans(EnPara_1S[1:4,])
}

colMeans(Res_ALL)
#save(Res_ALL,  file = "E:/Study-2024/09.RJ_2stage_resample_Oct24/1.Code_Oct24/L96_N100m40_1.RData")


data <- data.frame(
  Group = factor(rep(c("MSEs- New Method", "MSEs- pMCMC_0","MSEs- pMCMC_1",
                       "MSPEs- New Method", "MSPEs- pMCMC_0", "MSPEs- pMCMC_1"), each = 10)),
  Values = c(Res_ALL[1:10,1],Res_ALL[1:10,2],Res_ALL[1:10,3],
             Res_ALL[1:10,4],Res_ALL[1:10,5],Res_ALL[1:10,6])
)

# Create boxplot
ggplot(data, aes(x = Group, y = Values, fill = Group)) +
  geom_boxplot() +
  labs( 
    x = "Methods", 
    y = "MSEs/MSPEs", 
    fill = "Methods") +
  theme_minimal(base_size = 16) +  # Set base font size to 16, equivalent to 2x
  theme(
    axis.title = element_text(size = 20),  # Font size for y and x axis titles
    axis.text = element_text(size = 16),   # Font size for axis ticks
    legend.title = element_text(size = 20),  # Legend title font size
    legend.text = element_text(size = 16)    # Legend text font size
  )



# Below want to write parameter identification comparison --------------------------------------------------------------
# Directly add a list above
#L96-m=40



R_n040 <- c(0.02003933, 0.01837651, 0.02019963,
            1.26148704, 1.57694994, 0.55160206)
R_n070_2<-c(0.01983714, 0.01997749, 0.01648219,
            1.06771461, 1.58412866, 0.50036950)
R_n040 <- (R_n040*20+R_n040_2*17)/37
R_n070 <- c(0.01872140, 0.01809806, 0.02179405,
            0.56160288, 1.58547240, 0.65656756)#20
R_n070_2<-c(0.01769418, 0.02362979, 0.01712082,
            0.50114131, 1.58326663, 0.48713679)
R_n070 <- (R_n070*20+R_n070_2*10)/30
R_n100 <- c(0.01716069, 0.01934158, 0.01920814,
            0.55239000, 1.58790338, 0.62415856)#20
#c(0.01935273 0.01878655 0.01640051 0.10184189 1.62414479 0.54760980)



# Parameters F=1,6 simulation times 35/40
#load("E:/Study-2024/09.RJ_2stage_resample_Oct24/1.Code_Oct24/Results_Lorenz/nr_L96_N100m40.RData")

R_n040 <- c(0.01915066, 0.01978335, 0.01793502,
            0.04425217, 1.59046655, 0.72284056)
R_n070 <- c(0.01704360, 0.02181829, 0.01955970,
            0.23614419, 1.58672480, 0.60479998)
R_n100 <- c(0.01651126, 0.01974257, 0.01972677,
            0.04441996, 1.59339060, 0.55431425)

Res_MSEs <- data.frame(R_n040,R_n070,R_n100)

line_types <- c("solid","solid","solid",
                "dashed","dashed","dashed")
colors <- c('blue','green','red',
            'blue','green','red')
par(mar=c(5,5,5,2))
Vec.N <- c(40,70,100)
plot(Vec.N, Res_MSEs[4,], type = "l", 
     col = colors[4], lty = line_types[4],
     #ylim = range(Res_MSEs), xlab = "Index", ylab = "Values",
     ylim = c(0,2),xlab = "Size of the Ensemble", ylab = "MSPEs",
     cex.axis = 2,cex.lab = 2,cex.main=2,
     main = "MSPEs at t = 1100",lwd=3)
for (i in 5:6) {
  lines(Vec.N, Res_MSEs[i,], col = colors[i], 
        lty = line_types[i],lwd=3)}
legend("topright", legend = c('New Method',paste("pMCMC, k=", 0:1)),
       col = colors, lty = line_types, title = "MSPEs",
       cex = 1,           # Enlarge legend
       lwd = 3,             # Line width
       seg.len = 0.6,       # Line length
       y.intersp = 0.5,     # Increase vertical spacing between entries
       x.intersp = 0.5,     # Adjust horizontal distance between symbols and text
       text.width =8)

plot(Vec.N, Res_MSEs[1,], type = "l", 
     col = colors[1], lty = line_types[1],
     #ylim = range(Res_MSEs), xlab = "Index", ylab = "Values",
     ylim = c(0.01,0.03),xlab = "Size of the Ensemble", ylab = "MSEs",
     cex.axis = 2,cex.lab = 2,cex.main=2,
     main = "MSEs at t = 1000",lwd=3)
for (i in 2:3) {
  lines(Vec.N, Res_MSEs[i,], col = colors[i], 
        lty = line_types[i],lwd=3)}
legend("topright", legend = c('New Method',paste("pMCMC, k=", 0:1)),
       col = colors, lty = line_types, title = "MSEs",
       cex = 1,           # Enlarge legend
       lwd = 3,             # Line width
       seg.len = 0.6,       # Line length
       y.intersp = 0.5,     # Increase vertical spacing between entries
       x.intersp = 0.5,     # Adjust horizontal distance between symbols and text
       text.width =8)