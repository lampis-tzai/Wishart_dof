library(CholWishart)
library(abind)
library(rstan)
library(readxl)
library(writexl)
library(coda)
library(parallel)
library(R.utils)

## Maximum Likelihood approach

sum_det_list <- function(cov_list){
  sum_dlist <- sum(unlist(lapply(cov_list, function(x){ determinant(x,logarithm = TRUE)$modulus[1]})))
  return(sum_dlist)
}

log_v_estimation <- function(dof,M,p,cov_list){
  sum_list <- Reduce('+', cov_list)
  
  v <- sum_list/(dof*M)
  return(determinant(v,logarithm = TRUE)$modulus[1])
}


v_estimation <- function(dof,M,p,cov_list){
  sum_list <- 0
  for (mx in cov_list){
    sum_list <- sum_list + mx
  }
  
  v <-  sum_list/(dof*M)
  return(v)
}


degrees_of_freedom <- function(df,cov_data, p){
  k <- 1:p
  M <- length(cov_data)
  result <- sum(digamma((df-k+1)/2)) - (sum_det_list(cov_data) - M*(p*log(2) + log_v_estimation(df,M,p,cov_data)))/M
  return(result)
}

#### Bisection Method

bisection_method <- function(f, L, U, num = 1000, tol = 1e-5, ...) {
  # If the signs of the function at the evaluated points, a and b, stop the function and return message.
  if (L>U){
    stop('L must be smaller than U')
  } else if (!(f(L,...) < 0) && (f(U,...) > 0)) {
    #stop('signs of f(a) and f(b) must differ')
    return(L)
  } else if ((f(L, ...) > 0) && (f(U, ...) < 0)) {
    #stop('signs of f(a) and f(b) must differ')
    return(L)
  }
  hist <- c()
  for (i in 1:num) {
    C <- (L + U) / 2 # Calculate midpoint
    hist <- c(hist,C)
    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the 
    # function and return the root.
    if ((f(C, ...) == 0) || ((U - L) / 2) < tol) {
      return(hist)
    }
    
    # If another iteration is required, 
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    ifelse(sign(f(C,...)) == sign(f(L,...)), 
           L <- C,
           U <- C)
  }
  # If the max number of iterations is reached and no root has been found, 
  # return message and end function.
  print('Bisection Method failed')
  return(hist)
}


#### Newton-Raphson


second_derivative <-function(df,p,m){
  k<-1:p
  sec_der <- m*sum(trigamma((df-k+1)/2))
  return(sec_der)
}


NR_wishart<-function(model_data,p,initial ,tol = 1e-5) {
  
  m <- length(model_data)
  df <- p
  histn <- c(initial)
  t<-0
  while (t==0) {
    df_new <-df - degrees_of_freedom(df,model_data, p)/second_derivative(df,p,m)
    histn<-c(histn,df_new)
    if (abs(df-df_new)<tol) t<-10
    df<-df_new
  }
  return(histn)
}

## Bayesian approach


n_pdf <- function(model_data,V,p,dof,a,b,prior = 'gamma'){
  
  if (dof<=p-1){return(-Inf)}
  
  K=length(model_data)
  mult_x = sum(unlist(lapply(model_data, function(x){ determinant(x,logarithm = TRUE)$modulus[1]})))
  
  if (prior == "uniform"){
    n_prior = log(ifelse((dof>0 & dof<10000),1,0))
  } else if (prior == "exponential") {
    n_prior = log(dof) -b*dof
  }else if (prior == "gamma") {
    n_prior = (a-1)*log(dof) -b*dof
  }else if (prior == "inverse_gamma") {
    n_prior = (-a-1)*log(dof) -b/dof
  }else if (prior == "log_normal"){
    n_prior = -log(dof) - log(dof)/b + (log(dof)*a)/b
  }
  
  fx = (dof/2)*mult_x - K*(((dof*p)/2)*log(2) + (dof/2)*determinant(V,logarithm = TRUE)$modulus[1] +
                             lmvgamma(dof/2,p)) + n_prior
  
  return(fx)
}

#Random walk metropolis
RWM_for_dof <- function(model_data,V,p,dof,a,b,delta=1, prior = 'gamma'){
  
  old_n_pdf = n_pdf(dof=dof, model_data = model_data,V = V,p=p,a = a,b=b,prior = prior)
  
  y<- runif(1, dof-delta,dof+delta)
  accept<- exp(n_pdf(dof=y, model_data = model_data,V = V,p=p,a = a,b=b, prior = prior) - old_n_pdf)
  accept<- min(c(accept,1))
  u<-runif(1,0,1)
  if (u<accept){
    z <- y
  }else {
    z <- dof
  }
  
  
  return(z)
}

#Slice sampling metropolis
slice_sample <- function(model_data,V,p,dof,a,b,w=1, prior = 'gamma') {
  z = 0
  while (z<p){ 
    # Draw a vertical level (y) from the exponential distribution (equivalent to log(y) from uniform)
    y_level <- n_pdf(dof=dof, model_data = model_data,V = V,p=p,a = a,b=b, prior = prior) - rexp(1)
    
    # Find an interval [x_l, x_r] around the current sample
    x_l <- dof - runif(1) * w
    x_r <- x_l + w
    
    # Step out the interval until the log PDF is below the level
    while (n_pdf(dof = x_l, model_data = model_data,V = V,p=p,a = a,b=b, prior = prior) > y_level) x_l <- x_l - w
    while (n_pdf(dof = x_r, model_data = model_data,V = V,p=p,a = a,b=b, prior = prior) > y_level) x_r <- x_r + w
    
    # Propose a new sample within the interval and accept it if it's above the level
    repeat {
      x_proposed <- runif(1, x_l, x_r)
      if (n_pdf(dof = x_proposed, model_data = model_data,V = V,p=p,a = a,b=b, prior = prior) > y_level) {
        break
      } else {
        # Shrink the interval
        if (x_proposed < dof) {
          x_l <- x_proposed
        } else {
          x_r <- x_proposed
        }
      }
    }
    
    z<-x_proposed
  }
  
  return(z)
}

#Hameltonian monte carlo stan

n_pdf_gradient <- function(model_data,V,p,dof,a,b, prior){
  
  if (prior == "uniform"){
    n_prior <- log(ifelse((dof>0 & dof<10000),1,0))
  } else if (prior == "exponential") {
    n_prior <- 1/dof -b
  }else if (prior == "gamma") {
    n_prior <- (a-1)*(1/dof) -b
  }else if (prior == "inverse_gamma") {
    n_prior <- (-a-1)*(1/dof) -b
  }else if (prior == "log_normal"){
    n_prior <- -(1/dof) - (1/dof)*(1/b) + (1/dof)*(a/b)
  }
  
  kp=1:p
  K=length(model_data)
  mult_x = sum(unlist(lapply(model_data, function(x){ determinant(x,logarithm = TRUE)$modulus[1]})))
  
  fx = mult_x/2 - (K/2)*( p*log(2) + determinant(V,logarithm = TRUE)$modulus[1] +
                            sum(digamma((dof-kp+1)/2))) + n_prior
  
  return(fx)
}


HMC = function (model_data,V,p,dof,a,b, prior, epsilon=0.1, L=10){
  
  current_U <- -n_pdf(dof = dof, model_data = model_data,V = V,p=p,a = a,b=b, prior = prior)
  current_grad_U <- -n_pdf_gradient(dof = dof, model_data = model_data,V = V,p=p,a = a,b=b, prior = prior)
  
  momentum <- rnorm(length(dof)) # independent standard normal variates
  new_dof <- dof
  new_momentum <- momentum
  # Alternate full steps for position and momentum
  for (i in 1:L) {
    new_momentum <- new_momentum - epsilon * current_grad_U / 2
    new_dof <- new_dof + epsilon * new_momentum
    current_grad_U <- -n_pdf_gradient(dof = new_dof, model_data = model_data,V = V,p=p,a = a,b=b, prior = prior)
    new_momentum <- new_momentum - epsilon * current_grad_U / 2
  }
  
  # Proposed potential
  proposed_U <- -n_pdf(dof = new_dof, model_data = model_data,V = V,p=p,a = a,b=b, prior = prior)
  # Metropolis acceptance step
  if (runif(1) < exp(current_U - proposed_U + sum(momentum^2)/2 - sum(new_momentum^2)/2)){
    z = new_dof
  } else {
    z = dof
  }
  
  
  return (z)
}


metropolis_within_gibbs <- function(model_data,n_mcmc = 3000,
                                    burn_in = 1000, metropolis_alg = "RwM" , 
                                    prior = "gamma", def_time = 300){
  sum_xi = 0
  for (md in model_data){
    sum_xi = sum_xi + md
  }
  
  M = length(model_data)
  p = dim(model_data[[1]])[1]
  
  m = p
  U = diag(0.0001,p)
  
  if (prior == "log_normal"){
    b = 5 # var of log normal
    a = 0 # mean of log normal
  } else {
    b = 0.0001 #rate
    a = 0.0001 #shape
  }
  
  n_post = c(p)
  V_post = list()
  V_post[[1]] = diag(1,p)
  
  start.time <- Sys.time()
  time.taken <- 0
  
  if (metropolis_alg == "RwM"){
    for (i in 2:n_mcmc){
      if (time.taken<def_time){
        V_post[[i]] = rInvWishart(1,((n_post[i-1]*M) + m),(sum_xi + U))[,,1]
        n_post = c(n_post,RWM_for_dof(model_data=model_data,V=V_post[[i]], p=p, 
                                      dof=n_post[i-1], a=a, b=b, prior = prior))
      }
      end.time <- Sys.time()
      time.taken <- as.numeric(difftime(end.time,  start.time, units='secs'))
    }
  } else if (metropolis_alg == "Slice"){
    for (i in 2:n_mcmc){
      if (time.taken<def_time){
        V_post[[i]] = rInvWishart(1,((n_post[i-1]*M) + m),(sum_xi + U))[,,1]
        n_post = c(n_post,slice_sample(model_data=model_data,V=V_post[[i]], p=p, 
                                       dof=n_post[i-1], a=a, b=b, prior = prior))
      }
      end.time <- Sys.time()
      time.taken <- as.numeric(difftime(end.time,  start.time, units='secs'))
    }
  } else if (metropolis_alg == "HMC"){
    for (i in 2:n_mcmc){
      if (time.taken<def_time){
        V_post[[i]] = rInvWishart(1,((n_post[i-1]*M) + m),(sum_xi + U))[,,1]
        n_post = c(n_post,HMC(model_data=model_data,V=V_post[[i]], p=p, 
                              dof=n_post[i-1], a=a, b=b, prior = prior))
      }
      end.time <- Sys.time()
      time.taken <- as.numeric(difftime(end.time,  start.time, units='secs'))
    }
  }
  
  iter<- length(n_post)
  burn_in <- ceiling(iter/2)
  return(list(n_post[seq(burn_in,iter)],V_post[burn_in:iter]))
}


stan_model_lognormal <- stan_model(file = "Stan_models/stan_wishart_inv_wish_lognormal.stan", model_name = "whole_stan_model_lognormal")

all_modelling_nuts <- function(model_data,n_mcmc = 3000,
                               burn_in = 1000, prior = "gamma", def_time = 300){
  
  p = dim(model_data[[1]])[1]
  nv = p
  U = diag(0.001,p)
  M = length(model_data)
  
  if (prior == "log_normal"){
    b = 5 # var of log normal
    a = 0 # mean of log normal
  } else {
    b = 0.0001 #rate
    a = 0.0001 #shape
  }
  
  init_function <- function() {
    list(V = diag(1,p), n=p+1)
  }
  
  # Number of chains
  n_chains <- 1
  
  # Generate a list of initial values for each chain
  init_values <- lapply(1:n_chains, function(id) init_function())
  
  
  stan_data = list(
    "P" = p,
    "M" = M,
    "y" = model_data,
    "U"= U,
    "nv" = nv,
    "a"= as.numeric(a),
    "b"= as.numeric(b)
  )
  
  file_title = paste0("rstan_samples/",as.numeric(Sys.time()),"_obs = ",M,"_dim = ", p,"_samples.csv")
  
  if (prior == "uniform"){
    
    fit <- withTimeout(sampling(stan_model_unif, data = stan_data, iter = 3000,
                                chains = n_chains, init = init_values,refresh=0,
                                sample_file = file_title, save_warmup = TRUE), 
                       timeout=def_time)
    
  } else if (prior == "exponential") {
    
    fit <- withTimeout(sampling(stan_model_exp, data = stan_data, iter = 3000, 
                                chains = n_chains, init = init_values,refresh=0,
                                sample_file = file_title, save_warmup = TRUE), 
                       timeout=def_time)
  }else if (prior == "gamma") {
    fit <- withTimeout(sampling(stan_model_gamma, data = stan_data, iter = 3000, 
                                chains = n_chains, init = init_values,refresh=0,
                                sample_file = file_title, save_warmup = TRUE), 
                       timeout=def_time)
    
    
  }else if (prior == "inverse_gamma") {
    fit <- withTimeout(sampling(stan_model_invgamma, data = stan_data, iter = 3000, 
                                chains = n_chains, init = init_values,refresh=0,
                                sample_file = file_title, save_warmup = TRUE), 
                       timeout=def_time)
    
  }else if (prior == "log_normal"){
    fit <- withTimeout(sampling(stan_model_lognormal, data = stan_data, iter = 3000, 
                                chains = n_chains, init = init_values,refresh=0,
                                sample_file = file_title, save_warmup = TRUE), 
                       timeout=def_time)
  }
  
  
  post_samples <- read.csv(file_title, comment.char = "#")
  
  n_post = post_samples[,8]
  iter<- length(n_post)
  V_post = array(as.matrix(post_samples[,9:ncol(post_samples)]),dim = c(iter,p,p))
  
  burn_in <- ceiling(iter/2)
  
  #unlink(file_title, recursive = TRUE)
  
  #because the file is saved with one digit less than the extract function (some discrepencies)
  
  return(list(n_post[burn_in:iter],V_post[burn_in:iter,,]))
}



all_methods <- function(model_data,d, return_samples = ''){
  
  obs = length(model_data)
  
  df_bisection_hist <- bisection_method(degrees_of_freedom, L = d, U = 1000,num = 1000, tol = 1e-05,model_data, d)
  df_bisection <- df_bisection_hist[length(df_bisection_hist)]
  
  NR_hist <- NR_wishart(model_data =model_data, p=d, initial = df_bisection)
  df_NR <- NR_hist[length(NR_hist)]
  V_est_NR <- v_estimation(df_NR,obs,d,model_data)
  print('Maximum Likelihood')
  
  
  start.time <- Sys.time()
  RWM <- metropolis_within_gibbs(model_data, metropolis_alg = "RwM", prior = "gamma")
  end.time <- Sys.time()
  RWM_time.taken <- as.numeric(difftime(end.time,  start.time, units='mins'))
  
  n_post_RWM <- RWM[[1]]
  V_post_RWM <- RWM[[2]]
  if(return_samples=="RWM"){
    return(list(n_post_RWM,V_post_RWM))
  }
  ess_df_RWM <- effectiveSize(n_post_RWM)[[1]]
  geweke_df_RWM <- abs(geweke.diag(n_post_RWM)[[1]])
  n_post_RWM_est <- median(n_post_RWM)
  V_post_RWM_est <- apply(abind(V_post_RWM, along = 3),c(1,2),mean)
  mcmc_iter_RWM <- length(n_post_RWM)*2 -2 
  mcse_RWM <- sqrt(var(n_post_RWM) / geweke_df_RWM)
  print('RWM')
  
  start.time <- Sys.time()
  SS <- metropolis_within_gibbs(model_data, metropolis_alg = "Slice", prior = "inverse_gamma")
  end.time <- Sys.time()
  SS_time.taken <- as.numeric(difftime(end.time,  start.time, units='mins'))
  
  n_post_SS <- SS[[1]]
  V_post_SS <- SS[[2]]
  
  if(return_samples=="SS"){
    return(list(n_post_SS,V_post_SS))
  }
  
  ess_df_SS <- effectiveSize(n_post_SS)[[1]]
  geweke_df_SS <- abs(geweke.diag(n_post_SS)[[1]])
  n_post_SS_est <- median(n_post_SS)
  V_post_SS_est <- apply(abind(V_post_SS, along = 3),c(1,2),mean)
  mcmc_iter_SS <- length(n_post_SS)*2 -2 
  mcse_SS <- sqrt(var(n_post_SS) / geweke_df_SS)
  print('Slice sampling')
  
  start.time <- Sys.time()
  hmc_gibbs <- metropolis_within_gibbs(model_data, metropolis_alg = "HMC", prior = "log_normal")
  end.time <- Sys.time()
  hmc_gibbs_time.taken <- as.numeric(difftime(end.time,  start.time, units='mins'))
  
  n_post_hmc_gibbs <- hmc_gibbs[[1]]
  V_post_hmc_gibbs <- hmc_gibbs[[2]]
  if(return_samples=="HMC"){
    return(list(n_post_hmc_gibbs,V_post_hmc_gibbs))
  }
  ess_df_hmc_gibbs <- effectiveSize(n_post_hmc_gibbs)[[1]]
  geweke_df_hmc_gibbs <- abs(geweke.diag(n_post_hmc_gibbs)[[1]])
  n_post_hmc_est <- median(n_post_hmc_gibbs)
  V_post_hmc_est <- apply(abind(V_post_RWM, along = 3),c(1,2),mean)
  mcmc_iter_hmc_gibbs <- length(n_post_hmc_gibbs)*2 -2 
  mcse_hmc_gibbs <- sqrt(var(n_post_hmc_gibbs) / geweke_df_hmc_gibbs)
  print('HMC')
  
  start.time <- Sys.time()
  stan_all <- all_modelling_nuts(model_data = model_data, prior = "log_normal")
  end.time <- Sys.time()
  stan_all.taken <-as.numeric(difftime(end.time,  start.time, units='mins'))
  n_post_stan <- stan_all[[1]]
  V_post_stan <- stan_all[[2]]
  ess_df_stan <- effectiveSize(as.vector(n_post_stan))[[1]]
  
  geweke_df_stan <- abs(geweke.diag(as.vector(n_post_stan))[[1]])
  
  n_post_stan_est <- median(n_post_stan)
  V_post_stan_est <- apply(V_post_stan,c(2,3),mean)
  mcmc_iter_stan <- length(n_post_stan)*2 -2 
  mcse_stan <- sqrt(var(n_post_stan) / geweke_df_stan)
  print('stan')
  
  return ( list(df_NR, V_est_NR,
                RWM_time.taken, n_post_RWM_est,V_post_RWM_est,ess_df_RWM,geweke_df_RWM, mcmc_iter_RWM, mcse_RWM,
                SS_time.taken,n_post_SS_est,V_post_SS_est,ess_df_SS,geweke_df_SS, mcmc_iter_SS, mcse_SS,
                hmc_gibbs_time.taken, n_post_hmc_est,V_post_hmc_est,ess_df_hmc_gibbs,geweke_df_hmc_gibbs, mcmc_iter_hmc_gibbs, mcse_hmc_gibbs, 
                stan_all.taken, n_post_stan_est,V_post_stan_est,ess_df_stan,geweke_df_stan, mcmc_iter_stan, mcse_stan
  ))
  
}

set.seed(8)        
#### basketball data

basketball = as.data.frame(read.csv('real_datasets/2022-2023 NBA Player Stats - Regular.csv',
                                    header = TRUE, sep = ";"))

basketball = basketball[,c('Tm','TRB','AST','STL','BLK','TOV','PF','PTS')]

basketball_model_data = list()
i=1
for (team in unique(basketball$Tm)){
  basketball_model_data[[i]] = cov(basketball[basketball$Tm==team,2:ncol(basketball)])
  i = i+1
}

basketball_results = all_methods(basketball_model_data,d=(ncol(basketball)-1))
basketball_results

basketball_V_PE = array(0,dim = c(5,5))
pos = c(2,5,12,19,26)
for(i in 1:5){
  for(j in 1:5){
    upper_sigma_i = basketball_results[[pos[i]]][upper.tri(basketball_results[[pos[i]]], diag = T)]
    upper_sigma_j = basketball_results[[pos[j]]][upper.tri(basketball_results[[pos[j]]], diag = T)]
    basketball_V_PE[i,j] = mean(abs((upper_sigma_i - upper_sigma_j)/upper_sigma_i))
  }
}


#### air quality data


air = as.data.frame(read.csv('real_datasets/AQI and Lat Long of Countries.csv',
                             header = TRUE, sep = ","))

air = air[, (names(air) %in% c('Country','AQI.Value','Ozone.AQI.Value'))]

#alps countries
air_alps = air[air$Country %in% c("Italy" ,"France","Switzerland" ,"Austria", "Germany","Slovenia"),]

air_model_data = list()
i=1
for (my in unique(air_alps$Country)){
  if (sum(air_alps$Country==my)>=2){
    air_model_data[[i]] = cov(air_alps[air_alps$Country==my,2:3])
    i=i+1
  }
}

air_results = all_methods(air_model_data,d=2)
air_results

air_V_PE = array(0,dim = c(5,5))
pos = c(2,5,12,19,26)
for(i in 1:5){
  for(j in 1:5){
    upper_sigma_i = air_results[[pos[i]]][upper.tri(air_results[[pos[i]]], diag = T)]
    upper_sigma_j = air_results[[pos[j]]][upper.tri(air_results[[pos[j]]], diag = T)]
    air_V_PE[i,j] = mean(abs((upper_sigma_i - upper_sigma_j)/upper_sigma_i))
  }
}



#### hand data

hand_data <- as.data.frame(read_excel("real_datasets/handwriting_data_bootstrap.xlsx"))
hand_data['new_id'] = paste0(hand_data$Writer,'_',hand_data$Lettre)

hand_model_data = list()
i=1
for (t in unique(hand_data$new_id)){
  spdf = hand_data[hand_data$new_id ==t,2:21]
  spdf = spdf[complete.cases(spdf),]
  if (nrow(spdf)>2){
    hand_model_data[[i]] = cov(spdf)
    i=i+1
  }
}

hand_results = all_methods(hand_model_data,d=20)
hand_results


hand_V_PE = array(0,dim = c(5,5))
pos = c(2,5,12,19,26)
for(i in 1:5){
  for(j in 1:5){
    upper_sigma_i = hand_results[[pos[i]]][upper.tri(hand_results[[pos[i]]], diag = T)]
    upper_sigma_j = hand_results[[pos[j]]][upper.tri(hand_results[[pos[j]]], diag = T)]
    hand_V_PE[i,j] = mean(abs((upper_sigma_i - upper_sigma_j)/upper_sigma_i))
  }
}

round((basketball_V_PE + air_V_PE + hand_V_PE)/3,3)




#####################PLOTs##########################

library(bayesplot)
library(ggplot2)

file_list <- list.files(path = 'rstan_samples/', pattern = "\\.csv$", full.names = TRUE)
post_samples <- read.csv(file_list[[1]], comment.char = "#")
mcmc_iter = nrow(post_samples)
burn_in = ceiling(mcmc_iter/2)+1
# Create a 3D array from the vector
posterior <- array(post_samples[burn_in:mcmc_iter,8], dim = c(length(post_samples[burn_in:mcmc_iter,8]), 1, 1))
dimnames(posterior)[[3]] <- c("n")



# hand_results = all_methods(hand_model_data,d=20,return_samples = 'HMC')
# mcmc_iter = length(hand_results[[1]])
# burn_in = ceiling(mcmc_iter/2)+1
# # Create a 3D array from the vector
# posterior <- array(hand_results[[1]][burn_in:mcmc_iter], dim = c(length(hand_results[[1]][burn_in:mcmc_iter]), 1, 1))
# dimnames(posterior)[[3]] <- c("n")


acf(posterior[,,'n'])#[,'chain:1'])


color_scheme_set("red")
mcmc_intervals(posterior, pars = c("n"))



color_scheme_set("blue")
dist_n <-mcmc_areas(
  posterior, 
  pars = c("n"),
  prob = 0.95, # 95% intervals
  prob_outer = 0.99, # 99%
  point_est = "median")+
  labs(
    title = "Posterior distributions of the degrees of freedom ",
    subtitle = "for the NBA data-set")+
  theme(
    axis.text = element_text(size = 20),       # Change size of axis labels
    axis.title = element_text(size = 25),      # Change size of axis titles
    plot.title = element_text(size = 30, face = "bold"), # Change size of plot title
    plot.subtitle = element_text(size = 27, face = "bold")
  )



dist_n
jpeg("plots/post_dist_n_nba_time.jpg",width=3920, height=2000, res=300)
print(dist_n)
dev.off()



color_scheme_set("mix-blue-red")
trace_n <- mcmc_trace(posterior, pars = c("n"))+
  labs(
    title = "Trace plot for the degrees of freedom ",
    subtitle = "for the NBA data-set")+
  theme(
    axis.text = element_text(size = 20),       # Change size of axis labels
    axis.title = element_text(size = 25),      # Change size of axis titles
    plot.title = element_text(size = 30, face = "bold"), # Change size of plot title
    plot.subtitle = element_text(size = 27, face = "bold")
  )

trace_n
jpeg("plots/trace_n_nba_time.jpg",width=3920, height=2000, res=300)
print(trace_n)
dev.off()

p = 7
V_post = as.matrix(post_samples[burn_in:mcmc_iter,9:ncol(post_samples)])
V_post = lapply(1:dim(V_post)[1], function(i) array(V_post[i,], dim=c(p,p)))


#V_post = hand_results[[2]]


biggest_eigen <- unlist(lapply(V_post, function(x) eigen(x)$values[1]))
smallest_eigen <- unlist(lapply(V_post, function(x) eigen(x)$values[p]))


library(plot3D)

##  Create cuts:
biggest_eigen_c <- cut(biggest_eigen, 20)
smallest_eigen_c <- cut(smallest_eigen, 20)

##  Calculate joint counts at cut levels:
z <- table(biggest_eigen_c, smallest_eigen_c)

##  Plot as a 3D histogram:
jpeg("plots/3d_hist_v_nba_time.jpg",width=3920, height=2000, res=300)
par(cex.axis = 1, cex.lab = 1, cex.main = 1.5)
plot_V_3d <- hist3D(z=z, border="black", main = "3D histogram of the eigenvalues of the scale matrix V \n for the NBA data-set",
                    xlab = "Largest eigenvalue", ylab = "Smalest eigenvalue",
                    zlab = "Frequency")
dev.off()
