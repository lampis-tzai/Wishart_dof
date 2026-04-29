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
  
  if (dof<0){return(-Inf)}
  
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
RWM_for_dof <- function(model_data,V,p,dof,a,b,delta=3, prior = 'gamma'){
  z = 0
  old_n_pdf = n_pdf(dof=dof, model_data = model_data,V = V,p=p,a = a,b=b,prior = prior)
  while (z<p){ 
    y<- runif(1, dof-delta,dof+delta)
    if (y>=p){
      accept<- exp(n_pdf(dof=y, model_data = model_data,V = V,p=p,a = a,b=b, prior = prior) - old_n_pdf)
      accept<- min(c(accept,1))
      u<-runif(1,0,1)
      if (u<accept)z<-y
    } else {
      return(dof)
    }
    
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
  z = 0
  while (z<p-1){
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
    }
    else {
      return(dof)
    }
  }
  return (z)
}


metropolis_within_gibbs <- function(model_data,n_mcmc = 3000,
                                    metropolis_alg = "RwM" , 
                                    prior = "gamma",
                                    def_time = 60){
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


#all model in nuts algorithm
stan_model_unif <- stan_model(file = "Stan_models/stan_wishart_inv_wish_unif.stan", model_name = "whole_stan_model_unif")
stan_model_exp <- stan_model(file = "Stan_models/stan_wishart_inv_wish_exp.stan", model_name = "whole_stan_model_exp")
stan_model_gamma <- stan_model(file = "Stan_models/stan_wishart_inv_wish_gamma.stan", model_name = "whole_stan_model_gamma")
stan_model_invgamma <- stan_model(file = "Stan_models/stan_wishart_inv_wish_invgamma.stan", model_name = "whole_stan_model_invgamma")
stan_model_lognormal <- stan_model(file = "Stan_models/stan_wishart_inv_wish_lognormal.stan", model_name = "whole_stan_model_lognormal")

all_modelling_nuts <- function(model_data,n_mcmc = 3000, 
                               prior = "gamma", def_time = 60){
  
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
    list(V = diag(1,p), n=p+10)
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

obs_list <- c(5,10,30,50,100)
d_list_low <- c(2,5,10)
d_list_high <- c(30,50)
dof_list <- c(3,6,11,31,51,101)

#write_xlsx(data.frame(),"simulated_experiments_results_time_iter.xlsx")


experiments_def <- function(obs, d_list,dof_list){
  df_all_list <- list()
  
  mcmc_algos <- c("RwM", "Slice", "HMC")
  priors_dist <- c("uniform","exponential","gamma","inverse_gamma","log_normal")
  
  counter <- 1
  for (d in d_list[d_list<=obs]){
    for (dof in dof_list[dof_list>d]){
      for (iter in 1:2){
        A <- matrix(runif(d^2,-10,10), ncol=d) 
        Sigma <- t(A) %*% A
        upper_sigma <- Sigma[upper.tri(Sigma, diag = T)]
        
        model_data <- list()
        for(i in 1:obs){model_data[[i]] <- rWishart(1,dof,Sigma)[,,1]}
        
        
        df_bisection_hist <- bisection_method(degrees_of_freedom, L = d, U = 1000,num = 1000, tol = 1e-05,model_data, d)
        df_bisection <- df_bisection_hist[length(df_bisection_hist)]
        
        NR_hist <- NR_wishart(model_data =model_data, p=d, initial = df_bisection)
        df_NR <- NR_hist[length(NR_hist)]
        V_est_NR <- v_estimation(df_NR,obs,d,model_data)
        PE_df_Newton_Raphson <- abs((dof-df_NR)/dof)
        PE_V_Newton_Raphson <- mean(abs((V_est_NR[upper.tri(V_est_NR, diag = T)] - upper_sigma)/upper_sigma))
        
        df_all_iter <- data.frame(obs = obs, d = d, dof = dof, 
                                  method = 'Newton_Raphson', prior = '',
                                  iter = 0,PE_dof = PE_df_Newton_Raphson,
                                  PE_V = PE_V_Newton_Raphson,
                                  ESS = 0, MCSE = 0 , Geweke = 0)
        
        for (alg in mcmc_algos){
          for (pr in priors_dist){
            tryCatch( {

            MwG <- metropolis_within_gibbs(model_data, metropolis_alg = alg,
                                           prior = pr, def_time = 60)#30)#
            n_post <- MwG[[1]]
            V_post <- MwG[[2]]
            ess_df_mwg <- effectiveSize(n_post)[[1]]
            geweke_df_mwg <- abs(geweke.diag(n_post)[[1]])
            n_post_est <- median(n_post)
            V_post_est <- apply(abind(V_post, along = 3),c(1,2),mean)
            PE_df_mwg <- abs((dof-n_post_est)/dof)
            PE_V_mwg <-mean(abs((V_post_est[upper.tri(V_post_est, diag = T)] - upper_sigma)/upper_sigma))
            print(paste0('MwG_',alg))
            mcmc_iter <- length(n_post)
            mcse <- sqrt(var(n_post) / ess_df_mwg)
            df_all_iter <- rbind(df_all_iter, list(obs,d, dof, paste0('MwG_',alg),
                                                   pr, mcmc_iter, PE_df_mwg,
                                                   PE_V_mwg, ess_df_mwg, mcse,
                                                   geweke_df_mwg))

            }, error = function(e) {df_all_iter <- rbind(df_all_iter,
                                                         list(obs,d, dof, paste0('MwG_',alg),
                                                              pr, 0, NA,
                                                              NA, NA, NA,
                                                              NA))})
          }
        }
        
        for (pr in priors_dist){
          tryCatch( { 
            
              stan_all <- all_modelling_nuts(model_data = model_data, prior = pr, 
                                                     def_time = 60)#30)#
              n_post <- as.numeric(stan_all[[1]])
              V_post <- array(as.numeric(stan_all[[2]]), dim = dim(stan_all[[2]]))
              ess_df_stan <- effectiveSize(as.vector(n_post))[[1]]
              
              geweke_df_stan <- abs(geweke.diag(as.vector(n_post))[[1]])
              
              n_post_est <- median(n_post)
              V_post_est <- apply(V_post,c(2,3),mean)
              PE_df_stan <- abs((dof-n_post_est)/dof)
              PE_V_stan <-mean(abs((V_post_est[upper.tri(V_post_est, diag = T)] - upper_sigma)/upper_sigma))
              print('stan_all_model')
              mcmc_iter <- length(n_post)
              mcse <- sqrt(var(n_post) / ess_df_stan)
              df_all_iter <- rbind(df_all_iter, list(obs,d, dof, 'stan_all_model',
                                             pr, mcmc_iter, PE_df_stan,
                                             PE_V_stan, ess_df_stan, mcse, 
                                             geweke_df_stan))
              
          }, error = function(e) {df_all_iter <- rbind(df_all_iter, 
                                                       list(obs,d, dof, 'stan_all_model',
                                                            pr, 0, NA,
                                                            NA, NA, NA, 
                                                            NA))})
              
          
      
        }
        
        
        df_all_list[[counter]] <- df_all_iter
        counter <- counter+1
        
        print(paste0("obs:",obs,", d:",d,", dof:",dof,", iter:",iter))
        
      }
    }
    # df_all_append <- read_excel("simulated_experiments_results_time_iter.xlsx")
    # df_all_save <- do.call("rbind", df_all_list)
    # df_all_append <- rbind(df_all_append,df_all_save)
    # df_all_append <- df_all_append[!duplicated(df_all_append), ]
    # write_xlsx(df_all_append,"simulated_experiments_results_time_iter.xlsx")
  }
  df_all <- do.call("rbind", df_all_list)
  
  return(df_all)
}

example_df = experiments_def(20,5,15)
example_df

# detectCores()
# cl <- makeCluster(5,
#                   outfile="log.txt")
# clusterExport(cl,
#               list("sum_det_list","log_v_estimation","degrees_of_freedom",
#                    "v_estimation","bisection_method",
#                    "second_derivative","NR_wishart",
#                    "n_pdf","metropolis_within_gibbs",
#                    "read_excel","write_xlsx","abind",
#                    "rInvWishart","lmvgamma","RWM_for_dof",
#                    "slice_sample","n_pdf_gradient","HMC",
#                    "all_modelling_nuts","stan_model_unif",
#                    "stan_model_exp","stan_model_gamma",
#                    "stan_model_invgamma","stan_model_lognormal",
#                    "effectiveSize","geweke.diag","sampling","extract","withTimeout"),
#               envir=globalenv())
# 
# 
# 
# system.time({saves <- parLapply(cl, obs_list,
#                                 experiments_def,
#                                 d_list = d_list,
#                                 dof_list = dof_list)})
# 
# stopCluster(cl)
# 
# df_all_experiments <- do.call("rbind", saves)
# 
# 
# write_xlsx(df_all_experiments,"simulated_experiments_results_time.xlsx")
