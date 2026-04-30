# Wishart_dof
This is the public repository of the paper 'Bayesian Computation for Wishart Degrees of Freedom Estimation'

Abstract:

The main goal of this work is to efficiently estimate the posterior distribution of the degrees of freedom of the Wishart distribution, which is used to model covariance matrices, through various Markov Chain Monte Carlo (MCMC) sampling methods. It also explores different prior setups for the degrees of freedom of the Wishart distribution. Specifically, we implemented and evaluated several hybrid MCMC methods, including Random Walk Metropolis within Gibbs, Slice Sampling within Gibbs, and Hamiltonian Monte Carlo within Gibbs. These methods were compared against the No-U-Turn Sampler (NUTS) as implemented in Stan. The proposed hybrid methods efficiently utilize the conjugate form of the conditional distribution of the scale matrix (given the degrees of freedom), provided by the inverse Wishart distribution. All methods are illustrated using a thorough simulation study and three real-life datasets. These illustrations provide a controlled environment for assessing performance and practical validation. Each assessed computational method exhibits distinct strengths and weaknesses in terms of robustness, accuracy, and sampling efficiency. From the simulation study, we conclude that, for low-dimensional covariance matrices ($p\leq15$), NUTS  performed best. On the other hand, Slice Sampling within Gibbs was more effective for high-dimensional data ($p>15$), where the use of NUTS  becomes less feasible due to the substantial computational time it requires. The results from real datasets are consistent  with those from the simulation study, confirming our initial conclusions about the reliability of the methods.


Github structure:

Folders:
* 'Stan_models' consist the considered model for different prior setup.
* 'experimental_results' consist the simulation studies results with the descriptive R files reported in the paper.
* 'plots' the produced illustrations of the paper.
* 'real_datasets' contains the three real dataset of the paper.
* 'rstan_sample' help directory for the stan samples when you have a time constraint.

Files:
* 'real_dataset_experiments.R' contains all the methods and the three datasets results with predefined number of iterations.
* 'real_dataset_experiments_time.R' contains all the methods and the three datasets results with time constrainted.
* 'whole_process_experiments.R' contains all the methods and simulated experiments conducted code with simulated samples from Wishart distribution with predefined number of iterations.
* 'whole_process_experiments_time.R' contains all the methods and simulated experiments conducted code with simulated samples from Wishart distribution with time constrained.

  
