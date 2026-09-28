# Wishart_dof
This is the public repository of the paper 'Learning the Degrees of Freedom of a Wishart Likelihood from Positive-Definite Matrix Data: A Bayesian Approach'

Abstract:

In many applications of multivariate probabilistic modeling, the degrees of freedom parameter of the Wishart distribution is treated as fixed, although it is of great inferential interest. 
In this work, we consider Bayesian inference for the degrees of freedom of a Wishart likelihood when both the degrees of freedom and the scale matrix are unknown. 
In the latter conjugate structure of the model, the scale matrix is integrated out analytically, producing a collapsed approach with one-dimensional marginal posterior distribution for the degrees of freedom. 
The proposed approach is investigated through simulation studies examining posterior bias, uncertainty, and finite-sample behavior across different dimensions, numbers of observed positive-definite matrices, and degrees of freedom settings.
The collapsed approach is compared with a traditional Random-Walk Metropolis within Gibbs algorithm as an independent computational validation. 
The results show that the accuracy and concentration of posterior inference depend on the relationship between the number of observed matrices, the dimension, and the magnitude of the degrees of freedom, with greater uncertainty and bias arising in more demanding finite-sample settings. 
Finally, the methodology is illustrated using three real datasets, providing posterior inference for both the degrees of freedom and the scale matrix.
These applications demonstrate how the proposed Bayesian formulation can quantify uncertainty about the degrees of freedom when covariance-valued observations are modeled directly through a Wishart likelihood.


GitHub structure:

Folders:

* Simulation_experiments/: Contains the results for the simulation study.
* plots/: Contains the plots reported in the paper.
* real_datasets/: Contains the real datasets used in the paper.

Files:

* real_data_applications.R: Contains the code for Bayesian inference of the Wishart degrees of freedom and scale matrix for the real datasets.
* wishart_df_methods.R: Contains the implementations of the collapsed method, the Random-Walk Metropolis within Gibbs algorithm, and maximum likelihood estimation.
* wishart_df_simulations_parallel.R: Contains the code for the simulation study presented in the paper.




