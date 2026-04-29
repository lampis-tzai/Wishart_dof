data {
  int<lower=1> P;
  int<lower=1> M;
  cov_matrix[P] y[M];
  cov_matrix[P] U;
  real nv;
  real<lower=0> a;
  real<lower=0> b;
}


parameters {
real<lower=P-1> n; // Degrees of freedom
cholesky_factor_corr[P] L; // Cholesky factor of the scale matrix
}

model {
// Priors
n ~ gamma(a, b); // Gamma prior for degrees of freedom
L ~ lkj_corr_cholesky(0.0001); // LKJ prior for Cholesky factor of the scale matrix

// Likelihood
for (m in 1:M) {
y[m] ~ wishart_cholesky(n, L);
}
}