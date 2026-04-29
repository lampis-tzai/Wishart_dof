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
  real n;
  cov_matrix[P] V;
}

model {
  V ~ inv_wishart(nv, U);
  n ~ lognormal(a,sqrt(b));

  //Likelihood of the observed data
  for (m in 1:M) {
  y[m] ~ wishart(n, V);
  }
}

