# ============================================================
# wishart_df_methods.R
#
# Methods for Bayesian inference on the degrees of freedom of a
# Wishart likelihood with an unknown common scale matrix.
#
# Implements:
#   1. Maximum likelihood estimation (MLE)
#   2. Random-Walk Metropolis within Gibbs (RWM)
#   3. Collapsed Bayesian inference
#
# ============================================================

# -----------------------------
# Basic matrix/data utilities
# -----------------------------

as_matrix_list <- function(X) {
  if (is.list(X)) {
    X_list <- X
  } else if (length(dim(X)) == 3) {
    X_list <- lapply(seq_len(dim(X)[3]), function(i) X[, , i])
  } else if (is.matrix(X)) {
    X_list <- list(X)
  } else {
    stop("X must be a list of matrices or a 3D array.")
  }

  if (length(X_list) < 1L) {
    stop("X must contain at least one matrix.")
  }

  p <- nrow(X_list[[1]])
  if (ncol(X_list[[1]]) != p) {
    stop("Each X_i must be square.")
  }

  for (i in seq_along(X_list)) {
    Xi <- X_list[[i]]

    if (!is.matrix(Xi) || nrow(Xi) != p || ncol(Xi) != p) {
      stop("All observations must be square matrices of the same dimension.")
    }

    Xi <- (Xi + t(Xi)) / 2

    if (min(eigen(Xi, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      stop(sprintf("X[[%d]] is not positive definite.", i))
    }

    X_list[[i]] <- Xi
  }

  X_list
}

logdet_spd <- function(A) {
  A <- (A + t(A)) / 2
  R <- tryCatch(chol(A), error = function(e) NULL)

  if (is.null(R)) {
    stop("Matrix is not positive definite.")
  }

  2 * sum(log(diag(R)))
}

log_multigamma <- function(a, p) {
  if (a <= (p - 1) / 2) return(-Inf)

  p * (p - 1) / 4 * log(pi) +
    sum(lgamma(a - (0:(p - 1)) / 2))
}

# Vectorized version used for efficient collapsed-grid evaluation.
log_multigamma_vec <- function(a, p) {
  a <- as.numeric(a)

  if (any(a <= (p - 1) / 2)) {
    out <- rep(-Inf, length(a))
    ok <- a > (p - 1) / 2

    if (any(ok)) {
      mat <- outer(
        a[ok],
        (0:(p - 1)) / 2,
        FUN = "-"
      )

      out[ok] <-
        p * (p - 1) / 4 * log(pi) +
        rowSums(lgamma(mat))
    }

    return(out)
  }

  mat <- outer(
    a,
    (0:(p - 1)) / 2,
    FUN = "-"
  )

  p * (p - 1) / 4 * log(pi) +
    rowSums(lgamma(mat))
}

sample_inverse_wishart <- function(df, scale) {
  # If W ~ Wishart(scale^{-1}, df), then W^{-1} ~ IW(scale, df).
  scale <- (scale + t(scale)) / 2
  scale_inv <- solve(scale)

  W <- stats::rWishart(
    1,
    df = df,
    Sigma = scale_inv
  )[, , 1]

  W <- (W + t(W)) / 2
  V <- chol2inv(chol(W))

  (V + t(V)) / 2
}

prepare_wishart_data <- function(X) {
  X_list <- as_matrix_list(X)
  m <- length(X_list)
  p <- nrow(X_list[[1]])

  S <- Reduce("+", X_list)
  Xbar <- S / m

  logdet_X <- vapply(
    X_list,
    logdet_spd,
    numeric(1)
  )

  sum_logdet_X <- sum(logdet_X)
  mean_logdet_X <- mean(logdet_X)
  logdet_Xbar <- logdet_spd(Xbar)

  list(
    X = X_list,
    m = m,
    p = p,
    S = S,
    Xbar = Xbar,
    sum_logdet_X = sum_logdet_X,
    mean_logdet_X = mean_logdet_X,
    logdet_Xbar = logdet_Xbar
  )
}

# -----------------------------
# Random positive-definite V
# -----------------------------

generate_random_scale <- function(p, min = -10, max = 10) {
  B <- matrix(
    stats::runif(p * p, min = min, max = max),
    nrow = p,
    ncol = p
  )

  V <- crossprod(B)
  V <- V / mean(diag(V))
  V <- (V + t(V)) / 2

  # Very rare numerical protection for an exceptionally ill-conditioned draw.
  ev <- eigen(
    V,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  if (min(ev) <= 1e-10) {
    V <- V + diag(1e-6, p)
    V <- V / mean(diag(V))
  }

  V
}

# ============================================================
# 1. Maximum likelihood estimation
# ============================================================

wishart_profile_score <- function(n, prep) {
  p <- prep$p

  if (n <= p - 1) return(-Inf)

  j <- 1:p

  sum(digamma((n - j + 1) / 2)) -
    prep$mean_logdet_X +
    prep$logdet_Xbar -
    p * log(n / 2)
}

wishart_profile_score_derivative <- function(n, prep) {
  p <- prep$p
  j <- 1:p

  0.5 * sum(trigamma((n - j + 1) / 2)) -
    p / n
}

wishart_profile_loglik <- function(n, prep) {
  p <- prep$p
  m <- prep$m

  if (n <= p - 1) return(-Inf)

  ll <- 0.5 * (n - p - 1) * prep$sum_logdet_X
  ll <- ll - 0.5 * m * n * p * log(2)
  ll <- ll - m * log_multigamma(n / 2, p)
  ll <- ll - 0.5 * m * n *
    (prep$logdet_Xbar - p * log(n))
  ll <- ll - 0.5 * m * n * p

  ll
}

fit_mle <- function(X,
                    lower_eps = 1e-8,
                    upper_init = NULL,
                    upper_max = 1e6,
                    newton_tol = 1e-10,
                    newton_maxit = 25) {

  prep <- prepare_wishart_data(X)
  p <- prep$p
  lower <- p - 1 + lower_eps

  start_time <- proc.time()[["elapsed"]]

  if (is.null(upper_init)) {
    upper_init <- max(10 * p, 100)
  }

  upper <- upper_init

  # Increase the upper bracket until the score becomes positive.
  f_upper <- wishart_profile_score(upper, prep)

  while (
    is.finite(f_upper) &&
    f_upper <= 0 &&
    upper < upper_max
  ) {
    upper <- min(2 * upper, upper_max)
    f_upper <- wishart_profile_score(upper, prep)
  }

  boundary_solution <- FALSE

  if (is.finite(f_upper) && f_upper > 0) {

    # Bisection stage.
    lo <- lower
    hi <- upper

    for (iter in 1:200) {
      mid <- (lo + hi) / 2
      fmid <- wishart_profile_score(mid, prep)

      if (!is.finite(fmid)) {
        lo <- mid
      } else if (fmid > 0) {
        hi <- mid
      } else {
        lo <- mid
      }

      if ((hi - lo) < newton_tol * max(1, abs(mid))) {
        break
      }
    }

    n_hat <- (lo + hi) / 2

    # Newton-Raphson refinement.
    for (iter in 1:newton_maxit) {
      f <- wishart_profile_score(n_hat, prep)
      fp <- wishart_profile_score_derivative(n_hat, prep)

      if (!is.finite(f) || !is.finite(fp) || fp <= 0) {
        break
      }

      step <- f / fp
      candidate <- n_hat - step

      if (
        !is.finite(candidate) ||
        candidate <= lower ||
        candidate >= upper
      ) {
        break
      }

      n_hat <- candidate

      if (abs(step) < newton_tol * max(1, abs(n_hat))) {
        break
      }
    }

  } else {

    # Numerical fallback for unusual datasets.
    opt <- optimize(
      f = function(n) wishart_profile_loglik(n, prep),
      interval = c(lower + 1e-6, upper_max),
      maximum = TRUE
    )

    n_hat <- opt$maximum
    boundary_solution <- abs(n_hat - upper_max) <
      1e-5 * max(1, upper_max)
  }

  V_hat <- prep$Xbar / n_hat
  V_hat <- (V_hat + t(V_hat)) / 2

  runtime <- proc.time()[["elapsed"]] - start_time

  list(
    method = "MLE",
    n = list(
      estimate = n_hat,
      mean = n_hat,
      median = n_hat,
      mode = n_hat,
      sd = NA_real_,
      lower = NA_real_,
      upper = NA_real_
    ),
    V = list(
      estimate = V_hat,
      mean = V_hat
    ),
    runtime = runtime,
    diagnostics = list(
      boundary_solution = boundary_solution,
      score = wishart_profile_score(n_hat, prep),
      profile_loglik = wishart_profile_loglik(n_hat, prep)
    )
  )
}

# ============================================================
# 2. RWM within Gibbs
# ============================================================

# Broad proper Uniform prior used throughout the paper.
log_prior_uniform_n <- function(n, p, upper = 1e3) {
  lower <- p - 1

  if (n <= lower || n >= upper) {
    return(-Inf)
  }

  -log(upper - lower)
}

log_posterior_n_given_V <- function(
    n,
    prep,
    V,
    prior_upper = 1e3) {

  p <- prep$p
  m <- prep$m

  lp <- log_prior_uniform_n(
    n,
    p,
    prior_upper
  )

  if (!is.finite(lp)) {
    return(-Inf)
  }

  lp <- lp +
    0.5 * (n - p - 1) * prep$sum_logdet_X -
    0.5 * m * n * p * log(2) -
    0.5 * m * n * logdet_spd(V) -
    m * log_multigamma(n / 2, p)

  lp
}

rw_update_n <- function(
    n_current,
    prep,
    V,
    delta,
    prior_upper = 1e3) {

  p <- prep$p
  lower <- p - 1

  n_proposal <- n_current +
    stats::runif(1, -delta, delta)

  if (
    n_proposal <= lower ||
    n_proposal >= prior_upper
  ) {
    return(
      list(
        n = n_current,
        accepted = FALSE
      )
    )
  }

  lp_current <- log_posterior_n_given_V(
    n_current,
    prep,
    V,
    prior_upper
  )

  lp_proposal <- log_posterior_n_given_V(
    n_proposal,
    prep,
    V,
    prior_upper
  )

  log_alpha <- lp_proposal - lp_current

  if (
    log(stats::runif(1)) <
    min(0, log_alpha)
  ) {
    list(
      n = n_proposal,
      accepted = TRUE
    )
  } else {
    list(
      n = n_current,
      accepted = FALSE
    )
  }
}

summarize_mcmc_n <- function(draws) {
  qs <- stats::quantile(
    draws,
    probs = c(0.025, 0.5, 0.975),
    names = FALSE
  )

  dens <- stats::density(
    draws,
    n = 512
  )

  mode <- dens$x[which.max(dens$y)]

  list(
    # Primary point estimate = posterior mode.
    estimate = mode,
    mean = mean(draws),
    median = qs[2],
    mode = mode,
    sd = stats::sd(draws),
    lower = qs[1],
    upper = qs[3]
  )
}

fit_rwm <- function(
    X,
    iter = 4000,
    burn = 1000,
    delta = 2,
    prior_upper = 1e3,
    init_n = NULL,
    n_v = NULL,
    U = NULL,
    seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (iter <= burn) {
    stop("iter must be larger than burn.")
  }

  if (!is.finite(delta) || delta <= 0) {
    stop("delta must be a positive finite number.")
  }

  prep <- prepare_wishart_data(X)
  p <- prep$p
  m <- prep$m
  lower <- p - 1

  if (is.null(n_v)) {
    n_v <- p
  }

  if (is.null(U)) {
    U <- diag(1e-4, p)
  }

  # Fixed data-independent starting value.
  # No MLE is used to initialize or tune the RWM chain.
  if (is.null(init_n)) {
    init_n <- 4 * p
  }

  init_n <- min(
    max(init_n, lower + 1e-6),
    prior_upper - 1e-6
  )

  A <- prep$S + U

  n_chain <- numeric(iter)
  V_chain <- array(
    NA_real_,
    dim = c(p, p, iter)
  )

  n_chain[1] <- init_n

  V_chain[, , 1] <- sample_inverse_wishart(
    df = n_v + m * n_chain[1],
    scale = A
  )

  accepted <- 0L
  start_time <- proc.time()[["elapsed"]]

  for (b in 2:iter) {

    # Exact Gibbs update for V | n, X.
    V_chain[, , b] <- sample_inverse_wishart(
      df = n_v + m * n_chain[b - 1],
      scale = A
    )

    # RWM update for n | V, X.
    step <- rw_update_n(
      n_current = n_chain[b - 1],
      prep = prep,
      V = V_chain[, , b],
      delta = delta,
      prior_upper = prior_upper
    )

    n_chain[b] <- step$n
    accepted <- accepted +
      as.integer(step$accepted)
  }

  runtime <- proc.time()[["elapsed"]] - start_time

  keep <- seq.int(burn + 1, iter)

  n_post <- n_chain[keep]
  V_post <- V_chain[, , keep, drop = FALSE]

  V_mean <- apply(
    V_post,
    c(1, 2),
    mean
  )

  V_mean <- (V_mean + t(V_mean)) / 2

  list(
    method = "RWM",
    n = summarize_mcmc_n(n_post),
    V = list(
      estimate = V_mean,
      mean = V_mean,
      draws = dim(V_post)[3]
    ),
    n_draws = n_post,
    V_draws = V_post,
    runtime = runtime,
    diagnostics = list(
      acceptance_rate = accepted / (iter - 1),
      iter = iter,
      burn = burn,
      delta = delta,
      init_n = init_n
    )
  )
}

# ============================================================
# 3. Collapsed Bayesian inference
# ============================================================

log_collapsed_posterior_n <- function(
    n,
    prep,
    n_v = NULL,
    U = NULL,
    prior_upper = 1e3) {

  p <- prep$p
  m <- prep$m

  if (is.null(n_v)) {
    n_v <- p
  }

  if (is.null(U)) {
    U <- diag(1e-4, p)
  }

  lp <- log_prior_uniform_n(
    n,
    p,
    prior_upper
  )

  if (!is.finite(lp) || n <= p - 1) {
    return(-Inf)
  }

  A <- prep$S + U
  logdet_A <- logdet_spd(A)

  # The 2^{-mnp/2} term from the Wishart likelihood cancels
  # with the corresponding factor from the inverse-Wishart integral.
  lp <- lp +
    log_multigamma((n_v + m * n) / 2, p) -
    m * log_multigamma(n / 2, p) -
    0.5 * (n_v + m * n) * logdet_A +
    0.5 * (n - p - 1) * prep$sum_logdet_X

  lp
}

# Vectorized collapsed log-posterior for the full grid.
log_collapsed_posterior_grid <- function(
    n,
    prep,
    n_v,
    U,
    prior_upper) {

  p <- prep$p
  m <- prep$m
  lower <- p - 1

  if (any(n <= lower) || any(n >= prior_upper)) {
    stop("All grid values must lie inside the prior support.")
  }

  A <- prep$S + U
  logdet_A <- logdet_spd(A)

  log_prior <- -log(prior_upper - lower)

  log_prior +
    log_multigamma_vec(
      (n_v + m * n) / 2,
      p
    ) -
    m * log_multigamma_vec(
      n / 2,
      p
    ) -
    0.5 * (n_v + m * n) * logdet_A +
    0.5 * (n - p - 1) * prep$sum_logdet_X
}

trapezoid_weights <- function(x) {
  dx <- diff(x)

  if (length(dx) < 1L) {
    stop("At least two grid points are required.")
  }

  w <- numeric(length(x))

  w[1] <- dx[1] / 2
  w[length(x)] <- dx[length(dx)] / 2

  if (length(x) > 2) {
    w[2:(length(x) - 1)] <-
      (dx[-length(dx)] + dx[-1]) / 2
  }

  w
}

grid_quantile <- function(x, density_values, q) {
  if (q < 0 || q > 1) {
    stop("q must be between 0 and 1.")
  }

  w <- trapezoid_weights(x)
  mass <- density_values * w
  cdf <- cumsum(mass)
  cdf <- cdf / max(cdf)

  if (q <= cdf[1]) {
    return(x[1])
  }

  if (q >= cdf[length(cdf)]) {
    return(x[length(x)])
  }

  approx(
    cdf,
    x,
    xout = q,
    ties = "ordered"
  )$y
}

make_collapsed_grid <- function(
    prep,
    n_v,
    U,
    prior_upper,
    grid_size = 10000) {

  if (length(grid_size) != 1L ||
      !is.finite(grid_size) ||
      grid_size < 3) {
    stop("grid_size must be an integer greater than or equal to 3.")
  }

  grid_size <- as.integer(grid_size)

  p <- prep$p

  # Full support of the prior: p - 1 < n < prior_upper.
  lower <- p - 1 + 1e-7
  upper <- prior_upper - 1e-7

  if (upper <= lower) {
    stop("prior_upper must be greater than p - 1.")
  }

  # Fixed grid over the entire prior support.
  grid <- seq(
    lower,
    upper,
    length.out = grid_size
  )

  # Evaluate the collapsed log-posterior over the complete support.
  logpost <- log_collapsed_posterior_grid(
    n = grid,
    prep = prep,
    n_v = n_v,
    U = U,
    prior_upper = prior_upper
  )

  list(
    grid = grid,
    logpost = logpost
  )
}

fit_collapsed <- function(
    X,
    n_v = NULL,
    U = NULL,
    prior_upper = 1e3,
    grid_size = 10000,
    n_draws = 0,
    seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  prep <- prepare_wishart_data(X)
  p <- prep$p

  if (is.null(n_v)) {
    n_v <- p
  }

  if (is.null(U)) {
    U <- diag(1e-4, p)
  }

  start_time <- proc.time()[["elapsed"]]

  # Evaluate posterior over the full prior support.
  grid_obj <- make_collapsed_grid(
    prep = prep,
    n_v = n_v,
    U = U,
    prior_upper = prior_upper,
    grid_size = grid_size
  )

  n_grid <- grid_obj$grid
  logpost <- grid_obj$logpost

  # Numerical stabilization.
  max_lp <- max(logpost)
  unnorm <- exp(logpost - max_lp)

  # Numerical integration weights.
  weights <- trapezoid_weights(n_grid)
  mass <- unnorm * weights
  mass <- mass / sum(mass)

  # Posterior mean.
  n_mean <- sum(n_grid * mass)

  # Posterior standard deviation.
  n_sd <- sqrt(
    sum(
      (n_grid - n_mean)^2 * mass
    )
  )

  # Posterior median and credible interval.
  n_median <- grid_quantile(
    n_grid,
    unnorm,
    0.50
  )

  n_lower <- grid_quantile(
    n_grid,
    unnorm,
    0.025
  )

  n_upper <- grid_quantile(
    n_grid,
    unnorm,
    0.975
  )

  #locate maximum -> continuous local optimization
  mode_idx <- which.max(logpost)
  grid_step <- n_grid[2] - n_grid[1]

  left_idx <- max(
    1L,
    mode_idx - 2L
  )

  right_idx <- min(
    length(n_grid),
    mode_idx + 2L
  )

  mode_left <- max(
    p - 1 + 1e-8,
    n_grid[left_idx]
  )

  mode_right <- min(
    prior_upper - 1e-8,
    n_grid[right_idx]
  )

  if (mode_right <= mode_left) {
    n_mode <- n_grid[mode_idx]
  } else {

    mode_opt <- tryCatch(
      optimize(
        f = function(n) {
          log_collapsed_posterior_n(
            n = n,
            prep = prep,
            n_v = n_v,
            U = U,
            prior_upper = prior_upper
          )
        },
        interval = c(mode_left, mode_right),
        maximum = TRUE
      ),
      error = function(e) NULL
    )

    if (is.null(mode_opt) || !is.finite(mode_opt$maximum)) {
      n_mode <- n_grid[mode_idx]
    } else {
      n_mode <- mode_opt$maximum
    }
  }

  # Avoid an unused-variable note while retaining grid-step metadata.
  invisible(grid_step)

  # ----------------------------------------------------------
  # Posterior mean of V
  # ----------------------------------------------------------

  A <- prep$S + U

  denominator <-
    n_v + prep$m * n_grid - p - 1

  v_scalar_mean <- sum(
    (1 / denominator) * mass
  )

  V_mean <- A * v_scalar_mean
  V_mean <- (V_mean + t(V_mean)) / 2

  # ----------------------------------------------------------
  # Optional posterior draws
  # ----------------------------------------------------------

  n_post_draws <- NULL
  V_post_draws <- NULL

  if (n_draws > 0) {

    # Independent draws from the discretized marginal posterior of n.
    idx <- sample(
      seq_along(n_grid),
      size = n_draws,
      replace = TRUE,
      prob = mass
    )

    n_post_draws <- n_grid[idx]

    V_post_draws <- array(
      NA_real_,
      dim = c(p, p, n_draws)
    )

    for (b in seq_len(n_draws)) {
      V_post_draws[, , b] <-
        sample_inverse_wishart(
          df = n_v + prep$m * n_post_draws[b],
          scale = A
        )
    }
  }

  runtime <- (
    proc.time()[["elapsed"]] -
      start_time
  )

  list(
    method = "Collapsed",
    n = list(
      # Primary point estimate = posterior mode.
      estimate = n_mode,
      mean = n_mean,
      median = n_median,
      mode = n_mode,
      sd = n_sd,
      lower = n_lower,
      upper = n_upper
    ),
    V = list(
      estimate = V_mean,
      mean = V_mean
    ),
    n_draws = n_post_draws,
    V_draws = V_post_draws,
    runtime = runtime,
    diagnostics = list(
      grid_size = length(n_grid),
      grid_lower = min(n_grid),
      grid_upper = max(n_grid),
      grid_step = n_grid[2] - n_grid[1],
      posterior_at_upper = tail(unnorm, 1)
    )
  )
}






########################################
# optimize and integrate function
########################################



fit_collapsed_adaptive <- function(
    X,
    n_v = NULL,
    U = NULL,
    prior_upper = 1e3,
    grid_size = 10000,
    n_draws = 0,
    seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # ----------------------------------------------------------
  # Prepare data
  # ----------------------------------------------------------
  
  prep <- prepare_wishart_data(X)
  p <- prep$p
  
  if (is.null(n_v)) {
    n_v <- p
  }
  
  if (is.null(U)) {
    U <- diag(1e-4, p)
  }
  
  U <- (U + t(U)) / 2
  
  start_time <- proc.time()[["elapsed"]]
  
  # ----------------------------------------------------------
  # Prior support
  # ----------------------------------------------------------
  
  lower <- p - 1 + 1e-8
  upper <- prior_upper - 1e-8
  
  # ----------------------------------------------------------
  # Evaluate posterior on a grid
  #
  # The grid is retained for:
  #   1. locating the region of the posterior mode
  #   2. optional posterior draws
  #
  # Posterior moments and intervals below are obtained using
  # adaptive numerical integration.
  # ----------------------------------------------------------
  
  grid_obj <- make_collapsed_grid(
    prep = prep,
    n_v = n_v,
    U = U,
    prior_upper = prior_upper,
    grid_size = grid_size
  )
  
  n_grid <- grid_obj$grid
  logpost_grid <- grid_obj$logpost
  
  # ----------------------------------------------------------
  # Locate the maximum on the grid
  # ----------------------------------------------------------
  
  mode_idx <- which.max(logpost_grid)
  
  grid_step <- n_grid[2] - n_grid[1]
  
  left_idx <- max(
    1L,
    mode_idx - 2L
  )
  
  right_idx <- min(
    length(n_grid),
    mode_idx + 2L
  )
  
  mode_left <- max(
    lower,
    n_grid[left_idx]
  )
  
  mode_right <- min(
    upper,
    n_grid[right_idx]
  )
  
  # ----------------------------------------------------------
  # Continuous posterior mode using optimize()
  # ----------------------------------------------------------
  
  if (mode_right <= mode_left) {
    
    n_mode <- n_grid[mode_idx]
    
  } else {
    
    mode_opt <- tryCatch(
      
      optimize(
        f = function(n) {
          
          log_collapsed_posterior_n(
            n = n,
            prep = prep,
            n_v = n_v,
            U = U,
            prior_upper = prior_upper
          )
        },
        
        interval = c(
          mode_left,
          mode_right
        ),
        
        maximum = TRUE
      ),
      
      error = function(e) {
        NULL
      }
    )
    
    if (
      is.null(mode_opt) ||
      !is.finite(mode_opt$maximum)
    ) {
      
      n_mode <- n_grid[mode_idx]
      
    } else {
      
      n_mode <- mode_opt$maximum
    }
  }
  
  # ----------------------------------------------------------
  # Maximum log-posterior
  #
  # Used to stabilize the numerical integration.
  # ----------------------------------------------------------
  
  max_lp <- log_collapsed_posterior_n(
    n = n_mode,
    prep = prep,
    n_v = n_v,
    U = U,
    prior_upper = prior_upper
  )
  
  if (!is.finite(max_lp)) {
    stop("Could not obtain a finite posterior mode.")
  }
  
  # ----------------------------------------------------------
  # Scalar stabilized posterior kernel
  # ----------------------------------------------------------
  
  post_kernel_scalar <- function(n) {
    
    lp <- log_collapsed_posterior_n(
      n = n,
      prep = prep,
      n_v = n_v,
      U = U,
      prior_upper = prior_upper
    )
    
    if (!is.finite(lp)) {
      return(0)
    }
    
    exp(lp - max_lp)
  }
  
  # ----------------------------------------------------------
  # Vectorized wrapper for integrate()
  #
  # integrate() may evaluate the function at several points
  # simultaneously, so the scalar function is evaluated
  # point-by-point.
  # ----------------------------------------------------------
  
  adaptive_integral <- function(
    fun,
    a,
    b,
    subdivisions = 1000,
    rel.tol = 1e-8) {
    
    fun_vectorized <- function(x) {
      
      vapply(
        x,
        fun,
        numeric(1)
      )
    }
    
    integrate(
      f = fun_vectorized,
      lower = a,
      upper = b,
      subdivisions = subdivisions,
      rel.tol = rel.tol,
      abs.tol = 0
    )
  }
  
  # ----------------------------------------------------------
  # Normalizing constant
  # ----------------------------------------------------------
  
  Z_result <- adaptive_integral(
    fun = post_kernel_scalar,
    a = lower,
    b = upper
  )
  
  Z <- Z_result$value
  
  if (
    !is.finite(Z) ||
    Z <= 0
  ) {
    stop(
      "Numerical integration failed for the posterior normalizing constant."
    )
  }
  
  # ----------------------------------------------------------
  # Posterior mean
  # ----------------------------------------------------------
  
  mean_result <- adaptive_integral(
    fun = function(n) {
      n * post_kernel_scalar(n)
    },
    a = lower,
    b = upper
  )
  
  n_mean <- mean_result$value / Z
  
  # ----------------------------------------------------------
  # Posterior second moment
  # ----------------------------------------------------------
  
  second_result <- adaptive_integral(
    fun = function(n) {
      n^2 * post_kernel_scalar(n)
    },
    a = lower,
    b = upper
  )
  
  n_second <- second_result$value / Z
  
  # ----------------------------------------------------------
  # Posterior standard deviation
  # ----------------------------------------------------------
  
  n_variance <- max(
    n_second - n_mean^2,
    0
  )
  
  n_sd <- sqrt(n_variance)
  
  # ----------------------------------------------------------
  # Posterior CDF
  # ----------------------------------------------------------
  
  posterior_cdf <- function(q) {
    
    if (q <= lower) {
      return(0)
    }
    
    if (q >= upper) {
      return(1)
    }
    
    result <- adaptive_integral(
      fun = post_kernel_scalar,
      a = lower,
      b = q
    )
    
    cdf_value <- result$value / Z
    
    min(
      max(
        cdf_value,
        0
      ),
      1
    )
  }
  
  # ----------------------------------------------------------
  # Posterior median
  # ----------------------------------------------------------
  
  median_opt <- tryCatch(
    
    uniroot(
      function(q) {
        posterior_cdf(q) - 0.50
      },
      interval = c(
        lower,
        upper
      ),
      tol = 1e-8
    ),
    
    error = function(e) NULL
  )
  
  if (is.null(median_opt)) {
    n_median <- n_grid[
      which(
        cumsum(
          exp(logpost_grid - max(logpost_grid))
        ) /
          sum(
            exp(logpost_grid - max(logpost_grid))
          ) >= 0.50
      )[1]
    ]
  } else {
    n_median <- median_opt$root
  }
  
  # ----------------------------------------------------------
  # Lower 95% posterior quantile
  # ----------------------------------------------------------
  
  lower_opt <- tryCatch(
    
    uniroot(
      function(q) {
        posterior_cdf(q) - 0.025
      },
      interval = c(
        lower,
        upper
      ),
      tol = 1e-8
    ),
    
    error = function(e) NULL
  )
  
  if (is.null(lower_opt)) {
    
    n_lower <- grid_quantile(
      n_grid,
      exp(logpost_grid - max(logpost_grid)),
      0.025
    )
    
  } else {
    
    n_lower <- lower_opt$root
  }
  
  # ----------------------------------------------------------
  # Upper 95% posterior quantile
  # ----------------------------------------------------------
  
  upper_opt <- tryCatch(
    
    uniroot(
      function(q) {
        posterior_cdf(q) - 0.975
      },
      interval = c(
        lower,
        upper
      ),
      tol = 1e-8
    ),
    
    error = function(e) NULL
  )
  
  if (is.null(upper_opt)) {
    
    n_upper <- grid_quantile(
      n_grid,
      exp(logpost_grid - max(logpost_grid)),
      0.975
    )
    
  } else {
    
    n_upper <- upper_opt$root
  }
  
  # ----------------------------------------------------------
  # Posterior mean of V
  #
  # V | n, X ~ IW(A, n_v + mn)
  #
  # E[V | n, X] =
  # A / (n_v + mn - p - 1)
  #
  # Hence
  #
  # E[V | X] =
  # A E[
  #   1 / (n_v + mn - p - 1)
  #   | X
  # ]
  # ----------------------------------------------------------
  
  A <- prep$S + U
  
  v_scalar_result <- adaptive_integral(
    
    fun = function(n) {
      
      post_kernel_scalar(n) /
        (
          n_v +
            prep$m * n -
            p -
            1
        )
    },
    
    a = lower,
    b = upper
  )
  
  v_scalar_mean <-
    v_scalar_result$value / Z
  
  V_mean <- A * v_scalar_mean
  
  V_mean <- (
    V_mean + t(V_mean)
  ) / 2
  
  # ----------------------------------------------------------
  # Optional posterior draws
  #
  # These retain the original grid-based sampling mechanism.
  # ----------------------------------------------------------
  
  n_post_draws <- NULL
  
  V_post_draws <- NULL
  
  if (n_draws > 0) {
    
    # Stabilized grid posterior
    unnorm_grid <- exp(
      logpost_grid -
        max(logpost_grid)
    )
    
    # Trapezoidal weights
    weights <- trapezoid_weights(
      n_grid
    )
    
    # Approximate posterior mass on grid
    mass <- unnorm_grid * weights
    
    mass <- mass / sum(mass)
    
    # Draw n
    idx <- sample(
      seq_along(n_grid),
      size = n_draws,
      replace = TRUE,
      prob = mass
    )
    
    n_post_draws <- n_grid[idx]
    
    # Draw V | n, X
    V_post_draws <- array(
      NA_real_,
      dim = c(
        p,
        p,
        n_draws
      )
    )
    
    for (b in seq_len(n_draws)) {
      
      V_post_draws[, , b] <-
        
        sample_inverse_wishart(
          df =
            n_v +
            prep$m *
            n_post_draws[b],
          
          scale = A
        )
    }
  }
  
  # ----------------------------------------------------------
  # Runtime
  # ----------------------------------------------------------
  
  runtime <- (
    proc.time()[["elapsed"]] -
      start_time
  )
  
  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------
  
  list(
    
    method = "Collapsed",
    
    n = list(
      
      # Primary point estimate
      # = continuous posterior mode
      
      estimate = n_mode,
      
      mean = n_mean,
      
      median = n_median,
      
      mode = n_mode,
      
      sd = n_sd,
      
      lower = n_lower,
      
      upper = n_upper
    ),
    
    V = list(
      
      estimate = V_mean,
      
      mean = V_mean
    ),
    
    n_draws = n_post_draws,
    
    V_draws = V_post_draws,
    
    runtime = runtime,
    
    diagnostics = list(
      
      grid_size = length(n_grid),
      
      grid_lower = min(n_grid),
      
      grid_upper = max(n_grid),
      
      grid_step =
        n_grid[2] - n_grid[1],
      
      posterior_at_upper =
        tail(unnorm_grid, 1),
      
      integration = "adaptive",
      
      optimization = "continuous_local"
    )
  )
}