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

  if (length(X_list) < 1L) stop("X must contain at least one matrix.")

  p <- nrow(X_list[[1]])
  if (ncol(X_list[[1]]) != p) stop("Each X_i must be square.")

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
  if (is.null(R)) stop("Matrix is not positive definite.")
  2 * sum(log(diag(R)))
}

log_multigamma <- function(a, p) {
  if (a <= (p - 1) / 2) return(-Inf)
  p * (p - 1) / 4 * log(pi) +
    sum(lgamma(a - (0:(p - 1)) / 2))
}

sample_inverse_wishart <- function(df, scale) {
  # If W ~ Wishart(scale^{-1}, df), then W^{-1} ~ IW(scale, df).
  scale <- (scale + t(scale)) / 2
  scale_inv <- solve(scale)
  W <- stats::rWishart(1, df = df, Sigma = scale_inv)[, , 1]
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
  logdet_X <- vapply(X_list, logdet_spd, numeric(1))
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
  A <- matrix(stats::runif(p * p, min = min, max = max), nrow = p, ncol = p)
  V <- crossprod(A)
  V <- V / mean(diag(V))
  V <- (V + t(V)) / 2

  # Very rare numerical protection for an exceptionally ill-conditioned draw.
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
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

  0.5 * sum(trigamma((n - j + 1) / 2)) - p / n
}

wishart_profile_loglik <- function(n, prep) {
  p <- prep$p
  m <- prep$m

  if (n <= p - 1) return(-Inf)

  ll <- 0.5 * (n - p - 1) * prep$sum_logdet_X
  ll <- ll - 0.5 * m * n * p * log(2)
  ll <- ll - m * log_multigamma(n / 2, p)
  ll <- ll - 0.5 * m * n * (prep$logdet_Xbar - p * log(n))
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

  # The score tends to -infinity near p-1. Increase the upper
  # bracket until the score becomes positive.
  f_upper <- wishart_profile_score(upper, prep)
  while (is.finite(f_upper) && f_upper <= 0 && upper < upper_max) {
    upper <- min(2 * upper, upper_max)
    f_upper <- wishart_profile_score(upper, prep)
  }

  boundary_solution <- FALSE

  if (is.finite(f_upper) && f_upper > 0) {

    # Bisection stage.
    lo <- lower
    hi <- upper
    flo <- wishart_profile_score(lo, prep)
    fhi <- wishart_profile_score(hi, prep)

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

      if ((hi - lo) < newton_tol * max(1, abs(mid))) break
    }

    n_hat <- (lo + hi) / 2

    # Newton-Raphson refinement.
    for (iter in 1:newton_maxit) {
      f <- wishart_profile_score(n_hat, prep)
      fp <- wishart_profile_score_derivative(n_hat, prep)

      if (!is.finite(f) || !is.finite(fp) || fp <= 0) break

      step <- f / fp
      candidate <- n_hat - step

      if (!is.finite(candidate) || candidate <= lower || candidate >= upper) {
        break
      }

      n_hat <- candidate

      if (abs(step) < newton_tol * max(1, abs(n_hat))) break
    }

  } else {
    # Fallback: maximize the profile log-likelihood over a broad
    # interval. This protects the function against unusual datasets
    # for which the score does not cross zero numerically.
    opt <- optimize(
      f = function(n) wishart_profile_loglik(n, prep),
      interval = c(lower + 1e-6, upper_max),
      maximum = TRUE
    )
    n_hat <- opt$maximum
    boundary_solution <- abs(n_hat - upper_max) < 1e-5 * max(1, upper_max)
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

log_prior_uniform_n <- function(n, p, upper = 1e4) {
  lower <- p - 1
  if (n <= lower || n >= upper) return(-Inf)
  -log(upper - lower)
}

log_posterior_n_given_V <- function(n, prep, V,
                                   prior_upper = 1e4) {

  p <- prep$p
  m <- prep$m

  lp <- log_prior_uniform_n(n, p, prior_upper)
  if (!is.finite(lp)) return(-Inf)

  lp <- lp +
    0.5 * (n - p - 1) * prep$sum_logdet_X -
    0.5 * m * n * p * log(2) -
    0.5 * m * n * logdet_spd(V) -
    m * log_multigamma(n / 2, p)

  lp
}

rw_update_n <- function(n_current, prep, V,
                        delta,
                        prior_upper = 1e4) {

  p <- prep$p
  lower <- p - 1

  n_proposal <- n_current + stats::runif(1, -delta, delta)

  if (n_proposal <= lower || n_proposal >= prior_upper) {
    return(list(n = n_current, accepted = FALSE))
  }

  lp_current <- log_posterior_n_given_V(
    n_current, prep, V, prior_upper
  )
  lp_proposal <- log_posterior_n_given_V(
    n_proposal, prep, V, prior_upper
  )

  log_alpha <- lp_proposal - lp_current

  if (log(stats::runif(1)) < min(0, log_alpha)) {
    list(n = n_proposal, accepted = TRUE)
  } else {
    list(n = n_current, accepted = FALSE)
  }
}

summarize_mcmc_n <- function(draws) {
  qs <- stats::quantile(draws, probs = c(0.025, 0.5, 0.975), names = FALSE)
  dens <- stats::density(draws, n = 512)
  mode <- dens$x[which.max(dens$y)]

  list(
    estimate = mean(draws),
    mean = mean(draws),
    median = qs[2],
    mode = mode,
    sd = stats::sd(draws),
    lower = qs[1],
    upper = qs[3]
  )
}

fit_rwm <- function(X,
                    iter = 4000,
                    burn = 1000,
                    delta = NULL,
                    prior_upper = 1e3,
                    init_n = NULL,
                    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  if (iter <= burn) stop("iter must be larger than burn.")

  prep <- prepare_wishart_data(X)
  p <- prep$p
  m <- prep$m
  lower <- p - 1

  if (is.null(init_n)) {
    mle_fit <- fit_mle(X)
    init_n <- mle_fit$n$estimate
    if (!is.finite(init_n) || init_n <= lower) init_n <- p + 1
  }

  init_n <- min(max(init_n, lower + 1e-6), prior_upper - 1e-6)

  # A scale for the random walk that grows mildly with n.
  if (is.null(delta)) {
    delta <- max(1, 0.20 * init_n)
  }

  A <- prep$S + diag(1e-4, p)
  n_chain <- numeric(iter)
  V_chain <- array(NA_real_, dim = c(p, p, iter))

  n_chain[1] <- init_n
  V_chain[, , 1] <- sample_inverse_wishart(
    df = prep$m * n_chain[1] + p,
    scale = A
  )

  accepted <- 0L
  start_time <- proc.time()[["elapsed"]]

  for (b in 2:iter) {

    # Gibbs update for V | n, X
    V_chain[, , b] <- sample_inverse_wishart(
      df = prep$m * n_chain[b - 1] + p,
      scale = A
    )

    # Random-Walk Metropolis update for n | V, X
    step <- rw_update_n(
      n_current = n_chain[b - 1],
      prep = prep,
      V = V_chain[, , b],
      delta = delta,
      prior_upper = prior_upper
    )

    n_chain[b] <- step$n
    accepted <- accepted + as.integer(step$accepted)
  }

  runtime <- proc.time()[["elapsed"]] - start_time

  keep <- seq.int(burn + 1, iter)
  n_post <- n_chain[keep]
  V_post <- V_chain[, , keep, drop = FALSE]

  V_mean <- apply(V_post, c(1, 2), mean)
  V_mean <- (V_mean + t(V_mean)) / 2

  list(
    method = "RWM",
    n = summarize_mcmc_n(n_post),
    V = list(
      estimate = V_mean,
      mean = V_mean,
      draws = ncol(matrix(V_post, nrow = p * p))
    ),
    n_draws = n_post,
    V_draws = V_post,
    runtime = runtime,
    diagnostics = list(
      acceptance_rate = accepted / (iter - 1),
      iter = iter,
      burn = burn,
      delta = delta
    )
  )
}

# ============================================================
# 3. Collapsed Bayesian inference
# ============================================================

log_collapsed_posterior_n <- function(n, prep,
                                      n_v = NULL,
                                      U = NULL,
                                      prior_upper = 1e3) {

  p <- prep$p
  m <- prep$m

  if (is.null(n_v)) n_v <- p
  if (is.null(U)) U <- diag(1e-4, p)

  lp <- log_prior_uniform_n(n, p, prior_upper)
  if (!is.finite(lp) || n <= p - 1) return(-Inf)

  A <- prep$S + U
  logdet_A <- logdet_spd(A)

  # IMPORTANT: the 2^{-mnp/2} term from the likelihood cancels
  # exactly with the corresponding factor from the inverse-Wishart
  # integral. Therefore it does NOT appear in the collapsed posterior.
  lp <- lp +
    log_multigamma((n_v + m * n) / 2, p) -
    m * log_multigamma(n / 2, p) -
    0.5 * (n_v + m * n) * logdet_A +
    0.5 * (n - p - 1) * prep$sum_logdet_X

  lp
}

trapezoid_weights <- function(x) {
  dx <- diff(x)
  w <- numeric(length(x))
  w[1] <- dx[1] / 2
  w[length(x)] <- dx[length(dx)] / 2
  if (length(x) > 2) {
    w[2:(length(x) - 1)] <- (dx[-length(dx)] + dx[-1]) / 2
  }
  w
}

grid_quantile <- function(x, density_values, q) {
  w <- trapezoid_weights(x)
  mass <- density_values * w
  cdf <- cumsum(mass)
  cdf <- cdf / max(cdf)

  if (q <= cdf[1]) return(x[1])
  if (q >= cdf[length(cdf)]) return(x[length(x)])

  approx(cdf, x, xout = q, ties = "ordered")$y
}

make_collapsed_grid <- function(prep,
                                n_v,
                                U,
                                prior_upper,
                                grid_size,
                                n_upper_initial = NULL) {

  p <- prep$p
  lower <- p - 1 + 1e-7

  if (is.null(n_upper_initial)) {
    mle_try <- tryCatch(fit_mle(prep$X), error = function(e) NULL)
    n_guess <- if (!is.null(mle_try)) mle_try$n$estimate else 2 * p
    if (!is.finite(n_guess)) n_guess <- 2 * p

    n_upper_initial <- max(10 * p, 50, 2 * n_guess + 10)
  }

  upper <- min(prior_upper - 1e-7, n_upper_initial)

  repeat {
    grid <- seq(lower, upper, length.out = grid_size)
    lp <- vapply(
      grid,
      log_collapsed_posterior_n,
      numeric(1),
      prep = prep,
      n_v = n_v,
      U = U,
      prior_upper = prior_upper
    )

    mx <- max(lp)
    tail_lp <- max(lp[max(1, floor(0.95 * length(lp))):length(lp)])

    # If the posterior is visibly close to the upper boundary, expand.
    if ((tail_lp > mx - 25) &&
        upper < prior_upper - 1e-6) {

      new_upper <- min(prior_upper - 1e-7, 2 * upper)

      if (new_upper <= upper + 1e-8) break
      upper <- new_upper

    } else {
      break
    }
  }

  list(grid = grid, logpost = lp)
}

fit_collapsed <- function(X,
                          n_v = NULL,
                          U = NULL,
                          prior_upper = 1e3,
                          grid_size = 1500,
                          n_draws = 0,
                          seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  prep <- prepare_wishart_data(X)
  p <- prep$p

  if (is.null(n_v)) n_v <- p
  if (is.null(U)) U <- diag(1e-4, p)

  start_time <- proc.time()[["elapsed"]]

  grid_obj <- make_collapsed_grid(
    prep = prep,
    n_v = n_v,
    U = U,
    prior_upper = prior_upper,
    grid_size = grid_size
  )

  n_grid <- grid_obj$grid
  logpost <- grid_obj$logpost

  max_lp <- max(logpost)
  unnorm <- exp(logpost - max_lp)
  weights <- trapezoid_weights(n_grid)
  mass <- unnorm * weights
  mass <- mass / sum(mass)

  n_mean <- sum(n_grid * mass)
  n_sd <- sqrt(sum((n_grid - n_mean)^2 * mass))
  n_median <- grid_quantile(n_grid, unnorm, 0.5)
  n_lower <- grid_quantile(n_grid, unnorm, 0.025)
  n_upper <- grid_quantile(n_grid, unnorm, 0.975)
  n_mode <- n_grid[which.max(logpost)]

  A <- prep$S + U

  # E[V | X] = E_n[ E(V | n, X) ]
  # and E[V | n, X] = A / (n_v + m*n - p - 1).
  denominator <- n_v + prep$m * n_grid - p - 1
  v_scalar_mean <- sum((1 / denominator) * mass)

  V_mean <- A * v_scalar_mean
  V_mean <- (V_mean + t(V_mean)) / 2

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
      V_post_draws[, , b] <- sample_inverse_wishart(
        df = n_v + prep$m * n_post_draws[b],
        scale = A
      )
    }
  }

  runtime <- proc.time()[["elapsed"]] - start_time

  list(
    method = "Collapsed",
    n = list(
      estimate = n_mean,
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
      posterior_at_upper = tail(unnorm, 1)
    )
  )
}
