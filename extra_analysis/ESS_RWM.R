# ============================================================
# RWM ESS diagnostic simulation
#
# Purpose:
#   Assess RWM acceptance rate and effective sample size (ESS)
#   in representative easy, moderate, and difficult scenarios.
#
# The collapsed method is not run here because ESS is relevant
# only to the MCMC validation sampler.
#
# Assumes that fit_rwm(), prepare_wishart_data(), and related
# functions are already available from wishart_df_methods.R.
# ============================================================


# ------------------------------------------------------------
# Load methods
# ------------------------------------------------------------

source("wishart_df_methods.R")


# ------------------------------------------------------------
# Global settings
# ------------------------------------------------------------

set.seed(20260929)

R <- 100

rwm_iter <- 10000
rwm_burn <- 4000
rwm_delta <- 2

prior_upper <- 1e3

n_v <- NULL
U_scale <- 1e-4


# ------------------------------------------------------------
# Representative scenarios
# ------------------------------------------------------------

scenario_grid <- data.frame(
  
  scenario = c(
    "Easy",
    "Moderate",
    "Difficult"
  ),
  
  p = c(
    5,
    10,
    30
  ),
  
  m = c(
    200,
    50,
    10
  ),
  
  n_true = c(
    6,      # p + 1
    40,     # 4p
    240     # 8p
  ),
  
  stringsAsFactors = FALSE
)

scenario_grid$m_ratio <-
  scenario_grid$m /
  scenario_grid$p

scenario_grid$n_ratio <-
  scenario_grid$n_true /
  scenario_grid$p


# ------------------------------------------------------------
# Generate random non-spherical scale matrix
# ------------------------------------------------------------

generate_random_scale <- function(
    p,
    min = -10,
    max = 10) {
  
  B <- matrix(
    runif(
      p * p,
      min = min,
      max = max
    ),
    nrow = p,
    ncol = p
  )
  
  V <- crossprod(B)
  
  V <- V /
    mean(diag(V))
  
  V <- (
    V + t(V)
  ) / 2
  
  # Numerical protection.
  ev <- eigen(
    V,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  if (min(ev) <= 1e-10) {
    
    V <- V +
      diag(1e-6, p)
    
    V <- V /
      mean(diag(V))
  }
  
  V
}


# ------------------------------------------------------------
# Simulate Wishart observations
# ------------------------------------------------------------

simulate_wishart_dataset <- function(
    p,
    n0,
    m,
    V_true) {
  
  X_array <- stats::rWishart(
    n = m,
    df = n0,
    Sigma = V_true
  )
  
  lapply(
    seq_len(m),
    function(i) {
      X_array[, , i]
    }
  )
}


# ------------------------------------------------------------
# Safe mean
# ------------------------------------------------------------

mean_safe <- function(x) {
  
  if (all(is.na(x))) {
    return(NA_real_)
  }
  
  mean(
    x,
    na.rm = TRUE
  )
}


# ------------------------------------------------------------
# Safe median
# ------------------------------------------------------------

median_safe <- function(x) {
  
  if (all(is.na(x))) {
    return(NA_real_)
  }
  
  median(
    x,
    na.rm = TRUE
  )
}


# ============================================================
# Run ESS simulation
# ============================================================

all_results <- list()

result_id <- 0L


cat(
  "============================================================\n"
)

cat(
  "RWM ESS diagnostic simulation\n"
)

cat(
  "============================================================\n"
)

cat(
  "Scenarios       :", nrow(scenario_grid), "\n"
)

cat(
  "Replications    :", R, "\n"
)

cat(
  "Iterations      :", rwm_iter, "\n"
)

cat(
  "Burn-in         :", rwm_burn, "\n"
)

cat(
  "Proposal width  :", rwm_delta, "\n"
)

cat(
  "============================================================\n\n"
)


for (s in seq_len(nrow(scenario_grid))) {
  
  # ----------------------------------------------------------
  # Scenario parameters
  # ----------------------------------------------------------
  
  scenario_name <-
    scenario_grid$scenario[s]
  
  p <-
    scenario_grid$p[s]
  
  m <-
    scenario_grid$m[s]
  
  n0 <-
    scenario_grid$n_true[s]
  
  
  cat(
    sprintf(
      "\n------------------------------------------------------------\n"
    )
  )
  
  cat(
    sprintf(
      "Scenario: %s | p=%d | m=%d | n=%d | m/p=%.2f | n/p=%.2f\n",
      scenario_name,
      p,
      m,
      n0,
      m / p,
      n0 / p
    )
  )
  
  cat(
    "------------------------------------------------------------\n"
  )
  
  
  # ----------------------------------------------------------
  # One common V_true per scenario
  # ----------------------------------------------------------
  
  V_true <- generate_random_scale(
    p = p
  )
  
  
  # ----------------------------------------------------------
  # Replications
  # ----------------------------------------------------------
  
  for (r in seq_len(R)) {
    
    cat(
      sprintf(
        "  Replication %3d/%3d\r",
        r,
        R
      )
    )
    
    
    # --------------------------------------------------------
    # Generate dataset
    # --------------------------------------------------------
    
    X <- simulate_wishart_dataset(
      p = p,
      n0 = n0,
      m = m,
      V_true = V_true
    )
    
    
    # --------------------------------------------------------
    # RWM
    # --------------------------------------------------------
    
    fit <- fit_rwm(
      
      X = X,
      
      iter = rwm_iter,
      
      burn = rwm_burn,
      
      delta = rwm_delta,
      
      prior_upper = prior_upper,
      
      init_n = 4 * p,
      
      n_v = if (is.null(n_v)) {
        p
      } else {
        n_v
      },
      
      U = diag(
        U_scale,
        p
      ),
      
      seed = NULL
    )
    
    
    # --------------------------------------------------------
    # Store diagnostics
    # --------------------------------------------------------
    
    result_id <- result_id + 1L
    
    all_results[[result_id]] <-
      data.frame(
        
        scenario =
          scenario_name,
        
        p = p,
        
        m = m,
        
        n_true = n0,
        
        m_ratio =
          m / p,
        
        n_ratio =
          n0 / p,
        
        replicate =
          r,
        
        acceptance_rate =
          fit$diagnostics$acceptance_rate,
        
        effective_sample_size =
          fit$diagnostics$effective_sample_size,
        
        ess_fraction =
          fit$diagnostics$ess_fraction,
        
        post_burn_in_draws =
          fit$diagnostics$post_burn_in_draws,
        
        runtime =
          fit$runtime,
        
        n_mode =
          fit$n$mode,
        
        n_sd =
          fit$n$sd,
        
        n_lower =
          fit$n$lower,
        
        n_upper =
          fit$n$upper,
        
        stringsAsFactors = FALSE
      )
  }
  
  cat("\n")
}


# ------------------------------------------------------------
# Combine results
# ------------------------------------------------------------

rwm_ess_results <- do.call(
  rbind,
  all_results
)


# ------------------------------------------------------------
# Save replication-level results
# ------------------------------------------------------------

write.csv(
  rwm_ess_results,
  "RWM_ESS_diagnostic_replications.csv",
  row.names = FALSE
)


# ============================================================
# Summary by scenario
# ============================================================

summary_results <- do.call(
  
  rbind,
  
  lapply(
    split(
      rwm_ess_results,
      rwm_ess_results$scenario
    ),
    
    function(d) {
      
      data.frame(
        
        scenario =
          d$scenario[1],
        
        p =
          d$p[1],
        
        m =
          d$m[1],
        
        n_true =
          d$n_true[1],
        
        m_ratio =
          d$m_ratio[1],
        
        n_ratio =
          d$n_ratio[1],
        
        mean_acceptance_rate =
          mean(
            d$acceptance_rate,
            na.rm = TRUE
          ),
        
        sd_acceptance_rate =
          sd(
            d$acceptance_rate,
            na.rm = TRUE
          ),
        
        mean_ess =
          mean(
            d$effective_sample_size,
            na.rm = TRUE
          ),
        
        median_ess =
          median(
            d$effective_sample_size,
            na.rm = TRUE
          ),
        
        sd_ess =
          sd(
            d$effective_sample_size,
            na.rm = TRUE
          ),
        
        mean_ess_fraction =
          mean(
            d$ess_fraction,
            na.rm = TRUE
          ),
        
        mean_runtime =
          mean(
            d$runtime,
            na.rm = TRUE
          ),
        
        median_runtime =
          median(
            d$runtime,
            na.rm = TRUE
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)


# ------------------------------------------------------------
# Save summary
# ------------------------------------------------------------

write.csv(
  summary_results,
  "RWM_ESS_diagnostic_summary.csv",
  row.names = FALSE
)


# ============================================================
# Print summary
# ============================================================

cat(
  "\n\n============================================================\n"
)

cat(
  "RWM ESS diagnostic summary\n"
)

cat(
  "============================================================\n\n"
)

print(
  summary_results,
  row.names = FALSE
)
