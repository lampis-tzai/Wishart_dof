# ============================================================
# wishart_df_simulations_parallel.R
# ============================================================

rm(list = ls())

library(parallel)

source("wishart_df_methods.R")

# ============================================================
# CONFIG
# ============================================================

n_cores <- max(1L, parallel::detectCores() - 1L)

cat("Using", n_cores, "parallel cores.\n")

SEED <- 20260918

p_values <- c(3, 5, 10, 20, 30)
ratio_values <- c(1, 2, 4, 8)
m_values <- c(10, 20, 50, 100, 200)

# Use 3 for a test run; change to 1000 for the final run.
R_main <- 1000

rwm_iter <- 4000
rwm_burn <- 1000

prior_upper <- 1e3
n_v <- NULL
U_scale <- 1e-4

collapsed_grid_size <- 1500

output_results <- "wishart_df_simulation_results.csv"
output_raw <- "wishart_df_simulation_replications.csv"

SAVE_RAW <- TRUE

stopifnot(
  R_main >= 1,
  all(p_values > 0),
  all(m_values > 0),
  all(ratio_values > 0)
)

# ============================================================
# Helper functions
# ============================================================

simulate_wishart_dataset <- function(p, n0, m, V_true) {
  X_array <- stats::rWishart(
    n = m,
    df = n0,
    Sigma = V_true
  )

  lapply(seq_len(m), function(i) X_array[, , i])
}

stein_loss <- function(V_hat, V_true, V_true_inv = NULL) {
  p <- nrow(V_true)

  if (is.null(V_true_inv)) {
    V_true_inv <- solve(V_true)
  }

  V_hat <- (V_hat + t(V_hat)) / 2

  trace_term <- sum(diag(V_true_inv %*% V_hat))

  logdet_ratio <-
    logdet_spd(V_hat) -
    logdet_spd(V_true)

  loss <- trace_term - logdet_ratio - p

  max(loss, 0)
}

coverage_95 <- function(lower, upper, truth) {
  as.integer(
    is.finite(lower) &&
    is.finite(upper) &&
    truth >= lower &&
    truth <= upper
  )
}

median_safe <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

mean_safe <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

aggregate_method_results <- function(rep_df) {
  truth <- unique(rep_df$n_true)

  if (length(truth) != 1L) {
    stop("aggregate_method_results() received multiple true n values.")
  }

  coverage_value <- if (all(is.na(rep_df$coverage_95))) {
    NA_real_
  } else {
    mean(rep_df$coverage_95, na.rm = TRUE)
  }

  data.frame(
    n_rep = nrow(rep_df),

    bias =
      mean(rep_df$n_est - rep_df$n_true, na.rm = TRUE),

    relative_bias =
      mean(rep_df$n_est - rep_df$n_true, na.rm = TRUE) / truth,

    rmse =
      sqrt(
        mean(
          (rep_df$n_est - rep_df$n_true)^2,
          na.rm = TRUE
        )
      ),

    posterior_sd = mean_safe(rep_df$n_sd),

    coverage_95 = coverage_value,

    stein_loss = mean(rep_df$stein_loss, na.rm = TRUE),

    median_runtime = median_safe(rep_df$runtime),

    mean_runtime = mean(rep_df$runtime, na.rm = TRUE),

    mean_acceptance_rate = mean_safe(rep_df$acceptance_rate)
  )
}

scenario_n <- function(p, ratio) {
  if (ratio == 1) return(p + 1)
  ratio * p
}

# ============================================================
# Worker function
# ============================================================

run_one_replication <- function(
    r,
    p,
    n0,
    m,
    V_true,
    V_true_inv,
    scenario_name,
    n_v,
    U_scale,
    prior_upper,
    collapsed_grid_size,
    rwm_iter,
    rwm_burn,
    replication_seed) {

  set.seed(replication_seed)

  X <- simulate_wishart_dataset(
    p = p,
    n0 = n0,
    m = m,
    V_true = V_true
  )

  mle_fit <- fit_mle(X)

  collapsed_fit <- fit_collapsed(
    X = X,
    n_v = if (is.null(n_v)) p else n_v,
    U = diag(U_scale, p),
    prior_upper = prior_upper,
    grid_size = collapsed_grid_size,
    n_draws = 0
  )

  rwm_fit <- fit_rwm(
    X = X,
    iter = rwm_iter,
    burn = rwm_burn,
    prior_upper = prior_upper
  )

  rbind(
    data.frame(
      scenario = scenario_name,
      p = p,
      m = m,
      n_true = n0,
      n_ratio = n0 / p,
      m_ratio = m / p,
      replicate = r,
      method = "MLE",
      n_est = mle_fit$n$estimate,
      n_sd = NA_real_,
      n_lower = NA_real_,
      n_upper = NA_real_,
      coverage_95 = NA_real_,
      stein_loss = stein_loss(
        mle_fit$V$estimate,
        V_true,
        V_true_inv
      ),
      runtime = mle_fit$runtime,
      acceptance_rate = NA_real_,
      stringsAsFactors = FALSE
    ),

    data.frame(
      scenario = scenario_name,
      p = p,
      m = m,
      n_true = n0,
      n_ratio = n0 / p,
      m_ratio = m / p,
      replicate = r,
      method = "Collapsed",
      n_est = collapsed_fit$n$mode,
      n_sd = collapsed_fit$n$sd,
      n_lower = collapsed_fit$n$lower,
      n_upper = collapsed_fit$n$upper,
      coverage_95 = coverage_95(
        collapsed_fit$n$lower,
        collapsed_fit$n$upper,
        n0
      ),
      stein_loss = stein_loss(
        collapsed_fit$V$estimate,
        V_true,
        V_true_inv
      ),
      runtime = collapsed_fit$runtime,
      acceptance_rate = NA_real_,
      stringsAsFactors = FALSE
    ),

    data.frame(
      scenario = scenario_name,
      p = p,
      m = m,
      n_true = n0,
      n_ratio = n0 / p,
      m_ratio = m / p,
      replicate = r,
      method = "RWM",
      n_est = rwm_fit$n$mode,
      n_sd = rwm_fit$n$sd,
      n_lower = rwm_fit$n$lower,
      n_upper = rwm_fit$n$upper,
      coverage_95 = coverage_95(
        rwm_fit$n$lower,
        rwm_fit$n$upper,
        n0
      ),
      stein_loss = stein_loss(
        rwm_fit$V$estimate,
        V_true,
        V_true_inv
      ),
      runtime = rwm_fit$runtime,
      acceptance_rate = rwm_fit$diagnostics$acceptance_rate,
      stringsAsFactors = FALSE
    )
  )
}

# ============================================================
# Parallel worker wrapper
# ============================================================

run_replication_job <- function(
    job,
    p,
    n0,
    m,
    V_true,
    V_true_inv,
    scenario_name,
    n_v,
    U_scale,
    prior_upper,
    collapsed_grid_size,
    rwm_iter,
    rwm_burn) {

  run_one_replication(
    r = job$r,
    p = p,
    n0 = n0,
    m = m,
    V_true = V_true,
    V_true_inv = V_true_inv,
    scenario_name = scenario_name,
    n_v = n_v,
    U_scale = U_scale,
    prior_upper = prior_upper,
    collapsed_grid_size = collapsed_grid_size,
    rwm_iter = rwm_iter,
    rwm_burn = rwm_burn,
    replication_seed = job$seed
  )
}

# ============================================================
# Scenario grid
# ============================================================

set.seed(SEED)

scenario_grid <- expand.grid(
  p = p_values,
  ratio = ratio_values,
  m = m_values,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

n_scenarios <- nrow(scenario_grid)

all_summary <- list()
all_raw <- list()

summary_id <- 0L
raw_id <- 0L

cat("============================================================\n")
cat("Wishart degrees-of-freedom simulation\n")
cat("============================================================\n")
cat("Number of scenarios    :", n_scenarios, "\n")
cat("Replications/scenario  :", R_main, "\n")
cat("Methods                : MLE, Collapsed, RWM\n")
cat("RWM iterations         :", rwm_iter, "\n")
cat("RWM burn-in            :", rwm_burn, "\n\n")

# ============================================================
# Create cluster ONCE
# ============================================================

cl <- parallel::makeCluster(
  n_cores,
  type = "PSOCK"
)

methods_file <- normalizePath(
  "wishart_df_methods.R",
  mustWork = TRUE
)

parallel::clusterExport(
  cl,
  varlist = "methods_file",
  envir = environment()
)

parallel::clusterEvalQ(
  cl,
  source(methods_file)
)

parallel::clusterExport(
  cl,
  varlist = c(
    "simulate_wishart_dataset",
    "stein_loss",
    "coverage_95",
    "run_one_replication",
    "run_replication_job"
  ),
  envir = environment()
)

cat(
  "Parallel cluster initialized using",
  n_cores,
  "cores.\n\n"
)

# ============================================================
# Main scenario loop
# ============================================================

for (s in seq_len(n_scenarios)) {

  p <- scenario_grid$p[s]
  ratio <- scenario_grid$ratio[s]
  m <- scenario_grid$m[s]

  n0 <- scenario_n(p, ratio)

  scenario_name <- paste0(
    "p", p,
    "_m", m,
    "_n", n0
  )

  cat(
    sprintf(
      "[%3d/%3d] p=%d, m=%d, n=%d, n/p=%.2f, m/p=%.2f\n",
      s,
      n_scenarios,
      p,
      m,
      n0,
      n0 / p,
      m / p
    )
  )

  V_true <- generate_random_scale(p)
  V_true_inv <- solve(V_true)

  set.seed(SEED + s)

  replication_seeds <- sample.int(
    .Machine$integer.max,
    size = R_main,
    replace = FALSE
  )

  # Build an explicit job list. The seed travels with each job,
  # so workers do not need access to replication_seeds.
  replication_jobs <- lapply(
    seq_len(R_main),
    function(r) {
      list(
        r = r,
        seed = replication_seeds[r]
      )
    }
  )

  scenario_replications <- parallel::parLapplyLB(
    cl,
    X = replication_jobs,
    fun = run_replication_job,

    p = p,
    n0 = n0,
    m = m,

    V_true = V_true,
    V_true_inv = V_true_inv,

    scenario_name = scenario_name,

    n_v = n_v,
    U_scale = U_scale,
    prior_upper = prior_upper,

    collapsed_grid_size =
      collapsed_grid_size,

    rwm_iter = rwm_iter,
    rwm_burn = rwm_burn
  )

  scenario_rep_df <- do.call(
    rbind,
    scenario_replications
  )

  by_method <- split(
    scenario_rep_df,
    scenario_rep_df$method
  )

  for (method_name in names(by_method)) {

    d <- by_method[[method_name]]

    summary_id <- summary_id + 1L

    agg <- aggregate_method_results(d)

    all_summary[[summary_id]] <-
      data.frame(
        scenario = scenario_name,
        p = p,
        m = m,
        n_true = n0,
        n_ratio = n0 / p,
        m_ratio = m / p,
        method = method_name,
        agg,
        stringsAsFactors = FALSE
      )
  }

  if (SAVE_RAW) {
    raw_id <- raw_id + 1L
    all_raw[[raw_id]] <- scenario_rep_df
  }

  cat("    Scenario completed.\n\n")
}

# ============================================================
# Stop cluster ONCE
# ============================================================

parallel::stopCluster(cl)

cat("\nParallel cluster stopped.\n")

# ============================================================
# Save results
# ============================================================

results <- do.call(
  rbind,
  all_summary
)

write.csv(
  results,
  output_results,
  row.names = FALSE
)

if (SAVE_RAW) {

  raw_results <- do.call(
    rbind,
    all_raw
  )

  write.csv(
    raw_results,
    output_raw,
    row.names = FALSE
  )

} else {

  raw_results <- NULL
}

# ============================================================
# Sanity checks
# ============================================================

expected_summary_rows <- n_scenarios * 3

expected_raw_rows <- if (SAVE_RAW) {
  n_scenarios * R_main * 3
} else {
  NA_integer_
}

cat("\n============================================================\n")
cat("Sanity checks\n")
cat("============================================================\n")

cat(
  "Expected summary rows :",
  expected_summary_rows,
  "\n"
)

cat(
  "Actual summary rows   :",
  nrow(results),
  "\n"
)

if (SAVE_RAW) {

  cat(
    "Expected raw rows     :",
    expected_raw_rows,
    "\n"
  )

  cat(
    "Actual raw rows       :",
    nrow(raw_results),
    "\n"
  )
}

if (nrow(results) != expected_summary_rows) {
  stop(
    "Summary row count is incorrect. ",
    "Expected ",
    expected_summary_rows,
    ", got ",
    nrow(results),
    "."
  )
}

if (
  SAVE_RAW &&
  nrow(raw_results) != expected_raw_rows
) {
  stop(
    "Raw row count is incorrect. ",
    "Expected ",
    expected_raw_rows,
    ", got ",
    nrow(raw_results),
    "."
  )
}

cat("Simulation completed successfully.\n")




# ============================================================
# Plotting
# ============================================================
#
# The plots summarize the effect of m/p for each n/p regime.
# For readability, results are averaged over p within each
# (method, m/p, n/p) combination.
#
# The CSV files contain the complete scenario-level results.
# ============================================================

results <- read.csv("wishart_df_simulation_results.csv")

plot_data <- aggregate(
  cbind(
    bias,
    relative_bias,
    rmse,
    stein_loss,
    median_runtime
  ) ~ method + m_ratio + n_ratio,
  data = results,
  FUN = mean,
  na.rm = TRUE
)


make_panel_plot <- function(
    plot_data,
    yvar,
    ylab,
    filename,
    log_y = FALSE
) {
  
  methods <- c("RWM", "Collapsed")
  
  ## Appearance
  method_col <- c(
    RWM = "#2C7FB8",
    Collapsed = "#D95F02"
  )
  
  method_pch <- c(
    RWM = 16,
    Collapsed = 17
  )
  
  method_lty <- c(
    RWM = 1,
    Collapsed = 2
  )
  
  method_lwd <- c(
    RWM = 2.2,
    Collapsed = 2.2
  )
  
  n_ratios <- sort(unique(plot_data$n_ratio))
  
  png(
    filename,
    width = 2200,
    height = 1500,
    res = 200
  )
  
  ## 2 x 2 panels
  old_par <- par(
    mfrow = c(2, 2),
    mar = c(4.3, 4.5, 3.0, 1.2),
    oma = c(0, 0, 1.5, 0)
  )
  
  for (ratio_value in n_ratios) {
    
    d_ratio <- plot_data[
      plot_data$n_ratio == ratio_value,
      ,
      drop = FALSE
    ]
    
    d_ratio <- d_ratio[
      is.finite(d_ratio$m_ratio) &
        is.finite(d_ratio[[yvar]]),
      ,
      drop = FALSE
    ]
    
    d_ratio <- d_ratio[
      order(d_ratio$m_ratio),
      ,
      drop = FALSE
    ]
    
    ## ----------------------------
    ## Panel-specific y-axis
    ## ----------------------------
    
    y <- d_ratio[[yvar]]
    
    if (log_y) {
      
      y <- y[y > 0 & is.finite(y)]
      
      ylim <- range(y)
      
      ## Add some space
      ylim <- ylim * c(0.85, 1.15)
      
    } else {
      
      ylim <- range(
        y,
        finite = TRUE
      )
      
      ## Include zero for bias
      if (yvar == "bias") {
        ylim <- range(c(ylim, 0))
      }
      
      ## Add padding
      pad <- 0.08 * diff(ylim)
      
      if (pad == 0) {
        pad <- 0.1
      }
      
      ylim <- ylim + c(-pad, pad)
    }
    
    ## ----------------------------
    ## X-axis
    ## ----------------------------
    
    xlim <- range(
      d_ratio$m_ratio,
      finite = TRUE
    )
    
    ## Use logarithmic x-axis because m/p
    ## spans a wide range
    xlim <- xlim * c(0.9, 1.1)
    
    ## ----------------------------
    ## Empty plot
    ## ----------------------------
    
    plot(
      NA,
      NA,
      xlim = xlim,
      ylim = ylim,
      log = paste0(
        "x",
        if (log_y) "y" else ""
      ),
      xaxt = "n",
      xlab = "m / p",
      ylab = ylab,
      main = paste0(
        "n / p = ",
        formatC(
          ratio_value,
          format = "f",
          digits = 2
        )
      ),
      cex.main = 1.15,
      font.main = 2,
      cex.lab = 1.05,
      cex.axis = 0.9
    )
    
    ## Custom x-axis with one decimal place
    x_ticks <- axTicks(1)
    
    axis(
      1,
      at = x_ticks,
      labels = formatC(
        x_ticks,
        format = "f",
        digits = 1
      ),
      cex.axis = 0.9
    )
    
    ## Light grid
    grid(
      nx = 6,
      ny = 5,
      col = "grey90",
      lty = 1
    )
    
    ## Zero reference line for bias
    if (yvar == "bias" && !log_y) {
      abline(
        h = 0,
        lty = 3,
        lwd = 1.3,
        col = "grey30"
      )
    }
    
    ## ----------------------------
    ## Plot methods
    ## ----------------------------
    
    for (method_name in methods) {
      
      d <- d_ratio[
        d_ratio$method == method_name,
        ,
        drop = FALSE
      ]
      
      if (nrow(d) == 0) next
      
      d <- d[
        order(d$m_ratio),
        ,
        drop = FALSE
      ]
      
      lines(
        d$m_ratio,
        d[[yvar]],
        type = "o",
        pch = method_pch[method_name],
        lty = method_lty[method_name],
        lwd = method_lwd[method_name],
        col = method_col[method_name],
        cex = 1.0
      )
    }
    
    ## ----------------------------
    ## Legend only once
    ## ----------------------------
    
    if (ratio_value == max(n_ratios)) {
      
      legend(
        "topright",
        legend = methods,
        col = method_col,
        pch = method_pch,
        lty = method_lty,
        lwd = method_lwd,
        pt.cex = 1.0,
        bty = "n",
        cex = 0.95
      )
    }
  }
  
  ## Overall title
  mtext(
    paste0(
      "Simulation results: ",
      ylab
    ),
    side = 3,
    outer = TRUE,
    line = 0.2,
    font = 2,
    cex = 1.25
  )
  
  par(old_par)
  
  dev.off()
}


output_results <- "wishart_df_simulation_results.csv"
output_raw <- "wishart_df_simulation_replications.csv"

output_bias_plot <- "wishart_df_bias.png"
output_rmse_plot <- "wishart_df_rmse.png"
output_stein_plot <- "wishart_df_stein_loss.png"
output_runtime_plot <- "wishart_df_runtime.png"


# ------------------------------------------------------------
# Bias
# ------------------------------------------------------------

make_panel_plot(
  plot_data = plot_data,
  yvar = "bias",
  ylab = "Bias",
  filename = output_bias_plot
)


# ------------------------------------------------------------
# RMSE
# ------------------------------------------------------------

make_panel_plot(
  plot_data = plot_data,
  yvar = "rmse",
  ylab = "RMSE",
  filename = output_rmse_plot
)


# ------------------------------------------------------------
# Stein loss
# ------------------------------------------------------------

make_panel_plot(
  plot_data = plot_data,
  yvar = "stein_loss",
  ylab = "Mean Stein loss for V",
  filename = output_stein_plot
)


# ------------------------------------------------------------
# Runtime
# ------------------------------------------------------------

make_panel_plot(
  plot_data = plot_data,
  yvar = "median_runtime",
  ylab = "Median runtime per dataset (seconds)",
  filename = output_runtime_plot,
  log_y = TRUE
)


cat("\nOutput files:\n")
cat("  ", output_results, "\n")
if (SAVE_RAW) cat("  ", output_raw, "\n")
cat("  ", output_bias_plot, "\n")
cat("  ", output_rmse_plot, "\n")
cat("  ", output_stein_plot, "\n")
cat("  ", output_runtime_plot, "\n")



library(dplyr)

results <- read.csv("wishart_df_simulation_results.csv")

general_result <- results %>%
  group_by(method) %>%
  summarise(
    mean_bias          = mean(bias, na.rm = TRUE),
    sd_bias            = sd(bias, na.rm = TRUE),
    
    mean_relative_bias = mean(relative_bias, na.rm = TRUE),
    sd_relative_bias   = sd(relative_bias, na.rm = TRUE),
    
    mean_rmse          = mean(rmse, na.rm = TRUE),
    sd_rmse            = sd(rmse, na.rm = TRUE),
    
    mean_posterior_sd  = mean(posterior_sd, na.rm = TRUE),
    sd_posterior_sd    = sd(posterior_sd, na.rm = TRUE),
    
    mean_coverage      = mean(coverage_95, na.rm = TRUE),
    sd_coverage        = sd(coverage_95, na.rm = TRUE),
    
    mean_stein_loss    = mean(stein_loss, na.rm = TRUE),
    sd_stein_loss      = sd(stein_loss, na.rm = TRUE),
    
    mean_runtime       = mean(mean_runtime, na.rm = TRUE),
    sd_runtime         = sd(mean_runtime, na.rm = TRUE),
    
    mean_acceptance    = mean(mean_acceptance_rate, na.rm = TRUE),
    sd_acceptance      = sd(mean_acceptance_rate, na.rm = TRUE),
    
    .groups = "drop"
  )

as.data.frame(general_result)
