# ============================================================
# wishart_df_simulations.R
#
# Simulation study for:
#   1. Maximum likelihood estimation (MLE)
#   2. Collapsed Bayesian inference
#   3. Random-Walk Metropolis within Gibbs (RWM)
#
# All three methods are run on THE SAME simulated dataset in every
# replication and in every simulation scenario.
#
# Main design:
#   p       = {3, 5, 10, 20, 30}
#   n/p     = {1, 2, 4, 8}
#   m       = {10, 20, 50, 100, 200}
#
# For n/p = 1 we use n0 = p + 1, so that n0 > p - 1.
#
# The main quantity of interest is inference for n. V is also evaluated
# as a secondary matrix-estimation target using Stein's loss.
#
# Before running:
#   1. Put wishart_df_methods.R in the same working directory.
#   2. Adjust R_main and the RWM settings in CONFIG.
# ============================================================

rm(list = ls())

source("wishart_df_methods.R")

# ============================================================
# CONFIG
# ============================================================

library(parallel)

# Number of CPU cores to use.
# Leave one core free for the operating system.
n_cores <- max(1L, parallel::detectCores() - 1L)

cat("Using", n_cores, "parallel cores.\n")

SEED <- 20260918

# Main simulation design
p_values <- c(3, 5, 10, 20, 30)
ratio_values <- c(1, 2, 4, 8)       # n/p values
m_values <- c(10, 20, 50, 100, 200)

# Change this to 1000 for the final simulation.
R_main <- 1000

# RWM settings
rwm_iter <- 4000
rwm_burn <- 1000

# Bayesian specification
prior_upper <- 1e3
n_v <- NULL                         # NULL means n_v = p
U_scale <- 1e-4                     # U = U_scale * I_p

# Numerical grid for collapsed posterior
collapsed_grid_size <- 1500

# Output files
output_results <- "wishart_df_simulation_results.csv"
output_raw <- "wishart_df_simulation_replications.csv"

output_bias_plot <- "wishart_df_bias.png"
output_rmse_plot <- "wishart_df_rmse.png"
output_stein_plot <- "wishart_df_stein_loss.png"
output_runtime_plot <- "wishart_df_runtime.png"

SAVE_RAW <- TRUE

# ============================================================
# Checks
# ============================================================

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

  lapply(
    seq_len(m),
    function(i) X_array[, , i]
  )
}


stein_loss <- function(V_hat, V_true, V_true_inv = NULL) {

  p <- nrow(V_true)

  if (is.null(V_true_inv)) {
    V_true_inv <- solve(V_true)
  }

  # Ensure symmetry before the log determinant.
  V_hat <- (V_hat + t(V_hat)) / 2

  trace_term <- sum(diag(V_true_inv %*% V_hat))

  logdet_ratio <- logdet_spd(V_hat) -
    logdet_spd(V_true)

  loss <- trace_term - logdet_ratio - p

  # Numerical protection only.
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

  if (all(is.na(x))) {
    return(NA_real_)
  }

  median(x, na.rm = TRUE)
}


mean_safe <- function(x) {

  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}


aggregate_method_results <- function(rep_df) {

  truth <- unique(rep_df$n_true)

  if (length(truth) != 1L) {
    stop("aggregate_method_results() received multiple true n values.")
  }

  coverage_value <- if (
    all(is.na(rep_df$coverage_95))
  ) {
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

    posterior_sd =
      mean_safe(rep_df$n_sd),

    coverage_95 =
      coverage_value,

    stein_loss =
      mean(rep_df$stein_loss, na.rm = TRUE),

    median_runtime =
      median_safe(rep_df$runtime),

    mean_runtime =
      mean(rep_df$runtime, na.rm = TRUE),

    mean_acceptance_rate =
      mean_safe(rep_df$acceptance_rate)
  )
}


scenario_n <- function(p, ratio) {

  if (ratio == 1) {
    return(p + 1)
  }

  ratio * p
}


# ============================================================
# Main simulation
# ============================================================

set.seed(SEED)

scenario_grid <- expand.grid(
  p = p_values,
  ratio = ratio_values,
  m = m_values,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

all_summary <- list()
all_raw <- list()

summary_id <- 0L
raw_id <- 0L

n_scenarios <- nrow(scenario_grid)

cat("============================================================\n")
cat("Wishart degrees-of-freedom simulation\n")
cat("============================================================\n")
cat("Number of scenarios :", n_scenarios, "\n")
cat("Replications/scenario:", R_main, "\n")
cat("Methods             : MLE, Collapsed, RWM\n")
cat("RWM iterations      :", rwm_iter, "\n")
cat("RWM burn-in         :", rwm_burn, "\n\n")


for (s in seq_len(n_scenarios)) {

  # ----------------------------------------------------------
  # Scenario parameters
  # ----------------------------------------------------------

  p <- scenario_grid$p[s]
  ratio <- scenario_grid$ratio[s]
  m <- scenario_grid$m[s]

  # n/p = 1 is represented by n = p + 1 so that n > p - 1.
  n0 <- scenario_n(p, ratio)

  scenario_name <- paste0(
    "p", p,
    "_m", m,
    "_n", n0
  )

  cat(
    sprintf(
      "[%3d/%3d]  p=%d, m=%d, n=%d, n/p=%.2f, m/p=%.2f\n",
      s,
      n_scenarios,
      p,
      m,
      n0,
      n0 / p,
      m / p
    )
  )

  # ----------------------------------------------------------
  # One random non-spherical V_true per scenario
  # ----------------------------------------------------------
  #
  # The same V_true is used for all replications in this scenario.
  # This creates a controlled repeated-sampling experiment under a
  # common underlying scale structure.
  # ----------------------------------------------------------

  V_true <- generate_random_scale(p)
  V_true_inv <- solve(V_true)

  scenario_replications <- vector(
    "list",
    R_main
  )


  # ==========================================================
  # Replications
  # ==========================================================

  for (r in seq_len(R_main)) {

    # --------------------------------------------------------
    # Generate ONE dataset
    # --------------------------------------------------------
    #
    # The exact same X is analyzed by MLE, Collapsed, and RWM.
    # This makes the method comparison paired.
    # --------------------------------------------------------

    X <- simulate_wishart_dataset(
      p = p,
      n0 = n0,
      m = m,
      V_true = V_true
    )


    # --------------------------------------------------------
    # 1. MLE
    # --------------------------------------------------------

    mle_fit <- fit_mle(X)


    # --------------------------------------------------------
    # 2. Collapsed Bayesian method
    # --------------------------------------------------------

    collapsed_fit <- fit_collapsed(
      X = X,
      n_v = if (is.null(n_v)) p else n_v,
      U = diag(U_scale, p),
      prior_upper = prior_upper,
      grid_size = collapsed_grid_size,
      n_draws = 0
    )


    # --------------------------------------------------------
    # 3. RWM within Gibbs
    # --------------------------------------------------------

    rwm_fit <- fit_rwm(
      X = X,
      iter = rwm_iter,
      burn = rwm_burn,
      prior_upper = prior_upper
    )


    # --------------------------------------------------------
    # Store results from all three methods
    # --------------------------------------------------------

    scenario_replications[[r]] <- rbind(

      # =========================
      # MLE
      # =========================

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


      # =========================
      # Collapsed
      # =========================

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


      # =========================
      # RWM
      # =========================

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

        acceptance_rate =
          rwm_fit$diagnostics$acceptance_rate,

        stringsAsFactors = FALSE
      )
    )
  }


  # ==========================================================
  # Combine replications for this scenario
  # ==========================================================

  scenario_rep_df <- do.call(
    rbind,
    scenario_replications
  )


  # ==========================================================
  # Scenario summaries
  # ==========================================================

  by_method <- split(
    scenario_rep_df,
    scenario_rep_df$method
  )


  for (method_name in names(by_method)) {

    d <- by_method[[method_name]]

    summary_id <- summary_id + 1L

    agg <- aggregate_method_results(d)

    all_summary[[summary_id]] <- data.frame(
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


  # ==========================================================
  # Save replication-level results
  # ==========================================================

  if (SAVE_RAW) {

    raw_id <- raw_id + 1L

    all_raw[[raw_id]] <- scenario_rep_df
  }
}


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

if (SAVE_RAW &&
    nrow(raw_results) != expected_raw_rows) {

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

  png(
    filename,
    width = 1800,
    height = 1200,
    res = 180
  )

  old_par <- par(
    mfrow = c(2, 2),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 2, 0)
  )

  methods <- c(
    "MLE",
    "RWM",
    "Collapsed"
  )

  pch_values <- c(
    MLE = 15,
    RWM = 16,
    Collapsed = 17
  )

  lty_values <- c(
    MLE = 1,
    RWM = 2,
    Collapsed = 3
  )

  n_ratios <- sort(
    unique(plot_data$n_ratio)
  )

  y_all <- plot_data[[yvar]]

  if (log_y) {
    positive_y <- y_all[is.finite(y_all) & y_all > 0]
    ylim <- range(positive_y)
  } else {
    ylim <- range(
      y_all,
      finite = TRUE
    )

    if (yvar %in% c("bias")) {
      ylim <- range(
        c(ylim, 0),
        finite = TRUE
      )
    }
  }

  for (ratio_value in n_ratios) {

    d_ratio <- plot_data[
      plot_data$n_ratio == ratio_value,
      ,
      drop = FALSE
    ]

    xlim <- range(
      d_ratio$m_ratio
    )

    plot(
      NA,
      xlim = xlim,
      ylim = ylim,
      log = if (log_y) "y" else "",
      xlab = "m / p",
      ylab = ylab,
      main = paste0(
        "n / p = ",
        ratio_value
      )
    )

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
        pch = pch_values[method_name],
        lty = lty_values[method_name]
      )
    }

    legend(
      "topright",
      legend = methods,
      pch = pch_values,
      lty = lty_values,
      bty = "n"
    )
  }

  par(old_par)

  mtext(
    paste0(
      "Simulation summary: ",
      ylab
    ),
    outer = TRUE,
    line = 0.2,
    font = 2,
    cex = 1.2
  )

  dev.off()
}


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


results <- read.csv("wishart_df_simulation_replications.csv")

diagnostic <- results %>%
  tidyr::pivot_wider(
    id_cols = c(scenario, p, m, n_true, n_ratio, m_ratio, replicate),
    names_from = method,
    values_from = n_est
  ) %>%
  mutate(
    difference = RWM - Collapsed,
    abs_difference = abs(difference)
  )

summary(diagnostic$difference)
mean(abs(diagnostic$difference))
quantile(diagnostic$difference,
         c(0.025, 0.5, 0.975))



diagnostic2 <- diagnostic %>%
  left_join(
    results %>%
      select(
        scenario, p, m, n_true, n_ratio, m_ratio,
        replicate, method, acceptance_rate
      ) %>%
      filter(method == "RWM"),
    by = c(
      "scenario", "p", "m", "n_true",
      "n_ratio", "m_ratio", "replicate"
    )
  )

cor(
  abs(diagnostic2$difference),
  diagnostic2$acceptance_rate,
  use = "complete.obs"
)


diagnostic2 %>%
  arrange(desc(abs_difference)) %>%
  select(
    p, m, n_true, n_ratio, m_ratio,
    difference, abs_difference, acceptance_rate
  ) %>%
  head(20)


diagnostic %>%
  mutate(
    relative_difference =
      (RWM - Collapsed) / Collapsed,
    abs_relative_difference =
      abs(relative_difference)
  ) %>%
  summarise(
    mean_abs_relative_difference =
      mean(abs_relative_difference, na.rm = TRUE),
    median_abs_relative_difference =
      median(abs_relative_difference, na.rm = TRUE),
    q95 =
      quantile(
        abs_relative_difference,
        0.95,
        na.rm = TRUE
      )
  )



sd_compare <- diagnostic %>%
  left_join(
    results %>%
      filter(method %in% c("Collapsed", "RWM")) %>%
      select(
        scenario, p, m, n_true, n_ratio, m_ratio,
        replicate, method, n_sd
      ) %>%
      tidyr::pivot_wider(
        names_from = method,
        values_from = n_sd,
        names_prefix = "sd_"
      ),
    by = c(
      "scenario", "p", "m", "n_true",
      "n_ratio", "m_ratio", "replicate"
    )
  ) %>%
  mutate(
    sd_ratio = sd_RWM / sd_Collapsed
  )

summary(sd_compare$sd_ratio)
quantile(
  sd_compare$sd_ratio,
  c(0.025, 0.5, 0.975),
  na.rm = TRUE
)
