# ============================================================
# real_data_applications.R
#
# Real-data applications for:
#   1. MLE
#   2. Collapsed Bayesian inference
#   3. RWM within Gibbs validation
#
# The collapsed Bayesian method is the primary method.
# The posterior mode is the primary point estimate for n.
# The posterior mean is reported as an additional summary.
# The posterior mean is used for V.
#
# Required file:
#   wishart_df_methods.R
# ============================================================

rm(list = ls())

source("wishart_df_methods.R")


library(ggplot2)
library(plot3D)
library(dplyr)


# ============================================================
# Configuration
# ============================================================

prior_upper <- 1e3
grid_size <- 3000
n_draws_V <- 5000

# RWM settings for the reported application table
rwm_iter <- 10000
rwm_burn <- 3000

U_scale <- 1e-4

# ============================================================
# Helper: run the three methods for one real dataset
# ============================================================

run_real_dataset <- function(
    model_data,
    dataset_name,
    n_v = NULL) {
  
  p <- nrow(model_data[[1]])
  
  if (is.null(n_v)) {
    n_v <- p
  }
  
  U <- diag(U_scale, p)
  
  cat("\n============================================================\n")
  cat(dataset_name, "\n")
  cat("============================================================\n")
  
  # ----------------------------------------------------------
  # MLE
  # ----------------------------------------------------------
  
  cat("  MLE...\n")
  
  mle_fit <- fit_mle(
    X = model_data
  )
  
  # ----------------------------------------------------------
  # Collapsed Bayesian
  # ----------------------------------------------------------
  
  cat("  Collapsed...\n")
  
  collapsed_fit <- fit_collapsed(
    X = model_data,
    n_v = n_v,
    U = U,
    prior_upper = prior_upper,
    grid_size = grid_size,
    n_draws = n_draws_V
  )
  
  # ----------------------------------------------------------
  # RWM validation
  # ----------------------------------------------------------
  
  cat("  RWM...\n")
  
  rwm_fit <- fit_rwm(
    X = model_data,
    iter = rwm_iter,
    burn = rwm_burn,
    prior_upper = prior_upper
  )
  
  # ----------------------------------------------------------
  # Requested reporting table
  # ----------------------------------------------------------
  
  summary_table <- rbind(
    
    data.frame(
      dataset = dataset_name,
      method = "MLE",
      n_mode = mle_fit$n$estimate,
      n_mean = mle_fit$n$estimate,
      n_median = mle_fit$n$estimate,
      n_sd = NA_real_,
      n_lower_95 = NA_real_,
      n_upper_95 = NA_real_,
      runtime_seconds = mle_fit$runtime,
      acceptance_rate = NA_real_,
      stringsAsFactors = FALSE
    ),
    
    data.frame(
      dataset = dataset_name,
      method = "Collapsed",
      n_mode = collapsed_fit$n$mode,
      n_mean = collapsed_fit$n$mean,
      n_median = collapsed_fit$n$median,
      n_sd = collapsed_fit$n$sd,
      n_lower_95 = collapsed_fit$n$lower,
      n_upper_95 = collapsed_fit$n$upper,
      runtime_seconds = collapsed_fit$runtime,
      acceptance_rate = NA_real_,
      stringsAsFactors = FALSE
    ),
    
    data.frame(
      dataset = dataset_name,
      method = "RWM",
      n_mode = rwm_fit$n$mode,
      n_mean = rwm_fit$n$mean,
      n_median = rwm_fit$n$median,
      n_sd = rwm_fit$n$sd,
      n_lower_95 = rwm_fit$n$lower,
      n_upper_95 = rwm_fit$n$upper,
      runtime_seconds = rwm_fit$runtime,
      acceptance_rate = rwm_fit$diagnostics$acceptance_rate,
      stringsAsFactors = FALSE
    )
  )
  
  list(
    dataset = dataset_name,
    summary = summary_table,
    collapsed = collapsed_fit,
    mle = mle_fit,
    rwm = rwm_fit
  )
}


# ============================================================
# Build NBA data
# ============================================================

basketball <- as.data.frame(
  read.csv(
    "../real_datasets/2022-2023 NBA Player Stats - Regular.csv",
    header = TRUE,
    sep = ";"
  )
)

basketball <- basketball[
  ,
  c(
    "Tm",
    "TRB",
    "AST",
    "STL",
    "BLK",
    "TOV",
    "PF",
    "PTS"
  )
]

basketball_model_data <- list()

i <- 1L

for (team in unique(basketball$Tm)) {
  
  team_data <- basketball[
    basketball$Tm == team,
    2:ncol(basketball)
  ]
  
  basketball_model_data[[i]] <- cov(team_data)
  
  i <- i + 1L
}


# ============================================================
# Build air-quality data
# ============================================================

air <- as.data.frame(
  read.csv(
    "../real_datasets/AQI and Lat Long of Countries.csv",
    header = TRUE,
    sep = ","
  )
)

air <- air[
  ,
  names(air) %in%
    c(
      "Country",
      "AQI.Value",
      "Ozone.AQI.Value"
    )
]

air_alps <- air[
  air$Country %in%
    c(
      "Italy",
      "France",
      "Switzerland",
      "Austria",
      "Germany",
      "Slovenia"
    ),
]

air_model_data <- list()

i <- 1L

for (country in unique(air_alps$Country)) {
  
  country_data <- air_alps[
    air_alps$Country == country,
    c(
      "AQI.Value",
      "Ozone.AQI.Value"
    )
  ]
  
  if (nrow(country_data) >= 2) {
    
    air_model_data[[i]] <- cov(country_data)
    
    i <- i + 1L
  }
}


# ============================================================
# Build handwriting data
# ============================================================

library(readxl)

hand_data <- as.data.frame(
  read_excel(
    "../real_datasets/DB_loop_handwriting_ls.xlsx"
  )
)

hand_data[,1:9] = scale(hand_data[,1:9])

# hand_data$new_id <- paste0(
#   hand_data$writer_id,
#   "_",
#   hand_data$character
# )

hand_model_data <- list()

i <- 1L

for (id in unique(hand_data$writer_id)) {
  
  spdf <- hand_data[
    hand_data$writer_id == id,
    1:9
  ]
  
  spdf <- spdf[
    complete.cases(spdf),
    ,
    drop = FALSE
  ]
  
  if (nrow(spdf) > 2) {
    
    hand_model_data[[i]] <- cov(spdf)
    
    i <- i + 1L
  }
}


# ============================================================
# Run all three real-data analyses
# ============================================================

nba_results <- run_real_dataset(
  model_data = basketball_model_data,
  dataset_name = "NBA 2022--2023",
  n_v = 7
)

air_results <- run_real_dataset(
  model_data = air_model_data,
  dataset_name = "Alpine air quality",
  n_v = 2
)

hand_results <- run_real_dataset(
  model_data = hand_model_data,
  dataset_name = "Handwriting",
  n_v = 9
)


# ============================================================
# Combine exactly the requested estimated parameters
# ============================================================

real_data_results <- rbind(
  nba_results$summary,
  air_results$summary,
  hand_results$summary
)

# Round only for the displayed/saved reporting table.
real_data_results_report <- real_data_results

numeric_columns <- c(
  "n_mode",
  "n_mean",
  "n_median",
  "n_sd",
  "n_lower_95",
  "n_upper_95",
  "runtime_seconds",
  "acceptance_rate"
)

real_data_results_report[
  numeric_columns
] <- lapply(
  real_data_results_report[numeric_columns],
  function(x) round(x, 4)
)


# ============================================================
# Save CSV
# ============================================================

#write.csv(
#  real_data_results_report,
#  file.path(
#    out_dir,
#    "real_data_estimated_parameters.csv"
#  ),
#  row.names = FALSE
#)


# ============================================================
# Posterior density of n from the collapsed posterior
# ============================================================

# Reconstruct and plot the actual collapsed posterior density
# for each dataset. No MCMC samples are used here.
plot_collapsed_n <- function(
    fit,
    model_data,
    dataset_title,
    file_name) {
  
  prep <- prepare_wishart_data(model_data)
  
  p <- prep$p
  n_v <- p
  U <- diag(U_scale, p)
  
  grid_object <- make_collapsed_grid(
    prep = prep,
    n_v = n_v,
    U = U,
    prior_upper = prior_upper,
    grid_size = grid_size
  )
  
  n_grid <- grid_object$grid
  logpost <- grid_object$logpost
  
  # Stable normalization of the posterior density.
  density_values <- exp(logpost - max(logpost))
  
  weights <- trapezoid_weights(n_grid)
  
  normalization <- sum(density_values * weights)
  density_values <- density_values / normalization
  
  posterior_df <- data.frame(
    n = n_grid,
    density = density_values
  )
  
  # ----------------------------------------------------------
  # Calculate posterior summaries directly from the density.
  # ----------------------------------------------------------
  
  cdf_values <- cumsum(density_values * weights)
  cdf_values <- cdf_values / max(cdf_values)
  
  posterior_quantile <- function(prob) {
    approx(
      x = cdf_values,
      y = n_grid,
      xout = prob,
      ties = "ordered",
      rule = 2
    )$y
  }
  
  q001 <- posterior_quantile(0.001)
  q999 <- posterior_quantile(0.999)
  
  mode_value <- fit$n$mode
  mean_value <- fit$n$mean
  median_value <- fit$n$median
  
  # ----------------------------------------------------------
  # Use a data-adaptive x-axis.
  # The old version used the full prior range up to 1000,
  # which compressed the posterior into a narrow spike.
  # ----------------------------------------------------------
  
  plot_width <- q999 - q001
  
  if (!is.finite(plot_width) || plot_width <= 0) {
    plot_width <- max(1, abs(q999), abs(q001)) * 0.1
  }
  
  x_lower <- max(
    p - 1,
    q001 - 0.08 * plot_width
  )
  
  x_upper <- min(
    prior_upper,
    q999 + 0.08 * plot_width
  )
  
  # Keep all three point summaries visible even in very concentrated cases.
  summary_min <- min(mode_value, median_value, mean_value)
  summary_max <- max(mode_value, median_value, mean_value)
  
  x_lower <- min(x_lower, summary_min - 0.02 * plot_width)
  x_upper <- max(x_upper, summary_max + 0.02 * plot_width)
  
  x_lower <- max(p - 1, x_lower)
  x_upper <- min(prior_upper, x_upper)
  
  # Restrict the displayed curve to the informative posterior region.
  plot_df <- posterior_df[
    posterior_df$n >= x_lower & posterior_df$n <= x_upper,
    ,
    drop = FALSE
  ]
  
  line_data <- data.frame(
    n = c(mode_value, median_value, mean_value),
    label = c("Posterior mode", "Posterior median", "Posterior mean"),
    linetype = c("solid", "dashed", "dotted")
  )
  
  p_plot <- ggplot(
    plot_df,
    aes(
      x = n,
      y = density
    )
  ) +
    
    geom_ribbon(
      aes(
        ymin = 0,
        ymax = density
      ),
      alpha = 0.15
    ) +
    
    geom_line(
      linewidth = 1
    ) +
    
    geom_vline(
      data = line_data,
      aes(
        xintercept = n,
        linetype = label
      ),
      linewidth = 0.8,
      show.legend = TRUE
    ) +
    
    scale_linetype_manual(
      values = c(
        "Posterior mode" = "solid",
        "Posterior median" = "dashed",
        "Posterior mean" = "dotted"
      )
    ) +
    
    coord_cartesian(
      xlim = c(x_lower, x_upper),
      expand = FALSE
    ) +
    
    labs(
      title = "Collapsed posterior distribution of the degrees of freedom",
      subtitle = dataset_title,
      x = "Degrees of freedom n",
      y = "Posterior density",
      linetype = NULL
    ) +
    
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(size = 20),
      plot.subtitle = element_text(size = 16),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13)
    )
  
  ggsave(
    filename = file_name,
    plot = p_plot,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  p_plot
}

# ============================================================
# Generate posterior plots for all three datasets
# ============================================================

plot_collapsed_n(
  fit = nba_results$collapsed,
  model_data = basketball_model_data,
  dataset_title = "NBA 2022--2023",
  file_name = file.path(
    "plots/posterior_n_NBA.png"
  )
)

plot_collapsed_n(
  fit = air_results$collapsed,
  model_data = air_model_data,
  dataset_title = "Alpine air quality",
  file_name = file.path(
    "plots/posterior_n_air_quality.png"
  )
)

plot_collapsed_n(
  fit = hand_results$collapsed,
  model_data = hand_model_data,
  dataset_title = "Handwriting",
  file_name = file.path(
    "plots/posterior_n_handwriting.png"
  )
)


# ============================================================
# 3D histogram of extreme eigenvalues of V
# ============================================================
plot_V_eigenvalues <- function(
    fit,
    dataset_title,
    file_name) {
  
  V_post <- fit$V_draws
  
  if (is.null(V_post)) {
    stop(
      "V posterior draws are missing. ",
      "Set n_draws_V > 0."
    )
  }
  
  n_V <- dim(V_post)[3]
  
  largest_eigen <- numeric(n_V)
  smallest_eigen <- numeric(n_V)
  
  for (i in seq_len(n_V)) {
    
    ev <- eigen(
      V_post[, , i],
      symmetric = TRUE,
      only.values = TRUE
    )$values
    
    largest_eigen[i] <- max(ev)
    smallest_eigen[i] <- min(ev)
  }
  
  keep <- is.finite(largest_eigen) &
    is.finite(smallest_eigen)
  
  largest_eigen <- largest_eigen[keep]
  smallest_eigen <- smallest_eigen[keep]
  
  n_bins <- 15
  
  breaks_biggest <- seq(
    min(largest_eigen),
    max(largest_eigen),
    length.out = n_bins + 1
  )
  
  breaks_smallest <- seq(
    min(smallest_eigen),
    max(smallest_eigen),
    length.out = n_bins + 1
  )
  
  biggest_cut <- cut(
    largest_eigen,
    breaks = breaks_biggest,
    include.lowest = TRUE
  )
  
  smallest_cut <- cut(
    smallest_eigen,
    breaks = breaks_smallest,
    include.lowest = TRUE
  )
  
  z <- table(
    biggest_cut,
    smallest_cut
  )
  
  x <- (
    breaks_biggest[-1] +
      breaks_biggest[-length(breaks_biggest)]
  ) / 2
  
  y <- (
    breaks_smallest[-1] +
      breaks_smallest[-length(breaks_smallest)]
  ) / 2
  
  stopifnot(
    length(x) == nrow(z),
    length(y) == ncol(z)
  )
  
  jpeg(
    filename = file_name,
    width = 3600,
    height = 3000,
    res = 300
  )
  
  hist3D(
    x = x,
    y = y,
    z = z,
    border = "black",
    lwd = 0.6,
    
    main = paste0(
      "Posterior distribution of extreme eigenvalues of V\n",
      dataset_title
    ),
    
    xlab = "\n\nLargest eigenvalue",
    ylab = "\n\nSmallest eigenvalue",
    zlab = "Frequency",
    
    cex.main = 2.0,
    cex.lab = 1.6,
    cex.axis = 1.15,
    
    theta = 35,
    phi = 25,
    expand = 1,
    ticktype = "detailed"
  )
  
  dev.off()
}


plot_V_eigenvalues(
  fit = nba_results$collapsed,
  dataset_title = "NBA 2022--2023",
  file_name = file.path(
    "plots/eigenvalues_V_NBA_3D.jpg"
  )
)

plot_V_eigenvalues(
  fit = air_results$collapsed,
  dataset_title = "Alpine air quality",
  file_name = file.path(
    "plots/eigenvalues_V_air_quality_3D.jpg"
  )
)

plot_V_eigenvalues(
  fit = hand_results$collapsed,
  dataset_title = "Handwriting",
  file_name = file.path(
    "plots/eigenvalues_V_handwriting_3D.jpg"
  )
)


# ============================================================
# Save complete analysis objects
# ============================================================

#saveRDS(
#  list(
#    NBA = nba_results,
#    AirQuality = air_results,
#    Handwriting = hand_results
#  ),
#  file = file.path(
#    out_dir,
#    "real_data_results_complete.rds"
#  )
#)


# ============================================================
# Print the final table
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "FINAL REAL-DATA ESTIMATES\n"
)

cat(
  "============================================================\n\n"
)

print(
  real_data_results_report
)

