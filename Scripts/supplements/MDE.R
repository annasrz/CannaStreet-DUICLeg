# ==================================================================================================================================================================
#==================================================================================================================================================================
# ==================================================================================================================================================================

# PROJECT TITLE:  CANNASTREET - Analysis of Manuscript "Short-term effects of cannabis legalisation in Germany on driving under the influence of cannabis"
# SUPPLEMENT: Minimum Detectable Effect (MDE) for DiD
# CODE AUTHOR:    Anna Schranz
# DATE STARTED:   251125

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# clean workspace
rm(list=ls())

packages <- c("tidyverse", "kableExtra")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Load packages
invisible(lapply(packages, library, character.only = TRUE))


# FUNCTION: Simulation-based MDE calculation
# ----------------------------------
calculate_mde_simulation <- function(n_treat_pre, n_treat_post,
                                     n_control_pre, n_control_post,
                                     p_treat_pre, p_control_pre,
                                     time_trend_control,  
                                     time_trend_treat,  
                                     alpha = 0.05, 
                                     power = 0.80,
                                     n_sim = 100,
                                     effect_grid = seq(0, 0.2, by = 0.005)) {
  
  cat("\n=== Starting MDE Simulation ===\n")
  cat(sprintf("Target power: %.0f%%\n", power * 100))
  cat(sprintf("Alpha level: %.3f\n", alpha))
  cat(sprintf("Number of simulations per effect size: %d\n", n_sim))
  cat(sprintf("Sample sizes - Treatment: n_pre=%d, n_post=%d\n", n_treat_pre, n_treat_post))
  cat(sprintf("Sample sizes - Control: n_pre=%d, n_post=%d\n", n_control_pre, n_control_post))
  cat(sprintf("Baseline prevalences - Treatment: %.3f, Control: %.3f\n", p_treat_pre, p_control_pre))
  cat(sprintf("Assumed time trend (Control): %.3f (%.1f pp)\n", time_trend_control, time_trend_control * 100))
  cat(sprintf("Assumed time trend (Treatment): %.3f (%.1f pp)\n\n", time_trend_treat, time_trend_treat * 100))
  
  power_results <- numeric(length(effect_grid))
  mean_coef_log_or <- numeric(length(effect_grid))
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = length(effect_grid), style = 3)
  
  for (i in seq_along(effect_grid)) {
    effect <- effect_grid[i]
    significant_count <- 0
    coef_storage <- numeric(n_sim)
    
    for (sim in 1:n_sim) {
      # simulate data under the assumption of a true intervention effect
      # Treatment group (Germany): time trend + intervention effect
      y_treat_pre <- rbinom(n_treat_pre, 1, p_treat_pre)
      p_treat_post <- min(p_treat_pre + time_trend_treat + effect, 1.0)  
      y_treat_post <- rbinom(n_treat_post, 1, p_treat_post)
      
      # Control group (Austria) - only time trend, no intervention effect
      y_control_pre <- rbinom(n_control_pre, 1, p_control_pre)
      p_control_post <- min(p_control_pre + time_trend_control, 1.0)  
      y_control_post <- rbinom(n_control_post, 1, p_control_post)
      
      # dataset in DiD-format
      sim_data <- data.frame(
        y = c(y_treat_pre, y_treat_post, y_control_pre, y_control_post),
        time = rep(c(0, 1, 0, 1), 
                   c(n_treat_pre, n_treat_post, n_control_pre, n_control_post)),
        treat = rep(c(1, 1, 0, 0), 
                    c(n_treat_pre, n_treat_post, n_control_pre, n_control_post))
      )
      
      # logistic regression DiD model
      tryCatch({
        model <- glm(y ~ time * treat, data = sim_data, family = binomial())
        
        # check significance of interaction term (DiD effect)
        p_value <- summary(model)$coefficients["time:treat", "Pr(>|z|)"]
        
        coef_storage[sim] <- coef(model)["time:treat"]
        
        if (!is.na(p_value) && p_value < alpha) {
          significant_count <- significant_count + 1
        }
      }, error = function(e) {
        coef_storage[sim] <- NA
        # if model does not converge, count as non-significant
      })
    }
    
    power_results[i] <- significant_count / n_sim
    mean_coef_log_or[i] <- mean(coef_storage, na.rm = TRUE)
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # find MDE: smallest effect size with power >= target power
  mde_idx <- which(power_results >= power)[1]
  
  if (is.na(mde_idx)) {
    warning("Desired power not achieved in effect grid. Consider expanding the grid.")
    mde <- NA
    mde_or_model <- NA
    mde_log_or_model <- NA
    achieved_power <- max(power_results)
    mde_idx <- which.max(power_results)
  } else {
    mde <- effect_grid[mde_idx]
    achieved_power <- power_results[mde_idx]
    
    mde_log_or_model <- mean_coef_log_or[mde_idx]
    mde_or_model <- exp(mde_log_or_model)
  }
  
  p_treat_post_mde <- p_treat_pre + time_trend_treat + mde
  p_control_post_mde <- p_control_pre + time_trend_control
  
  cat("\n\n=== Results ===\n")
  cat(sprintf("MDE (absolute): %.4f (%.2f percentage points)\n", mde, mde * 100))
  cat(sprintf("MDE (relative to baseline): %.2f%%\n", (mde / p_treat_pre) * 100))
  cat(sprintf("MDE (Odds Ratio - from model): %.3f\n", mde_or_model))
  cat(sprintf("Achieved power at MDE: %.3f\n", achieved_power))
  
  return(list(
    mde_absolute = mde,
    mde_percentage_points = mde * 100,
    mde_relative = mde / p_treat_pre,
    mde_relative_percent = (mde / p_treat_pre) * 100,
    mde_or = mde_or_model,
    mde_log_or = mde_log_or_model,
    achieved_power = achieved_power,
    effect_grid = effect_grid,
    power_curve = power_results,
    mean_log_or_curve = mean_coef_log_or, 
    mde_index = mde_idx,
    baseline_prev = p_treat_pre,
    baseline_prev_control = p_control_pre,
    time_trend_control = time_trend_control,
    time_trend_treat = time_trend_treat,
    alpha = alpha,
    target_power = power,
    n_sim = n_sim,
    sample_sizes = list(
      n_treat_pre = n_treat_pre,
      n_treat_post = n_treat_post,
      n_control_pre = n_control_pre,
      n_control_post = n_control_post
    )
  ))
}
# ==================================================================================================
# SETUP
# ==================================================================================================

# Set paths
dataexport <- TRUE
DATE <- format(Sys.Date(), "%Y%m%d")
folder_path_plots <- "Output/figures/" 
folder_path_tables <- "Output/tables/"

# ==================================================================================================
# ANALYSIS 1: Cannabis Use (12-month prevalence)
# ==================================================================================================

cat("\n\n")
cat("================================================================================\n")
cat("  ANALYSIS 1: Cannabis Use (12M prevalence)\n")
cat("================================================================================\n")

# Sample sizes
n_DE_t0 <- 6670
n_DE_t1 <- 9692
n_AT_t0 <- 2132
n_AT_t1 <- 2102

# Baseline prevalences 
p_DE_t0_canuse <- 0.121  # 12.1%
p_AT_t0_canuse <- 0.094  # 9.4%
p_AT_t1_canuse <- 0.096  # 9.6%

# Time trend - PARALLELE TRENDS assumption (Standard DiD)
time_trend <- p_AT_t1_canuse - p_AT_t0_canuse  # +0.2 pp

time_trend_control <- time_trend  # AT: +0.2 pp
time_trend_treat <- time_trend    # DE: +0.2 pp (parallel!)

cat(sprintf("Assumed timetrend (both countries): %.4f (%.2f pp)\n", time_trend, time_trend * 100))
cat(sprintf("Counterfactual GER-prevalence (without intervention): %.1f%%\n", 
            (p_DE_t0_canuse + time_trend_treat) * 100))

# MDE Simulation
set.seed(1234)
mde_cannabis <- calculate_mde_simulation(
  n_treat_pre = n_DE_t0,
  n_treat_post = n_DE_t1,
  n_control_pre = n_AT_t0,
  n_control_post = n_AT_t1,
  p_treat_pre = p_DE_t0_canuse,
  p_control_pre = p_AT_t0_canuse,
  time_trend_control = time_trend_control,  
  time_trend_treat = time_trend_treat,
  alpha = 0.05,
  power = 0.80,
  n_sim = 500,
  effect_grid = seq(0, 0.1, by = 0.002)
)

# ==================================================================================================
# ANALYSIS 2: DUIC (12-month prevalence)
# ==================================================================================================

cat("\n\n")
cat("================================================================================\n")
cat("  ANALYSIS 2: DUIC (12M prevalence)\n")
cat("================================================================================\n")

# Sample sizes
n_DE_duic_t0 <- 393
n_DE_duic_t1 <- 589
n_AT_duic_t0 <- 86
n_AT_duic_t1 <- 92

# Baseline prevalences
p_DE_duic_t0 <- 0.285  
p_AT_duic_t0 <- 0.128  
p_AT_duic_t1 <- 0.163 

# Time trend - PARALLELE TRENDS assumption (Standard DiD)
time_trend_duic <- p_AT_duic_t1 - p_AT_duic_t0  # +3.5 pp

time_trend_control_duic <- time_trend_duic  # AT: +3.5 pp
time_trend_treat_duic <- time_trend_duic    # DE: +3.5 pp (parallel!)

cat(sprintf("Assumed timetrend (both countries): %.4f (%.2f pp)\n", time_trend_duic, time_trend_duic * 100))
cat(sprintf("Counterfactual GER-prevalence (without intervention): %.1f%%\n", 
            (p_DE_duic_t0 + time_trend_treat_duic) * 100))

# MDE Simulation
set.seed(1234)
mde_duic <- calculate_mde_simulation(
  n_treat_pre = n_DE_duic_t0,
  n_treat_post = n_DE_duic_t1,
  n_control_pre = n_AT_duic_t0,
  n_control_post = n_AT_duic_t1,
  p_treat_pre = p_DE_duic_t0,
  p_control_pre = p_AT_duic_t0,
  time_trend_control = time_trend_control_duic,  
  time_trend_treat = time_trend_treat_duic,
  alpha = 0.05,
  power = 0.80,
  n_sim = 500,
  effect_grid = seq(0, 0.4, by = 0.002) 
)

# ==================================================================================================
# POWER CURVE PLOTS
# ==================================================================================================

cat("\n\n")
cat("================================================================================\n")
cat("  Creating Power Curve Plots\n")
cat("================================================================================\n")

# Function to create power curve plot
create_power_curve_plot <- function(mde_result, outcome_name, observed_or = NULL) {
  
  # Prepare data
  plot_data <- data.frame(
    effect_pp = mde_result$effect_grid * 100,
    effect_or = exp(mde_result$mean_log_or_curve),
    power = mde_result$power_curve * 100
  )
  
  # Create plot with percentage points on x-axis
  p1 <- ggplot(plot_data, aes(x = effect_pp, y = power)) +
    geom_line(linewidth = 1.2, color = "#2C3E50") +
    geom_hline(yintercept = 80, linetype = "dashed", color = "#E74C3C", linewidth = 1) +
    geom_vline(xintercept = mde_result$mde_percentage_points, 
               linetype = "dashed", color = "#3498DB", linewidth = 1) +
    annotate("text", x = mde_result$mde_percentage_points, 
             y = 5, label = sprintf("MDE = %.2f pp", mde_result$mde_percentage_points),
             hjust = -0.1, color = "#3498DB", size = 4, fontface = "bold") +
    annotate("text", x = max(plot_data$effect_pp) * 0.5, 
             y = 85, label = "Target Power = 80%",
             color = "#E74C3C", size = 4, fontface = "bold") +
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
    labs(
      title = paste0("Power Curve: ", outcome_name),
      subtitle = sprintf("n = %d/%d (DE), %d/%d (AT) | %d simulations per effect size",
                         mde_result$sample_sizes$n_treat_pre,
                         mde_result$sample_sizes$n_treat_post,
                         mde_result$sample_sizes$n_control_pre,
                         mde_result$sample_sizes$n_control_post,
                         mde_result$n_sim),
      x = "Effect Size (DiD-effect in percentage points)",
      y = "Statistical Power (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90")
    )
  
  # Create plot with OR on x-axis
  p2 <- ggplot(plot_data, aes(x = effect_or, y = power)) +
    geom_line(linewidth = 1.2, color = "#2C3E50") +
    geom_hline(yintercept = 80, linetype = "dashed", color = "#E74C3C", linewidth = 1) +
    geom_vline(xintercept = mde_result$mde_or, 
               linetype = "dashed", color = "#3498DB", linewidth = 1) +
    geom_vline(xintercept = 1, linetype = "dotted", color = "gray50", linewidth = 0.8) +
    annotate("text", x = mde_result$mde_or, 
             y = 5, label = sprintf("MDE (OR) = %.2f", mde_result$mde_or),
             hjust = -0.1, color = "#3498DB", size = 4, fontface = "bold") +
    scale_x_continuous(limits = c(0, max(plot_data$effect_or) * 1.1)) +
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
    labs(
      title = paste0("Power Curve: ", outcome_name),
      subtitle = "Effect size: DiD effect, expressed as Odds Ratio; alpha = 0.05",
      x = "Effect Size (Odds Ratio)",
      y = "Statistical Power (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90")
    )
  
  # Add observed OR if provided
  if (!is.null(observed_or)) {
    p2 <- p2 + 
      geom_vline(xintercept = observed_or, 
                 linetype = "dotted", color = "#27AE60", linewidth = 1) +
      annotate("text", x = observed_or, 
               y = 95, label = sprintf("Observed\nOR = %.2f", observed_or),
               hjust = 1.1, color = "#27AE60", size = 3.5, fontface = "bold")
  }
  
  return(list(plot_pp = p1, plot_or = p2))
}

# Create plots for Cannabis Use
plots_cannabis <- create_power_curve_plot(
  mde_cannabis, 
  "Cannabis Use (12-month prevalence)",
  observed_or = 1.18
)

# Create plots for DUIC
plots_duic <- create_power_curve_plot(
  mde_duic, 
  "DUIC (12-month prevalence)",
  observed_or = 0.68
)

# Save plots
if (dataexport) {
  # Cannabis Use - pp scale
  ggsave(
    filename = paste0(folder_path_plots, "power_curve_cannabis_pp_", DATE, ".png"),
    plot = plots_cannabis$plot_pp,
    width = 10, height = 6, dpi = 300
  )
  
  
  ggsave(
    filename = paste0(folder_path_plots, "power_curve_cannabis_OR_", DATE, ".png"),
    plot = plots_cannabis$plot_or,
    width = 10, height = 6, dpi = 300
  )
  
  
  ggsave(
    filename = paste0(folder_path_plots, "power_curve_DUIC_pp_", DATE, ".png"),
    plot = plots_duic$plot_pp,
    width = 10, height = 6, dpi = 300
  )
  
  
  ggsave(
    filename = paste0(folder_path_plots, "power_curve_DUIC_OR_", DATE, ".png"),
    plot = plots_duic$plot_or,
    width = 10, height = 6, dpi = 300
  )
  
  cat("Power curve plots saved!\n")
}

# Display plots
print(plots_cannabis$plot_pp)
print(plots_cannabis$plot_or)
print(plots_duic$plot_pp)
print(plots_duic$plot_or)

# ==================================================================================================
# SUMMARY TABLE
# ==================================================================================================

mde_summary <- data.frame(
  Analysis = c("Cannabis Use (12M)", "DUIC (12M)"),
  
  # Sample sizes
  n_DE_t0 = c(n_DE_t0, n_DE_duic_t0),
  n_DE_t1 = c(n_DE_t1, n_DE_duic_t1),
  n_AT_t0 = c(n_AT_t0, n_AT_duic_t0),
  n_AT_t1 = c(n_AT_t1, n_AT_duic_t1),
  
  # Baseline prevalences
  Baseline_DE = c(p_DE_t0_canuse, p_DE_duic_t0) * 100,
  Baseline_AT = c(p_AT_t0_canuse, p_AT_duic_t0) * 100,
  
  # Time trends - PARALLEL 
  Time_Trend_pp = c(time_trend, time_trend_duic) * 100,
  Counterfactual_DE = c(p_DE_t0_canuse + time_trend_treat, 
                        p_DE_duic_t0 + time_trend_treat_duic) * 100,
  
  # MDE results
  MDE_pp = c(mde_cannabis$mde_percentage_points, mde_duic$mde_percentage_points),
  MDE_relative_pct = c(mde_cannabis$mde_relative_percent, mde_duic$mde_relative_percent),
  MDE_OR = c(mde_cannabis$mde_or, mde_duic$mde_or), 
  
  # Power settings
  Achieved_Power = c(mde_cannabis$achieved_power, mde_duic$achieved_power) * 100,
  Alpha = c(mde_cannabis$alpha, mde_duic$alpha),
  N_Simulations = c(mde_cannabis$n_sim, mde_duic$n_sim)
)

# Format table
mde_summary_formatted <- mde_summary %>%
  mutate(
    across(c(Baseline_DE, Baseline_AT, Time_Trend_pp, Counterfactual_DE,
             MDE_pp, MDE_relative_pct, Achieved_Power), 
           ~round(.x, 2)),
    across(c(MDE_OR), ~round(.x, 3)), 
    across(c(Alpha), ~round(.x, 3))
  )

print(mde_summary_formatted)

# Export as CSV
if (dataexport) {
  write.csv(
    mde_summary_formatted, 
    paste0(folder_path_tables, "MDE_summary_parallel_trends_", DATE, ".csv"),
    row.names = FALSE
  )
}

# Export as formatted table (kable)
if (dataexport) {
  mde_table <- mde_summary_formatted %>%
    kbl(
      caption = "Minimum Detectable Effect (MDE) Analysis",
      col.names = c("Analysis", "DE t₀", "DE t₁", "AT t₀", "AT t₁", 
                    "Baseline DE (%)", "Baseline AT (%)", 
                    "Trend (pp)", "Counterfact. DE (%)",
                    "MDE (pp)", "MDE (% rel.)", "MDE (OR)", 
                    "Power (%)", "Alpha", "N Sims"),
      align = c("l", rep("r", 14))  
    ) %>%
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE
    ) %>%
    add_header_above(c(" " = 1, "Sample Sizes" = 4, "Baseline Prev." = 2, 
                       "Time Trends" = 2, "MDE Results" = 3, "Settings" = 3))  
  
  save_kable(
    mde_table,
    file = paste0(folder_path_tables, "MDE_summary_table_", DATE, ".html")
  )
}

cat("\n\n=== MDE Analysis Complete ===\n")
cat("Results saved to:", folder_path_tables, "\n")