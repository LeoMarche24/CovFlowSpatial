library(ggplot2)
library(dplyr)
library(tidyr)

#### General settings ####

third <- rgb(63/255, 128/255, 97/255)
col5 <- rgb(204/255, 0/255, 0/255)

domain <- "Sardinia"
plot_size <- 7

model_colors <- c(
  Network = third,
  Euclidean = col5)

plot_theme <- theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = plot_size*1.7),
    axis.title = element_text(size = plot_size*1.7),
    axis.text = element_text(size = plot_size*1.3),
    strip.text = element_text(size = plot_size*1.5)
  )

#### Data loading ####

load(paste0(
  "Data/Cross_Validation_Rolling_Blocks_SST_", domain, ".RData"))

results$model <- factor(
  results$model, levels = c("Network", "Euclidean"))
point_results$model <- factor(
  point_results$model, levels = c("Network", "Euclidean"))
years <- sort(unique(results$year))

hole_levels <- sort(hole_centers$hole_id)

results$hole <- factor(
  results$hole_id,
  levels = hole_levels,
  labels = paste0("H", hole_levels))
metric_results <- results %>%
  pivot_longer(
    cols = c(MSE, CRPS),
    names_to = "metric",
    values_to = "loss")

metric_results$metric <- factor(
  metric_results$metric,
  levels = c("MSE", "CRPS"))

#### Overall scores ####

p_model <- ggplot(
  metric_results,
  aes(x = model, y = loss, fill = model, color = model)) +
  geom_boxplot(alpha = 0.55, outlier.shape = NA) +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = model_colors) +
  scale_color_manual(values = model_colors, guide = "none") +
  labs(x = NULL, y = "Loss (lower is better)") +
  plot_theme

print(p_model)

#### Paired model comparisons ####

paired_results <- results %>%
  select(year, hole_id, model, MSE, CRPS) %>%
  pivot_wider(
    names_from = model,
    values_from = c(MSE, CRPS))

#### Numerical summaries ####

# Average spatial blocks within each year before estimating uncertainty.
# The t intervals below assume independent, approximately normal yearly means;
# they do not account for dependence between years from rolling estimation.
confidence_level <- 0.95

yearly_scores <- results %>%
  group_by(model, year) %>%
  summarise(
    MSE = mean(MSE),
    CRPS = mean(CRPS),
    .groups = "drop")

model_summary <- yearly_scores %>%
  group_by(model) %>%
  summarise(
    n_years = n(),
    MSE_se = sd(MSE) / sqrt(n_years),
    MSE = mean(MSE),
    CRPS_se = sd(CRPS) / sqrt(n_years),
    CRPS = mean(CRPS),
    .groups = "drop") %>%
  mutate(
    CRPS_margin = qt((1 + confidence_level) / 2, df = pmax(n_years - 1, 1)) * CRPS_se,
    CRPS_lower = CRPS - CRPS_margin,
    CRPS_upper = CRPS + CRPS_margin,
    MSE_margin = qt((1 + confidence_level) / 2, df = pmax(n_years - 1, 1)) * MSE_se,
    MSE_lower = MSE - MSE_margin,
    MSE_upper = MSE + MSE_margin) %>%
  select(-c(MSE_margin,CRPS_margin))

p_mse_ci <- ggplot(model_summary, aes(x = model, y = MSE, color = model)) +
  geom_errorbar(aes(ymin = MSE_lower, ymax = MSE_upper), width = 0.15) +
  geom_point(size = 3) +
  scale_color_manual(values = model_colors) +
  labs(x = NULL, y = "Mean MSE",
       subtitle = paste0(100 * confidence_level, "% t intervals across yearly means")) +
  plot_theme

print(p_mse_ci)

paired_summary <- data.frame(
  metric = c("MSE", "CRPS"),
  mean_network_minus_euclidean = c(
    mean(paired_results$MSE_Network-paired_results$MSE_Euclidean),
    mean(paired_results$CRPS_Network-paired_results$CRPS_Euclidean)),
  network_wins = c(
    sum(paired_results$MSE_Network < paired_results$MSE_Euclidean),
    sum(paired_results$CRPS_Network < paired_results$CRPS_Euclidean)),
  comparisons = nrow(paired_results)
)

print(model_summary)
print(paired_summary)
print(hole_centers)
