library(ggplot2)
library(dplyr)
library(tidyr)

#### General settings ####

third <- rgb(63/255, 128/255, 97/255)
col5 <- rgb(204/255, 0/255, 0/255)

domain <- "Sardinia"
plot_size <- 7
output <- "Plots/Cross_Validation/"
dir.create(output, recursive = T, showWarnings = F)

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
ggsave(
  filename = paste0(output, "Scores_By_Model.pdf"),
  plot = p_model,
  width = 12,
  height = 8,
  units = "in")

#### Paired model comparisons ####

paired_results <- results %>%
  select(year, hole_id, model, MSE, CRPS) %>%
  pivot_wider(
    names_from = model,
    values_from = c(MSE, CRPS))

#### Numerical summaries ####

model_summary <- results %>%
  group_by(model) %>%
  summarise(
    MSE = mean(MSE),
    CRPS = mean(CRPS),
    .groups = "drop")

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
