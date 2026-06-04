library(reshape2)
library(ggplot2)
library(tidyverse)

source("functions_observations.R")

#### General settings ####

col1 <- "cyan"
col2 <- "darkblue"
third <- rgb(63/255, 128/255, 97/255)
fourth <- rgb(252/255, 233/255, 3/255)
grey <- "grey80"
custom_colorscale <- list(
  c(0, 1),
  c(col2, col1)
)

years <- 2006:2022
domain <- "Sardinia"

plot_size <- 7
line_size <- 2

####RCP 4.5 ####

load(paste0("HPC/Export/output_data/", domain, "/Covariance_Data_45_", domain, ".RData"))

# Residuals analysis

plot(years, unlist(lapply(res_total, mean)), pch = 16, col = third)
plot(bias, pch = 16, col = third)

ind <- 15
hist(res_total[[ind]], breaks = 20, main = paste("Year", as.character(years[ind])), 
     xlab = "Residuals", col = grey)

plot(years, sill, pch = 16, col = third)
plot(years, range, pch = 16, col = third)

plot(fuvs)

mean_covariance_45 <- data.frame(dist = round(colMeans(distances_45, na.rm = T)), 
                                gamma = colMeans(covariances_45), 
                                np = 1)

covariances_45 <- melt(covariances_45)
colnames(covariances_45) <- c("Row", "Column", "Value")
covariances_45$Distance <- 
  round(colMeans(distances_45, na.rm = T)[covariances_45$Column])

ggplot(covariances_45, aes(x = factor(Distance), y = Value, group = Row)) + 
  geom_line(color = col2) +  # Set color for the lines
  labs(x = "Distance", y = "Values") +
  theme_minimal()

covariances_45 <- 
  covariances_45[which(!is.na(covariances_45$Distance)),]
mean_covariance_45 <- 
  mean_covariance_45[which(!is.na(mean_covariance_45$dist)) ,]

quant_df <- covariances_45 %>%
  group_by(Distance) %>%
  summarise(
    q05 = quantile(Value, 0.2),
    q25 = quantile(Value, 0.4),
    q50 = quantile(Value, 0.50),
    q75 = quantile(Value, 0.6),
    q95 = quantile(Value, 0.8)
  )

p <- ggplot(quant_df, aes(x = as.numeric(Distance))) +
  geom_ribbon(aes(ymin = q25, ymax = q75),
              fill = third, alpha = 0.6) +
  geom_ribbon(aes(ymin = q05, ymax = q95),
              fill = third, alpha = 0.3) +
  geom_line(data = mean_covariance_45,
            aes(x = dist, y = gamma),
            color = "black", linewidth = 1.2) +
  labs(x = "Distance (km)", y = expression(hat(C)(h))) +
  scale_x_continuous(breaks = seq(0, max(quant_df$Distance), length = 5),
            labels = floor(seq(0, max(quant_df$Distance)/1000, length = 5))) +
  ylim(c(-0.5, 1.1)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = plot_size*3),
    axis.title = element_text(size = plot_size*3),
    legend.text = element_text(size = plot_size*3),
    legend.title = element_text(size = plot_size*3),
    strip.text = element_text(size = plot_size*3)
  )
p
ggsave(
  filename = paste0("HPC/Export/output_data/", domain, "/Plots/Real Data 45/Covariance_45.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

mean_covariance_45_eucl <- data.frame(dist = colMeans(distances_45_eucl), 
                                     gamma = colMeans(variograms_45_eucl))

variograms_45_eucl <- melt(variograms_45_eucl)
colnames(variograms_45_eucl) <- c("Row", "Column", "Value")
variograms_45_eucl$Distance <- 
  colMeans(distances_45_eucl)[variograms_45_eucl$Column]

ggplot(variograms_45_eucl, aes(x = factor(Distance), y = Value, group = Row)) + 
  geom_line(color = col2) +  # Set color for the lines
  labs(x = "Distance", y = "Values") +
  theme_minimal()

quant_df <- variograms_45_eucl %>%
  group_by(Distance) %>%
  summarise(
    q05 = quantile(Value, 0.05),
    q25 = quantile(Value, 0.25),
    q50 = quantile(Value, 0.50),
    q75 = quantile(Value, 0.75),
    q95 = quantile(Value, 0.95)
  )

p <- ggplot(quant_df, aes(x = as.numeric(Distance))) +
  geom_ribbon(aes(ymin = q25, ymax = q75),
              fill = third, alpha = 0.6) +     # Central band
  geom_ribbon(aes(ymin = q05, ymax = q95),
              fill = third, alpha = 0.3) +     # External band
  geom_line(data = mean_covariance_45_eucl,
            aes(x = dist, y = gamma),
            color = "black", linewidth = 1.2) + 
  labs(x = "Distance (km)", y = "Semi-variance") +
  scale_x_continuous(breaks = seq(0, max(quant_df$Distance), length = 5),
                     labels = floor(seq(0, max(quant_df$Distance)/1000, length = 5))) +
  coord_cartesian(ylim = c(-0.5, 2)) + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = plot_size*3),
    axis.title = element_text(size = plot_size*3),
    legend.text = element_text(size = plot_size*3),
    legend.title = element_text(size = plot_size*3),
    strip.text = element_text(size = plot_size*3)
  )
p
ggsave(
  filename = paste0("HPC/Export/output_data/", domain, "/Plots/Real Data 45/Variogram_45_eucl.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

sill45 <- mean(fuvs)

range45 <- fit_range(mean_covariance_45, exponential_covariance,
                     max(mean_covariance_45$dist)/3, sill45)

res <- fit_covariance(mean_covariance_45_eucl, exponential_kernel,
        c(max(mean_covariance_45_eucl$gamma), max(mean_covariance_45_eucl$dist)/3))

sill45_eucl <- res[1]
range45_eucl <- res[2]
bias45 <- mean(bias)

#### Save results ####

save(sill45, range45, sill45_eucl, range45_eucl,
     bias45,
     file = paste0("HPC/Export/input_data/Covariance_Data_45_", domain, ".RData"))





