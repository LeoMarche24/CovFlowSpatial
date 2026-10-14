library(sp)
library(gstat)
library(ncdf4)
library(ggplot2)
library(spatstat)
library(igraph)
library(sf)
library(progress)
library(parallel)

source("functions_network.R")
source("functions_observations.R")

col1 <- "#3182BD"
col2 <- "#E6550D"
third <- rgb(63/255, 128/255, 97/255)
fourth <- "#333399"
col5 <- rgb(204/255, 0/255, 0/255)
grey <- "grey80"

domain <- "Sardinia"
load(paste0("Data/Data_", domain, ".RData"))
load(paste0("Data/bbox_", domain, ".RData"))
df <- df[which(df$years == 2022) ,]

df$value <- df$temperature
inx <- which(!is.na(df$east) & !is.na(df$north) & !is.na(df$value))
B <- length(inx)
results <- compute_matrices(df = df)
dist <- results[[1]]
PI <- results[[2]]
inx <- results[[3]]
lines <- build_network_lines(df, dist)

PI <- PI[inx, inx]
dist <- dist[inx, inx]

sources <- which(colSums(!is.na(PI)) == 0)
outlets <- which(rowSums(PI, na.rm=T) < 0.99)

B <- dim(PI)[1]
P <- PI
B <- dim(P)[1]
for(i in 1:B)
{
  for(j in 1:B)
  {
    if(i!=j & !is.na(P[i,j]))
    {
      P[i,j] <- PI[i,j]/sqrt(sum(PI[,j], na.rm = T))
    }
  }
}

print("start composing distance object")
distances <- compose_distances(dist, P)

U <- evaluate_U(PI)

#### Sensibility study ####

sill <- 1
B <- length(distances)
R <- 5
ex1 <- 5e4
ex2 <- 2e5
range <- seq(ex1, ex2, length = R)
K <- 100
bins <- c(10, 13, 15, 17, 20)
model <- "Exponential"
output <- paste0("Plots/Sensibility_Estimator/", model, "/", domain, "/")
dir.create(output, recursive = T, showWarnings = F)

covariances <- lapply(1:R, function(x)
  lapply(seq_along(bins), function(y) vector("list", K)))
sills <- array(0, dim = c(K, R, length(bins)))
ranges <- array(0, dim = c(K, R, length(bins)))
MSE <- array(0, dim = c(K, R, length(bins)))

xx <- seq(0, max(range)*3, length.out = 200)

simulation_distribution <- do.call(rbind, lapply(1:K, function(i)
{
  set.seed(i)
  rnorm(B,0,1)
}))

for (r in 1:R)
{
  cov <- build_covariance_exponential(sill, range[r], dist, P, U)

  if (!all(eigen(cov)$values>1e-15))
  {
    cov <- diag(rep(sill, B))
    print("Covariance matrix is not positive definite, using diagonal matrix instead.")
  }

  chol_cov <- chol(cov)

  for (iter in 1:K)
  {
    print(iter)
    simulation <- t(chol_cov) %*% simulation_distribution[iter,]

    for (b in seq_along(bins))
    {
      covariance <- evaluate_covariance_penalization(
        simulation, distances, bins[b]-1, U, return_fuv = T)
      sill_est <- covariance[[1]]
      covariance <- covariance[[2]]
      covariances[[r]][[b]][[iter]] <- covariance
      sills[iter, r, b] <- sill_est

      range_est <- fit_range(covariance, exponential_covariance,
                             max(covariance$dist)/3, sill_est)
      ranges[iter, r, b] <- range_est

      covariance_est <- exponential_covariance(c(sill_est, range_est), xx)
      covariance_true <- exponential_covariance(c(sill, range[r]), xx)
      MSE[iter, r, b] <- mean((covariance_est-covariance_true)^2)
    }
  }
}

results <- expand.grid(
  iteration = 1:K,
  range = range,
  bins = bins
)
results$sill_estimated <- c(sills)
results$range_estimated <- c(ranges)
results$relative_range_error <-
  100*(results$range_estimated-results$range)/results$range
results$MSE <- c(MSE)

plot_size <- 5
line_size <- 2
bin_colors <- c(col1, col2, third, col5, fourth)

true_ranges <- data.frame(
  range = range,
  range_estimated = range/1000
)

p <- ggplot(results,
            aes(x = factor(range), y = range_estimated/1000,
                fill = factor(bins), color = factor(bins))) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(data = true_ranges,
             aes(x = factor(range), y = range_estimated),
             inherit.aes = F, color = "black", size = 3) +
  scale_fill_manual(values = bin_colors) +
  scale_color_manual(values = bin_colors, guide = "none") +
  scale_x_discrete(labels = function(x)
    as.character(as.integer(as.numeric(x)/1000))) +
  labs(x = "True range (km)", y = "Estimated range (km)", fill = "Bins") +
  theme_minimal() +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = plot_size*3),
    legend.text        = element_text(size = plot_size*3),
    axis.title.x       = element_text(size = plot_size*3),
    axis.title.y       = element_text(size = plot_size*3),
    axis.text.x        = element_text(size = plot_size*3),
    axis.text.y        = element_text(size = plot_size*3)
  )
p
ggsave(
  filename = paste0(output, "Estimated_Range.pdf"),
  plot     = p,
  width    = plot_size*1.5,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

p <- ggplot(results,
            aes(x = factor(range), y = relative_range_error,
                fill = factor(bins), color = factor(bins))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  scale_fill_manual(values = bin_colors) +
  scale_color_manual(values = bin_colors, guide = "none") +
  scale_x_discrete(labels = function(x)
    as.character(as.integer(as.numeric(x)/1000))) +
  labs(x = "True range (km)", y = "Relative range error (%)", fill = "Bins") +
  theme_minimal() +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = plot_size*3),
    legend.text        = element_text(size = plot_size*3),
    axis.title.x       = element_text(size = plot_size*3),
    axis.title.y       = element_text(size = plot_size*3),
    axis.text.x        = element_text(size = plot_size*3),
    axis.text.y        = element_text(size = plot_size*3)
  )
p
ggsave(
  filename = paste0(output, "Relative_Error.pdf"),
  plot     = p,
  width    = plot_size*1.5,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

p <- ggplot(results,
            aes(x = factor(range), y = MSE,
                fill = factor(bins), color = factor(bins))) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  scale_fill_manual(values = bin_colors) +
  scale_color_manual(values = bin_colors, guide = "none") +
  scale_x_discrete(labels = function(x)
    as.character(as.integer(as.numeric(x)/1000))) +
  labs(x = "True range (km)", y = "MSE", fill = "Bins") +
  ylim(c(0,0.3)) +
  theme_minimal() +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = plot_size*3),
    legend.text        = element_text(size = plot_size*3),
    axis.title.x       = element_text(size = plot_size*3),
    axis.title.y       = element_text(size = plot_size*3),
    axis.text.x        = element_text(size = plot_size*3),
    axis.text.y        = element_text(size = plot_size*3)
  )
p
ggsave(
  filename = paste0(output, "MSE_Covariance.pdf"),
  plot     = p,
  width    = plot_size*1.5,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

curve_results <- do.call(rbind, lapply(1:R, function(r)
  do.call(rbind, lapply(seq_along(bins), function(b)
    do.call(rbind, lapply(1:K, function(iter) {
      covariance <- covariances[[r]][[b]][[iter]]
      covariance$iteration <- iter
      covariance$range <- range[r]/1000
      covariance$bins <- bins[b]
      covariance$dist <- covariance$dist/1000
      covariance
    }))
  ))
))

true_covariances <- do.call(rbind, lapply(range, function(x)
  data.frame(
    range = x/1000,
    dist = xx/1000,
    gamma = exponential_covariance(c(sill, x), xx)
  )
))

p <- ggplot(curve_results,
            aes(x = dist, y = gamma, group = iteration)) +
  geom_line(color = grey, alpha = 0.4) +
  geom_line(data = true_covariances,
            aes(x = dist, y = gamma, group = 1),
            inherit.aes = F, color = third, linewidth = line_size) +
  facet_grid(bins ~ range, labeller = label_both) +
  labs(x = "Network distance (km)", y = "C") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = plot_size*3),
    axis.title.y = element_text(size = plot_size*3),
    strip.text = element_text(size = plot_size*2),
    axis.text.x = element_text(size = plot_size*2),
    axis.text.y = element_text(size = plot_size*2)
  )
p
