## Libraries ##

library(sp)
library(gstat)
library(ncdf4)
library(ggplot2)
library(spatstat)
library(igraph)
library(sf)
library(progress)
library(tidyverse)
library(plotly)
library(reshape2)
library(geosphere)

## Functions ##

source("functions_network.R")
source("functions_observations.R")
col1 <- "#3182BD"
col2 <- "#E6550D"
third <- rgb(63/255, 128/255, 97/255)
fourth <- "#333399"
col5 <- rgb(204/255, 0/255, 0/255)
grey <- "grey80"
custom_colorscale <- list(
  c(0, 1),
  c(col2, col1)
)
plot_size <- 7

#### Data loading ####

domain <- "Sardinia"
load(paste0("Data/Covariance_Data_45_", domain, ".RData"))
load(paste0("Data/Projections_", domain, ".RData"))
data_real <- df

#### Hot spot identification ####

year <- 2050
rcp <- 'rcp45'
output <- paste0("Plots/Projections/")

df <- data[which(data$year==year & data$RCP == rcp), ]

ggplot(df[which(!is.na(df$east) & !is.na(df$north) & !is.na(df$value)) ,], 
       aes(x=lon, y=lat, fill=value)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)+
  theme_minimal()

results <- compute_matrices(df = df)
dist <- results[[1]]
PI <- results[[2]]
inx <- results[[3]]
PI <- PI[inx, inx]
dist <- dist[inx, inx]
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

U <- evaluate_U(PI)

# Create a contour plot #

distances_eucl <- as.matrix(dist(df[inx, c("longitude", "latitude")], method = 'euclidean'))
cov_eucl <- exponential_covariance(c(1, max(distances_eucl)*10), distances_eucl)
load(paste0("Data/range_heatmap_", domain, ".RData"))

cov <- build_covariance_exponential(1, range, dist, P, U)
if (domain == "Sardinia") {
  points_high <- c(1800, 1521)
} else {
  points_high <- c(1,30)
}
points_high_inx <- sapply(points_high, FUN = function(x) which(inx == x))

df$covariances <- NA
df$covariances[inx] <- 0
df$covariances_euclidean1 <- NA
df$covariances_euclidean2 <- NA
df$covariances_euclidean1[inx] <- 0
df$covariances_euclidean2[inx] <- 0
df$covariances[inx] <- rowSums(sapply(points_high_inx, FUN = function(x) cov[x,]))
df$covariances_euclidean1[inx] <- sapply(points_high_inx[1], 
                                         FUN = function(x) cov_eucl[x,])
if(length(points_high_inx) > 1) {
  df$covariances_euclidean2[inx] <- sapply(points_high_inx[2], 
                                           FUN = function(x) cov_eucl[x,])
} else {
  df$covariances_euclidean2[inx] <- 0 
}

cov <- build_covariance_exponential(sill45, range45, dist, P, U)

cov_chol <- chol(cov)

simulation <- NULL
K <- 500

for (i in 1:K)
{
  set.seed(i)
  wn <- rnorm(B, 0, 1)
  simulation <- rbind(simulation,t(t(cov_chol) %*% wn))
}

simulation <- simulation + 
  matrix(rep(df$value[inx], times = K), nrow = K, ncol = length(inx), byrow = T) + mean(bias45)
std <- apply(simulation, 2, sd)

threshold <- 25:30

for (index in seq_along(threshold))
{
  df_temp <- df[inx, 1:4]
  
  over <- apply(simulation, MARGIN = 1, function(x) ifelse(x>threshold[index], 1, 0))
  mins <- NULL
  predictions <- (df$value[inx] - threshold[index])/std
  for (i in 1:K)
    mins <- c(mins, min(predictions[which(over[,i]==1)]))
  thre_temp <- quantile(mins, probs = 0.05)
  df_temp$hot_plus <- NA
  df_temp$hot_plus <- ifelse(predictions>thre_temp, 1, 0)
  
  under <- apply(simulation, MARGIN = 1, function(x) ifelse(x<threshold[index], 1, 0))
  maxs <- NULL
  predictions <- (df$value[inx] - threshold[index])/std
  for (i in 1:K)
    maxs <- c(maxs, max(predictions[which(under[,i]==1)]))
  thre_temp <- quantile(maxs, probs = 0.95)
  df_temp$hot_minus <- NA
  df_temp$hot_minus <- ifelse(predictions>=thre_temp, 1, 0)
  
  df_temp$hot_plot <- as.factor(rowSums(df_temp[, c("hot_plus", "hot_minus")]))
  
  col_name <- paste0("exceedance_", threshold[index])
  df[[col_name]] <- NA
  df[[col_name]][inx] <- df_temp$hot_plot
  
}

# Joint exceedance probability

coords <- c(10.3, 42.7)
coords_data <- df[inx,1:2]
n <- dim(coords_data)[1]
coordinates_index <- which.min(apply((coords_data - 
                                        matrix(as.numeric(rep(coords, times = n)), nrow = n, byrow = TRUE))^2, 
                                     MARGIN = 1, sum))

radius <- c(0, 10, 15, 20, 30, 50)*1e3

distances <- as.matrix(dist(df[inx, c("longitude", "latitude")]))

find_prob <- function(t_val, r, simulation)
{
  to_check <- which(distances[coordinates_index ,] <= r)
  if (length(to_check) > 1)
  {
    n <- dim(simulation)[1]
    sums <- rowSums(apply(simulation[,to_check], MARGIN = 1, function(x) 
      c(any(x>t_val), all(x>t_val))))
  }
  else
  {
    n <- dim(simulation)[1]
    sums <- rowSums(sapply(simulation[,to_check], function(x) 
      c(any(x>t_val), all(x>t_val))))
  }
  return(sums/n)
}

probs <- lapply(radius, function(x) sapply(threshold, FUN = function(y) 
  find_prob(y, x, simulation)))

exceedance_probs <- lapply(seq_along(probs[-1]), function(i) {
  temp <- NULL
  temp$value <- (c(t(probs[[i]])))
  temp$group <- rep(1:2, each = length(threshold))
  temp$level <- rep(threshold, times = 2)
  temp$plot_id <- rep(radius[i], each = length(threshold)*2)
  as.data.frame(temp)
}) %>%
  bind_rows()

base_line <- NULL
base_line$value <- c(t(probs[[1]]))
base_line$group <- rep(1:2, each = length(threshold))
base_line$level <- rep(threshold, times = 2)
base_line$plot_id <- rep(0, each = length(threshold)*2)
base_line <- as.data.frame(base_line)

ribbon_df <- exceedance_probs %>%
  pivot_wider(
    names_from    = group,
    values_from   = value,
    names_prefix = "g"
  ) %>%
  rename(
    ymin    = g1,
    ymax    = g2,
    radius  = plot_id,
    level   = level
  )

p <- ggplot() +
  geom_ribbon(
    data = ribbon_df,
    aes(
      x      = level,
      ymin   = ymin,
      ymax   = ymax,
      group = factor(radius),
      alpha = radius
    ),
    fill = third
  ) +
  geom_line(
    data = base_line,
    aes(x = level, y = value),
    color      = "black",
    linewidth = 1.5
  ) +
  scale_alpha_continuous(
    name  = "Radius",
    range = c(1.2, 0.2) 
  ) +
  labs(
    x = "Levels",
    y = "Exceedance probability"
  ) +
  theme_minimal()
p
ggsave(
  filename = paste0(output, "Exceedance_Probabilities_", as.character(year),
                    "_", as.character(rcp), ".pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

#### Euclidean framework ####

distances_eucl <- as.matrix(dist(df[inx,c("longitude", "latitude")], method = 'euclidean'))
cov <- exponential_covariance(c(sill45_eucl, range45_eucl), distances_eucl)

cov_chol <- chol(cov)

simulation <- NULL
K <- 500

for (i in 1:K)
{
  set.seed(i)
  wn <- rnorm(B, 0, 1)
  simulation <- rbind(simulation,t(t(cov_chol) %*% wn))
}

simulation <- simulation + 
  matrix(rep(df$value[inx], times = K), nrow = K, ncol = length(inx), byrow = T) + mean(bias45)
std <- apply(simulation, 2, sd)

for (index in seq_along(threshold))
{
  df_temp <- df[inx, 1:4]
  
  over <- apply(simulation, MARGIN = 1, function(x) ifelse(x>threshold[index], 1, 0))
  mins <- NULL
  predictions <- (df$value[inx] - threshold[index])/std
  for (i in 1:K)
    mins <- c(mins, min(predictions[which(over[,i]==1)]))
  thre_temp <- quantile(mins, probs = 0.05)
  df_temp$hot_plus <- NA
  df_temp$hot_plus <- ifelse(predictions>thre_temp, 1, 0)
  
  under <- apply(simulation, MARGIN = 1, function(x) ifelse(x<threshold[index], 1, 0))
  maxs <- NULL
  predictions <- (df$value[inx] - threshold[index])/std
  for (i in 1:K)
    maxs <- c(maxs, max(predictions[which(under[,i]==1)]))
  thre_temp <- quantile(maxs, probs = 0.95)
  df_temp$hot_minus <- NA
  df_temp$hot_minus <- ifelse(predictions>=thre_temp, 1, 0)
  
  df_temp$hot_plot <- as.factor(rowSums(df_temp[, c("hot_plus", "hot_minus")]))
  
  col_name <- paste0("exceedance_", threshold[index], "_Eucl")
  df[[col_name]] <- NA
  df[[col_name]][inx] <- df_temp$hot_plot
  
}

# Joint exceedance probability (Euclidean)

coords <- c(10.3, 42.7)
coords_data <- df[inx,1:2]
n <- dim(coords_data)[1]
coordinates_index <- which.min(apply((coords_data - 
                                        matrix(as.numeric(rep(coords, times = n)), nrow = n, byrow = TRUE))^2, 
                                     MARGIN = 1, sum))

radius <- c(0, 10, 15, 20, 30, 50)*1e3

find_prob <- function(t_val, r, simulation)
{
  to_check <- which(distances[coordinates_index ,] <= r)
  if (length(to_check) > 1)
  {
    n <- dim(simulation)[1]
    sums <- rowSums(apply(simulation[,to_check], MARGIN = 1, function(x) 
      c(any(x>t_val), all(x>t_val))))
  }
  else
  {
    n <- dim(simulation)[1]
    sums <- rowSums(sapply(simulation[,to_check], function(x) 
      c(any(x>t_val), all(x>t_val))))
  }
  return(sums/n)
}

probs <- lapply(radius, function(x) sapply(threshold, FUN = function(y) 
  find_prob(y, x, simulation)))

exceedance_probs <- lapply(seq_along(probs[-1]), function(i) {
  temp <- NULL
  temp$value <- (c(t(probs[[i]])))
  temp$group <- rep(1:2, each = length(threshold))
  temp$level <- rep(threshold, times = 2)
  temp$plot_id <- rep(radius[i], each = length(threshold)*2)
  as.data.frame(temp)
}) %>%
  bind_rows()

base_line <- NULL
base_line$value <- c(t(probs[[1]]))
base_line$group <- rep(1:2, each = length(threshold))
base_line$level <- rep(threshold, times = 2)
base_line$plot_id <- rep(0, each = length(threshold)*2)
base_line <- as.data.frame(base_line)

ribbon_df <- exceedance_probs %>%
  pivot_wider(
    names_from    = group,
    values_from   = value,
    names_prefix = "g"
  ) %>%
  rename(
    ymin    = g1,
    ymax    = g2,
    radius  = plot_id,
    level   = level
  )

# --- PLOT FINALE: Exceedance Probability RIBBON (Euclidean Framework) ---
p <- ggplot() +
  geom_ribbon(
    data = ribbon_df,
    aes(
      x      = level,
      ymin   = ymin,
      ymax   = ymax,
      group = factor(radius),
      alpha = radius
    ),
    fill = third
  ) +
  geom_line(
    data = base_line,
    aes(x = level, y = value),
    color      = "black",
    linewidth = 1.5
  ) +
  scale_alpha_continuous(
    name  = "Radius",
    range = c(1.2, 0.2) 
  ) +
  labs(
    x = "Levels",
    y = "Exceedance probability"
  ) +
  theme_minimal()
p
ggsave(
  filename = paste0(output, "Exceedance_Probabilities_", as.character(year),
                    "_", as.character(rcp), "_Eucl.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

save(file = paste0("Data/data_projection_", domain, "_", year, ".RData"), df)

