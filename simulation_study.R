library(sp)
library(gstat)
library(ncdf4)
library(ggplot2)
library(spatstat)
library(igraph)
library(sf)
library(progress)
library(ggmap)
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

library(ggmap)
register_stadiamaps("424ec34b-4e14-4acc-8aa9-9e8ab4f4fbc9", write = FALSE)
map <- get_stadiamap(bbox, zoom = 8, maptype = "alidade_smooth")

ggplot(df[which(!is.na(df$east) & !is.na(df$north) & !is.na(df$temperature)) ,], 
       aes(x=lon, y=lat, fill=temperature)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1) +
  theme_minimal()

ggmap(map, darken = c(.356,"white")) + 
  geom_tile(
    data = df[which(!is.na(df$east) & !is.na(df$north) & !is.na(df$temperature)), ],
    mapping = aes(x = lon, y = lat, fill = temperature)
  ) + 
  scale_fill_gradient(low = col2, high = col1) +
  theme_minimal()

df$value <- df$temperature
inx <- which(!is.na(df$east) & !is.na(df$north) & !is.na(df$value))
B <- length(inx)
results <- compute_matrices(df = df)
lines <- results[[1]]
dist <- results[[2]]
PI <- results[[3]]
inx <- results[[4]]

PI <- PI[inx, inx]
dist <- dist[inx, inx]

sources <- which(colSums(!is.na(PI)) == 0)
outlets <- which(rowSums(PI, na.rm=T) < 0.99)

g <- graph_from_adjacency_matrix(ifelse(is.na(dist), 0, 1), mode = "directed")
has_cycles <- !is_acyclic(g)
print("The built network has a cycle: ")
print(has_cycles)
plot_size <- 5
p <- plot_lin_net(lines, df[, 1:2], plot_size, map, inx, sources, outlets)
p

df_eucl <- df[inx,]
coordinates(df_eucl) <- c("lon", "lat")

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

distances <- initialize(dist, P)
updated <- TRUE
print("start composing distance object")
iter <- 1
while (updated)
{
  print(paste0("iteration ", iter))
  iter <- iter + 1
  updated <- FALSE
  distances <- lapply(distances, function(p) update(p,dist, P))
}

U <- evaluate_U(PI)

#### Evaluate positive definiteness of our model ####

sill <- 1
range <- 2e6
cov <- build_covariances(sill, range, distances, U, spherical_covariance)
image(cov)
eigen(cov)$values

#### Exponential model ####

sill <- 1
range <- 2e6
cov <- build_covariance_exponential(sill, range, dist, P, U)
image(cov)
eigen(cov)$values

#### Simulation study ####

values <- df$temperature[inx]
sill <- 1
B <- length(distances)
R <- 5
ex1 <- 5e5
ex2 <- 2e6
range <- seq(ex1, ex2, length = R)
K <- 100
l <- 15
model <- "Exponential"
covariances <- lapply(range, function(x) vector('list', K))
covariances_eucl <- lapply(range, function(x) vector('list', K))
sills <- matrix(0, nrow = K, ncol = R)
ranges <- matrix(0, nrow = K, ncol = R)
sills_eucl <- matrix(0, nrow = K, ncol = R)
ranges_eucl <- matrix(0, nrow = K, ncol = R)
frobenius_norm <- matrix(0, nrow = K, ncol = R)
frobenius_norm_eucl <- matrix(0, nrow = K, ncol = R)
kl_divergence <- matrix(0, nrow = K, ncol = R)
kl_divergence_eucl <- matrix(0, nrow = K, ncol = R)
errors <- vector("list", length = R)
errors_eucl <- vector("list", length = R)
MSE <- matrix(0, nrow = K, ncol = R)
MSE_eucl <- matrix(0, nrow = K, ncol = R)

coords <- df_eucl@coords
lat_mean <- mean(coords[,2])
km_per_deg_lon <- 111 * cos(lat_mean * pi / 180)
km_per_deg_lat <- 111

# Trasformo le coordinate in "pseudo-km"
coords_km <- coords
coords_km[,1] <- coords[,1] * km_per_deg_lon
coords_km[,2] <- coords[,2] * km_per_deg_lat

dist_eucl <- as.matrix(dist(coords_km))

simulation_distribution <- do.call(rbind, lapply(1:K, function(i)
{
  set.seed(i)
  rnorm(B,0,1)
}))

for (r in 1:R)
{
  if (model == "Exponential")
  {
    cov <- build_covariance_exponential(sill, range[r], dist, P, U)
  }
  else
  {
    cov <- build_covariances(sill, range[r], distances, U, fun = spherical_covariance)
  }
  
  errors_temp <- matrix(NA, nrow = K, ncol = B)
  errors_temp_eucl <- matrix(NA, nrow = K, ncol = B)

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
    set.seed(iter*r)
    test <- sample(seq_along(inx), length(inx)/3)
    train <- setdiff(1:B, test)
    
    covariance <- evaluate_covariance_penalization_train(simulation,
                                                         distances, train, l-1, U, return_fuv = T)
    sill_est <- covariance[[1]]
    covariance <- covariance[[2]]
    covariances[[r]][[iter]] <- covariance
    sills[iter, r] <- sill_est
    range_est <- fit_range(covariance, exponential_covariance,
                           max(covariance$dist)/3, sill_est)
    ranges[iter, r] <- range_est
    
    if(model == "Exponential")
    {
      covariance_est <- build_covariance_exponential(sill_est, range_est, dist, P, U)
    }
    else
    {
      covariance_est <- build_covariances(sill_est, range_est, distances, U,
                                          fun = spherical_covariance)
    }
    
    if(!all(eigen(covariance_est)$values>1e-15))
    {
      covariance_est <- diag(rep(sill_est, B))
    }
    
    frobenius_norm[iter, r] <- sqrt(sum((covariance_est - cov)^2))
    kl_divergence[iter, r] <- 0.5 * (sum(diag(solve(cov) %*% covariance_est)) -
                                       nrow(covariance_est) + log(det(cov) / det(covariance_est)))
    
    predictions <- evaluate_test(simulation[train], covariance_est, train, test)
    
    errors_temp[iter, test] <- simulation[test] - predictions$preds
    MSE[iter, r] <- mean((predictions$preds - simulation[test])^2)
    
    # Eucliean framework
    
    cov_eucl <- calc_covariance(df_eucl[train,]@coords, simulation[train])
    
    covariances_eucl[[r]][[iter]] <- cov_eucl
    
    variogram_eucl <- calc_variogram(df_eucl[train,]@coords, simulation[train])
    
    params <- fit_covariance(variogram_eucl, exponential_kernel,
                            c(max(variogram_eucl$gamma), max(variogram_eucl$dist/3)))
    sill_eucl <- params[1]
    sills_eucl[iter, r] <- sill_eucl
    range_eucl <- params[2]
    ranges_eucl[iter, r] <- range_eucl
    
    if (params[2] == 0)
    {
      covariance_est_eucl <- diag(rep(sill_eucl, B))
    }
    else
    {
      covariance_est_eucl <- exponential_covariance(params,
                                                    as.matrix(dist_eucl))
    }
    
    frobenius_norm_eucl[iter, r] <- sqrt(sum((covariance_est_eucl - cov)^2))
    kl_divergence_eucl[iter, r] <- 0.5 * (sum(diag(solve(cov) %*% covariance_est_eucl)) -
                                            nrow(covariance_est_eucl) + log(det(cov) / det(covariance_est_eucl)))
    
    predictions_eucl <- evaluate_test(simulation[train],
                                      covariance_est_eucl, train, test)
    
    errors_temp_eucl[iter, test] <- simulation[test] - predictions_eucl$preds
    MSE_eucl[iter, r] <- mean((predictions_eucl$preds - simulation[test])^2)
    
  }
  errors[[r]] <- errors_temp
  errors_eucl[[r]] <- errors_temp_eucl

}

plot_size <- 5
line_size <- 2

data <- data.frame(
  value = c(rowMeans(sills), rowMeans(sills_eucl)),
  group = factor(c(rep("New framework", K),
                   rep("Euclidean framework", K)))
)

p <- ggplot(data, aes(x = value, fill = group)) +
  geom_density(aes(y = after_stat(density)), alpha = 0.5, position = "identity") +
  geom_vline(xintercept = sill, linetype = "dashed", color = "black") +
  facet_wrap(~ group, ncol = 1, scales = "fixed") +
  scale_fill_manual(values = c(third, col5)) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  xlim(0, 2 * sill) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = plot_size*3, face = "bold"),
    axis.text.y = element_text(size = plot_size*3),
    axis.text.x = element_text(size = plot_size*3)
  )

p

index <- 1

empirical_df <- do.call(
  rbind,
  lapply(seq_along(covariances[[index]]), function(i) {
    df <- covariances[[index]][[i]]
    df$group <- factor(i)
    df
  })
)

xx <- seq(
  0,
  max(sapply(covariances[[index]], function(x) max(x$dist))),
  length.out = 200
)
yy <- exponential_covariance(c(sill, range[index]), xx)
theoretical_df <- data.frame(dist = xx, gamma = yy)

p <- ggplot() +
  geom_line(
    data = empirical_df,
    aes(x = dist, y = gamma, group = group),
    colour = "grey"
  ) +
  geom_line(
    data = theoretical_df,
    aes(x = dist, y = gamma),
    colour = third,
    linewidth = line_size,
    linetype = "dashed"
  ) +
  scale_x_continuous(breaks = seq(0, max(theoretical_df$dist), length = 5),
                     labels = floor(seq(0, max(theoretical_df$dist)/1000, length = 5))) +
  coord_cartesian(ylim = c(-sill, 2 * sill)) +
  labs(
    x = "Network Distance",
    y = "C",
    title = NULL
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = plot_size*4),
    axis.title.y = element_text(size = plot_size*4),
    plot.title = element_text(size = plot_size*4),
    strip.text = element_text(size = plot_size*4),
    axis.text.x = element_text(size = plot_size*4),
    axis.text.y = element_text(size = plot_size*4)
  )

p

xx <- seq(
  0,
  max(sapply(covariances[[index]], function(x) max(x$dist))),
  length.out = l
)
empirical_df <- do.call(rbind,
                        lapply(1:K, function(x)
                        {
                          df <- data.frame(dist = xx, gamma =
                                             exponential_covariance(c(sills[x,index], ranges[x,index]), xx))
                          df$group <- factor(x)
                          df}))
xx <- seq(
  0,
  max(sapply(covariances[[index]], function(x) max(x$dist))),
  length.out = 200
)
yy <- exponential_covariance(c(sill, range[index]), xx)
theoretical_df <- data.frame(dist = xx, gamma = yy)

p <- ggplot() +
  geom_line(
    data = empirical_df,
    aes(x = dist, y = gamma, group = group),
    colour = "grey"
  ) +
  geom_line(
    data = theoretical_df,
    aes(x = dist, y = gamma),
    colour = third,
    linewidth = line_size,
    linetype = "dashed"
  ) +
  coord_cartesian(ylim = c(-sill, 2 * sill)) +
  labs(
    x = "Network Distance",
    y = "C",
    title = NULL
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = plot_size*4),
    axis.title.y = element_text(size = plot_size*4),
    plot.title = element_text(size = plot_size*4),
    strip.text = element_text(size = plot_size*4),
    axis.text.x = element_text(size = plot_size*4),
    axis.text.y = element_text(size = plot_size*4)
  )

p

empirical_df <- do.call(
  rbind,
  lapply(seq_along(covariances_eucl[[index]]), function(i) {
    df <- covariances_eucl[[index]][[i]]
    df$group <- factor(i)
    df
  })
)

p <- ggplot(empirical_df, aes(x = dist, y = gamma, group = group)) +
  geom_line(colour = grey) +                              # linee grigie
  coord_cartesian(ylim = c(0, 2 * sill)) +                 # y ∈ [0, 2*sill]
  labs(
    x = "Distance",
    y = expression(C),
    title = NULL
  ) +
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = plot_size*4),
    axis.title.y = element_text(size = plot_size*4),
    plot.title = element_text(size = plot_size*4),
    strip.text = element_text(size = plot_size*4),
    axis.text.x = element_text(size = plot_size*4),
    axis.text.y = element_text(size = plot_size*4)
  )
p

frobenius <- data.frame(
  error = c(frobenius_norm, frobenius_norm_eucl),
  range = rep(rep(range, each = K),2),
  process = factor(rep(c("Covariance function",
                         "Euclidean"),
                       each = K*R))
)

p <- ggplot(frobenius, aes(x = factor(range), y = error, fill = process, color = process)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(
    values = c(third, col5),
    labels = c(expression(hat(Sigma)), expression(hat(Sigma)[eucl]))
  ) +
  scale_color_manual(values = c(third, col5), guide = 'none') +
  labs(x = "Range (km in network distance)", y = "Difference in Frobenius norm", fill = "Framework") +
  scale_x_discrete(labels = function(x) as.character(as.integer(as.numeric(x) / 1000))) +
  ylim(0,quantile(frobenius$error, 0.95, na.rm = T)) +
  theme_minimal()+
  theme(
    legend.position    = "bottom",
    legend.key         = element_rect(),
    legend.title       = element_blank(),
    legend.text        = element_text(size = plot_size*4),
    axis.title.x       = element_text(size = plot_size*4),
    axis.title.y       = element_blank(),
    axis.text.x        = element_text(size = plot_size*4),
    axis.text.y        = element_text(size = plot_size*4),
    legend.key.width  = unit(plot_size*0.1, "cm"),
    legend.key.height = unit(plot_size*0.1,   "cm")
    
  )
p

kl <- data.frame(
  error = c(kl_divergence, kl_divergence_eucl),
  range = rep(rep(range, each = K),2),
  process = factor(rep(c("Covariance function",
                         "Euclidean"),
                       each = K*R))
)

p <- ggplot(kl, aes(x = factor(range), y = error, fill = process, color = process)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(
    values = c(third, col5),
    labels = c(expression(hat(Sigma)), expression(hat(Sigma)[eucl]))
  ) +
  scale_color_manual(values = c(third, col2), guide = 'none') +
  labs(x = "Range (km in network distance)", y = "Difference in KL norm", fill = "Framework") +
  scale_x_discrete(labels = function(x) as.character(as.integer(as.numeric(x) / 1000))) +
  ylim(0,quantile(kl$error, 0.95, na.rm = T)) +
  theme_minimal() +
  theme(
    legend.position    = "bottom",
    legend.key         = element_rect(),
    legend.title       = element_blank(),
    legend.text        = element_text(size = plot_size*4),
    axis.title.x       = element_text(size = plot_size*4),
    axis.title.y       = element_blank(),
    axis.text.x        = element_text(size = plot_size*4),
    axis.text.y        = element_text(size = plot_size*4),
    legend.key.width  = unit(plot_size*0.1, "cm"),
    legend.key.height = unit(plot_size*0.1,   "cm")
    
  )
p

MSE_temp <- data.frame(
  error = c(sqrt(MSE), sqrt(MSE_eucl)),
  range = rep(rep(range, each = K),2),
  process = factor(rep(c("Covariance function", "Euclidean"), each = K*R))
)

p <- ggplot(MSE_temp, aes(x = factor(range), y = error, fill = process, color = process)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(
    values = c(third, col5),
    labels = c(expression(hat(Sigma)), expression(hat(Sigma)[eucl]))
  ) +
  scale_color_manual(values = c(third, col5), guide = 'none') +
  labs(x = "Range (km in network distance)", y = "MSE", fill = "Framework") +
  scale_x_discrete(labels = function(x) as.character(as.integer(as.numeric(x) / 1000))) +
  ylim(0, max(MSE_temp$error)) +
  theme_minimal() +
  theme(
    legend.position    = "bottom",
    legend.key         = element_rect(),
    legend.title       = element_blank(),
    legend.text        = element_text(size = plot_size*3),
    axis.title.x       = element_text(size = plot_size*3),
    axis.title.y       = element_blank(),
    axis.text.x        = element_text(size = plot_size*3),
    axis.text.y        = element_text(size = plot_size*3),
    legend.key.width  = unit(plot_size*0.1, "cm"),
    legend.key.height = unit(plot_size*0.1,   "cm")
    
  )
p

