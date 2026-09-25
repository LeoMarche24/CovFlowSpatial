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

domain <- "Temp"
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

df_eucl <- df[inx,]
coordinates(df_eucl) <- c("longitude", "latitude")

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

#### Evaluate positive definiteness of our model ####

sill <- 1
range <- 2e6
cov <- build_covariances(sill, range, distances, U = U, fun = spherical_covariance)
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
ex1 <- 5e4
ex2 <- 2e5
range <- seq(ex1, ex2, length = R)
K <- 100
l <- 10
model <- "Exponential"
output <- paste0("Plots/Simulation/", model, "/", domain, "/")
dir.create(output, recursive = T, showWarnings = F)
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
dist_eucl <- as.matrix(dist(coords))

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
    
    covariance_eucl_temp <- calc_covariance(df_eucl[train,]@coords, simulation[train])
    
    covariances_eucl[[r]][[iter]] <- covariance_eucl_temp
    
    sill_eucl <- covariance_eucl_temp$gamma[1]
    
    range_eucl <- fit_range(covariance_eucl_temp, exponential_covariance,
                            max(covariance_eucl_temp$dist/3), sill_eucl)
    sills_eucl[iter, r] <- sill_eucl
    ranges_eucl[iter, r] <- range_eucl
    
    if (range_eucl == 0)
    {
      covariance_est_eucl <- diag(rep(sill_eucl, B))
    }
    else
    {
      covariance_est_eucl <- exponential_covariance(c(sill_eucl, range_eucl),
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
ggsave(
  filename = paste0(output, "Sills.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

for (index in seq_along(range))
{

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
ggsave(
  filename = paste0(output, "Covariance_New", as.character(index), ".pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

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
ggsave(
  filename = paste0(output, "Covariance_New_Fitted", as.character(index), ".pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

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
ggsave(
  filename = paste0(output, "Covariances_Eucl", as.character(index), ".pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)
}

MSE_estimate <- data.frame(MSE = c(unlist(lapply(1:R, function(x)
  unlist(lapply(covariances[[x]], function(y)
    sum((y$gamma - exponential_covariance(c(sill, range[x]), y$dist))^2))))),
  unlist(lapply(1:R, function(x)
    sapply(1:K, function(y) sum((rep(sills_eucl[y,x],
              length(covariances[[x]][[y]]$dist)) -
     exponential_kernel(c(sill, range[x]), covariances[[x]][[y]]$dist))^2))))),
  range = rep(rep(range, each = K), 2),
  process = c(rep("Covariance funct.", K*R),
              rep("Variance est.", K*R)))

expr_labels <- expression(
  paste("||", hat(C)(h) - C(h), "||"),
  paste("||", hat(C)[eucl](h) - C(h), "||")
)

p <- ggplot(MSE_estimate, aes(x = factor(range), y = MSE, fill = process,
                           color = process)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = c(third, col5), labels = expr_labels) +
  scale_color_manual(values = c(third, col5), guide = 'none') +
  labs(x = "Range (km in network distance)", y = "MSE of the estimator", fill = "Framework") +
  scale_x_discrete(labels = function(x) as.character(as.integer(as.numeric(x) / 1000))) +
  guides(
    fill = guide_legend(label = label_parsed)
  ) +
  ylim(0,quantile(MSE_estimate$MSE, 0.95, na.rm = T)) +
  theme_minimal() +
  theme(
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = plot_size*3),
    legend.spacing   = unit(plot_size * 0.5, "cm"),
    axis.title.x       = element_text(size = plot_size*3),
    axis.title.y       = element_blank(),
    axis.text.x        = element_text(size = plot_size*3),
    axis.text.y        = element_text(size = plot_size*3),
    legend.key.width  = unit(plot_size*0.1, "cm"),
    legend.key.height = unit(plot_size*0.1,   "cm")
  )

ggsave(
  filename = paste0(output, "MSE_Estimate.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

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
ggsave(
  filename = paste0(output, "Frobenius.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

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
ggsave(
  filename = paste0(output, "KL.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

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
ggsave(
  filename = paste0(output, "MSE.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)
