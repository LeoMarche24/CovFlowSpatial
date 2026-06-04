# Libraries
library(sp)
library(igraph)
library(reshape2)
library(dplyr)
library(tidyr)

library(ggplot2)
library(labeling)
library(farver)
library(parallel)

#### Load functions ####

source(file.path(getwd(), "Functions.R"))

#### Load necessary data ####

col1 <- "cyan"
col2 <- "darkblue"
third <- rgb(63/255, 128/255, 97/255)
fourth <- rgb(173/255, 38/255, 36/255)
grey <- "grey80"

domain <- "Sardinia"

load(paste0("input_data/Covariance_Data_45_", domain, ".RData"))
load(paste0("input_data/Projections_", domain, ".RData"))

data_proj <- data[which(data$year == 2022), ]

lat_max <- max(data_proj$lat)
lat_min <- min(data_proj$lat)
lon_max <- max(data_proj$lon)
lon_min <- min(data_proj$lon)

lat_range <- lat_max - lat_min
lon_range <- abs(lon_max - lon_min)

lat_mean_rad <- ((lat_max + lat_min) / 2) * (pi / 180)

aspect_ratio <- lat_range / (lon_range * cos(lat_mean_rad))

plot_size <- 7
line_size <- 2
point_size <- 5

output_plot <- paste0("output_data/", domain, "/Plots/Projections/")
output_data <- paste0("output_data/", domain, "/")

threshold <- seq(25,30, by=1)

B <- length(unique(data$lon)) * length(unique(data$lat))

K <- 500
simulated_distribution <- do.call(rbind, lapply(1:K, function(i)
  {
  set.seed(i)
  rnorm(B,0,1)
}))

find_prob <- function(distances, threshold, r, simulation, coordinates_index)
{
  to_check <- which(distances[coordinates_index ,] <= r)
  if (length(to_check) > 1)
  {
    n <- dim(simulation)[1]
    sums <- rowSums(apply(simulation[,to_check], MARGIN = 1, function(x)
      c(any(x>threshold), all(x>threshold))))
  }
  else
  {
    n <- dim(simulation)[1]
    sums <- rowSums(sapply(simulation[,to_check], function(x)
      c(any(x>threshold), all(x>threshold))))
  }
  return(sums/n)
}


make_hotspots <- function(df, bias, cov, inx, eucl = FALSE)
{
  B <- length(inx)
  
  cov_chol <- chol(cov)

  simulation <- simulated_distribution[,1:B]
  simulation <- simulation %*% cov_chol
  
  simulation <- simulation + 
    matrix(rep(df$value[inx], times = K), nrow = K, ncol = length(inx), byrow = T) + bias
  
  des <- ifelse(eucl, "_Eucl", "")

  # Exceedance region

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

    col_name <- paste0("exceedance_", threshold[index], des)
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
  
  coords_data <- df[inx, c("longitude", "latitude")]
  distances <- as.matrix(dist(coords_data))

  probs <- lapply(radius, function(x) sapply(threshold, FUN = function(y)
    find_prob(distances, y, x, simulation, coordinates_index)))
  
  exceedance_probs <- lapply(seq_along(probs[-1]), function(i) {
    temp <- NULL
    temp$value <- c(t(probs[[i]]))
    temp$group <- rep(1:2, each = length(threshold))
    temp$level <- rep(threshold, times = 2)
    temp$plot_id <- rep(radius[i], each = length(threshold)*2)/1000
    as.data.frame(temp)
  }) %>%
    bind_rows()
  
  base_line <- NULL
  base_line$value <- probs[[1]][1,]
  base_line$group <- rep(1:2, each = length(threshold))
  base_line$level <- rep(threshold, times = 2)
  base_line$plot_id <- rep(0, each = length(threshold)*2)
  base_line <- as.data.frame(base_line)

  ribbon_df <- exceedance_probs %>%
    pivot_wider(
      names_from   = group,
      values_from  = value,
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
    x     = level,
    ymin  = ymin,
    ymax  = ymax,
    group = factor(radius),
    alpha = radius
  ),
  fill = third
  ) +
  scale_alpha_continuous(
  name  = "Radius",
  range = c(0.8, 0.2) # più piccolo il radius, più opaco
  ) +
  geom_line(
    data = base_line,
    aes(
      x   = level,
      y   = value,
      group = group,
      linetype = factor(group)
    ),
    color = "black",
    size = 1
    ) +
    labs(
      x = "Levels",
      y = "Probability"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = plot_size*4),
      axis.text.y = element_text(size = plot_size*4),
      axis.title = element_text(size = plot_size*4),
      strip.text = element_text(size = plot_size*4),
      legend.position = 'none'
    )

  ggsave(
    filename = paste0(output_plot, "Exceedance_Probabilities_", as.character(year), 
                      "_", as.character(rcp), des, ".pdf"),
    plot     = p,
    width    = plot_size,
    height   = plot_size * aspect_ratio,
    units    = "in",
    dpi      = 92,
    limitsize = FALSE
  )
  return(df)
}

make_analysis <- function(year, rcp)
{
  df <- data[which(data$year==year & data$RCP == rcp), c(1:5,8:9)]
  results <- compute_matrices(df = df)
  dist <- results[[1]]
  PI <- results[[2]]
  inx <- results[[3]]
  B <- length(inx)
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

  distances_eucl <- as.matrix(dist(df[inx, 6:7]))
  cov_eucl <- exponential_covariance(c(1, max(distances_eucl)*10), distances_eucl)
  load(paste0("input_data/range_heatmap_", domain, ".RData"))

  U <- evaluate_U(PI)
  cov <- build_covariance_exponential(1, range, dist, P, U)
  {
    if (domain == "Sardinia")
      points_high <- c(1800, 1521)
    else
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
  df$covariances_euclidean2[inx] <- sapply(points_high_inx[2],
                                           FUN = function(x) cov_eucl[x,])
  
  ## RCP Hotspot analysis ##
  
  sill <- ifelse(rcp=='rcp45', sill45, sill85)
  range <- ifelse(rcp=='rcp45', range45, range85)
  bias <- ifelse(rcp == 'rcp45', bias45, bias85)

  cov <- build_covariance_exponential(sill, range, dist, P, U)

  df <- make_hotspots(df, bias, cov, inx)
  
  #### Euclidean framework ####
  
  sill <- ifelse(rcp=='rcp45', sill45_eucl, sill85_eucl)
  range <- ifelse(rcp=='rcp45', range45_eucl, range85_eucl)

  cov <- as.matrix(exponential_covariance(c(sill, range),
                                          dist(df[inx,6:7]))) + diag(sill, B)
  
  df <- make_hotspots(df, bias, cov, inx, eucl = TRUE)
  
  save(df, file = paste0(output_data, "data_projection_", domain, "_", as.character(year), ".RData"))
}

#### RCP 45 - Year 2050 ####

year <- 2050
rcp <- 'rcp45'

make_analysis(year, rcp)

cat("RCP 45 - Year 2050 done\n")