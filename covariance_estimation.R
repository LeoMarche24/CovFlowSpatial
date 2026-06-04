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

domain <- "Temp"

load(paste0("Data/Data_", domain, ".RData"))
load(paste0("Data/Projections_", domain, ".RData"))

#### Extraction of the residuals ####

#### rcp45 ####

plot_size <- 7
line_size <- 2
point_size <- 4

years <- unique(df$years)
len <- length(years)
res_total <- vector('list', len)
sill <- rep(0, len)
range <- rep(0, len)
sill_eucl <- rep(0, len)
range_eucl <- rep(0, len)

fuvs <- rep(0, len)
bias45 <- rep(0,len)

covariances_45 <- NULL
distances_45 <- NULL
variograms_45_eucl <- NULL
distances_45_eucl <- NULL

l <- 15

B <- length(unique(df$lon)) * length(unique(df$lat))

for(y in 1:len)
{
  print(y)
  data_proj <- data[which(data$year==years[y] & data$RCP == 'rcp45'), c(1:5,8:9)]
  data_real <- df[which(df$year==years[y]), c(1:5,7:8)]
  names(data_real)[5] <- "value"
  results <- compute_matrices(df = data_real)
  dist <- results[[2]]
  PI <- results[[3]]
  inx <- results[[4]]
  data_real <- data_real[inx,]
  dist <- dist[inx, inx]
  PI <- PI[inx, inx]
  
  if(all(is.na(PI)))
    Error("All transition probabilities are NA")
  
  coords <- data_real[, 1:2]
  
  # I want a vector that for each location, find the nearest location based on data's locations
  
  coords_data <- data_proj[,1:2]
  n <- dim(coords_data)[1]
  
  nearest <- apply(coords, MARGIN = 1, function(x) 
    which.min(apply((coords_data - 
                       matrix(as.numeric(rep(x, times = n)), nrow = n, byrow = TRUE))^2, 
                    MARGIN = 1, sum)))
  
  inx2 <- which(!is.na(data_proj$value[nearest]))
  data_real <- data_real[inx2,]
  projections <- data_proj$value[nearest][inx2]
  bias45[y] <- mean(data_real$value - projections)
  
  projections_adj <- projections + bias45[y]
  
  res <- data_real$value - projections_adj
  res_total[[y]] <- res
  dist <- dist[inx2, inx2]
  PI <- PI[inx2, inx2]
  
  P <- PI
  B <- dim(P)[1]
  
  for(i in 1:B)
  {
    for(j in 1:B)
    {
      if(i!=j & !is.na(P[i,j]))
      {
        P[i,j] <- PI[i,j]/sqrt(sum(PI[,j], na.rm=T))
      }
    }
  }
  distances <- initialize(dist, P)
  updated <- TRUE
  print("start composing distance object")
  while (updated)
  {
    updated <- FALSE
    distances <- lapply(distances, function(p) update(p,dist, P))
  }
  
  U <- evaluate_U(PI)
  
  covariance <- evaluate_covariance_penalization(res, distances,
                                                 l = l-1, U, return_fuv = T)
  
  fuv <- covariance[[1]]
  covariance <- covariance[[2]]
  fuvs[y] <- fuv
  
  covariances_45 <- rbind(covariances_45,
                          c(covariance$gamma, rep(0, l-length(covariance$gamma))))
  distances_45 <- rbind(distances_45,
                        c(covariance$dist, rep(NA, l-length(covariance$dist))))
  
  sill[y] <- mean(fuvs[1:y])
  
  mean_covariance_temp <-
    data.frame(dist = round(colMeans(distances_45, na.rm = T)),
               gamma = colMeans(covariances_45),
               np = 1)
  mean_covariance_temp <- na.omit(mean_covariance_temp)
  
  range[y] <- fit_range(mean_covariance_temp, exponential_covariance,
                        max(mean_covariance_temp$dist)/3, sill[y])

  variogram <- calc_variogram(data_real[,c('longitude', 'latitude')], res)
  
  variograms_45_eucl <- rbind(variograms_45_eucl,
                              variogram$gamma)
  distances_45_eucl <- rbind(distances_45_eucl,
                                  variogram$dist)
  
  mean_variogram_temp <-
    data.frame(dist = round(colMeans(distances_45_eucl, na.rm = T)),
               gamma = colMeans(variograms_45_eucl),
               np = 1)
  mean_variogram_temp <- na.omit(mean_variogram_temp)
  params_temp <- fit_covariance(mean_variogram_temp, exponential_kernel,
                                c(max(variogram$gamma), max(variogram$dist)/3))
  
  sill_eucl[y] <- params_temp[1]
  range_eucl[y] <- params_temp[2]
  
}

#### Assess the final estimation ####

plot(years, unlist(lapply(res_total, mean)), pch = 16, col = third)
plot(bias45, pch = 16, col = third)

mean_covariance_45 <- data.frame(dist = round(colMeans(distances_45, na.rm = T)), 
                                 gamma = colMeans(covariances_45), 
                                 np = 1)

covariances_45 <- melt(covariances_45)
colnames(covariances_45) <- c("Row", "Column", "Value")
covariances_45$Distance <- 
  round(colMeans(distances_45, na.rm = T)[covariances_45$Column])

ggplot(covariances_45, aes(x = factor(Distance), y = Value, group = Row)) + 
  geom_line(color = "black") +  # Set color for the lines
  labs(x = "Distance", y = "Values") +
  theme_minimal()

covariances_45 <- 
  covariances_45[which(!is.na(covariances_45$Distance)),]

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
              fill = third, alpha = 0.6) +     # Central band
  geom_ribbon(aes(ymin = q05, ymax = q95),
              fill = third, alpha = 0.3) +     # External band
  geom_line(data = mean_covariance_45,
            aes(x = dist, y = gamma),
            color = "black", linewidth = 1.2) +     # variogramma medio
  labs(x = "Distance (meters)", y = "Semi-variogram") +
  ylim(c(-0.5, 1.1)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = plot_size*2),
    axis.title = element_text(size = plot_size*2),
    legend.text = element_text(size = plot_size*2),
    legend.title = element_text(size = plot_size*2),
    strip.text = element_text(size = plot_size*2)
  )
p

mean_covariance_45_eucl <- data.frame(dist = colMeans(distances_45_eucl), 
                                     gamma = colMeans(variograms_45_eucl))

variograms_45_eucl <- melt(variograms_45_eucl)
colnames(variograms_45_eucl) <- c("Row", "Column", "Value")
variograms_45_eucl$Distance <- 
  colMeans(distances_45_eucl)[variograms_45_eucl$Column]

ggplot(variograms_45_eucl, aes(x = factor(Distance), y = Value, group = Row)) + 
  geom_line(color = "black") +  # Set color for the lines
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
            color = "black", linewidth = 1.2) +     # variogramma medio
  labs(x = "Distance (meters)", y = "Semi-variance") +
  coord_cartesian(ylim = c(0, 2)) +   # ZOOM senza eliminare i ribbon
  theme_minimal() +
  theme(
    axis.text = element_text(size = plot_size*2),
    axis.title = element_text(size = plot_size*2),
    legend.text = element_text(size = plot_size*2),
    legend.title = element_text(size = plot_size*2),
    strip.text = element_text(size = plot_size*2)
  )
p

sill45 <- mean(fuvs)

range45 <- fit_range(mean_covariance_45, exponential_covariance,
                     max(mean_covariance_45$dist)/3, sill45)

mean_variogram_temp$np <- 1
res <- fit_covariance(mean_variogram_temp, exponential_kernel,
                     c(max(mean_variogram_temp$gamma), max(mean_variogram_temp$dist)/3))

sill45_eucl <- res[1]
range45_eucl <- res[2]

#### Saving ####

save(sill45, range45, sill45_eucl, range45_eucl, bias45,
     file = paste0("Data/Covariance_Data_45_", domain, ".RData"))

