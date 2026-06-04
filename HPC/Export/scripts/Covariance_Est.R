# Libraries

library(sp)
library(igraph)
library(reshape2)

library(ggplot2)
library(labeling)
library(farver)
library(parallel)

#### Functions loading ####

source(file.path(getwd(), "Functions.R"))

#### Load necessary data ####

col1 <- "cyan"
col2 <- "darkblue"
third <- rgb(63/255, 128/255, 97/255)
fourth <- rgb(173/255, 38/255, 36/255)
col5 <- rgb(204/255, 0/255, 0/255)
grey <- "grey80"

domain <- "Sardinia"

load(paste0("input_data/Data_", domain, ".RData"))
load(paste0("input_data/Projections_", domain, ".RData"))
data_proj <- df[which(df$years == 2022), ]

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
point_size <- 4

l <- 15

#### rcp45 ####

output <- paste0("output_data/", domain, "/Plots/Real Data 45/")
delta <- 10

years <- unique(df$years)
len <- length(years)
res_total <- vector('list', len)
sill <- rep(0, len)
range <- rep(0, len)
sill_eucl <- rep(0, len)
range_eucl <- rep(0, len)

fuvs <- rep(0, len)
bias <- rep(0,len)

covariances_45 <- NULL
distances_45 <- NULL
variograms_45_eucl <- NULL
distances_45_eucl <- NULL

K <- 200

# Create the simulated vectors

B <- length(unique(df$lon)) * length(unique(df$lat))

for(y in 1:len)
{
  print(y)
  data_proj <- data[which(data$year==years[y] & data$RCP == 'rcp45'), c(1:5,8:9)]
  data_real <- df[which(df$year==years[y]), c(1:5,7:8)]
  names(data_real)[5] <- "value"
  results <- compute_matrices(df = data_real)
  dist <- results[[1]]
  PI <- results[[2]]
  inx <- results[[3]]
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
  bias[y] <- mean(data_real$value - projections)
  
  projections_adj <- projections + bias[y]
  
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
  n_cores <- as.numeric(Sys.getenv("PBS_NP"))
  if (is.na(n_cores) || n_cores < 1) {
    n_cores <- as.numeric(Sys.getenv("NCPUS"))
  }
  if (is.na(n_cores) || n_cores < 1) {
    n_cores <- 1
  }
  cat("Uso", n_cores, "core\n")
  cl <- makeForkCluster(n_cores)
  clusterExport(cl, c("dist", "P", "update"))
  updated <- TRUE

  while (updated) {
  results <- parLapply(cl, seq_along(distances), function(i, distances_now) {
    p <- distances_now[[i]]
    update(p, dist, P)
  }, distances_now = distances)

    distances <- lapply(results, `[[`, "p")
    updated <- any(vapply(results, function(x) x$updated, logical(1)))
  }

  stopCluster(cl)

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

save(sill, range, bias,
     sill_eucl, range_eucl, res_total, fuvs,
     covariances_45, distances_45,
     variograms_45_eucl, distances_45_eucl,
     file = paste0("output_data/", domain, "/Covariance_Data_45_", domain, ".RData"))

cat("Covariance estimation for rcp45 completed\n")