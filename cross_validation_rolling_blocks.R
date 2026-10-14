## Functions ##

source("functions_network.R")
source("functions_observations.R")

#### Data loading ####

domain <- "Sardinia"

load(paste0("Data/Data_", domain, ".RData"))
load(paste0("Data/Projections_", domain, ".RData"))

#### Cross-validation settings ####

years_total <- sort(unique(df$years))
rolling_window <- 10

l <- 15
n_holes <- 10
hole_radius <- 120e3
hole_separation <- hole_radius

#### Estimation of the yearly covariance functions ####

year_data <- vector("list", length(years_total))

fuvs <- rep(0, length(years_total))
covariances <- matrix(NA, nrow = length(years_total), ncol = l)
distances <- matrix(NA, nrow = length(years_total), ncol = l)
variograms_eucl <- matrix(NA, nrow = length(years_total), ncol = l+1)
distances_eucl <- matrix(NA, nrow = length(years_total), ncol = l+1)

for(y in 1:length(years_total))
{
  print(years_total[y])

  data_proj <- data[which(data$year == years_total[y] &
                            data$RCP == "rcp45"), c(1:5, 8:9)]
  data_real <- df[which(df$years == years_total[y]), c(1:5, 7:8)]
  names(data_real)[5] <- "value"

  results <- compute_matrices(df = data_real)
  dist <- results[[1]]
  PI <- results[[2]]
  inx <- results[[3]]

  data_real <- data_real[inx,]
  dist <- dist[inx, inx]
  PI <- PI[inx, inx]

  coords <- data_real[,1:2]
  coords_data <- data_proj[,1:2]
  n <- dim(coords_data)[1]

  nearest <- apply(coords, MARGIN = 1, function(x)
    which.min(apply((coords_data-
      matrix(as.numeric(rep(x, times = n)), nrow = n, byrow = TRUE))^2,
      MARGIN = 1, sum)))

  inx2 <- which(!is.na(data_proj$value[nearest]))
  data_real <- data_real[inx2,]
  projections <- data_proj$value[nearest][inx2]
  dist <- dist[inx2, inx2]
  PI <- PI[inx2, inx2]

  residuals <- data_real$value-projections
  P <- PI
  B <- dim(P)[1]

  for(i in 1:B)
  {
    for(j in 1:B)
    {
      if(i != j & !is.na(P[i,j]))
      {
        P[i,j] <- PI[i,j]/sqrt(sum(PI[,j], na.rm = T))
      }
    }
  }

  print("start composing distance object")
  paths <- compose_distances(dist, P)
  U <- evaluate_U(PI)

  covariance <- evaluate_covariance_penalization(
    residuals, paths, l = l-1, U, return_fuv = T)

  fuvs[y] <- covariance[[1]]
  covariance <- covariance[[2]]

  covariances[y,1:nrow(covariance)] <- covariance$gamma
  distances[y,1:nrow(covariance)] <- covariance$dist

  variogram <- calc_variogram(
    data_real[,c("longitude", "latitude")], residuals, bins = l)

  variograms_eucl[y,] <- variogram$gamma
  distances_eucl[y,] <- variogram$dist

  coords <- as.matrix(data_real[,c("longitude", "latitude")])

  year_data[[y]] <- list(
    data_real = data_real,
    projections = projections,
    residuals = residuals,
    dist = dist,
    PI = PI,
    P = P,
    U = U,
    coords = coords,
    distances_eucl = as.matrix(dist(coords))
  )
}

#### Spatially balanced validation holes ####

reference <- year_data[[rolling_window]]
points_per_hole <- rowSums(reference$distances_eucl <= hole_radius)
hole_candidates <- which(points_per_hole >= 40 & points_per_hole <= 60)

if(length(hole_candidates) < n_holes)
  stop("Not enough geometrically admissible validation holes.")

candidate_coords <- reference$coords[hole_candidates,,drop = F]
candidate_distances <- as.matrix(stats::dist(candidate_coords))
domain_center <- colMeans(candidate_coords)

# Start near the center of the domain, then repeatedly select the candidate
# farthest from its nearest selected center. No current or model information
# is used in the selection.
selected_candidates <- which.min(rowSums(
  (candidate_coords-matrix(
    domain_center, nrow(candidate_coords), 2, byrow = T))^2))

while(length(selected_candidates) < n_holes)
{
  remaining <- setdiff(1:nrow(candidate_coords), selected_candidates)
  distance_to_selected <- candidate_distances[
    remaining, selected_candidates, drop = F]
  nearest_selected <- apply(distance_to_selected, 1, min)
  selected_candidates <- c(
    selected_candidates, remaining[which.max(nearest_selected)])
}

hole_indices <- hole_candidates[selected_candidates]
selected_distances <- reference$distances_eucl[
  hole_indices, hole_indices, drop = F]

if(any(selected_distances[upper.tri(selected_distances)] < hole_separation))
  stop("The selected validation holes do not satisfy the separation rule.")

# Give the holes a stable geographical order from south to north.
hole_indices <- hole_indices[order(
  reference$data_real$lat[hole_indices],
  reference$data_real$lon[hole_indices])]

hole_centers <- data.frame(
  hole_id = 1:n_holes,
  lon = reference$data_real$lon[hole_indices],
  lat = reference$data_real$lat[hole_indices],
  longitude = reference$coords[hole_indices,1],
  latitude = reference$coords[hole_indices,2],
  reference_points = points_per_hole[hole_indices]
)

#### Rolling cross-validation ####

results <- NULL
point_results <- NULL

for(y in (rolling_window+1):length(years_total))
{
  history <- (y-rolling_window):(y-1)

  print(paste("predicting", years_total[y]))

  mean_covariance <- data.frame(
    dist = round(colMeans(distances[history,,drop = F], na.rm = T)),
    gamma = colMeans(covariances[history,,drop = F], na.rm = T),
    np = 1)
  mean_covariance <- na.omit(mean_covariance)

  sill <- mean(fuvs[history])
  range <- fit_range(
    mean_covariance, exponential_covariance,
    max(mean_covariance$dist)/3, sill)

  mean_variogram <- data.frame(
    dist = colMeans(distances_eucl[history,,drop = F], na.rm = T),
    gamma = colMeans(variograms_eucl[history,,drop = F], na.rm = T),
    np = 1)
  mean_variogram <- na.omit(mean_variogram)

  params_eucl <- fit_covariance(
    mean_variogram, exponential_kernel,
    c(max(mean_variogram$gamma), max(mean_variogram$dist)/3))

  current <- year_data[[y]]
  B <- length(current$residuals)

  covariance_network <- build_covariance_exponential(
    sill, range, current$dist, current$P, current$U)
  covariance_euclidean <- exponential_covariance(
    params_eucl, current$distances_eucl)

  for(h in 1:n_holes)
  {
    center <- which.min(
      (current$coords[,1]-hole_centers$longitude[h])^2+
        (current$coords[,2]-hole_centers$latitude[h])^2)

    distance_center <- sqrt(rowSums((current$coords-
      matrix(current$coords[center,], nrow = B, ncol = 2, byrow = T))^2))

    test <- which(distance_center <= hole_radius)
    train <- setdiff(1:B, test)

    predictions_network <- evaluate_test(
      current$residuals[train], covariance_network, train, test)
    predictions_euclidean <- evaluate_test(
      current$residuals[train], covariance_euclidean, train, test)

    for(model in c("Network", "Euclidean"))
    {
      if(model == "Network")
        predictions <- predictions_network
      if(model == "Euclidean")
        predictions <- predictions_euclidean

      residual_predictions <- as.numeric(predictions$preds)
      prediction_variance <- as.numeric(predictions$vars)
      temperature_predictions <-
        current$projections[test]+residual_predictions
      prediction_sd <- sqrt(prediction_variance)

      errors <- temperature_predictions-current$data_real$value[test]
      z <- (current$data_real$value[test]-temperature_predictions)/prediction_sd
      crps <- prediction_sd*(z*(2*pnorm(z)-1)+2*dnorm(z)-1/sqrt(pi))

      results <- rbind(results,
        data.frame(
          year = years_total[y],
          hole_id = h,
          model = model,
          center_lon = current$data_real$lon[center],
          center_lat = current$data_real$lat[center],
          n_test = length(test),
          MSE = mean(errors^2),
          CRPS = mean(crps)
        ))

      point_results <- rbind(point_results,
        data.frame(
          year = years_total[y],
          hole_id = h,
          model = model,
          distance_boundary_km = (hole_radius-distance_center[test])/1000,
          squared_error = errors^2,
          CRPS = crps
        ))
    }
  }
}

save(hole_centers, results, point_results,
     file = paste0(
       "Data/Cross_Validation_Rolling_Blocks_SST_", domain, ".RData"))
