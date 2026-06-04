####
# Covariances' models
####

linear_covariance <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  ifelse(h <= range, sill * (1 - h / range), 0)
}
spherical_covariance <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  ifelse(h <= range, sill * (1 - (1.5 * (h / range)) + (0.5 * (h / range)^3)), 0)
}
exponential_covariance <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  return(sill * exp(-h/range))
}
mariah_covariance <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  return(sill * ((log(h/range + 1))/(h/range)))
}

exponential_kernel <- function(params, h)
{
  range <- params[1]
  return(exp(-h/range))
}

linear_kernel <- function(params, h)
{
  range <- params[1]
  ifelse(h <= range, (1 - h / range), 0)
}

spherical_kernel <- function(params, h)
{
  range <- params[1]
  ifelse(h <= range, (1 - (1.5 * (h / range)) + (0.5 * (h / range)^3)), 0)
}

mariah_kernel <- function(params, h)
{
  range <- params[1]
  return(((log(h/range + 1))/(h/range)))
}

####
## Kriging functions #
####

evaluate_covariance_penalization <- function(ma, paths, l, U,
                                             cutoff = NULL, return_fuv = F)
{
  if (is.null(cutoff))
  {
    cutoff <- max(unlist(lapply(paths, function(x) max(unlist(lapply(x$lengths,
                                                                     function(y) max(y[,1])))))))
  }
  B <- length(paths)
  lags <- seq(0, ifelse(is.null(cutoff), maxs, cutoff), length = l)
  vars <- rep(0, (B*(B-1))/2)
  omega <- matrix(0, nrow = (B*(B-1))/2, ncol = l)
  dists <- vector("list", l)
  np <- rep(0, l)
  connected <- rep(0, (B*(B-1)/2))
  index <- 1
  for (i in 1:(B-1))
  {
    for (j in (i+1):B)
    {
      v <- ((ma[i] - ma[j]) ^ 2)/2
      vars[index] <- v
      a <- paths[[i]]$lengths[[j]]
      b <- paths[[j]]$lengths[[i]]
      if (!is.null(a))
      {
        mat <- a
        interval <- findInterval(mat[,1], lags)
        interval_unique <- unique(interval)
        for (k in seq_along(interval_unique))
        {
          omega[index, interval_unique[k]] <- omega[index, interval_unique[k]] +
            sum(mat[interval==interval_unique[k],2]) * (U[j,i] / sqrt(U[i,i]*U[j,j]))
          np[interval_unique[k]] <- np[interval_unique[k]] + 1
          dists[[interval_unique[k]]] <- rbind(dists[[interval_unique[k]]], mean(mat[interval==interval_unique[k],1]))
        }
        connected[index] <- 1
      }
      if (!is.null(b))
      {
        mat <- b
        interval <- findInterval(mat[,1], lags)
        interval_unique <- unique(interval)
        for (k in seq_along(interval_unique))
        {
          omega[index, interval_unique[k]] <- omega[index, interval_unique[k]] +
            sum(mat[interval==interval_unique[k],2]) * (U[i,j] / sqrt(U[i,i]*U[j,j]))
          np[interval_unique[k]] <- np[interval_unique[k]] + 1
          dists[[interval_unique[k]]] <- rbind(dists[[interval_unique[k]]], mean(mat[interval==interval_unique[k],1]))
        }
        connected[index] <- 1
      }
      index <- index + 1
    }
  }
  
  fuv <- mean(vars[which(!connected)])
  vars <- vars[connected==1]
  omega <- omega[connected==1,]
  
  cols <- which(colSums(omega) == 0)
  
  if (length(cols))
  {
    omega <- omega[,-cols]
    aux <- t(omega)%*%omega
    q <- t(omega)%*%(fuv - vars)
    nor <- max(abs(q))/fuv
    mins <- min(2*diag(abs(aux)) - rowSums(abs(aux)))
    lambda <- nor - mins
    covariance <- c(solve(aux +
                            lambda*diag(1,l-length(cols)))%*%q,
                    rep(0, length(cols)))
  }
  else
  {
    aux <- t(omega)%*%omega
    q <- t(omega)%*%(fuv - vars)
    nor <- max(abs(q))/fuv
    mins <- min(2*diag(abs(aux)) - rowSums(abs(aux)))
    lambda <- nor - mins
    covariance <- c(solve(aux +
                            lambda*diag(1,l-length(cols)))%*%q,
                    rep(0, length(cols)))
  }
  
  covariance <- c(fuv, covariance)
  
  dists <- c(0, unlist(lapply(dists, mean, na.rm=T)))
  weights <- np/dists
  if (any(is.na(weights)) || any(is.na(covariance)))
  {
    covariance <- data.frame(dist = dists,
                             gamma = covariance,
                             np = weights)[-unique(c(which(is.na(weights)),
                                                     which(is.na(covariance)))) ,]
  }
  else
  {
    covariance <- data.frame(dist = dists,
                             gamma = covariance,
                             np = weights)
  }
  if(return_fuv)
  {
    return(list(fuv, covariance))
  }
  else
    return(covariance)
}

evaluate_covariance_penalization_train <- function(ma, paths, train, l, U,
                                                   cutoff = NULL, return_fuv = F)
{
  if (is.null(cutoff))
  {
    cutoff <- max(unlist(lapply(paths, function(x) max(unlist(lapply(x$lengths,
               function(y) max(y[,1])))))))
  }
  B <- length(train)
  lags <- seq(0, ifelse(is.null(cutoff), maxs, cutoff), length = l)
  vars <- rep(0, (B*(B-1))/2)
  omega <- matrix(0, nrow = (B*(B-1))/2, ncol = l)
  dists <- vector("list", l)
  np <- rep(0, l)
  connected <- rep(0, (B*(B-1)/2))
  index <- 1
  for (in1 in 1:(B-1))
  {
    for (in2 in (in1+1):B)
    {
      i <- train[in1]
      j <- train[in2]
      v <- ((ma[i] - ma[j]) ^ 2)/2
      vars[index] <- v
      a <- paths[[i]]$lengths[[j]]
      b <- paths[[j]]$lengths[[i]]
      if (!is.null(a))
      {
        mat <- a
        interval <- findInterval(mat[,1], lags)
        interval_unique <- unique(interval)
        for (k in seq_along(interval_unique))
        {
          omega[index, interval_unique[k]] <- omega[index, interval_unique[k]] +
            sum(mat[interval==interval_unique[k],2]) * (U[j,i] / sqrt(U[i,i]*U[j,j]))
          np[interval_unique[k]] <- np[interval_unique[k]] + 1
          dists[[interval_unique[k]]] <- rbind(dists[[interval_unique[k]]], mean(mat[interval==interval_unique[k],1]))
        }
        connected[index] <- 1
      }
      if (!is.null(b))
      {
        mat <- b
        interval <- findInterval(mat[,1], lags)
        interval_unique <- unique(interval)
        for (k in seq_along(interval_unique))
        {
          omega[index, interval_unique[k]] <- omega[index, interval_unique[k]] +
            sum(mat[interval==interval_unique[k],2]) * (U[j,i] / sqrt(U[i,i]*U[j,j]))
          np[interval_unique[k]] <- np[interval_unique[k]] + 1
          dists[[interval_unique[k]]] <- rbind(dists[[interval_unique[k]]], mean(mat[interval==interval_unique[k],1]))
        }
        connected[index] <- 1
      }
      index <- index + 1
    }
  }
  
  fuv <- mean(vars[which(rowSums(omega) == 0)])
  vars <- vars[connected==1]
  omega <- omega[connected==1,]
  
  cols <- which(colSums(omega) == 0)
  
  if (length(cols))
  {
    omega <- omega[,-cols]
    aux <- t(omega)%*%omega
    q <- t(omega)%*%(fuv - vars)
    nor <- max(abs(q))/fuv
    mins <- min(2*diag(abs(aux)) - rowSums(abs(aux)))
    lambda <- nor - mins
    covariance <- c(solve(aux +
                            lambda*diag(1,l-length(cols)))%*%q,
                    rep(0, length(cols)))
  }
  else
  {
    aux <- t(omega)%*%omega
    q <- t(omega)%*%(fuv - vars)
    nor <- max(abs(q))/fuv
    mins <- min(2*diag(abs(aux)) - rowSums(abs(aux)))
    lambda <- nor - mins
    covariance <- c(solve(aux +
                            lambda*diag(1,l-length(cols)))%*%q,
                    rep(0, length(cols)))
  }
  
  covariance <- c(fuv, covariance)
  
  dists <- c(0, unlist(lapply(dists, mean, na.rm=T)))
  weights <- np/dists
  if (any(is.na(weights)) || any(is.na(covariance)))
  {
    covariance <- data.frame(dist = dists,
                             gamma = covariance,
                             np = weights)[-unique(c(which(is.na(weights)),
                                                     which(is.na(covariance)))) ,]
  }
  else
  {
    covariance <- data.frame(dist = dists,
                             gamma = covariance,
                             np = weights)
  }
  if(return_fuv)
  {
    return(list(fuv, covariance))
  }
  else
    return(covariance)
}

####
# Simulation functions
####

build_covariances <- function(sill, range, paths, U, fun)
{
  B <- length(paths)
  covs1 <- matrix(0, B, B)
  for (i in 1:(B-1))
  {
    if (!is.null(paths[[i]]$lengths[[i]]))
    {
      sums <- 0  
      mat <- paths[[i]]$lengths[[i]]
      l <- dim(mat)[1]
      for (k in 1:l) { 
        sums <- sums + mat[k,2]*fun(c(sill, range), mat[k,1])
      }
      covs1[i,i] <- sums
    }
    for (j in (i+1):B)
    {
      sums_1 <- sums_2 <- 0
      if (!is.null(paths[[i]]$lengths[[j]]))
      {
        mat <- paths[[i]]$lengths[[j]]
        l <- dim(mat)[1]
        for (k in 1:l)
        {
          sums_1 <- sums_1 + mat[k,2]*fun(c(sill, range), mat[k,1])
        }
        covs1[j,i] <- sums_1
      }
      if (!is.null(paths[[j]]$lengths[[i]]))
      {
        mat <- paths[[j]]$lengths[[i]]
        l <- dim(mat)[1]
        for (k in 1:l)
        {
          sums_2 <- sums_2 + mat[k,2]*fun(c(sill, range), mat[k,1])
        }
        covs1[i,j] <- sums_2
      }
    }
  }
  covs1 <- covs1 * U
  D <- diag(sqrt(1/diag(U)))
  covs1 <- D %*% covs1 %*% D
  covs1 <- covs1+t(covs1)+diag(sill,B)
  return(covs1)
}

build_covariance_exponential <- function(sill, range, dist, P, U)
{
  P[is.na(P)] <- 0
  dist[is.na(dist)] <- 0
  B <- dim(P)[1]
  cov <- solve(diag(1,B) - P*exp(-dist/range))
  D <- matrix(1/diag(cov), nrow = B, ncol = B)
  cov <- cov*D
  cov <- cov * U
  D <- diag(sqrt(1/diag(U)))
  cov <- D %*% cov %*% D
  cov <- cov + t(cov)
  diag(cov) <- 1
  return(sill*cov)
}

create_process_covariance <- function(seed, cov)
{
  B <- dim(cov)[1]
  set.seed(seed)
  wn <- rnorm(B, 0, 1)
  L <- chol(cov)
  x <- t(L) %*% wn
  return(x)
}

fit_covariance <- function(covariance, fun, initials)
{
  # Extract distances and empirical semivariances
  h <- covariance$dist
  c_empirical <- covariance$gamma
  
  # Objective function to minimize (sum of squared errors)
  obj_fun <- function(params) {
    c_model <- fun(params, h)
    sum((c_empirical - c_model)^2)
  }
  
  # Fit using optim
  fit <- optim(par = initials, fn = obj_fun, method = "L-BFGS-B",
               lower = c(0.1, 0.1))  # set sensible lower bounds
  
  # Return fit details
  fit$par
}

fit_range <- function(empirical, fun, initial, param_s)
{
  # Extract distances and empirical semivariances
  h <- empirical$dist
  gamma_empirical <- empirical$gamma
  
  # Objective function to minimize (sum of squared errors)
  obj_fun <- function(param) {
    gamma_model <- fun(c(param_s, param), h)
    sum((gamma_empirical - gamma_model)^2)
  }
  
  # Fit using optim
  fit <- optim(par = initial, fn = obj_fun, method = "Brent",
               lower = 0.1, upper = max(h))  # set sensible lower bounds
  
  # Return fit details
  fit$par
}

evaluate_test <- function(values, covs, train, test)
{
  sigma <- covs[train,train]
  D <- length(test)
  B <- length(train)
  ones <- rep(1, B)
  predictions <- matrix(0, 1, D)
  vars <- matrix(0, 1, D)
  for (i in 1:D)
  {
    cc <- covs[test[i], train]
    if (qr(sigma)$rank < ncol(sigma))
    {
      sigma <- sigma + diag(rep(1e-6, nrow(sigma)))
    }
    lam <- solve(sigma, cbind(cc, ones))
    ViX <- lam[,-1]
    skwts <- lam[,1]
    beta <- solve(t(ones) %*% ViX, t(ViX) %*% values)
    predictions[i] <- sum(beta) + t(skwts) %*% (values - sum(beta))
    
    Q <- t(1) - t(ViX) %*% cc
    vars[i] <- covs[1,1] - apply(cbind(cc*skwts), 2, sum) + 
      apply(Q * solve(t(ones) %*% ViX, Q), 2, sum)  
  }
  
  return(data.frame(preds = t(predictions), vars = t(vars)))
}

calc_covariance <- function(coords, values, bins = 10) {
  lat_mean <- mean(coords[,2])
  km_per_deg_lon <- 111 * cos(lat_mean * pi / 180)
  km_per_deg_lat <- 111
  
  coords_km <- coords
  coords_km[,1] <- coords[,1] * km_per_deg_lon
  coords_km[,2] <- coords[,2] * km_per_deg_lat
  
  dists <- as.matrix(dist(coords_km))
  
  # Center values
  z <- values - mean(values)
  cross_prods <- outer(z, z, "*")
  
  max_dist <- max(dists)
  bin_breaks <- seq(0, max_dist, length.out = bins + 1)
  cov_data <- data.frame(dist = numeric(bins),
                         cov  = numeric(bins),
                         np   = numeric(bins))
  
  for (i in 1:bins) {
    indices <- which(dists >= bin_breaks[i] & dists < bin_breaks[i+1], arr.ind = TRUE)
    if (nrow(indices) > 0) {
      valid <- indices[indices[,1] < indices[,2], , drop = FALSE]
      if (nrow(valid) > 0) {
        cov_data$dist[i] <- (bin_breaks[i] + bin_breaks[i+1]) / 2
        cov_data$gamma[i]  <- mean(cross_prods[valid])
        cov_data$np[i]   <- nrow(valid)
      }
    }
  }
  cov_data$gamma[1] <- var(values)
  return(cov_data)
}

calc_variogram <- function(coords, values, bins = 10) {
  lat_mean <- mean(coords[,2])
  km_per_deg_lon <- 111 * cos(lat_mean * pi / 180)
  km_per_deg_lat <- 111
  
  coords_km <- coords
  coords_km[,1] <- coords[,1] * km_per_deg_lon
  coords_km[,2] <- coords[,2] * km_per_deg_lat
  
  dists <- as.matrix(dist(coords_km))
  
  # Center values
  diffs <- outer(values, values, "-")
  sq_diffs <- 0.5 * (diffs^2)
  
  max_dist <- max(dists)
  bin_breaks <- seq(0, max_dist, length.out = bins + 1)
  vario_data <- data.frame(dist = numeric(bins),
                           gamma = numeric(bins),
                           np    = numeric(bins))
  
  for (i in 1:bins) {
    indices <- which(dists >= bin_breaks[i] & dists < bin_breaks[i+1], arr.ind = TRUE)
    if (nrow(indices) > 0) {
      valid <- indices[indices[,1] < indices[,2], , drop = FALSE]
      if (nrow(valid) > 0) {
        vario_data$dist[i] <- (bin_breaks[i] + bin_breaks[i+1]) / 2
        vario_data$gamma[i] <- mean(sq_diffs[valid]) / 2
        vario_data$np[i]   <- nrow(valid)
      }
    }
  }
  return(vario_data)
}
