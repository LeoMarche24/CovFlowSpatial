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
  scaled_distance <- h/range
  correlation <- ifelse(scaled_distance == 0, 1,
                        log1p(scaled_distance)/scaled_distance)
  return(sill * correlation)
}

exponential_kernel <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  return(sill * (1 - exp(-h/range)))
}

linear_kernel <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  sill * (1 - ifelse(h <= range, (1 - h / range), 0))
}

spherical_kernel <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  correlation <- ifelse(
    h <= range,
    1 - (1.5 * (h / range)) + (0.5 * (h / range)^3),
    0
  )
  sill * (1 - correlation)
}

mariah_kernel <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  scaled_distance <- h/range
  correlation <- ifelse(scaled_distance == 0, 1,
                        log1p(scaled_distance)/scaled_distance)
  return(sill * (1 - correlation))
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
  omega <- omega[connected==1, , drop = FALSE]
  
  cols <- which(colSums(omega) == 0)
  
  active_cols <- setdiff(seq_len(l), cols)
  omega <- omega[, active_cols, drop = FALSE]
  aux <- t(omega)%*%omega
  q <- t(omega)%*%(fuv - vars)
  nor <- max(abs(q))/fuv
  mins <- min(2*diag(abs(aux)) - rowSums(abs(aux)))
  lambda <- max(0, nor - mins)
  estimated <- solve(aux + lambda*diag(1, ncol(omega)))%*%q
  covariance <- numeric(l)
  covariance[active_cols] <- estimated
  
  covariance <- c(fuv, covariance)
  
  dists <- c(0, unlist(lapply(dists, mean, na.rm=T)))
  weights <- c(0, np/dists[-1])
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
            sum(mat[interval==interval_unique[k],2]) * (U[i,j] / sqrt(U[i,i]*U[j,j]))
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
  omega <- omega[connected==1, , drop = FALSE]
  
  cols <- which(colSums(omega) == 0)
  
  active_cols <- setdiff(seq_len(l), cols)
  omega <- omega[, active_cols, drop = FALSE]
  aux <- t(omega)%*%omega
  q <- t(omega)%*%(fuv - vars)
  nor <- max(abs(q))/fuv
  mins <- min(2*diag(abs(aux)) - rowSums(abs(aux)))
  lambda <- max(0, nor - mins)
  estimated <- solve(aux + lambda*diag(1, ncol(omega)))%*%q
  covariance <- numeric(l)
  covariance[active_cols] <- estimated
  
  covariance <- c(fuv, covariance)
  
  dists <- c(0, unlist(lapply(dists, mean, na.rm=T)))
  weights <- c(0, np/dists[-1])
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

evaluate_covariance_penalization_train_correctsill <- function(ma, paths, train,
                                                               l, U, cutoff = NULL,
                                                               sill)
{
  if (is.null(cutoff))
  {
    cutoff <- max(unlist(lapply(paths, function(x) max(unlist(lapply(x$lengths,
      function(y) max(y[,1])))))))
  }
  B <- length(train)
  lags <- seq(0, cutoff, length = l)
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
      vars[index] <- ((ma[i] - ma[j]) ^ 2)/2
      a <- paths[[i]]$lengths[[j]]
      b <- paths[[j]]$lengths[[i]]

      if (!is.null(a))
      {
        interval <- findInterval(a[,1], lags)
        for (k in unique(interval))
        {
          omega[index, k] <- omega[index, k] +
            sum(a[interval == k, 2]) * (U[j,i] / sqrt(U[i,i]*U[j,j]))
          np[k] <- np[k] + 1
          dists[[k]] <- rbind(dists[[k]], mean(a[interval == k, 1]))
        }
        connected[index] <- 1
      }

      if (!is.null(b))
      {
        interval <- findInterval(b[,1], lags)
        for (k in unique(interval))
        {
          omega[index, k] <- omega[index, k] +
            sum(b[interval == k, 2]) * (U[i,j] / sqrt(U[i,i]*U[j,j]))
          np[k] <- np[k] + 1
          dists[[k]] <- rbind(dists[[k]], mean(b[interval == k, 1]))
        }
        connected[index] <- 1
      }
      index <- index + 1
    }
  }

  vars <- vars[connected == 1]
  omega <- omega[connected == 1, , drop = FALSE]
  cols <- which(colSums(omega) == 0)

  active_cols <- setdiff(seq_len(l), cols)
  omega <- omega[, active_cols, drop = FALSE]

  aux <- t(omega) %*% omega
  q <- t(omega) %*% (sill - vars)
  nor <- max(abs(q))/sill
  mins <- min(2*diag(abs(aux)) - rowSums(abs(aux)))
  lambda <- max(0, nor - mins)
  estimated <- solve(aux + lambda*diag(1, ncol(omega))) %*% q
  covariance_bins <- numeric(l)
  covariance_bins[active_cols] <- estimated
  covariance <- c(sill, covariance_bins)

  dists <- c(0, unlist(lapply(dists, mean, na.rm = TRUE)))
  weights <- c(0, np/dists[-1])
  covariance <- data.frame(dist = dists, gamma = covariance, np = weights)
  covariance <- covariance[complete.cases(covariance[, c("dist", "gamma")]), ]
  covariance
}

####
# Simulation functions
####

build_covariances <- function(sill, range, paths, paths2 = NULL, U, fun, mat = NULL)
{
  B <- length(paths)
  covs1 <- matrix(0, B, B)
  for (i in 1:(B-1))
  {
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
  valid <- is.finite(h) & is.finite(c_empirical)
  if ("np" %in% names(covariance))
    valid <- valid & covariance$np > 0
  h <- h[valid]
  c_empirical <- c_empirical[valid]

  if (!length(h))
    stop("No non-empty distance bins are available for covariance fitting.")
  
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
  valid <- is.finite(h) & is.finite(gamma_empirical)
  if ("np" %in% names(empirical))
    valid <- valid & empirical$np > 0
  h <- h[valid]
  gamma_empirical <- gamma_empirical[valid]

  if (!length(h))
    stop("No non-empty distance bins are available for range fitting.")
  
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

calc_covariance <- function(coords, values, bins = 15) {
  coords <- as.matrix(coords)
  valid_observations <- complete.cases(coords) & is.finite(values)
  coords <- coords[valid_observations, , drop = FALSE]
  values <- values[valid_observations]

  if (length(values) < 2)
    stop("At least two complete observations are required.")

  dists <- as.matrix(dist(coords))
  upper <- upper.tri(dists)
  pair_distances <- dists[upper]
  z <- values - mean(values)
  pair_covariances <- outer(z, z, "*")[upper]
  positive <- pair_distances > 0

  if (!any(positive))
    stop("At least two distinct spatial locations are required.")

  pair_distances <- pair_distances[positive]
  pair_covariances <- pair_covariances[positive]
  bin_breaks <- seq(0, max(pair_distances), length.out = bins + 1)
  bin_index <- findInterval(
    pair_distances, bin_breaks, rightmost.closed = TRUE, all.inside = TRUE
  )
  cov_data <- data.frame(
    dist = c(0, rep(NA_real_, bins)),
    gamma = c(var(values), rep(NA_real_, bins)),
    np = c(length(values), integer(bins))
  )

  for (i in seq_len(bins)) {
    selected <- bin_index == i
    if (any(selected)) {
      cov_data$dist[i + 1] <- mean(pair_distances[selected])
      cov_data$gamma[i + 1] <- mean(pair_covariances[selected])
      cov_data$np[i + 1] <- sum(selected)
    }
  }

  cov_data
}

calc_variogram <- function(coords, values, bins = 15) {
  coords <- as.matrix(coords)
  valid_observations <- complete.cases(coords) & is.finite(values)
  coords <- coords[valid_observations, , drop = FALSE]
  values <- values[valid_observations]

  if (length(values) < 2)
    stop("At least two complete observations are required.")

  dists <- as.matrix(dist(coords))
  upper <- upper.tri(dists)
  pair_distances <- dists[upper]
  pair_semivariances <- (0.5 * outer(values, values, "-")^2)[upper]
  positive <- pair_distances > 0

  if (!any(positive))
    stop("At least two distinct spatial locations are required.")

  pair_distances <- pair_distances[positive]
  pair_semivariances <- pair_semivariances[positive]
  bin_breaks <- seq(0, max(pair_distances), length.out = bins + 1)
  bin_index <- findInterval(
    pair_distances, bin_breaks, rightmost.closed = TRUE, all.inside = TRUE
  )
  vario_data <- data.frame(
    dist = c(0, rep(NA_real_, bins)),
    gamma = c(0, rep(NA_real_, bins)),
    np = c(length(values), integer(bins))
  )

  for (i in seq_len(bins)) {
    selected <- bin_index == i
    if (any(selected)) {
      vario_data$dist[i + 1] <- mean(pair_distances[selected])
      vario_data$gamma[i + 1] <- mean(pair_semivariances[selected])
      vario_data$np[i + 1] <- sum(selected)
    }
  }

  vario_data
}
