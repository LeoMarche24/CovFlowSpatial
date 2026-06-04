#### Moving average functions ####

exponential_covariance <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  return(sill * exp(-h/range))
}

spherical_covariance <- function(params, h)
{
  sill <- params[1]
  range <- params[2]
  cov <- ifelse(h < range, sill * (1 - 1.5*(h/range) + 0.5*(h/range)^3), 0)
  return(cov)
}

exponential_kernel <- function(params, h)
{
  range <- params[2]
  return(exp(-h/range))
}

spherical_kernel <- function(params, h)
{
  range <- params[2]
  kernel <- ifelse(h < range, (1 - 1.5*(h/range) + 0.5*(h/range)^3), 0)
  return(kernel)
}

#### Construction of the network ####

normalize_column <- function(col)
{
  non_na_col <- col[!is.na(col)]
  if(length(non_na_col) < 2)
    normalized_col <- non_na_col
  else
    normalized_col <- non_na_col/sum(non_na_col)
  col[!is.na(col)] <- normalized_col
  return(col)
}

generate_line_id <- function(coords)
{
  paste(coords[, 1], coords[, 2])
}

compute_matrices <- function(df)
{
  coord_sub <- data.frame(longitude = df$longitude, latitude = df$latitude)
  nor <- df$north
  eas <- df$east
  nor[is.na(nor)] <- 0
  eas[is.na(eas)] <- 0
  valid <- which(!is.na(nor) & !is.na(eas) & !is.na(df$value) & ((abs(nor)>0) | (abs(eas)>0)))
  dist_eucl <- as.matrix(dist(as.matrix(coord_sub), method = 'euclidean'))
  B <- dim(df)[1]
  dist_mat <- matrix(NA, B,B)
  PI <- matrix(NA, B, B)
  
  extremes <- c(max(df$lon, na.rm = TRUE), min(df$lon, na.rm = TRUE),
                max(df$lat, na.rm = TRUE), min(df$lat, na.rm = TRUE))
  
  for (i in 1:B)
  {
    if (! (i %in% valid))
      next
    vel <- c(nor[i], eas[i])
    if(is.na(vel[1]) & is.na(vel[2]))
      next
    vel <- ifelse(is.na(vel), 0, vel)
    angle <- atan2(vel[1], vel[2])
    lon0 <- df$lon[i]
    lat0 <- df$lat[i]
    aux <- (lon0 %in% extremes) + (lat0 %in% extremes)
    # I take the nearest and remove that at 0 distance (the point itself)
    
    if(aux == 0)
    {
      nearest <- order(dist_eucl[i,])[2:9]
    }
    else if(aux == 1)
    {
      nearest <- order(dist_eucl[i,])[2:6]
    }
    else if(aux == 2)
    {
      nearest <- order(dist_eucl[i,])[2:4]
    }
    else
    {
      nearest <- NULL
    }
    
    if (is.null(nearest))
      next
    
    directions <- cbind(unlist(lapply(nearest, function(j)
    {
      lon1 <- df$lon[j]
      lat1 <- df$lat[j]
      angle_j <- atan2(lat1 - lat0, lon1 - lon0)
    })), nearest)
    directions <- directions[!duplicated(directions[, 1]), ]
    if (is.null(nrow(directions)))
    {
      directions <- rbind(directions, directions)
    }
    if (any(directions[,1] == pi))
    {
      row <- which(directions[,1] == pi)
      directions <- rbind(directions, c(-pi, directions[row,2]))
    }
    directions <- directions[order(abs(directions[,1] - angle))[1:2],]
    tol <- 1e-12
    mask_sq <- (cos(directions[, 1]) * vel[2] >= -tol) &
      (sin(directions[, 1]) * vel[1] >= -tol)
    
    if(nrow(directions) == 0)
      next
    
    D <- do.call(cbind, lapply(directions[,1], function(d) {
      c(sin(d), cos(d))
    }))
    
    r <- det(D)
    if(abs(r) < 1e-10)
    {
      directions <- directions[order(abs(directions[,1] - angle))[1],]
      is_east <- abs(directions[1]) < pi/4 || abs(directions[1]) > 3*pi/4
      if(is_east & mask_sq[1] & (directions[2] %in% valid))
      {
        dist_mat[i, directions[2]] <- dist_eucl[i, directions[2]]
        PI[i, directions[2]] <- abs(vel[2])/sqrt(sum(vel^2))
      }
      if(!is_east & mask_sq[1] & (directions[2] %in% valid))
      {
        dist_mat[i, directions[2]] <- dist_eucl[i, directions[2]]
        PI[i, directions[2]] <- abs(vel[1])/sqrt(sum(vel^2))
      }
      next
    }
    
    ab <- solve(D) %*% matrix(c(vel[1],vel[2]), ncol=1)
    
    if(ab[1]>0 & mask_sq[1] & (directions[1,2] %in% valid))
    {
      dist_mat[i, directions[1,2]] <- dist_eucl[i, directions[1,2]]
      PI[i, directions[1,2]] <- ab[1]/sqrt(sum(ab^2))
    }
    if(ab[2]>0 & mask_sq[2] & (directions[2,2] %in% valid))
    {
      dist_mat[i, directions[2,2]] <- dist_eucl[i, directions[2,2]]
      PI[i, directions[2,2]] <- ab[2]/sqrt(sum(ab^2))
    }
  }
  
  sources <- which(colSums(!is.na(dist_mat)) == 0)
  outlets <- which(rowSums(!is.na(dist_mat)) < 2)
  
  singletons <- which(colSums(!is.na(dist_mat)) == 0 & rowSums(!is.na(dist_mat)) == 0)
  
  valid <- setdiff(1:B, singletons)
  PI <- t(apply(PI, MARGIN =  1, normalize_column))
  return(list(dist_mat, PI, valid))
}

#### Distances Object ####

initialize <- function(dist, prob)
{
  B <- dim(dist)[1]
  distances <- vector("list", length = B)
  
  for (i in 1:B)
  {
    temp <- vector("list", length = B)
    distances[[i]] <- list(lengths = temp, update = NULL, Id = i, visited = temp)
  }
  
  for (i in 1:B)
  {
    for (j in 1:B)
    {
      if ((!is.na(dist[j, i])) & (!is.na(prob[j,i])))
      {
        distances[[i]]$lengths[[j]] <- 
          rbind(distances[[i]]$lengths[[j]], c(dist[j, i], prob[j, i]))
        distances[[i]]$update <- rbind(distances[[i]]$update, 
                                       cbind(j, nrow(distances[[i]]$lengths[[j]])))
        distances[[i]]$visited[[j]] <- list(c(j,i))
      }
    }
  }
  return(distances)
}

update <- function(p, dist, prob) 
{
  list_update <- NULL
  updated_flag <- FALSE
  
  if (!is.null(p$update) && nrow(p$update) > 0)
  {
    for (row in seq_len(nrow(p$update)))
    {
      i <- p$update[row, 1]
      j <- p$update[row, 2]
      inx <- which(!is.na(dist[, i]))
      
      if (length(inx) > 0)
      {
        for (k in seq_along(inx))
        {
          new_dist <- dist[inx[k], i] + p$lengths[[i]][j, 1]
          
          step <- prob[inx[k], i]
          if (is.na(step))
            new_PI <- 0
          else
            new_PI <- p$lengths[[i]][j, 2] * step
          
          if (new_PI > 1e-5 && !(inx[k] %in% p$visited[[i]][[j]]))
          {
            p$lengths[[inx[k]]] <- rbind(
              p$lengths[[inx[k]]],
              c(new_dist, new_PI)
            )
            
            l <- length(p$visited[[inx[k]]])
            if (is.null(l)) l <- 0
            p$visited[[inx[k]]][[l+1]] <- c(p$visited[[i]][[j]], inx[k])
            
            list_update <- rbind(list_update,
                                 cbind(inx[k],
                                       nrow(p$lengths[[inx[k]]])))
            
            updated_flag <- TRUE
          }
        }
      }
    }
    
    if (!is.null(list_update)) {
      list_update <- unique(list_update)
    }
    
    if (is.null(list_update) || nrow(list_update) == 0) {
      p$update <- NULL
    } else {
      p$update <- list_update
      updated_flag <- TRUE
    }
  }
  
  return(list(p = p, updated = updated_flag))
}

build_covariances <- function(sill, range, paths, paths2 = NULL, U, fun, mat)
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
        covs1[j,i] <- sums_1 * (U[j,i] * sqrt(U[i,i]*U[j,j]))
      }
      if (!is.null(paths[[j]]$lengths[[i]]))
      {
        mat <- paths[[j]]$lengths[[i]]
        l <- dim(mat)[1]
        for (k in 1:l)
        {
          sums_2 <- sums_2 + mat[k,2]*fun(c(sill, range), mat[k,1])
        }
        covs1[i,j] <- sums_2 * (U[i,j] * sqrt(U[i,i]*U[j,j]) )
      }
    }
  }
  covs1 <- covs1 * U
  D <- matrix(sqrt(1/diag(U)), nrow = B, ncol = B)
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

compute_probs <- function(G, i, j, PI_aux)
{
  G_AA <- G[c(i,j), c(i,j)]
  b <- solve(G_AA)%*%c(1,1)
  G_xa <- G[-c(i,j), c(i,j)]
  h_xA <- G_xa %*% b
  g_x <- 1-h_xA
  pr_ij <- sum(PI_aux[j, -c(i,j)] * g_x)
  pr_ji <- sum(PI_aux[i, -c(i,j)] * g_x)
  return(list(i = i, j = j, p_ij = pr_ij, p_ji = pr_ji))
}

evaluate_U <- function(PI)
{
  PI_aux <- PI
  PI_aux[which(is.na(PI_aux))] <- 0
  PI_aux <- cbind(PI_aux, 1 - rowSums(PI_aux))
  PI_aux <- rbind(PI_aux, c(rep(0, nrow(PI_aux)), 1))
  B <- nrow(PI_aux)
  G <- graph_from_adjacency_matrix(PI_aux > 0, mode = "directed")
  scc <- components(G, mode = "strong")
  recurrent <- c()
  for (comp_id in unique(scc$membership)) {
    states <- which(scc$membership == comp_id)
    subP <- PI_aux[states, , drop = FALSE]
    if (all(colSums(subP[, -states, drop = FALSE]) == 0)) {
      recurrent <- c(recurrent, states)
    }
  }
  recurrent <- setdiff(recurrent, B)
  if (length(recurrent))
    warning("Some recurrent classes are absorbing.")
  
  B <- B - 1
  PI_aux <- PI_aux[1:B, 1:B]
  nrc <- setdiff(1:B, recurrent)
  if(length(recurrent) > 0)
  {
    PI_aux <- PI_aux[nrc, nrc]
  }
  B <- nrow(PI_aux)
  G <- solve(diag(B) - PI_aux)
  not_self <- 1/diag(G)
  
  pair_inx <- which(upper.tri(matrix(0, B, B)), arr.ind = TRUE)
  n_cores <- as.numeric(Sys.getenv("PBS_NP"))
  if (is.na(n_cores) || n_cores < 1) {
    n_cores <- as.numeric(Sys.getenv("NCPUS"))
  }
  if (is.na(n_cores) || n_cores < 1) {
    n_cores <- 1
  }
  cat("Uso", n_cores, "core\n")
  
  cl <- makeCluster(n_cores)
  
  clusterExport(cl, c("pair_inx", "G", "PI_aux", "compute_probs"), envir = environment())
  
  res_list <- parLapply(cl, seq_len(nrow(pair_inx)), function(k) {
    compute_probs(G, pair_inx[k, 1], pair_inx[k, 2], PI_aux)
  })
  
  stopCluster(cl)
  G_aux <- matrix(0, nrow = B, ncol = B)
  for (r in res_list) {
    if (!is.null(r)) {
      i   <- r[["i"]]
      j  <- r[["j"]]
      pij <- r[["p_ij"]]
      pji <- r[["p_ji"]]
      
      if (!is.na(pij) && !is.na(pji)) {
        G_aux[i, j] <- pij
        G_aux[j, i] <- pji
      }
    }
  }
  
  if (length(recurrent) > 0)
  {
    B <- nrow(PI)
  }
  U <- diag(B)
  U[nrc, nrc] <- G_aux
  diag(U)[nrc] <- not_self
  
  return(U)
}

#### Covariance estimation ####

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
  
  fuv <- mean(vars[!connected])
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

calc_covariance <- function(coords, values, bins = 15)
{
  
  dists <- as.matrix(dist(coords))
  
  # Center values
  z <- values - mean(values)
  cross_prods <- outer(z, z, "*")
  
  max_dist <- max(dists)
  bin_breaks <- seq(0, max_dist, length.out = bins + 1)
  cov_data <- data.frame(dist = numeric(bins),
                         gamma  = numeric(bins),
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

evaluate_covariance_penalization_train_correctsill <- function(ma, paths, train,
                          l, U, cutoff = NULL, sill)
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

  vars <- vars[connected==1]
  omega <- omega[connected==1,]

  cols <- which(colSums(omega) == 0)

  if (length(cols))
  {
    omega <- omega[,-cols]
    aux <- t(omega)%*%omega
    q <- t(omega)%*%(sill - vars)
    nor <- max(abs(q))/sill
    mins <- min(2*diag(abs(aux)) - rowSums(abs(aux)))
    lambda <- nor - mins
    covariance <- c(solve(aux +
                                 lambda*diag(1,l-length(cols)))%*%q,
                         rep(0, length(cols)))
  }
  else
  {
    aux <- t(omega)%*%omega
    q <- t(omega)%*%(sill - vars)
    nor <- max(abs(q))/sill
    mins <- min(2*diag(abs(aux)) - rowSums(abs(aux)))
    lambda <- nor - mins
    covariance <- c(solve(aux +
                                 lambda*diag(1,l-length(cols)))%*%q,
                         rep(0, length(cols)))
  }

  covariance <- c(sill, covariance)

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
  return(covariance)

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

calc_variogram <- function(coords, values, bins = 15)
{
  dists <- as.matrix(dist(coords))

  # Center values
  diffs <- outer(values, values, "-")
  sq_diffs <- 0.5 * (diffs^2)

  max_dist <- max(dists)
  bin_breaks <- seq(0, max_dist, length.out = bins + 1)
  variogram_data <- data.frame(dist = numeric(bins),
                              gamma = numeric(bins),
                              np   = numeric(bins))

  for (i in 1:bins) {
    indices <- which(dists >= bin_breaks[i] & dists < bin_breaks[i+1], arr.ind = TRUE)
    if (nrow(indices) > 0) {
      valid <- indices[indices[,1] < indices[,2], , drop = FALSE]
      if (nrow(valid) > 0) {
        variogram_data$dist[i] <- (bin_breaks[i] + bin_breaks[i+1]) / 2
        variogram_data$gamma[i]  <- mean(sq_diffs[valid])
        variogram_data$np[i]   <- nrow(valid)
      }
    }
  }
  return(variogram_data)
}