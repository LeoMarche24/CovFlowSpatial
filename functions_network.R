##################################
##### CREATION OF THE NETWORK ####
##################################

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

# The following function works if the grid is regular - COPERNICUS data
compute_matrices <- function(df)
{
  nor <- df$north
  eas <- df$east
  nor[is.na(nor)] <- 0
  eas[is.na(eas)] <- 0
  valid <- which(!is.na(nor) & !is.na(eas) & !is.na(df$value) & ((abs(nor)>0) | (abs(eas)>0)))
  coord_sub <- data.frame(longitude = df$longitude, latitude = df$latitude)
  dist_eucl <- as.matrix(dist(as.matrix(coord_sub), method = 'euclidean'))
  coord_sub <- data.frame(longitude = df$lon, latitude = df$lat)
  B <- dim(df)[1]
  dist_mat <- matrix(NA, B,B)
  mat <- matrix(F, nrow=B, ncol=B)
  PI <- matrix(NA, B, B)
  lines <- list()
  
  lines <- list()
  for (i in 1:B) {
    line_id <- generate_line_id(coord_sub[i, ])
    lines[[i]] <- sp::Lines(slinelist = list(), ID = line_id)
  }
  
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
    temp <- list()
    line <- NULL
    line2 <- NULL
    
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
        line <- sp::Line(rbind(coord_sub[i,], coord_sub[directions[2] ,]))
      }
      if(!is_east & mask_sq[1] & (directions[2] %in% valid))
      {
        dist_mat[i, directions[2]] <- dist_eucl[i, directions[2]]
        PI[i, directions[2]] <- abs(vel[1])/sqrt(sum(vel^2))
        line2 <- sp::Line(rbind(coord_sub[i,], coord_sub[directions[2] ,]))
      }
      id <- paste(coord_sub[i,1], coord_sub[i,2])
      inx <- which(sapply(lines, function(line) attr(line, "ID") == id))
      temp <- list(line, line2)
      if (!is.null(temp[[1]]) || !is.null(temp[[2]]))
      {
        lines[[inx]] <- sp::Lines(temp[!sapply(temp, is.null)], ID = id)
      }
      
      next
    }
    
    ab <- solve(D) %*% matrix(c(vel[1],vel[2]), ncol=1)
    
    if(ab[1]>0 & mask_sq[1] & (directions[1,2] %in% valid))
    {
      dist_mat[i, directions[1,2]] <- dist_eucl[i, directions[1,2]]
      PI[i, directions[1,2]] <- ab[1]/sqrt(sum(ab^2))
      line <- sp::Line(rbind(coord_sub[i,], coord_sub[directions[1,2] ,]))
    }
    if(ab[2]>0 & mask_sq[2] & (directions[2,2] %in% valid))
    {
      dist_mat[i, directions[2,2]] <- dist_eucl[i, directions[2,2]]
      PI[i, directions[2,2]] <- ab[2]/sqrt(sum(ab^2))
      line2 <- sp::Line(rbind(coord_sub[i,], coord_sub[directions[2,2] ,]))
    }
    id <- paste(coord_sub[i,1], coord_sub[i,2])
    inx <- which(sapply(lines, function(line) attr(line, "ID") == id))
    temp <- list(line, line2)
    if (!is.null(temp[[1]]) || !is.null(temp[[2]]))
    {
      lines[[inx]] <- sp::Lines(temp[!sapply(temp, is.null)], ID = id)
    }
  }
  
  sources <- which(colSums(!is.na(dist_mat)) == 0)
  outlets <- which(rowSums(!is.na(dist_mat)) < 2)
  
  singletons <- which(colSums(!is.na(dist_mat)) == 0 & rowSums(!is.na(dist_mat)) == 0)
  
  valid <- setdiff(1:B, singletons)
  PI <- t(apply(PI, MARGIN =  1, normalize_column))
  
  return(list(lines, dist_mat, PI, valid))
}

#############################
#### PLOT OF THE NETWORK ####
#############################

plot_lin_net <- function(lines, coord_sub, plot_size, map, inx, sources, outlets)
{
  B <- dim(coord_sub)[1]
  for (i in B:1)
  {
    if (!length(lines[[i]]@Lines))
      lines[[i]] <- NULL
  }
  lin_net <- SpatialLines(LinesList = lines)
  lin_net_sf <- st_as_sf(lin_net, crs = 4326)
  lin_net_df <- st_coordinates(lin_net_sf) |> as.data.frame()
  names(lin_net_df)[1:2] <- c("lon", "lat")
  lin_net_df$group <- interaction(lin_net_df$L1, lin_net_df$L2, drop = TRUE)
  
  lines_coords <- NULL
  for (i in 1:length(lin_net@lines))
  {
    for (j in 1:length(lin_net@lines[[i]]@Lines))
    {
      new_cord <- as.vector(lin_net@lines[[i]]@Lines[[j]]@coords)
      if (length(new_cord) == 4)
      {
        lines_coords <- rbind(lines_coords, new_cord)
      }
      else
      {
        print(i,j)
        errorCondition("Problem")
      }
    }
  }
  lines_coords <- data.frame(lines_coords)
  names(lines_coords) <- c("x0", "x1", "y0", "y1")
  midpoints_x <- (lines_coords[, "x0"] + lines_coords[, "x1"]) / 2
  midpoints_y <- (lines_coords[, "y0"] + lines_coords[, "y1"]) / 2
  arrow_data <- data.frame(
    x0 = lines_coords$x0,
    y0 = lines_coords$y0,
    x1 = midpoints_x,
    y1 = midpoints_y
  )
  
  coord_points <- as.matrix(coord_sub)
  coord_points <- coord_points[inx, ]
  
  # Plot using ggplot2
  ggmap(map, darken = c(.356,"white")) +
    geom_path(data = lin_net_df,
              aes(x = lon, y = lat, group = group),
              color = grey, linewidth = 3) +
    geom_point(data=data.frame(lon=coord_sub[, 1], lat=coord_sub[, 2]), 
               aes(x=lon, y=lat), col=grey, size = 3, alpha = 0.4) +
    geom_point(data=data.frame(x=coord_points[, 1], y=coord_points[, 2]), 
               aes(x=x, y=y), col="black", size = 3) +
    # Punti sources
    geom_point(data=data.frame(x=coord_points[sources, 1], y=coord_points[sources, 2]), 
               aes(x=x, y=y), shape = 21, color = "black", fill = "gold", size = 13) +
    # Punti outlets in col1
    geom_point(data=data.frame(x=coord_points[outlets, 1], y=coord_points[outlets, 2]), 
               aes(x=x, y=y), shape = 21, color = "black", fill = col2, size = 13) +
    geom_segment(
      data = arrow_data,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      arrow = arrow(length = unit(0.15, "inches"), type = "closed"),
      color = fourth, linewidth = 1, alpha = 0.8
    ) +
    theme_minimal() +
    labs(x = "Longitude")  +
    labs(y = "Latitude")  +
    theme(
      axis.text = element_text(size = plot_size*3),
      axis.title = element_text(size = plot_size*3),
      strip.text = element_text(size = plot_size*3)
    )
  
}

###########################################
#### CREATION OF THE DISTANCES OBJECTS ####
###########################################

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
  if (length(nrow(p$update)))
  {
    for (row in 1:nrow(p$update))
    {
      i <- p$update[row, 1]
      j <- p$update[row, 2]
      inx <- which(!is.na(dist[, i]))
      if (length(inx))
      {
        for (k in 1:length(inx))
        {
          new_dist <- dist[inx[k],i]+p$lengths[[i]][j, 1]
          step <- prob[inx[k],i]
          if(is.na(step))
            new_PI <- 0
          else
            new_PI <- p$lengths[[i]][j, 2]*step
          if (new_PI > 1e-3 & !(inx[k] %in% p$visited[[i]][[j]]))
          {
            p$lengths[[inx[k]]] <- rbind(p$lengths[[inx[k]]], 
                                         c(new_dist, new_PI))
            l <- length(p$visited[[inx[k]]])
            if(is.null(l))
              l <- 0
            p$visited[[inx[k]]][[l+1]] <- c(p$visited[[i]][[j]], inx[k])
            list_update <- rbind(list_update, cbind(inx[k], nrow(p$lengths[[inx[k]]])))
          }
        }
        updated <<- TRUE
      }
    }
    p$update <- list_update
  }
  return(p)
}

compute_probs <- function(G, i, j, PI_aux)
{
  not_self <- 1/diag(G)
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
  G <- solve(diag(1, B) - PI_aux)
  not_self <- 1/diag(G)
  
  pair_inx <- which(upper.tri(matrix(0, B, B)), arr.ind = TRUE)
  res_list <- lapply(seq_len(nrow(pair_inx)), function(k) {
    compute_probs(G, pair_inx[k, 1], pair_inx[k, 2], PI_aux)
  })
  
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
