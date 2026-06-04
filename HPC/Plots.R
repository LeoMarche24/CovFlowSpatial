library(sp)
library(gstat)
library(ncdf4)
library(ggplot2)
library(spatstat)
library(igraph)
library(sf)
library(progress)
library(ggmap)
library(ggforce)

source("functions_network.R")
source("functions_observations.R")

col1 <- "#3182BD"
col2 <- "#E6550D"
third <- rgb(63/255, 128/255, 97/255)
fourth <- "#333399"
col5 <- "#C51B8A"
grey <- "grey80"

domain <- "Sardinia"
load(paste0("Data/Data_", domain, ".RData"))
df <- df[which(df$years == 2022) ,]
register_stadiamaps("424ec34b-4e14-4acc-8aa9-9e8ab4f4fbc9", write = FALSE)
load(paste0("Data/bbox_", domain, ".Rdata"))
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
plot_size <- 20
p <- plot_lin_net(lines, df[, 1:2], plot_size, map, inx, sources, outlets)
p
lat_max <- max(df$lat)
lat_min <- min(df$lat)
lon_max <- max(df$lon)
lon_min <- min(df$lon)

lat_range <- lat_max - lat_min
lon_range <- abs(lon_max - lon_min)

lat_mean_rad <- ((lat_max + lat_min) / 2) * (pi / 180)

aspect_ratio <- lat_range / (lon_range * cos(lat_mean_rad))
output <- paste0("Plots/Exploratory/")

ggsave(
  filename = paste0(output, "Linear_Network_Simulation.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size * aspect_ratio,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

#### Plot the hot spot and contour ####

rcp <- 'rcp45'
year <- 2050
plot_size <- 25
register_stadiamaps("424ec34b-4e14-4acc-8aa9-9e8ab4f4fbc9", write = FALSE)
fourth <- "#333399"

load(paste0("HPC/Export/output_data/", domain, "/data_projection_", domain, "_2050.RData"))
output <- paste0("HPC/Export/output_data/", domain, "/Plots/Projections/")

lat_max <- max(df$lat)
lat_min <- min(df$lat)
lon_max <- max(df$lon)
lon_min <- min(df$lon)

bbox <- c(left = lon_min, bottom = lat_min, right = lon_max, top = lat_max)
map <- get_stadiamap(bbox, zoom = 8, maptype = "alidade_smooth")

lat_range <- lat_max - lat_min
lon_range <- abs(lon_max - lon_min)

lat_mean_rad <- ((lat_max + lat_min) / 2) * (pi / 180)

aspect_ratio <- lat_range / (lon_range * cos(lat_mean_rad))

points <- which(df$covariances %in% (sort(df$covariances, decreasing = T)[1:2]))
inx <- which(!is.na(df$covariances))

p <- ggmap(map, darken = c(.356,"white")) +
  geom_tile(data = df[inx ,], 
            mapping = aes(x = lon, y = lat, fill=covariances)) +
  geom_segment(data = df[inx ,], 
               aes(x = lon, y = lat, xend = lon + east, yend = lat + north),
               arrow = arrow(length = unit(0.3, "cm")), color = fourth, alpha=0.4, 
               linewidth = 2) +
  geom_point(data = df[points, ], aes(x = lon, y = lat), color = "black", 
             size = plot_size*0.7) +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient(low = grey, high = third, name = "Temperature") + 
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = plot_size*5),
    axis.title = element_text(size = plot_size*5)
  )
p
ggsave(
  filename = paste0(output, "Contour_New.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size * aspect_ratio,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

p <- ggmap(map, darken = c(.356,"white")) +
  geom_tile(data = df[inx ,], 
            mapping = aes(x = lon, y = lat, fill=covariances_euclidean1)) +
  geom_segment(data = df[inx ,], 
               aes(x = lon, y = lat, xend = lon + east, yend = lat + north),
               arrow = arrow(length = unit(0.3, "cm")), color = fourth, alpha = 0.4,
               linewidth = 2) +
  geom_point(data = df[points[2], ], aes(x = lon, y = lat), color = "black", size = plot_size*0.5) +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient(low = grey, high = third, name = "Temperature") + 
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = plot_size*5),
    axis.title = element_text(size = plot_size*5)
  )
p
ggsave(
  filename = paste0(output, "Contour_Eucl1.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size * aspect_ratio,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

p <- ggmap(map, darken = c(.356,"white")) +
  geom_tile(data = df[inx ,], 
            mapping = aes(x = lon, y = lat, fill=covariances_euclidean2)) +
  geom_segment(data = df[inx,], 
               aes(x = lon, y = lat, xend = lon + east, yend = lat + north),
               arrow = arrow(length = unit(0.3, "cm")), color = fourth, alpha=0.4, 
               linewidth = 2) +
  geom_point(data = df[points[1], ], aes(x = lon, y = lat), color = "black", size = plot_size*0.5) +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient(low = grey, high = third, name = "Temperature") + 
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = plot_size*5),
    axis.title = element_text(size = plot_size*5)
  )
p
ggsave(
  filename = paste0(output, "Contour_Eucl2.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size * aspect_ratio,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

threshold <- 25:30

for(index in threshold)
{
  col_name <- paste0("exceedance_", index)
  df_temp <- df[,1:4]
  df_temp$hot_plot <- df[[col_name]]
  df_temp$hot_plot <- factor(df_temp$hot_plot, levels = c(1, 2, 3))
  
  p <- ggmap(map, darken = c(.356, "white")) +
    geom_tile(data = df_temp[inx,],
              mapping = aes(x = lon, y = lat, fill=hot_plot), alpha = 0.7) +
    geom_segment(data = df_temp[inx,],
                 aes(x = lon, y = lat, xend = lon + east, yend = lat + north),
                 arrow = arrow(length = unit(0.1, "cm")), color = fourth, alpha=0.8,
                 linewidth = 2) +
    scale_fill_manual(
      values = c(
        "1" = grey,    # color for hot_plot == 1
        "2" = third,   # color for hot_plot == 2
        "3" = col2
      )) +
    labs(
      x = "Longitude",          # Change X-axis label
      y = "Latitude",           # Change Y-axis label
      fill = NULL
    )+
    theme_minimal() +
    theme(
      legend.position = 'none',
      axis.title.x = element_text(size = plot_size*5), # Change X-axis label font size
      axis.title.y = element_text(size = plot_size*5), # Change Y-axis label font size
      strip.text = element_text(size = plot_size*5),
      axis.text.x = element_text(size = plot_size*5),
      axis.text.y = element_text(size = plot_size*5)
    )
  
  ggsave(
    filename = paste0(output, "Hot_spots_VF_", index, ".pdf"),
    plot     = p,
    width    = plot_size,
    height   = plot_size * aspect_ratio,
    units    = "in",
    dpi      = 92,
    limitsize = FALSE
  )
}

for(index in threshold)
{
  col_name <- paste0("exceedance_", index, "_Eucl")
  df_temp <- df[,1:4]
  df_temp$hot_plot <- df[[col_name]]
  df_temp$hot_plot <- factor(df_temp$hot_plot, levels = c(1, 2, 3))
  
  p <- ggmap(map, darken = c(.356, "white")) +
    geom_tile(data = df_temp[inx,],
              mapping = aes(x = lon, y = lat, fill=hot_plot), alpha = 0.7) +
    geom_segment(data = df_temp[inx,],
                 aes(x = lon, y = lat, xend = lon + east, yend = lat + north),
                 arrow = arrow(length = unit(0.1, "cm")), color = fourth, alpha=0.8,
                 linewidth = 2) +
    scale_fill_manual(
      values = c(
        "1" = grey,    # color for hot_plot == 1
        "2" = third,   # color for hot_plot == 2
        "3" = col2
      )) +
    labs(
      x = "Longitude",          # Change X-axis label
      y = "Latitude",           # Change Y-axis label
      fill = NULL
    )+
    theme_minimal() +
    theme(
      legend.position = 'none',
      axis.title.x = element_text(size = plot_size*5), # Change X-axis label font size
      axis.title.y = element_text(size = plot_size*5), # Change Y-axis label font size
      strip.text = element_text(size = plot_size*5),
      axis.text.x = element_text(size = plot_size*5),
      axis.text.y = element_text(size = plot_size*5)
    )
  
  
  ggsave(
    filename = paste0(output, "Hot_spots_VF_", index, "_Eucl.pdf"),
    plot     = p,
    width    = plot_size,
    height   = plot_size * aspect_ratio,
    units    = "in",
    dpi      = 92,
    limitsize = FALSE
  )
}

