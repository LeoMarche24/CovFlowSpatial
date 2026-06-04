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
register_stadiamaps("424ec34b-4e14-4acc-8aa9-9e8ab4f4fbc9", write = FALSE)
load(paste0("Data/bbox_", domain, ".Rdata"))
map <- get_stadiamap(bbox, zoom = 8, maptype = "alidade_smooth")

rcp <- 'rcp45'
plot_size <- 7

load(paste0("Data/data_projection_", domain, "_2050.RData"))

lat_max <- max(df$lat)
lat_min <- min(df$lat)
lon_max <- max(df$lon)
lon_min <- min(df$lon)

lat_range <- lat_max - lat_min
lon_range <- abs(lon_max - lon_min)

lat_mean_rad <- ((lat_max + lat_min) / 2) * (pi / 180)

aspect_ratio <- lat_range / (lon_range * cos(lat_mean_rad))

points <- which(df$covariances %in% (sort(df$covariances, decreasing = T)[1:2]))
inx <- which(!is.na(df$covariances))

p <- ggmap(map, darken = c(.356,"white")) +
  geom_tile(data = df[inx ,], 
            mapping = aes(x = lon, y = lat, fill=covariances)) +
  geom_segment(data = df[which(!is.na(df$value)) ,], 
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

index <- 25

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

p


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


p

