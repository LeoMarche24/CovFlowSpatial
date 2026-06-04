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
output <- "Plots/Exploratory/"

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

ggsave(
  filename = paste0(output, "Linear_Network_Simulation.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size * aspect_ratio,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

#### Definition of the network ####

inx <- which(!is.na(df$east) & !is.na(df$north) & !is.na(df$temperature))
set.seed(2405)
point <- sample(inx, size = 1)
plot_size <- 5

vel <- c(df$north[point], df$east[point])
angle <- atan2(vel[1], vel[2])
lon0 <- df$lon[point]
lat0 <- df$lat[point]
nearest <- order(as.matrix(dist(df[,1:2]))[point,])[2:9]

p <- ggplot() +
  geom_point(data=data.frame(x=df[c(point, nearest), 1], y=df[c(point, nearest), 2], 
                             label = c("a", "c", "i", "e", "g", "b", "d", "h", "f")), aes(x=x, y=y), col='black') + 
  geom_text(data=data.frame(x=df[c(point, nearest), 1], y=df[c(point, nearest), 2], 
                            label = c("a", "c", "i", "e", "g", "b", "d", "h", "f"), 
                            hjust = c(1.2,1.2,-.2,1.2,1.2,-.2,1.2,-.2,1.2),
                            vjust = c(1,-.2,1,1,1,-.2,-.2,1,1)), 
            aes(x=x, y=y, label = label, hjust=hjust, vjust=vjust), size = 10) + 
  geom_segment(
    data = df[point ,],
    aes(x = lon, y = lat, xend = lon + east*2, yend = lat + north*2),
    arrow = arrow(length = unit(.3, "inches"), type = "closed"),
    color = col1, linewidth = 1) + 
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank()
  )
p
ggsave(
  filename = paste0(output, "8points.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

directions <- cbind(unlist(lapply(nearest, function(j)
{
  lon1 <- df$lon[j]
  lat1 <- df$lat[j]
  angle_j <- atan2(lat1 - lat0, lon1 - lon0)
})), nearest)
directions <- directions[!duplicated(directions[, 1]), ]
directions <- directions[order(abs(directions[,1] - angle))[1:2],]

D <- do.call(cbind, lapply(directions[,1], function(d) {
  c(sin(d), cos(d))
}))

ab <- solve(D) %*% matrix(c(vel[1],vel[2]), ncol=1)

coord_1 <- as.numeric(df[point,2:1] + ab[1] * D[,1] *2)
x1 <- coord_1[2]
y1 <- coord_1[1]

coord_2 <- as.numeric(df[point,2:1] + ab[2] * D[,2] *2)
x2 <- coord_2[2]
y2 <- coord_2[1]

arrow_data <- data.frame(
  x0 = rep(df$lon[point], 2),
  y0 = rep(df$lat[point], 2),
  x1 = c(x1, x2),
  y1 = c(y1,y2)
)
parallelogram_data <- data.frame(
  x0 = c(x1,x2),
  y0 = c(y1,y2),
  x1 = rep(with(df[point,], lon+east*2, 2)),
  y1 = rep(with(df[point,], lat+north*2, 2))
)
x1 <- df$lon[directions[1,2]]
y1 <- df$lat[directions[1,2]]
x2 <- df$lon[directions[2,2]]
y2 <- df$lat[directions[2,2]]
segment_data <- data.frame(
  x0 = rep(df$lon[point], 2),
  y0 = rep(df$lat[point], 2),
  x1 = c(x1, x2),
  y1 = c(y1,y2)
)

# Plot using ggplot2
p <- ggplot() +
  geom_point(data=data.frame(x=df[c(point, nearest), 1], y=df[c(point, nearest), 2], 
                             label = c("a", "c", "i", "e", "g", "b", "d", "h", "f")), aes(x=x, y=y), col='black') + 
  geom_text(data=data.frame(x=df[c(point, nearest), 1], y=df[c(point, nearest), 2], 
                            label = c("a", "c", "i", "e", "g", "b", "d", "h", "f"), 
                            hjust = c(1.2,1.2,-.2,1.2,1.2,-.2,1.2,-.2,1.2),
                            vjust = c(1,-.2,1,1,1,-.2,-.2,1,1)), 
            aes(x=x, y=y, label = label, hjust=hjust, vjust=vjust), size = 10) + 
  geom_segment(
    data = segment_data,
    aes(x = x0, y = y0, xend = x1, yend = y1),
    color = grey, linewidth = 2) + 
  geom_segment(
    data = arrow_data,
    aes(x = x0, y = y0, xend = x1, yend = y1),
    arrow = arrow(length = unit(0.2, "inches"), type = "closed"),
    color = fourth, linewidth = 1) +
  # Linee tratteggiate con parallelogramma
  geom_segment(
    data = parallelogram_data,
    aes(x = x0, y = y0, xend = x1, yend = y1),
    linetype = "dashed", color = grey, linewidth = 2) +
  geom_segment(
    data = df[point ,],
    aes(x = lon, y = lat, xend = lon + east*2, yend = lat + north*2),
    arrow = arrow(length = unit(.3, "inches"), type = "closed"),
    color = col1, linewidth = 1) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank()
  )
p
ggsave(
  filename = paste0(output, "Edges.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

#### Plot the data ####

load(paste0("Data/Projections_", domain, ".RData"))
df <- data
rm(data)
df <- df[which(df$year == 2050 & df$RCP == "rcp45") ,]

lat_max <- max(df$lat)
lat_min <- min(df$lat)
lon_max <- max(df$lon)
lon_min <- min(df$lon)

lat_range <- lat_max - lat_min
lon_range <- abs(lon_max - lon_min)

lat_mean_rad <- ((lat_max + lat_min) / 2) * (pi / 180)

aspect_ratio <- lat_range / (lon_range * cos(lat_mean_rad))

plot_size <- 25
axis_title_pt <- plot_size * 2
axis_text_pt  <- plot_size * 2

p <- ggmap(map, darken = c(.356,"white")) +
  geom_tile(data = df[which(!is.na(df$east) & !is.na(df$north) & !is.na(df$value)) ,], 
            mapping = aes(x = lon, y = lat, fill=value), alpha = 0.8) +
  geom_segment(data = df[which(!is.na(df$value)) ,], 
               aes(x = lon, y = lat, xend = lon + east, yend = lat + north),
               arrow = arrow(length = unit(0.3, "cm")), color = fourth, alpha=0.8, 
               linewidth = 2) +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient(low = col1, high = col2, name = "Temperature") + 
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = plot_size*3),
    axis.title = element_text(size = plot_size*3)
  )
p
ggsave(
  filename = paste0(output, "Whole_Picture.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size * aspect_ratio,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

results <- compute_matrices(df = df)
lines <- results[[1]]
dist <- results[[2]]
PI <- results[[3]]
inx <- results[[4]]

B <- length(inx)
PI <- PI[inx, inx]

g <- graph_from_adjacency_matrix(ifelse(is.na(PI), 0, 1), mode = "directed")
has_cycles <- !is_acyclic(g)
print("The built network has a cycle: ")
print(has_cycles)

sources <- which(colSums(PI, na.rm = T) == 0)
outlets <- which(rowSums(PI, na.rm = T) < 0.99)

p <- plot_lin_net(lines, df[, 1:2], plot_size, map, inx, sources, outlets)
p

ggsave(
  filename = paste0(output, "Linear_Network.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size * aspect_ratio,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)

#### Plot the joint exceedance region

coords <- c(10.2, 42.8)

radii_km <- c(10, 15, 20, 35, 50)
radii_deg <- radii_km / 111  # Approximation: 1 degree ~ 111 km

circles <- data.frame(
  lon = coords[1],
  lat = coords[2],
  r = radii_deg,
  radius = radii_km
)

center <- data.frame(
  lon = coords[1],
  lat = coords[2]
)

plot_size <- 7

p <- ggmap(map, darken = c(.356,"white")) +
  geom_circle(
    data = circles,
    aes(x0 = lon, y0 = lat, r = r, alpha = radius),
    fill = third,
    color = NA
  ) +
  geom_point(
    data = center,
    aes(x = lon, y = lat),
    color = "black", size = 2
  ) +
  scale_alpha_continuous(
    name = "Radius (km)",
    range = c(1.2, 0.2),   # stesso range che hai usato
    guide = "legend"
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = plot_size*4),
    axis.title = element_text(size = plot_size*4)
  )
p

ggsave(
  filename = paste0(output, "Joint_Exceedance_Region.pdf"),
  plot     = p,
  width    = plot_size,
  height   = plot_size * aspect_ratio,
  units    = "in",
  dpi      = 92,
  limitsize = FALSE
)
