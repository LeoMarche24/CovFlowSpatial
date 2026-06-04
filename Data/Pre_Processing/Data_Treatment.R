library(ggplot2)
library(sp)

col1 <- "#3182BD"
col2 <- "#E6550D"
third <- rgb(63/255, 128/255, 97/255)
fourth <- "#333399"
col5 <- rgb(204/255, 0/255, 0/255)
grey <- "grey80"

#### Sardinia domain ####

df <- read.csv('Data/Pre_Processing/RealData_August.csv')
df_temp <- df[which(df$years == 2022) ,]
domain <- "Sardinia"

ggplot(df_temp[which(!is.na(df_temp$east) & !is.na(df_temp$north) & !is.na(df_temp$temperature)) ,],
       aes(x=lon, y=lat, fill=temperature)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

df <- df[which(df$lat>38 & df$lon<12 & df$lat<44 & df$lon>7) ,]
bbox <- c(left = 7, bottom = 38, right = 12, top = 44)
df_temp <- df[which(df$years == 2022) ,]

ggplot(df_temp[which(!is.na(df_temp$east) & !is.na(df_temp$north) & !is.na(df_temp$temperature)) ,],
       aes(x=lon, y=lat, fill=temperature)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

data <- read.csv("Data/Pre_Processing/Projections_August.csv")
data <- data[which(data$lat>38 & data$lon<12 & data$lat<44 & data$lon>7) ,]

## Load in the data for high performance computing platform ##

coords <- df[,1:2]
coordinates(coords) <- c("lon", "lat")
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")
coords <- spTransform(coords, CRS("EPSG:3035"))
df$longitude <- coords@coords[,1]
df$latitude <- coords@coords[,2]
coords <- data[,1:2]
coordinates(coords) <- c("lon", "lat")
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")
coords <- spTransform(coords, CRS("EPSG:3035"))
data$longitude <- coords@coords[,1]
data$latitude <- coords@coords[,2]

save(data, file = paste0("HPC/Export/input_data/Projections_", domain, ".RData"))
save(df, file = paste("HPC/Export/input_data/Data_", domain, ".RData", sep = ""))
save(data, file = paste0("Data/Projections_", domain, ".RData"))
save(df, file = paste0("Data/Data_", domain, ".RData"))
save(bbox, file = paste0("Data/bbox_", domain, ".RData"))


#### Temp domain ####

# A smaller domain to test the code before running the full model

df <- read.csv('Data/Pre_Processing/RealData_August.csv')
df_temp <- df[which(df$years == 2022) ,]
domain <- "Temp"

ggplot(df_temp[which(!is.na(df_temp$east) & !is.na(df_temp$north) & !is.na(df_temp$temperature)) ,],
       aes(x=lon, y=lat, fill=temperature)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

df <- df[which(df$lat>42 & df$lon<9 & df$lat<44 & df$lon>7) ,]
bbox <- c(left = 7, bottom = 42, right = 9, top = 44)
df_temp <- df[which(df$years == 2022) ,]

ggplot(df_temp[which(!is.na(df_temp$east) & !is.na(df_temp$north) & !is.na(df_temp$temperature)) ,],
       aes(x=lon, y=lat, fill=temperature)) + 
  geom_raster() + 
  scale_fill_gradient(low = col2, high = col1)

data <- read.csv("Data/Pre_Processing/Projections_August.csv")
data <- data[which(data$lat>42 & data$lon<9 & data$lat<44 & data$lon>7) ,]

## Load in the data for high performance computing platform ##

coords <- df[,1:2]
coordinates(coords) <- c("lon", "lat")
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")
coords <- spTransform(coords, CRS("EPSG:3035"))
df$longitude <- coords@coords[,1]
df$latitude <- coords@coords[,2]
coords <- data[,1:2]
coordinates(coords) <- c("lon", "lat")
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")
coords <- spTransform(coords, CRS("EPSG:3035"))
data$longitude <- coords@coords[,1]
data$latitude <- coords@coords[,2]

save(data, file = paste0("HPC/Export/input_data/Projections_", domain, ".RData"))
save(df, file = paste("HPC/Export/input_data/Data_", domain, ".RData", sep = ""))
save(data, file = paste0("Data/Projections_", domain, ".RData"))
save(df, file = paste0("Data/Data_", domain, ".RData"))
save(bbox, file = paste0("Data/bbox_", domain, ".RData"))
