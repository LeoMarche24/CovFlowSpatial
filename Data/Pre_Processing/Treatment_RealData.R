library(ncdf4)

dataset <- nc_open('Data/Pre_Processing/data.nc')
latitude <- ncvar_get(dataset, "latitude")
longitude <- ncvar_get(dataset, "longitude")
temperature <- ncvar_get(dataset, "to")
time <- ncvar_get(dataset, "time")
uo_clean <- ncvar_get(dataset, "ugo")
vo_clean <- ncvar_get(dataset, "vgo")
nc_close(dataset)
aug <- seq(8,length(time), by = 12)
years <- 2006:2022
y <- 1

lon_sub <- longitude[which(longitude>5 & longitude<22)]
lat_sub <- latitude[which(latitude>35)]

coords <- expand.grid(lon_sub, lat_sub)
B <- dim(coords)[1]

lon_idx <- match(coords[, 1], longitude)
lat_idx <- match(coords[, 2], latitude)

df_aug <- data.frame()

for (y in 1:length(years))
{
  uo_df <- data.frame(lon = rep(longitude, each = length(latitude)),
                      lat = rep(latitude, times = length(longitude)),
                      value = c(t(uo_clean[,,aug[y]])))
  vo_df <- data.frame(lon = rep(longitude, each = length(latitude)),
                      lat = rep(latitude, times = length(longitude)),
                      value = c(t(vo_clean[,,aug[y]])))
  tem_df <- data.frame(lon = rep(longitude, each = length(latitude)),
                       lat = rep(latitude, times = length(longitude)),
                       value = c(t(temperature[,,aug[y]])))
  
  rows <- match(paste(longitude[lon_idx], latitude[lat_idx]), paste(uo_df$lon, uo_df$lat))
  
  temp <- tem_df$value[rows]
  east_values <- uo_df$value[rows]
  north_values <- vo_df$value[rows]
  df <- data.frame(
    lon = coords[, 1],
    lat = coords[, 2],
    east = east_values,
    north = north_values,
    temperature = temp,
    years = rep(years[y], dim(coords[,1]))
  )
  df_aug <- rbind(df_aug, df)
}

write.csv(df_aug, file = "Data/Pre_Processing/RealData_August.csv", row.names = F)
