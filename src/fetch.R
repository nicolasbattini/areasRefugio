# LIBRERÍAS ----
library(dplyr)
library(sf)
library(terra)
library(raster)

source("funciones.R")


# Remuestreo del área de estudio ----



# LOOP

sampled_points <- spatSample(r_weights, size = 1000, method = "weights",
                             na.rm = TRUE, as.points = TRUE) %>%
  st_as_sf() %>%
  st_transform(crs = utm_crs)

plot(r_weights)
points(sampled_points$geometry)

dir_viento <- terra::extract(r_viento_dir, sampled_points)
vel_viento <- terra::extract(r_viento_vel, sampled_points)

olas <- vector(mode = "numeric", length = length(vel_viento$ID))

for (i in 1:length(vel_viento$ID)) {
  fetchi <- fetch_en_mar(sampled_points$geometry[i], dir_viento[i,2], tierra_region)
  Hsi <- Hs(fetchi, vel_viento[i,2])
  olas[i] <- Hsi
}

sampled_points_df <- as.data.frame(cbind(st_coordinates(sampled_points), olas))

sampled_points_sp <- st_as_sf(sampled_points_df, coords = c("X", "Y"), crs = utm_crs)

#modelo <- fit.variogram(variogram(olas ~ 1, sampled_points_sp), vgm(2, c("Exp", "Sph"), 500000, 1))
modelo <- gstat(formula = olas ~ 1, data = sampled_points_sp, 
                model = fit.variogram(variogram(olas ~ 1, sampled_points_sp), vgm(2, c("Exp", "Sph"), 500000, 1)))

r_olas <- interpolate(r_template, modelo, debug.level = 0) %>%
  rast()%>%
  terra::mask(mask_viento)

plot(r_viento_vel)
plot(r_olas)

olas_umbral <- 1.5
r_index <- r_viento_vel>=viento_umbral/2 & r_olas<=olas_umbral
plot(r_index)

# Polygonize: each cell becomes a polygon (or dissolve = TRUE groups equal values)
spoly <- as.polygons(r_index, dissolve = TRUE)   # returns SpatialPolygonsDataFrame
plot(spoly)

# Convert to sf
sf_poly <- st_as_sf(as.polygons(r_index, dissolve = TRUE))
plot(sf_poly)
sf_poly_1 <- sf_poly[sf_poly$layer == 1, ]
plot(sf_poly_1)

sf_poly_1 <- sf_poly_1 |>
  #sf::st_as_sfc() |>
  sf::st_as_sf() |>
  sf::st_transform(crs = 4326)




df_olas <- as.data.frame(na.omit(r_olas))
df_viento_vel <- as.data.frame(na.omit(r_viento_vel))

plot(df_olas$var1.pred, df_viento_vel$layer)

index <- df_viento_vel/df_olas

plot(index)

min(index)
max(index)
