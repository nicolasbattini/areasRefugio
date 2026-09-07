# LIBRERÍAS ----
library(dplyr)
library(sf)
library(terra)
library(raster)

source("funciones.R")


# Remuestreo del área de estudio ----
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

modelo <- gstat(formula = olas ~ 1, data = sampled_points_sp, 
                model = fit.variogram(variogram(olas ~ 1, sampled_points_sp), vgm(2, c("Exp", "Sph"), 500000, 1)))

r_olas <- interpolate(r_template, modelo, debug.level = 0) %>%
  rast()%>%
  terra::mask(mask_viento)

plot(r_viento_vel)
plot(r_olas)

# Asignar un valor umbral a la altura de olas 
olas_umbral <- 1.5
r_areas_refugio <- r_viento_vel>=viento_umbral/2 & r_olas<=olas_umbral
plot(r_areas_refugio)

# Poligonizar el área de refugio para extraer datos de pesca.
spoly <- as.polygons(r_areas_refugio, dissolve = TRUE)   # devuelve un SpatialPolygonsDataFrame
plot(spoly)

# Convert to sf
sf_poly <- st_as_sf(as.polygons(r_index, dissolve = TRUE))
sf_poly <- sf_poly[sf_poly$layer == 1, ]
plot(sf_poly)

sf_poly <- sf_poly |>
  #sf::st_as_sfc() |>
  sf::st_as_sf() |>
  sf::st_transform(crs = 4326)

plot(index)

min(index)
max(index)
