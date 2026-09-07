# LIBRERÍAS ----

library(lubridate)
library(tidyverse)
library(dplyr)
library(sf)
library(gstat)
library(terra)
library(rnaturalearth)
library(raster)
library(gfwr)
library(ropenmeteo)
library(knitr)

source("funciones.R")

## Altura de olas ----
resol_olas <- 5000

r_template <- raster(extent(bbox_region), resolution = resol_olas, crs = CRS(paste0("EPSG:", utm_crs)))

analizar_periodo <- function(x, olas_umbral = 1.5) { 
  
  inicio <- as.character(periodos_ventosos[x,]$fecha_inicio)
  fin <- as.character(periodos_ventosos[x,]$fecha_fin)
  
  print(paste0("Analizando periodo: ", inicio, " a ", fin, " x=", x))
  
  r_viento <- interpolar_viento_periodo(
    df = viento_periodo,
    inicio = inicio,
    fin = fin,
    capa_raster = r_template,
    utm = utm_crs)
  
  mask_viento <- rasterize(vect(mar_proj), rast(r_template))
  
  r_viento_dir <- terra::mask(rast(r_viento$direction), mask = mask_viento)
  r_viento_vel <- terra::mask(rast(r_viento$speed), mask = mask_viento)
  r_viento_th <- r_viento_vel >= viento_umbral/2
  
  ruta_mapas <- paste0("data/chubut/mapas/")
  writeRaster(r_viento_vel, paste0(ruta_mapas, inicio, "_", fin, "_r_viento_vel.tiff"), overwrite = T) 
  writeRaster(r_viento_dir, paste0(ruta_mapas, inicio, "_", fin, "_r_viento_dir.tiff"), overwrite = T)
  writeRaster(r_viento_th, paste0(ruta_mapas, inicio, "_", fin, "_r_viento_th.tiff"), overwrite = T)
  
  sampled_points <- spatSample(r_weights, size = 1000, method = "weights",
                               na.rm = TRUE, as.points = TRUE) %>%
    st_as_sf() %>%
    st_transform(crs = utm_crs)
  
  # plot(r_weights)
  # points(sampled_points$geometry)
  
  print(paste0("Generando mapa de olas"))
  
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
  
  r_areas_refugio <- r_viento_vel>=viento_umbral/2 & r_olas<=olas_umbral
  
  writeRaster(r_olas, paste0(ruta_mapas, inicio, "_", fin, "_r_olas.tiff"), overwrite = T)
  writeRaster(r_areas_refugio, paste0(ruta_mapas, inicio, "_", fin, "_r_areas_refugio.tiff"), overwrite = T)
  
  if(any(values(r_areas_refugio) == 1, na.rm = TRUE)){
    # Polygonize: each cell becomes a polygon (or dissolve = TRUE groups equal values)
    
    # Convert to sf
    sf_poly <- st_as_sf(as.polygons(r_areas_refugio, dissolve = TRUE))
    sf_poly <- sf_poly[sf_poly$layer == 1, ] |>
      sf::st_as_sf() |>
      sf::st_transform(crs = 4326)
    
    print(paste0("Se detectaron áreas de refugio, extrayendo esfuerzo de pesca"))
    
    esfuerzoPesca <- gfw_ais_fishing_hours(spatial_resolution = "HIGH",
                                           temporal_resolution = "DAILY",
                                           start_date = inicio,
                                           end_date = fin,
                                           region = sf_poly,
                                           region_source = "USER_SHAPEFILE",
                                           print_request = TRUE)
    
    nombre <- paste0("data/chubut/esfuerzo/", inicio, "_", fin, "_esfuerzo.csv")
    write.csv2(esfuerzoPesca, nombre)
  } else {
    
    print(paste0("No hay áreas de refugio en el periodo: ", inicio, " a ", fin, " x=", x))
    
  }
  
}

# Loop para extraer todos los valores del período.

length(periodos_ventosos$fecha_inicio) # 407

for (i in 1:100) {
  
  try(analizar_periodo(i))
  
}


