# MAPAS

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


# Mapas Base ----

# Extracción de polígonos de países para la región de interés

tierra_union <- ne_countries(scale = 10, country = c("Argentina", "Chile", "Uruguay"), returnclass = "sf") %>%
  st_union()  # unir todas las geometrías (incluye islas)

# Proyectar a un CRS métrico (ej. ARGENTINA)

utm_crs <- 5345  # UTM Sur
utm_crs <- 5347  # UTM Norte
utm <- 4326
tierra <- st_transform(tierra_union, utm_crs)
# plot(tierra)

# Definir bounding box de interés (por ejemplo, costa central)

# Región de Buenos Aires
bbox_region <- st_bbox(c(xmin = 5000000, xmax = 6000000, ymin = 5200000, ymax = 6200000), crs = st_crs(utm_crs))
# Mapa de Patagonia
# bbox_region <- st_bbox(c(xmin = 3100000, xmax = 4000000, ymin = 3800000, ymax = 5800000), crs = st_crs(utm_crs))
# Mapa de Chubut
bbox_region <- st_bbox(c(xmin = 3300000, xmax = 4100000, ymin = 4700000, ymax = 5600000), crs = st_crs(utm_crs))
# Mapa del Golfo Nuevo
# bbox_region <- st_bbox(c(xmin = 3500000, xmax = 3750000, ymin = 5200000, ymax = 5400000), crs = st_crs(utm_crs))
# Mapa de Santa Cruz y TdF
bbox_region <- st_bbox(c(xmin = 3200000, xmax = 4000000, ymin = 3800000, ymax = 4900000), crs = st_crs(utm_crs))

# Recortar tierra al bounding box
tierra_region <- st_crop(tierra, bbox_region)
plot(tierra_region)

# Crear polígono del mar como el bounding box proyectado menos la tierra
bbox_proj <- st_as_sfc(st_bbox(bbox_region))  # bounding box del área
mar_proj <- st_difference(bbox_proj, st_union(tierra_region))

## Centroides ----
resol <- 150000 # Para extraccion de datos climáticos | Resolución en metros
r_template_150000 <- raster(extent(bbox_region), resolution = resol, crs = utm_crs)

# Obtener todos los centroides de las celdas (como matriz de coordenadas)
centroides <- xyFromCell(r_template_150000, 1:ncell(r_template_150000)) %>%
  as.data.frame() %>%
  st_as_sf(coords = c("x", "y"), crs = utm_crs)

# Seleccionar solo los puntos que están dentro del mar (no sobre tierra)
puntos_mar <- st_intersection(centroides, mar_proj)

plot(st_geometry(mar_proj),
     main = "Puntos para la extracción de datos meteorológicos",
     col = "lightblue")
points(centroides$geometry, col = "blue")
points(puntos_mar$geometry, col = "red")
length(puntos_mar$geometry)

# Crear raster vacío con la extensión del bounding box proyectado
resol_fetch <- 5000
r_template <- raster(extent(bbox_region), resolution = resol_fetch, crs = utm_crs)

# Raster de la distancia a la costa
coast_vec <- vect(tierra_region)
dist_raster <- distance(project(rast(r_template), paste0("EPSG:", utm_crs)), coast_vec)
plot(dist_raster)

# Mapa para extraer sampled points
max_val <- global(dist_raster, "max", na.rm = TRUE) # Dist raster en archivo viento.R
r_weights <- 1-(dist_raster/max_val$max)
r_weights[r_weights == 1] <- NA

plot(r_weights)
