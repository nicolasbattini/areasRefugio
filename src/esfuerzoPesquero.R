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

# key <- gfw_auth()

esfuerzoPesca <- gfw_ais_fishing_hours(spatial_resolution = "HIGH",
                                       temporal_resolution = "DAILY",
                                       start_date = inicio,
                                       end_date = fin,
                                       region = sf_poly_1,
                                       region_source = "USER_SHAPEFILE",
                                       print_request = TRUE)
getwd()
nombre <- paste0("esfuerzo/bsas/", inicio, "_", fin)
write.csv2(esfuerzoPesca, "esfuerzo/")



## Calculo de esfuerzo total

esfuerzoPescaDF <- esfuerzoPesca %>%
  group_by(Lat, Lon, `Time Range`, `Vessel ID`) %>%
  summarise(esfuerzoTotal = sum(`Apparent Fishing Hours`))

length(esfuerzoPescaDF$`Vessel ID`)
length(unique(esfuerzoPescaDF$`Vessel ID`))

# Convertir a objeto sf (WGS84, EPSG:4326)
puntos_wgs84 <- st_as_sf(as.data.frame(esfuerzoPescaDF), coords = c("Lon", "Lat"), crs = 4326)

puntos_5345 <- st_transform(puntos_wgs84, crs = utm_crs)

# Rasterizar los puntos: asignar valor de "esfuerzo" a cada celda ----
# Usamos terra::rasterize. Si varios puntos caen en la misma celda,
# podemos promediarlos (fun = mean) o usar la primera ocurrencia.

r_esfuerzo <- rasterize(puntos_5345, r_template, field = "esfuerzoTotal", fun = sum, na.rm = TRUE)

# Ver el resultado
plot(r_esfuerzo)
plot(r_esfuerzo, r_viento$speed)

plot(r_esfuerzo,
     main = paste0("Esfuerzo (POSGAR Argentina 5) ", inicio, " | ", fin),
     col = hcl.colors(256, "viridis"))

tierra_proj <- st_transform(tierra, utm_crs)
plot(st_geometry(tierra_proj), add = TRUE, col = "gray", border = NA)


# Datos de los buques ----

# datosBuques <- get.vesselData("fd6dd9e42-2634-a6c2-f4a5-8efed6e523e6")

for(i in unique(esfuerzoPescaDF$`Vessel ID`)){
  
  buqueID <- i
  
  print(paste0("Analizando buque ", buqueID))
  
  if(buqueID %in% datosBuques$ID) {
    
    warning("Este buque ya esta en la base de datos")
    
    } else {
    
    datos <- get.vesselData(buqueID)
    
    datosBuques <- rbind(datosBuques, datos)
    
    Sys.sleep(sample(5:10, 1))
  }
}

# write.csv(datosBuques, "data/datosBuques.csv")


