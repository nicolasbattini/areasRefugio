# VIENTO ----

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

# Definir el período de estudio
inicio <- "2015-01-01"
fin <- "2025-12-31"
# inicio <- "2026-05-06"
# fin <- "2026-05-10"

# Transformar los puntos para poder extraer los datos desde la API de Open-Weather
puntos_mar_dec <- st_transform(puntos_mar, 4326)
coordenadas_viento <- st_coordinates(puntos_mar_dec)
plot(coordenadas_viento)

# Generar un data.frame para almacenar los datos de viento.
# datos_viento <- data.frame()
  
for (i in 1:nrow(coordenadas_viento)){
  
  lat <- coordenadas_viento[i,2]
  lon <- coordenadas_viento[i,1]
  
  datos_punto <- get_historical_weather(
    latitude = lat,
    longitude = lon,
    site_id = NULL,
    inicio,
    fin,
    variables = c("wind_speed_10m", "wind_direction_10m")
  )
  
  datos_punto <- cbind(lat, lon, datos_punto)
  datos_viento <- rbind(datos_viento, datos_punto)
}

# write.csv2(datos_viento, "C:/Users/Usuario/OneDrive/Documentos/Proyectos/Areas Refugio Buques/areasRefugioGFW/data/chubut/datos_viento_chubut.csv")
# write.csv2(datos_viento, "C:/Users/Usuario/OneDrive/Documentos/Proyectos/Areas Refugio Buques/areasRefugioGFW/data/bsas/datos_viento_bsas.csv", row.names = F)

# Datos viento ----

datos_viento <- read.csv2("C:/Users/Usuario/OneDrive/Documentos/Proyectos/Areas Refugio Buques/areasRefugioGFW/data/chubut/datos_viento_chubut.csv", header = T)
datos_viento <- read.csv("C:/Users/Usuario/OneDrive/Documentos/Proyectos/Areas Refugio Buques/areasRefugioGFW/data/bsas/datos_viento_bsas.csv", header = T)

# Temporales ----

# Filtrar los períodos con intensidad de viento mayor a 40 nudos (valores de tabla expresados en m/s)
viento_umbral <- 30*0.514444 # 30 nudos
# viento_umbral <- 40*1000/(60*60) # 40 km / h (1000 m/km, 1/60*60 h/seg)

datos_filtrados <- datos_viento %>%
  
  pivot_wider(id_cols = c(datetime, site_id, lat, lon),
              names_from = variable,
              values_from = prediction) %>%
  filter(wind_speed_10m >= viento_umbral)


chequeo <- datos_viento %>%
  pivot_wider(id_cols = c(datetime, site_id, lat, lon),
              names_from = variable,
              values_from = prediction) %>%
  mutate(fecha = as.Date(!!sym("datetime"))) %>%
  group_by(lat, lon) %>%
  summarise(N = n())

points(chequeo$lon, chequeo$lat, col = "red")

kable(head(datos_filtrados), caption = "Valores de viento mayor a 30 nudos para el área de estudio")

ss_datos_filtrados <- datos_viento %>%
  pivot_wider(id_cols = c(datetime, site_id, lat, lon),
              names_from = variable,
              values_from = prediction) %>%
  filter(datetime %in% datos_filtrados$datetime)

dias_ventosos <- datos_filtrados %>%
  dplyr::select(datetime) %>%
  mutate(fecha = as.Date(!!sym("datetime"))) %>%
  distinct(fecha) %>%
  arrange()

periodos_ventosos <- get_consecutive_periods(dias_ventosos, fecha = "fecha")
# write.csv2(periodos_ventosos, "data/bsas/periodos_ventosos.csv")
# write.csv2(periodos_ventosos, "data/chubut/periodos_ventosos.csv")

kable(head(periodos_ventosos), caption = "Períodos identificados con vientos fuertes")

ggplot(periodos_ventosos) +
  geom_bar(aes(as.numeric(duracion))) +
  xlab("Duracion (dias)") +
  ylab("Cantidad de eventos") +
  theme_light()

ggplot(periodos_ventosos) +
  geom_point(aes(fecha_inicio, duracion)) +
  xlab("Fecha") +
  ylab("Duracion (dias)") +
  theme_light()
