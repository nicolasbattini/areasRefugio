# FUNCIONES ----
promedio_viento_diario <- function(df, 
                                   col_sitio = "site_id",
                                   col_fecha = "datetime", 
                                   col_vel = "velocidad", 
                                   col_dir = "direccion") {
  
  # Convertir a radianes y calcular componentes u, v
  # Nota: Convención meteorológica estándar (0° = Norte, ángulo aumenta hacia el Este)
  # Si tu convención es diferente (0° = Este, antihorario), ajusta los senos/cosenos.
  df_proc <- df %>%
    group_by(col_sitio) %>%
    mutate(
      rad = !!sym(col_dir) * pi / 180,
      # Componentes vectoriales (viento hacia donde sopla, no desde)
      u = -!!sym(col_vel) * sin(rad),   # componente Este-Oeste
      v = -!!sym(col_vel) * cos(rad)    # componente Norte-Sur
    ) %>%
    # Excluir filas con velocidad = 0 (dirección indefinida)
    filter(!!sym(col_vel) > 0) %>%
    # Extraer fecha (sin hora)
    mutate(fecha = as.Date(!!sym(col_fecha)))
  
  # Promedio diario de componentes
  diario <- df_proc %>%
    group_by(fecha) %>%
    mutate(
      u_prom = mean(u, na.rm = TRUE),
      v_prom = mean(v, na.rm = TRUE),
      velocidad_prom_escalar = mean(!!sym(col_vel), na.rm = TRUE), # promedio aritmético
      n_horas = n(),   # número de horas con viento > 0 en el día
      .groups = "drop"
    )
  
  # Calcular dirección media a partir de (u_prom, v_prom)
  diario <- diario %>%
    mutate(
      dir_media_rad = atan2(u_prom, v_prom),   # resultado en [-pi, pi]
      dir_media_grados = dir_media_rad * 180 / pi,
      dir_media_grados = (dir_media_grados + 360) %% 360,  # pasar a [0,360)
      # Velocidad media vectorial (magnitud del viento resultante)
      vel_media_vectorial = sqrt(u_prom^2 + v_prom^2)
    ) %>%
    dplyr::select(site_id, fecha, lat, lon, n_horas, velocidad_prom_escalar, vel_media_vectorial, dir_media_grados)
  
  return(diario)
}

calcular_flecha <- function(lon, lat, direccion, velocidad, escala) {
  # Convertir dirección a radianes (convención matemática: 0° = Este, antihorario)
  # Nuestra dirección meteorológica: 0° = Norte, aumenta hacia el Este.
  # Para pasar a coordenadas cartesianas (dx, dy) donde:
  #   dx = velocidad * sin(ángulo_desde_Norte_hacia_Este)
  #   dy = velocidad * cos(ángulo_desde_Norte_hacia_Este)
  # Luego el destino será (lon + dx*escala, lat + dy*escala)
  rad <- direccion * pi / 180
  dx <- velocidad * sin(rad)      # desplazamiento en dirección Este (positivo)
  dy <- velocidad * cos(rad)      # desplazamiento en dirección Norte (positivo)
  xend <- lon + dx * escala
  yend <- lat + dy * escala
  return(data.frame(xend = xend, yend = yend))
}

get_consecutive_periods <- function(df, fecha = "fecha", dias_pre = 1, dias_post = 1) {
  # Ensure the date column is of class Date
  df <- df %>%
    mutate({{ fecha }} := as.Date(!!sym(fecha)))
  
  # Sort by date
  df <- df %>% arrange(!!sym(fecha))
  
  # Create a grouping variable: increment whenever gap > 1 day
  df <- df %>%
    mutate(gap = c(0, diff(!!sym(fecha)) > 1),
           group = cumsum(gap))
  
  # For each group, get first and last date
  result <- df %>%
    group_by(group) %>%
    summarise(fecha_inicio = min(!!sym(fecha)-dias_pre),
              fecha_fin = max(!!sym(fecha)+dias_post), 
              duracion = fecha_fin - fecha_inicio - 1,
              .groups = "drop") %>%
    dplyr::select(fecha_inicio, fecha_fin, duracion)
  
  return(result)
}

interpolar_viento_periodo <- function(df, # df <- viento_periodo,
                                      inicio = inicio, # inicio <- "2020-07-07"
                                      fin = fin, # fin <- "2025-09-27"
                                      col_sitio = "site_id",
                                      col_vel = "wind_speed_10m",
                                      col_dir = "wind_direction_10m",
                                      capa_raster,
                                      utm) {
  
  # Convertir a radianes y calcular componentes u, v
  # Nota: Convención meteorológica estándar (0° = Norte, ángulo aumenta hacia el Este)
  # Si tu convención es diferente (0° = Este, antihorario), ajusta los senos/cosenos.
  
  df_proc <- df %>%
    filter(fecha < fin & fecha > inicio) %>%
    group_by(!!sym(col_sitio)) %>%
    mutate(
      rad = !!sym(col_dir) * pi / 180,
      # Componentes vectoriales (viento hacia donde sopla, no desde)
      u = -!!sym(col_vel) * sin(rad),   # componente Este-Oeste
      v = -!!sym(col_vel) * cos(rad)    # componente Norte-Sur
    ) %>%
    
    # Excluir filas con velocidad = 0 (dirección indefinida)
    filter(!!sym(col_vel) > 0) %>%
    mutate(
      u_prom = mean(u, na.rm = TRUE),
      v_prom = mean(v, na.rm = TRUE),
      velocidad_prom_escalar = mean(!!sym(col_vel), na.rm = TRUE), # promedio aritmético
      n_horas = n(),   # número de horas con viento > 0 en el día
      .groups = "drop"
    ) %>%
    dplyr::select(site_id, lat, lon, n_horas, u_prom, v_prom, velocidad_prom_escalar, n_horas)
  
  df_proc <- unique(df_proc)
  
  # Convert points to SpatialPointsDataFrame for gstat
  df_proc_sp <- st_as_sf(df_proc, coords = c("lon", "lat"), crs = 4326) %>%
    st_transform(utm)
  
  # Define a template raster covering the extent of the points
  # r_template <- rast(capa_raster)
  
  # Setup the IDW model structure using gstat
  # 'value ~ 1' implies that the prediction relies strictly on spatial location
  # idw_model_u <- gstat(formula = u_prom ~ 1, data = df_proc_sp)
  # idw_model_v <- gstat(formula = v_prom ~ 1, data = df_proc_sp)
  
  # Try Kriging smoothing
  # Fit a variogram and run Kriging through terra's interpolate interface
   sf_points <- sf::st_as_sf(df_proc_sp)
   
   idw_model_u <- gstat(formula = u_prom ~ 1, data = df_proc_sp,
                        model = fit.variogram(variogram(u_prom ~ 1, sf_points), vgm(2, c("Exp", "Sph"), 500000, 1)))
   
   idw_model_v <- gstat(formula = v_prom ~ 1, data = df_proc_sp,
                        model = fit.variogram(variogram(v_prom ~ 1, sf_points), vgm(2, c("Exp", "Sph"), 500000, 1)))
  
  # Interpolate across the empty grid
  idw_raster_u <- interpolate(capa_raster, idw_model_u, debug.level = 0)
  idw_raster_v <- interpolate(capa_raster, idw_model_v, debug.level = 0)
  
  # 5. Back‑transform to speed and direction
  r_ws <- sqrt(idw_raster_u^2 + idw_raster_v^2)
  r_wd_rad <- atan2(idw_raster_u, idw_raster_v)
  
  coordinates(r_ws)
  
  df <- cbind(coordinates(r_ws),
              as.data.frame(r_ws),
              as.data.frame(r_wd_rad))
  
  names(df) <- c("x", "y", "speed", "angle")
  
  # Convert to meteorological direction (from which wind comes)
  r_wd_deg <- (270 - r_wd_rad * 180/pi) %% 360   # adjust for meteorological convention
  
  return(list(speed = r_ws,
              direction = r_wd_deg,
              df = df))
}

fetch_en_mar <- function(punto, direccion_viento_deg, pol_tierra, max_fetch = 500000) {
  # punto: sf POINT (en el mar)
  # direccion_viento: de dónde sopla el viento (0=N, 90=E)
  # pol_tierra: polígono de tierra (proyectado)
  # max_fetch: distancia máxima (metros) para limitar el rayo (ej. 500 km)
  
  # Dirección opuesta: hacia donde medimos el fetch (aguas arriba)
  #dir_fetch <- (direccion_viento_deg + 180) %% 360
  dir_fetch <- direccion_viento_deg
  rad <- dir_fetch * pi / 180
  dir_vec <- c(sin(rad), cos(rad))
  # punto <- resultados[60,]
  # plot(tierra_region)
  # points(punto)
  coords <- st_coordinates(punto)[1:2]
  end_x <- coords[1] + dir_vec[1] * max_fetch
  end_y <- coords[2] + dir_vec[2] * max_fetch
  rayo <- st_linestring(matrix(c(coords[1], coords[2], end_x, end_y), ncol=2, byrow=TRUE)) %>%
    st_sfc(crs = st_crs(pol_tierra))
  
  # Intersección del rayo con el polígono de tierra
  intersecciones <- st_intersection(rayo, tierra_region)
  
  if (length(intersecciones) == 0) {
    return(max_fetch)
  }
  
  # Distancias desde el punto a cada intersección
  distancias <- as.numeric(st_distance(punto, intersecciones))
  # Eliminar posibles intersecciones en el mismo punto (si el punto está sobre tierra - no debería)
  distancias_validas <- distancias[distancias > 0.001]
  if (length(distancias_validas) == 0) return(max_fetch)
  
  fetch <- min(distancias_validas)
  return(as.numeric(min(fetch, max_fetch)))
}

U <- 10  # m/s
g <- 9.81
Hs <- function(F, U){
  0.0016 * (U^2/g) * sqrt(g * F / U^2)  # opcional: dividir F/1000 para ajustar
} 
# Nota: la fórmula original daba valores muy altos; ajusta coeficiente según calibración.


# Buques -----

cleanFun <- function(htmlString) {
  return(gsub("<.*?>", "", htmlString))
}

get.vesselData <- function(buqueID){
  
  require(xml2)
  require(rvest)
  require(htmltools)
  require(Hmisc)
  
  info <- gfwr::gfw_vessel_info(search_type = 'id',
                                ids = as.character("2d757da57-7c47-a67d-7e55-ffb39691c259")
  )
  
  if(is.null(info)) {
    
    datos <- data.frame("ID" = buqueID,
                        "vesselOmi" = NA,
                        "vesselName" = NA,
                       "vesselType" = NA,
                       "vesselFlag" = NA,
                       "yearConstructed" = NA,
                       "vesselLength"= NA,
                       "vesselBeam" = NA,
                       "vesselGT" = NA
    )
    return(datos)
    
  } else {
    
    omi <- info$registryInfo$imo
    
    if(is.null(omi)) {
      
      mmsi <- info$selfReportedInfo$ssvid
      
      if(is.null(mmsi)) {
        
        datos <- data.frame("ID" = buqueID,
                            "vesselOmi" = NA,
                            "vesselName" = NA,
                            "vesselType" = NA,
                            "vesselFlag" = NA,
                            "yearConstructed" = NA,
                            "vesselLength"= NA,
                            "vesselBeam" = NA,
                            "vesselGT" = NA
        )
        
        return(datos)
        
      } else {
        
        if(length(mmsi)==1){
          mmsi <- mmsi
        } else {
          
          if(length(unique(omi)) == 1){
            mmsi <- mmsi[1]
          } else {
            mmsi <- mmsi[1]
            warning(paste0("Hay más de un identificador para el buque (SSVID:", mmsi, ")"))
          }
        }
        
        id <- as.numeric(mmsi)
        
        urlVF <- paste("https://www.vesselfinder.com/vessels/details/", id, sep = "")
        
        B <- try(read_html(urlVF), silent = T)
        
        if (inherits(B, "try-error")){
          
          warning("No se pudo realizar la búsqueda")
          
        } else {
          
          output <- html_elements(read_html(urlVF), "td.tpc2")
          
          datos <- data.frame("ID" = buqueID,
                              "vesselOmi" = as.numeric(cleanFun(output[1])),
                              "vesselName" = cleanFun(output[2]),
                              "vesselType" = cleanFun(output[3]),
                              "vesselFlag" = cleanFun(output[4]),
                              "yearConstructed" = as.numeric(cleanFun(output[5])),
                              "vesselLength"= as.numeric(cleanFun(output[6])),
                              "vesselBeam" = as.numeric(cleanFun(output[8])),
                              "vesselGT" = as.numeric(cleanFun(output[11]))
          )
          
          Sys.sleep(sample(5:10, 1))
          
          return(datos)
        
        }
      }
      
    } else {
      if(length(omi)==1){
        omi <- omi
      } else {
        if(length(unique(omi)) == 1){
          omi <- omi[1]
        } else {
          omi <- omi[1]
          warning(paste0("Hay más de un identificador para el buque (OMI:", omi, ")"))
        }
      }
      
      id <- as.numeric(omi)
      
      urlVF <- paste("https://www.vesselfinder.com/vessels/details/", id, sep = "")
      
      B <- try(read_html(urlVF), silent = T)
      
      if (inherits(B, "try-error")){
        
        warning("No se pudo realizar la búsqueda")
        
      } else {
        
        output <- html_elements(read_html(urlVF), "td.tpc2")
        
        datos <- data.frame("ID" = buqueID,
                            "vesselOmi" = as.numeric(cleanFun(output[1])),
                            "vesselName" = cleanFun(output[2]),
                            "vesselType" = cleanFun(output[3]),
                            "vesselFlag" = cleanFun(output[4]),
                            "yearConstructed" = as.numeric(cleanFun(output[5])),
                            "vesselLength"= as.numeric(cleanFun(output[6])),
                            "vesselBeam" = as.numeric(cleanFun(output[8])),
                            "vesselGT" = as.numeric(cleanFun(output[11]))
        )
        
        Sys.sleep(sample(5:10, 1))
        
        return(datos)
        
      }
    } 
  }
  
}


checkAndAdd.vesselData <- function(df, buqueID) {
  
  ids <- unique(df[,1])
  
  if(buqueID %in% ids) {
    warning("Este buque ya esta en la base de datos")
  } else {
    datos <- get.vesselData(buqueID)
    datos <- rbind(df, datos)
  }
  
}
