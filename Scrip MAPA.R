#Nombre:Scrip mapa, Proyecto R 
#Autor: Daniel Alonso-VIllarreal 
#Fecha: 19_octubre_25
#===============================
#libreria utilizada 

library(sf)
library(terra)
library(tidyverse) 
library(ggspatial) 
library(patchwork) 
library(scales)    
library(ggnewscale)
#=========================

#Parámetros estimados

#Alzheimer ====
#A nivel nacional: Un millón 300 mil personas padecen la enfermedad de Alzheimer en México, lo que representa entre el 60% y el 70% de los diagnósticos de demencia en adultos mayores de 65 años. En la Ciudad de México (estimado): Un reporte de la Fundación para Ancianos Concepción Béistegui, aunque de 2019, indicó que había 104,674 personas con algún tipo de demencia en la Ciudad de México. 
#====

#Parkison====
#A nivel nacional: La estimación de personas con Parkinson en México varía considerablemente según la fuente, con cifras recientes que van desde 125,000 hasta más de 500,000 casos. En la Ciudad de México: Al no haber cifras específicas para la ciudad, una estimación basada en la población de la Ciudad de México (aproximadamente 9.2 millones) y la tasa nacional de 50 casos nuevos por cada 100,000 habitantes anualmente, podría proporcionar una aproximación del riesgo. 
#====

#Huntington====
#A nivel nacional: Se estima que alrededor de 8,000 personas padecen la enfermedad de Huntington en México. En la Ciudad de México (estimado): No hay una cifra específica disponible para la Ciudad de México. A partir de la prevalencia nacional (1 caso por cada 10,000 habitantes) reportada en 2025 por la asociación Huntington México, y la población de la Ciudad de México (9.2 millones), se podría estimar que hay alrededor de 920 personas con esta enfermedad en la ciudad.
#====

#Parámetros Estimados
t_alzheimer_cdmx <- 104674
t_huntington_cdmx <- 920

poblacion_cdmx <- 9.2e6
poblacion_mexico <- 126e6
t_parkinson_cdmx <- round(125000 * (poblacion_cdmx / poblacion_mexico))

#Cargar polígonos de alcaldías (capa vectorial)

alcaldias_shp <- "poligonos_alcaldias_cdmx.shp" 

alcaldias <- st_read(alcaldias_shp) %>%
  st_transform(crs = 4326) # Reproyectar a WGS84 

# Cargar raster de población (capa raster)
raster_poblacion_mx <- "poblacion_mexico_2020.tif"

pop_rast_full <- rast(raster_poblacion_mx)

# extensión de las alcaldías, esto hace el procesamiento y ploteo más rápido

pop_rast_cdmx <- terra::crop(pop_rast_full, vect(alcaldias))

pop_rast_cdmx <- terra::mask(pop_rast_cdmx, vect(alcaldias))


# Calcular población total por alcaldía usando terra::extract
pop_extract <- terra::extract(pop_rast_cdmx, 
                              vect(alcaldias), 
                              fun = sum, 
                              na.rm = TRUE)

# Asignar población extraída al objeto 'alcaldias'

alcaldias$pop_tot <- pop_extract[, 2] 

# Estimar casos por alcaldía, donde "prop_pop" es la proporción de la población de la alcaldía vs. el total de poblacion de la CDMX

alcaldias <- alcaldias %>%
  mutate(prop_pop = pop_tot / sum(pop_tot, na.rm = TRUE),
    alzheimer_casos = round(prop_pop * t_alzheimer_cdmx),
    parkinson_casos = round(prop_pop * t_parkinson_cdmx),
    huntington_casos = round(prop_pop * t_huntington_cdmx))

#Prepararamos los datos para ggplot

# primero convertimos el RASTER de población a un dataframe
pop_df <- as.data.frame(pop_rast_cdmx, 
                        xy = TRUE, 
                        na.rm = TRUE)

#Renombraramos la columna de "población" 

names(pop_df)[3] <- "poblacion" 

# Limpiar valores de población

pop_df <- pop_df %>% filter(poblacion >= 0)

#Dado que trabajamos con un modelo de la enfermedad de huntington solo tomare los datos estimados de este sin embargo, genere las bases de datos por si se requiere modificar o "jugar" con los datos estimados 

# Mapa para Huntington ====

# Para hacer los otros mapas, solo cambia la variable en geom_sf y los títulos

mapa_huntington <- ggplot(data = alcaldias) +
  geom_sf(aes(fill = huntington_casos), 
          color = "white", 
          lwd = 0.3) +
  scale_fill_viridis_c(
    option = "plasma",
    labels = comma,
    name = "Casos Huntington \n(Estimado)") +
  labs( title = "Estimación de Casos de Huntington por Alcaldía, CDMX (2025)",
    subtitle = paste0("Total CDMX (estimado): ", 
                      format(t_huntington_cdmx, big.mark=",")),
    caption = "Elaboración propia. Fuente: Estimaciones de prevalencia y datos de WorldPop.") +
  annotation_scale(location = "bl", 
                   width_hint = 0.4, 
                   style = "ticks") +
  annotation_north_arrow(location = "tr", 
                         style = north_arrow_fancy_orienteering)+
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", 
                                  size = 16),
        axis.title = element_blank())

print(mapa_huntington)

ggsave("mapa_huntington_cdmx.png",
       path = "Figuras",
       plot = mapa_huntington,
       width = 8, 
       height = 8, 
       dpi = 300, 
       bg = "white")

#Crear zoom de escala en el mapa

#Identificar la alcaldía con más casos para el zoom(por probabilidad de dencidad es de poblacion con los datos que proporcionamos, teoricamente la fila 1 seria la de mayor casos teoricos de huntington)

alcaldia_max_casos <- alcaldias %>%
  arrange(desc(huntington_casos)) %>%
  slice(1) 

# Necesitamos el nombre de la columna que tiene el nombre de la alcaldía use 'NOM_ALCALD'. 

alcaldia_zoom_nombre <- alcaldia_max_casos$NOMGEO

#graficamos 

mapa_zoom <- ggplot() +
   geom_raster(data = pop_df, 
               aes(x = x, 
                   y = y, 
                   fill = poblacion)) +
  scale_fill_viridis_c(option = "cividis", 
                       trans = "log10", 
                       guide = "none") +
  geom_sf(data = alcaldia_max_casos, 
          aes(fill = huntington_casos), 
          color = "red", 
          lwd = 1) +
  scale_fill_viridis_c(
    option = "plasma", 
    trans = "sqrt",
    name = "Casos Estimados") +
# Aplicar el zoom usando las coordenadas de la alcaldía seleccionada
  coord_sf(
    xlim = st_bbox(alcaldia_max_casos)[c(1, 3)],
    ylim = st_bbox(alcaldia_max_casos)[c(2, 4)],
    expand = FALSE) +
  labs(
    title = paste(alcaldia_zoom_nombre),
    subtitle = paste("Población:", format(alcaldia_max_casos$pop_tot, big.mark=","),
                     "\nCasos (Huntington.):", format(alcaldia_max_casos$huntington_casos, big.mark=","))) +
  theme_void() +
  theme(
    plot.title = element_text( face = "bold",
                               size = 12,
                               hjust = 0.5),
    legend.position = "none", 
    panel.border = element_rect(color = "black", 
                                fill = NA, 
                                linewidth = 1))

# Imprimir el mapa de zoom para verlo
print(mapa_zoom)

#generamos un mapa convinado 

mapa_final <- mapa_huntington + mapa_zoom +
  plot_layout(
    widths = c(3, 1))

print(mapa_final)

ggsave("Mapa_huntington_CDMX.png", 
       path = "Figuras",
       plot = mapa_final, 
       width = 10, 
       height = 7, 
       units = "in", 
       bg="white",
       dpi = 320)  


#==========================================a







#Mapa Alzheimer====

mapa_principal <- ggplot() +
  geom_raster(data = pop_df, 
              aes(x = x, 
                  y = y, 
                  fill = poblacion)) +
  scale_fill_viridis_c(
    option = "cividis", 
    trans = "log10",
    labels = comma,
    name = "Densidad Pobl. \n(Raster)") +
  new_scale_fill() +
  geom_sf(data = alcaldias, 
          aes(fill = alzheimer_casos), 
          color = "white", 
          lwd = 0.3, 
          alpha = 0.7) +
  scale_fill_viridis_c(
    option = "magma", 
    trans = "sqrt",
    labels = comma,
    name = "Casos Alzheimer \n(Vector Estimado)") +
  labs(
    title = "Estimación de Casos de Alzheimer por Alcaldía, CDMX (2025)",
    subtitle = paste0("Total CDMX: ", format(t_alzheimer_cdmx, big.mark=","), " casos estimados, distribuidos por población (WorldPop 2020)"),
    caption = "Elaboración propia. Fuente: Estimaciones de prevalencia y datos de WorldPop."
  ) +
  annotation_scale(location = "bl", 
                   width_hint = 0.4, 
                   style = "ticks") +
    annotation_north_arrow(
    location = "tr", 
    which_north = "true",
    pad_x = unit(0.1, "in"), 
    pad_y = unit(0.2, "in"),
    style = north_arrow_fancy_orienteering) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", 
                              size = 16),
    plot.subtitle = element_text(size = 10, 
                                 colour = "grey30"),
    plot.caption = element_text(size = 8, 
                                colour = "grey50"),
    legend.title = element_text(size = 9, 
                                face = "bold"),
    legend.position = "right",
    axis.title = element_blank())

print(mapa_principal)

#Crear zoom de escala en el mapa

# Identificar la alcaldía con más casos para el zoom

alcaldia_max_casos <- alcaldias %>%
  arrange(desc(alzheimer_casos)) %>%
  slice(1) 

# Necesitamos el nombre de la columna que tiene el nombre de la alcaldía use 'NOM_ALCALD'. 
alcaldia_zoom_nombre <- alcaldia_max_casos$NOM_ALCALD 

mapa_zoom <- ggplot() +
  geom_raster(data = pop_df, 
              aes(x = x, 
                  y = y, 
                  fill = poblacion)) +
  scale_fill_viridis_c(option = "cividis", 
                       trans = "log10", 
                       guide = "none") +
  geom_sf(data = alcaldia_max_casos, 
          aes(fill = alzheimer_casos), 
          color = "red", 
          lwd = 1) +
  scale_fill_viridis_c(
    option = "magma", 
    trans = "sqrt",
    name = "Casos Estimados") +
  coord_sf(
    xlim = st_bbox(alcaldia_max_casos)[c(1, 3)],
    ylim = st_bbox(alcaldia_max_casos)[c(2, 4)],
    expand = FALSE) +
  labs(
    title = paste("Zoom:", alcaldia_zoom_nombre),
    subtitle = paste("Población:", format(alcaldia_max_casos$pop_tot, big.mark=","),
                     "\nCasos (Alz.):", format(alcaldia_max_casos$alzheimer_casos, big.mark=","))) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", 
                              size = 12),
    plot.subtitle = element_text(size = 9),
    legend.position = "none")


print(mapa_zoom)

#Combinar Mapas 

# Combinar ambos mapas: el principal (A) y el zoom (B)
mapa_final_combinado <- mapa_principal + mapa_zoom +
  plot_layout(widths = c(2.5, 1)) 


print(mapa_final_combinado)
#=======================





