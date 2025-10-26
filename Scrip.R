
#nombre: " Proyecto final clase R"
#realizo: Daniel Alonso Villarreal 
#fecha: 01 octubre_25
#==========================================
#libreria utilizada ====

library(tidyverse)
library(file2meco)
library(microeco)
library(ggplot2)
library(magrittr)
library(dplyr)
library(GUniFrac)
library(tidyr)
library(phyloseq)
library(patchwork)
library(ggpubr)

#=========================================

#Carga base de datos====
abund_file.path <- ("DATA/QZA/table-dada2.qza")
muestra.path <- ("DATA/sample-metadata_V2.txt") 
a_taxonomia.path <- ("DATA/QZA/taxonomy.qza") 
Arbol_filogenetico.path <- ("DATA/QZA/rooted-tree.qza")
seq_datos.path <- ("DATA/QZA/seqs-dada2.qza")


#==========================================

#aqui lo que estoy haciendo es tranformar los datos de microbiota obtenidos, a datos legibles en tidyverse, mediante la paqueteria "microeco". 

Datos_mb <- qiime2meco(
  feature_table =  "DATA/QZA/table-dada2.qza",
  taxonomy_table = "DATA/QZA/taxonomy.qza",
  sample_table = "DATA/sample-metadata_V2.txt",
  phylo_tree = "DATA/QZA/rooted-tree.qza",
  rep_fasta = "DATA/QZA/seqs-dada2.qza",
  auto_tidy = TRUE)

#ahora sembramos una semilla para que los datos sean reproducibles
set.seed(123)

#Mi seccion de Datos =======
#Los datos que usare son como tal una microtabla que en R hace referencia  a un tipo de objeto que se llama "microtable" por eso se instalo la libreria microeco, esta libreria en particular nos ayuda a analizar datos en metagenomica, en estas microtablas se almacenan datos  básicos como tablas de abundancia de muestras y metadatos.


#que datos usar====
#En mi base de datos tengo grupos que no me pertenecen y no puedo usar pero dado que el objeto Datos_mb es una microtabla no puedo usar la funcion summary, por lo que para saber cuantas filas tiene y donde se encuentran mis datos primedo debo usar la funcion "Datos_mb[["sample_table"]][["Grupo"]]; #yo se que mis datos son los ultimos 20 por lo que para visualizarlos uso. 
#====

Datos_mb[["sample_table"]][["Grupo"]]

tail(Datos_mb$sample_table,20)


#Inspeccion incial para determinar la estructura de los datos de abundancia microbiana ====
#aqui lo que estamos haciendo es una inspección inicial de las primeras 8 muestras para ver cómo están estructurados los datos de abundancia microbiana. Ya que en la tabla de abundancia de OTU (o ASV) las filas suelen ser las muestras y las columnas son los identificadores taxonómicos (OTU/ASV).
#====

head(Datos_mb$otu_table, 8)


# Verificando que los identificadores microbianos====
#ahora lo que estamos haciendo es verificando que los identificadores microbianos (OTU/ASV) tengan asignadas correctamente sus clasificaciones taxonómicas, ya que, tax tale es una matriz que contiene la información de clasificación biológica para cada OTU (o ASV) de un conjunto de datos. Las filas corresponden a los OTU/ASV y las columnas contienen sus clasificaciones (Reino, Filo, Clase, Orden, Familia, Género, Especie).
#====

head(Datos_mb$tax_table, 8)

#una vez determinados los datos a usar y verificamos si sample_table, otu_tabe y tax table son data_frame 

class(Datos_mb$otu_table)
class(Datos_mb$sample_table)
class(Datos_mb$tax_table)


#Filtrar los datos
#Aqui debemos filtrar grupos especificos====
#dado que la base de datos contiene datos que no quiero usar, debemos filtrar los grupos experimentales que nos interesan los cuales son "Control2, 3NP2, CGA,3NP_ CGA"
#====

#primero generamos una copia de seguridad 

Datos_mb2 <- unserialize(serialize(Datos_mb, NULL))

grupos_mantener <- c("Control2", "3NP2", "CGA", "3NP_CGA")

# Filtrar la tabla de muestras de mb2
Datos_mb2$sample_table <- Datos_mb2$sample_table %>%
  filter(Grupo %in% grupos_mantener)


Datos_mb2[["sample_table"]][["Grupo"]]

#para hacer que siempre respete el orden de las columnas 

orden_datos <- c("Control2", "3NP2", "CGA", "3NP_CGA")

#aqui cambiamos el nombre de las columnas 
n_nombres <- c("Control", "3NP", "CGA", "3NP_CGA")

Datos_mb2$sample_table$Grupo <- factor(Datos_mb2$sample_table$Grupo,levels = orden_datos,
                                       labels = n_nombres)

#verificamos el cambio de nombre 
Datos_mb2[["sample_table"]][["Grupo"]]

#Ordenar los datos dentro de Datos_mb====
#Ahora lo que debemos hacer es estandarizar y formatear los datos taxonómicos dentro del objeto Datos_mb, para que sean fáciles de usar y analizar. Para lo cual utilizaremos la funcion :tidy_taxonomy, ya que esta funcion lo que hace es  "ordenar"  los datos taxonómicos, evitando redundancias y poniendolos todos en un formato legible e identico para todos los datos. Ademas estos datos se anexaran a un nuevo objeto, pero guardando una copia de los datos originales, en este caso usaremos el "%<>%" (Pipe de Asignación Compuesta): Transfiere el objeto original, realiza la secuencia de operaciones y sobrescribe el objeto original con el resultado final.
#====

Datos_mb2$tax_table %>% tidy_taxonomy()

mb <- Datos_mb2

mb_copia <- mb

class(mb)

#Ahora filtraremos los taxa que nos interesan, en mi caso especifico es Archea y Bacteria. 

mb$tax_table <- mb$tax_table %>%
  filter(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")

# Ahora dentro de mi base de datos se encuentran OTUs con asignación taxonómica de “mitochondria”#==== 
#or “chloroplast”, por lo que debemos eliminar las líneas que contienen la palabra del taxón independientemente de los rangos taxonómicos e ignorar las mayúsculas y minúsculas en la "tax_table"
#=====

mb$filter_pollution(taxa = c("mitochondria", "chloroplast"))

#aqui para poder filtrar dos taxones necesito agruparlos en un vector para lo cual estoy utilizando la funcion C(), la cual agurpa datos en vectores, una vez aplicados hay que organizar de nuevo los datos en un formato tidy 

mb$tidy_dataset()

mb

# Con la función sample_sums() checamos el número de secuencias en cada muestra, en nuestro caso usaremos 20,000 para casa muestra 
mb$sample_sums() %>% range()

mb$rarefy_samples(sample.size = 20000)

#revisamos si modifico el numero de secuencia por muestra 
mb$sample_sums() %>% range

#calcular diversidad alpha, calculando el indice de diversidad filogenetica de Faith 

mb$cal_alphadiv(PD= FALSE)

class(mb$alpha_diversity)

#Guardar datos de diversidad alpha====    
#ahora guardademos por separado la base de datos de diversidad alfa 

mb$save_alphadiv(dirpath = "Salidas")

#Calcular la diversidad beta====
#aqui utilizaremos la funcion "UniFrac" la cual es una métrica de diversidad beta que incluye información filogenética. de modo que no solo se comparan las especies que se encuentran presentes, sino también qué tan relacionadas están entre sí, para poder calcular esto ahora debemos tener cargada la libreria "GUniFrac"
#====

mb$cal_betadiv(unifrac = TRUE)

class(mb$beta_diversity)

#guardar resultados de divercidad beta 
mb$save_betadiv(dirpath = "Salidas")

mb

# Genero una copia de seguridad 

mb_v2 <- unserialize(serialize(mb, NULL))

# Establecemos la abundancia relativa mayor a 0.001 (0.1%)

mb_v2$filter_taxa(rel_abund = 0.0001, freq = 0.1)

#en mi caso me interesan los 10 generos y familias mas abundantes, ya que es lo que usualmente se reporta en los articulos de microbiota, por lo que ahora debemos filtrar las 10 familias mas abundantes en nuestra base de datos 

f1 <- trans_abund$new (dataset= mb_v2, 
                       taxrank = "Family", 
                       ntaxa = 10)


f1$plot_bar(others_color = "grey50",
            facet = "Grupo",
            xtext_keep = FALSE, 
            legend_text_italic = TRUE,
            bar_full = TRUE,
            barwidth = 0.8)



# Personalización del formato de color y tipo de barra ====
f1$plot_bar(manual.palette <- c (
                        "#0000FF", "#FF0000", "#00FF00",
                       
                       "#808080","#ADD8E6","#FFFF00",
                       
                       "#FF00FF", "#FF8C00","#00FFFF",
                       
                       "#4B0082", "#808080", "#90EE90"), 
                       others_color = "grey70",                                                 bar_full = TRUE,
                      barwidth = 0.8,
                      facet = "Grupo",
                      xtext_keep = FALSE)
#====

#como se puede ver me da las familias taxonomicas por cada individuo dentro del grupo experimental por lo que debemos calcular los promedios de cada grupo experimental. 

#calcular promedios de grupos para familia 
t1_prom_family <- trans_abund$new(dataset = mb_v2,
                           taxrank = "Family",
                           ntaxa = 10,
                           groupmean = "Grupo",
                           group_morestats = TRUE)

#Graficamos 
g1 <- t1_prom_family$plot_bar(others_color = "grey70")

g1 <- g1 + 
  theme_classic() + 
  theme(
    axis.title.y = element_text(size = 20),    
    axis.text.x = element_text(size = 16),     
    axis.title.x = element_text(size = 16),
    legend.text = element_text(size = 15, face = "italic"),
    legend.title = element_text(size = 17),
    legend.position = "right")

print(g1)


ggsave("Barras_familia.png", 
       path = "Figuras",
       plot = g1, 
       width = 9, #ancho de la figura
       height = 7, #alto de la figura
       units = "in", #unidad del tamaño figura "pulgadas"
       bg="white",
       dpi = 320) #calidad de la figura 

#calcular promedios de grupos para genero

t1_prom_genero <- trans_abund$new(dataset = mb_v2,
                                  taxrank = "Genus",
                                  ntaxa = 10,
                                  groupmean = "Grupo",
                                  group_morestats = TRUE)

#Graficamos 

g2 <- t1_prom_genero$plot_bar(others_color = "grey70")

g2 <- g2 + 
  theme_classic() + 
  theme(
    axis.title.y = element_text(size = 20),    
    axis.text.x = element_text(size = 16),     
    axis.title.x = element_text(size = 16),
    legend.text = element_text(size = 15, face = "italic"),
    legend.title = element_text(size = 17),
    legend.position = "right")

print(g2)

ggsave("Barras_Genero.png", 
       path = "Figuras",
       plot = g2, 
       width = 9, 
       height = 7, 
       units = "in", 
       bg="white",
       dpi = 320)  

#hacemos un mapa de calor para familia ====

g3 <- t1_prom_family$plot_heatmap(xtext_keep = FALSE, 
                                    withmargin = FALSE, 
                                    plot_breaks = c(0.01, 
                                                    0.1,
                                                    1, 
                                                    10))+
  theme(axis.text.y = element_text( face = "italic",
                                    size = 13))+
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5, 
                                   size = 12.5
                                   ))
print(g3)

ggsave("Mapa_calor_Familia.png", 
       path = "Figuras",
       plot = g3, 
       width = 9, 
       height = 7, 
       units = "in", 
       bg="white",
       dpi = 320) 


#mapa de calor para genero====

g4 <- t1_prom_genero$plot_heatmap(xtext_keep = FALSE, 
                                  withmargin = FALSE, 
                                  plot_breaks = c(0.01, 
                                                  0.1, 
                                                  1, 
                                                  10))+
  theme(axis.text.y = element_text(face = "italic",
                                   size = 13)) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5, 
                                   size = 12.5))

print(g4)

ggsave("Mapa_calor_Genero.png", 
       path = "Figuras",
       plot = g4, 
       width = 9, 
       height = 7, 
       units = "in", 
       bg="white",
       dpi = 320) 

#====

#ahora vamos a crear un objeto para calcular la diversidad alpha

#utilizamos la funcion trans$alpha$new tomamos la base de datos que hemos usado mb_v2, y especificamos que lo haga solo por grupo 

t1_div_alpha <- trans_alpha$new(dataset = mb_v2,
                             group = "Grupo")

#revisamos si existe la tabla estadictica que contenga los valores de diversidad alfa calculados para cada muestra, junto con estadísticas descriptivas (como la media, mediana, desviación estándar) para cada nivel de la variable  "Grupo".

t1_div_alpha$data_stat


#ejecutamos una la prueba estadistica de wilcoxon para determinar si existen diferencias significativas 

t1_div_alpha$cal_diff(method = "wilcox")

print(t1_div_alpha$res_diff, max = 99999)

#graficamos diversidad alpha con indice Shannon 

g_div_alfa_1 <- t1_div_alpha$plot_alpha(measure = "Shannon",
                        add = "dotplot",
                        xtext_size = 15) + 
  ggtitle("Diversidad Alfa por grupos")+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 0.5,
                                  face = "bold"))

print(g_div_alfa_1)

#graficamos diversidad alpha con indice Simpson

g_div_alfa_2 <- t1_div_alpha$plot_alpha(measure = "Simpson",
                        add = "dotplot",
                        xtext_size = 15) + 
  ggtitle("Diversidad Alfa por grupos")+
  theme(plot.title = element_text(size = 18, 
                                  hjust = 0.5,
                                  face = "bold"))
print(g_div_alfa_2)

#Generamos una grafica compuesta 

g_div_alfa_t <- g_div_alfa_1 + g_div_alfa_2

print(g_div_alfa_t)

ggsave("Diversidad alpha.png", 
       path = "Figuras",
       plot = g_div_alfa_t, 
       width = 9, 
       height = 7, 
       units = "in", 
       bg="white",
       dpi = 320)


 #====

#calculamos divercidad beta====
# "bray" se refiere a la disimilitud de Bray-Curtis

t1_div_beta<- trans_beta$new(dataset = mb_v2, 
                      group = "Grupo", 
                      measure = "bray")
#ahora dado que tenemos una matriz con muchas dimenciones en la diversidad beta, utilizaremos la funcion "PCoA"(Principal Coordinates Analysis o Análisis de Coordenadas Principales) para ordenar la matriz y poder visualizarla de modo 2D, ademas este metodo es el metodo estandar para usar con matrices resultantes del indice Bray-Curtis.

t1_div_beta$cal_ordination(method = "PCoA")

class(t1_div_beta$res_ordination)

#graficamos 

PCoA_1 <- t1_div_beta$plot_ordination(plot_color = "Grupo",
                                      point_size = 5,
                                      plot_type = c("point", "ellipse"),
                                      ellipse_chull_fill = FALSE,
                                      ellipse_chull_alpha = 0.1) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 18),    
        axis.text.x = element_text(size = 12),     
        axis.title.x = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "right")  

print(PCoA_1)


#box plot para diversidad beta por grupos 

t1_div_beta$cal_group_distance(within_group = TRUE)

t1_div_beta$cal_group_distance_diff(method = "wilcox")

g1_div_beta <- t1_div_beta$plot_group_distance(add = "mean",) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 18),    
    axis.text.x = element_text(size = 12),     
    axis.title.x = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.position = "right")


print(g1_div_beta)


#anova diversidad beta
t1_div_beta$cal_manova(manova_all = TRUE)
t1_div_beta$res_manova
resultado_global_div_B <- t1_div_beta$res_manova

print(resultado_global_div_B)

t1_div_beta$cal_manova(manova_all = FALSE)
t1_div_beta$res_manova
resultado_porgrupo_div_B <- t1_div_beta$res_manova

print(resultado_porgrupo_div_B)

#taxas significatvos de genero por grupos 

t2_grupo<- trans_diff$new(dataset = mb_copia, 
                     method = "lefse", 
                     group = "Grupo", 
                     taxa_level = "Genus", 
                     alpha = 0.05, 
                     lefse_subgroup = NULL, 
                     p_adjust_method = "none")

#Graficamos 
g5 <- t2_grupo$plot_diff_bar(use_number = 1:10, width = 0.6)+
  xlab("Genero")+
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 20, 
                                  face = "bold"))

print(g5)


#taxas significatvos de familias por grupos 

t2_familia <- trans_diff$new(dataset = mb_copia, 
                     method = "lefse", 
                     group = "Grupo", 
                     taxa_level = "Family", 
                     alpha = 0.05, 
                     lefse_subgroup = NULL, 
                     p_adjust_method = "none")
  

g6 <- t2_familia$plot_diff_bar(use_number = 1:10, width = 0.6)+
  xlab("Familia")+
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 20, 
                                  face = "bold"))
 
print(g6)

#Generamos una grafica compuesta 

g7_sig <- g6+g5

print(g7_sig)

ggsave("comparacion abundancia.png",
       path = "Figuras",
       plot = g7_sig,
       width = 15, 
       height = 7, 
       units = "in", 
       bg="white",
       dpi = 320)



















