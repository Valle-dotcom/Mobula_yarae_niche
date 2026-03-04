---
title: "Ejercicio de modelado de nicho de Mobula yarae"
---

# Flujo de trabajo
<nav>
  <a href="#paso-1">I. Preparación de datos</a> |
  <a href="#paso-2">II. </a> |
  <a href="#paso-3">III. </a> |
  <a href="#paso-4">Paso</a> |
  <a href="#paso-5">Anexo: Instalación de ellipsenm</a> |
</nav>

---
## Introducción
El siguiente ejercicio está diseñado con la idea de que el estudiante pueda poner en práctica los pasos esenciales del flujo de trabajo para la 
realización de modelos de distribución de especies (MDE o SDM en inglés) y modelos de nicho ecológico (MNE o ENM en inglés) para organismos marinos. 
Este *workflow* está basado y adaptado a partir del manual de [Simoes et al. (2020)](https://doi.org/10.17161/bi.v15i2.13376), el cual recomiendo ampliamente 
como una guía práctica introductoria al tema. Además, recomiendo también para la limpieza de datos el manual de [Cobos et al. (2018)](https://doi.org/10.17161/bi.v13i0.7600) y el libro de 
[Peterson et al. (2011)](https://doi.org/10.1515/9781400840670?urlappend=%3Futm_source%3Dresearchgate.net%26utm_medium%3Darticle) como texto base para esta disciplina.

Este ejercicio tiene como objetivo realizar un SDM y un ENM (para aprender más sobre las diferencias conceptuales entre estos dos tipos de modelos recomiendo
[Soberón et al. 2017](https://doi.org/10.1016/j.rmb.2017.03.011)) de la mantarraya del Atlántico (*Mobula yarae*), una especie recién descrita por [Bucair et al. (2025)](https://doi.org/10.1007/s10641-025-01727-2) y de la cual aún no se conocen varios aspectos de su historia natural, incluyendo su área de distribución.
En este ejercicio usaremos [*ellipsenm* (Cobos et al. 2020)](https://github.com/marlonecobos/ellipsenm2) un paquete de R que utiliza envolturas elipsoidales para caracterizar el nicho de la especie [(Farber and Kadmon 2003)](https://doi.org/10.1016/S0304-3800(02)00327-7)





<a id="paso-1"></a>
## I. Preparación de datos

Debido a que *M. yarae* es una especie recién descrita, los datos presentes en [GBIF](https://www.gbif.org/es/)
aún no presentan registros de esta especie. Ya que *M. yarae* es simpátrica con *M. birostris* [(Bucair et al.,2025)](https://doi.org/10.1007/s10641-025-01727-2), es necesario volver a reclasificar los registros de *M. birostris* que se encuentran dentro del área de distribución de *M. yarae*. Si bien para esto, y según [Bucair et al.,2025](https://doi.org/10.1007/s10641-025-01727-2) es necesario una identificación genética; para fines de este ejercicio vamos a realizar un proceso de limpieza estandarizado para MDE/MNE con los registros provenientes de GBIF para *M. birostris*.

El archivo se encuentra dentro de la carpeta [input](https://github.com/Valle-dotcom/Mobula_yarae_niche/tree/main/input) con el nombre de *m_birostris.csv* Este es el [DOI](https://doi.org/10.15468/dl.jnav8h) con los metadatos de la base de datos descargada 

### 1) Datos de ocurrencia
#### 1. Preparación del Entorno y Carga de Datos 
Antes de iniciar cualquier, es necesario establecer un entorno de trabajo limpio y cargar todos los paquetes al inicio para asegurar que las funciones necesarias estén disponibles. 

```r
# --- 1. PREPARACIÓN DEL ENTORNO ---

# Carga de librerías para manejo de datos, limpieza y análisis espacial
library(dplyr)
library(CoordinateCleaner)
library(terra)
library(geodata)
library(spThin)
library(rnaturalearth)
library(rnaturalearthdata)

# Definición del directorio de trabajo e importación
setwd("ruta")
occurrences_raw <- read.csv("m_bevirostris.csv")
```

#### 2. Exploración Visual del Espacio Geográfico
Visualizar los datos crudos permite realizar un diagnóstico rápido para detectar errores comunes 

```r
# --- 2. EXPLORACIÓN VISUAL ---

# 1. Filtro de seguridad informático (remover NAs en coordenadas)
occ_plot <- occurrences_raw[!is.na(occurrences_raw$decimalLongitude) & 
                              !is.na(occurrences_raw$decimalLatitude), ]

# 2. Conversión a vector espacial (SpatVector)
occ_raw_vect <- vect(occ_plot, 
                     geom = c("decimalLongitude", "decimalLatitude"), 
                     crs = "EPSG:4326")

# 3. Carga del mapa base de forma segura (offline)
# Extraemos el mapa como objeto 'sf' y lo convertimos a 'SpatVector' para terra
world_sf <- ne_countries(scale = "medium", returnclass = "sf")
world_map <- vect(world_sf)

# 4. Visualización global (Ahora terra::plot() controlará el gráfico)
plot(world_map, col = "antiquewhite", border = "gray50",
     main = "Distribución Cruda de M. bevirostris (Sin Curaduría)",
     background = "aliceblue",
     mar = c(3, 3, 2, 2))

# Añadimos los puntos empíricos
plot(occ_raw_vect, col = "darkred", pch = 20, cex = 0.8, add = TRUE)

# Líneas de referencia para identificar errores de "Null Island"
abline(h = 0, v = 0, col = "blue", lty = 2)
```

#### 3. Limpieza Temática, Temporal y Espacial
Este bloque es el núcelo de la limpieza de datos. 

```r
# --- 3. LIMPIEZA RIGUROSA ---

# Selección de columnas críticas
occurrences <- occurrences_raw %>% 
  dplyr::select(gbifID, species, decimalLongitude, decimalLatitude, 
                countryCode, stateProvince, locality, year)

# Filtro temporal (1955-2010)
occurrences <- subset(occurrences, year >= 1955 & year <= 2010)

# Remoción de NAs y ceros en coordenadas
occurrences <- occurrences[!is.na(occurrences$decimalLongitude) | !is.na(occurrences$decimalLatitude), ]
occurrences <- occurrences[occurrences$decimalLongitude != 0 & occurrences$decimalLatitude != 0, ]
occurrences <- occurrences[!is.na(occurrences$year), ]

# Función de precisión decimal y filtro (retener > 2 decimales)
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
occurrences <- occurrences[sapply(occurrences$decimalLongitude, decimalplaces) > 2 &
                           sapply(occurrences$decimalLatitude, decimalplaces) > 2, ]

# Limpieza de duplicados con CoordinateCleaner
occurrences <- cc_dupl(occurrences, lon = "decimalLongitude", lat = "decimalLatitude", species = "species")
```

#### 4. Delimitación dentro de la Área M
Según el diagrama BAM, el algoritmo debe entrenarse exclusivamente dentro de las regiones que la especie ha podido explorar a lo largo de su historia evolutiva (el área M).

```r
# --- 4. RESTRICCIÓN AL ÁREA DE ESTUDIO (M) ---

# Vectorización de los datos limpios
occ_vect <- vect(occurrences, 
                 geom = c("decimalLongitude", "decimalLatitude"), 
                 crs = "EPSG:4326")

# Carga del polígono de estudio
zona_estudio <- vect("D:/lab/Clases/shapefiles/mi_zona_estudio.shp")

# Homologación de sistemas de referencia (CRS)
if (crs(occ_vect) != crs(zona_estudio)) {
  message("Reproyectando shapefile...")
  zona_estudio <- project(zona_estudio, crs(occ_vect))
}

# Recorte espacial (Clipping)
occ_M <- occ_vect[zona_estudio, ]

# Respaldo de datos restringidos a M
occ_finales_df <- as.data.frame(occ_M, geom = "XY")
write.csv(occ_finales_df, "m_bevirostris_dentro_M.csv", row.names = FALSE)

# Comprobación visual
plot(zona_estudio, main = "Registros en el Área Accesible (M)", col = "aliceblue")
plot(occ_M, col = "darkred", pch = 20, cex = 1.2, add = TRUE)
```

#### 5. Spatial thining
Los esfuerzos de muestreo no son aleatorios; tienden a concentrarse cerca de infraestructura humana. Se debe de corrigir este *spatial bias*.

```r
# --- 5. SPATIAL THINNING ---

# Aplicación de spThin sobre los datos recortados
# Nota para alumnos: verifiquen qué data.frame están introduciendo aquí.
p1 <- thin(loc.data = occ_finales_df, 
           lat.col = "y", long.col = "x", 
           spec.col = "species", 
           thin.par = 5, reps = 100, 
           locs.thinned.list.return = TRUE, 
           write.files = TRUE,
           max.files = 1,
           out.dir = "E:/lab/tesis/bd_5",
           out.base = "thinned",
           write.log.file = FALSE)

# Visualización y resumen del adelgazamiento espacial
plotThin(p1)
summaryThin(p1)
```
### 2) Variables ambientales
#### 1. Preparación del Entorno y Carga de Datos 

```r
# --- 1. PREPARACIÓN ---

install.packages(c("terra","predicts","corrplot","rJava"))

library(terra)
library(predicts)
#instalar Java, descargar maxent, system.file("java", package="predicts")
library(corrplot)
library(rJava)

#Objetos


# 1. Cargar variables

clima_M <- rast("C:/Users/jovac/Downloads/MARSPEC_5m/Variables/MARSPEC_recortado.tif")
print(names(clima_M)) # Comprobamos las 18 variables

# 2. Cargar las presencias


# 1. Leer el archivo CSV 
tabla_presencias <- read.csv("D:/lab/Clases/myarae_niche/m_yarae.csv")

print(names(tabla_presencias))

# 2. Convertir la tabla a un objeto espacial (SpatVector)

occ_thin <- vect(tabla_presencias, 
                 geom = c("longitude", "latitude"),  
                 crs = "EPSG:4326")

print(occ_thin)

# Opcional: Una comprobación visual rápida para asegurar que los puntos caen sobre el mapa
plot(clima_M[[1]], main = "Comprobación de Datos Cargados")
plot(occ_thin, add = TRUE, col = "red", pch = 16, cex = 0.5)

```
#### 2. Jacknnife en maxent
```r
# --- 2. JACKNIFE MAXENT---


# Definimos el vector de argumentos (args) 
mis_argumentos <- c(
  "autofeature=false",   
  "linear=true",         
  "quadratic=true",      
  "product=true",        
  "hinge=false",         
  "threshold=false",     
  "responsecurves=true", 
  "pictures=true",       
  "jackknife=true",      
  "outputformat=raw",    
  "outputfiletype=asc"   
)
library(terra)
library(predicts)


coords_presencia <- crds(occ_thin)

puntos_fondo <- spatSample(clima_M, size = 10000, method = "random", na.rm = TRUE, xy = TRUE)

coords_fondo <- puntos_fondo[, c("x", "y")]





# Ejecutamos el modelo 
modelo_exploratorio <- MaxEnt(x = clima_M, 
                              p = coords_presencia, 
                              a = coords_fondo,
                              args = mis_argumentos)

# Diagnóstico
plot(modelo_exploratorio, main = "Importancia de Variables (L, Q, P - Raw)")
modelo_exploratorio # Abre el HTML en el navegador
```
#### 2. Correlación de Pearson
```r
# --- 3. EXTRACCIÓN Y CORRELACIÓN DE PEARSON ---

# Extraer valores climáticos
valores_presencia <- extract(clima_M, occ_thin, ID = FALSE)
valores_presencia <- na.omit(valores_presencia)

# Matriz de correlación y visualización
matriz_cor <- cor(valores_presencia)
corrplot(matriz_cor, method = "number", type = "lower", 
         tl.col = "black", tl.cex = 0.8, number.cex = 0.7,
         main = "Matriz de Correlación de Pearson")
```



*CONTINUARÁ..........*

<a id="paso-2"></a>
## Paso 2 — Datos
- Descargar datos
- Limpiar datos
- Formato esperado

<a id="paso-3"></a>
## Paso 3 — Variables ambientales
- Fuente
- Recorte / resolución
- Preparación

<a id="paso-4"></a>
## Paso 4 — Modelado
- Algoritmo
- Parámetros
- Ejecución

<a id="paso-5"></a>
## Anexo — Instalación de ellipsenm
Vamos a utilizar el paquete de R [*ellipsenm* (Cobos et al. 2020)](https://github.com/marlonecobos/ellipsenm2) para realizar el modelo. La instalación de este paquete se realiza de la siguiente manera:

```r

# Installing and loading packages
if(!require(remotes)){
    install.packages("remotes")
}
if(!require(ellipsenm)){
    remotes::install_github("marlonecobos/ellipsenm2")
}
library(ellipsenm)

```

- Mapas
- Métricas
- Exportar productos

---

### Volver arriba
[↑ Ir al índice](#indice)
