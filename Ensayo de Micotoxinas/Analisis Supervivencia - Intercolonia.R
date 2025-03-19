#Librerias 

#install.packages("dplyr")  
#install.packages("tidyr")
#install.packages("survival")
#install.packages("bdsmatrix")
#install.packages("FSA")

library(dplyr)  # Para combinar los data frames de forma eficiente
library(tidyr)
library(survival)
library(bdsmatrix)
library(FSA)  # Para pruebas post hoc

#Análisis de Efectos de supervivencia Mixtos

#Efecto mixto - Cox


ensayo_10 <- read.csv(file = "Ensayos de Micotoxinas de B. bassiana en L. humile - Ensayo #10 - Para analisis R.csv", header = T, sep=",", skip=1)

df_ensayo_10 <- data.frame(tiempoControlSac = ensayo_10$tiempo_dias,
                           controlSac = ensayo_10$X.Muertas,
                           tiempoControlPGBQ = ensayo_10$tiempo_dias.1,
                           controlPGBQ = ensayo_10$X.Muertas.1,
                           tiempoBbQ150= ensayo_10$tiempo_dias.2,
                           BbQ150 = ensayo_10$X.Muertas.2,
                           tiempoBbQ300= ensayo_10$tiempo_dias.3,
                           BbQ300 = ensayo_10$X.Muertas.3,
                           tiempoBbQ450= ensayo_10$tiempo_dias.4,
                           BbQ450 = ensayo_10$X.Muertas.4,
                           tiempoBbQ600= ensayo_10$tiempo_dias.5,
                           BbQ600 = ensayo_10$X.Muertas.5,
                           tiempoBbQ750= ensayo_10$tiempo_dias.6,
                           BbQ750 = ensayo_10$X.Muertas.6)

ensayo_11 <- read.csv(file = "Ensayos de Micotoxinas de B. bassiana en L. humile - Ensayo #11 - Para analisis R.csv", header = T, sep=",", skip=1)

df_ensayo_11 <- data.frame(tiempoControlSac = ensayo_11$tiempo_dias,
                           controlSac = ensayo_11$X.Muertas,
                           tiempoControlPGBQ = ensayo_11$tiempo_dias.1,
                           controlPGBQ = ensayo_11$X.Muertas.1,
                           tiempoBbQ150= ensayo_11$tiempo_dias.2,
                           BbQ150 = ensayo_11$X.Muertas.2,
                           tiempoBbQ300= ensayo_11$tiempo_dias.3,
                           BbQ300 = ensayo_11$X.Muertas.3,
                           tiempoBbQ450= ensayo_11$tiempo_dias.4,
                           BbQ450 = ensayo_11$X.Muertas.4,
                           tiempoBbQ600= ensayo_11$tiempo_dias.5,
                           BbQ600 = ensayo_11$X.Muertas.5,
                           tiempoBbQ750= ensayo_11$tiempo_dias.6,
                           BbQ750 = ensayo_11$X.Muertas.6)

ensayo_12 <- read.csv(file = "Ensayos de Micotoxinas de B. bassiana en L. humile - Ensayo #12 - Para analisis R.csv", header = T, sep=",", skip=1)

df_ensayo_12 <- data.frame(tiempoControlSac = ensayo_12$tiempo_dias,
                           controlSac = ensayo_12$X.Muertas,
                           tiempoControlPGBQ = ensayo_12$tiempo_dias.1,
                           controlPGBQ = ensayo_12$X.Muertas.1,
                           tiempoBbQ150= ensayo_12$tiempo_dias.2,
                           BbQ150 = ensayo_12$X.Muertas.2,
                           tiempoBbQ300= ensayo_12$tiempo_dias.3,
                           BbQ300 = ensayo_12$X.Muertas.3,
                           tiempoBbQ450= ensayo_12$tiempo_dias.4,
                           BbQ450 = ensayo_12$X.Muertas.4,
                           tiempoBbQ600= ensayo_12$tiempo_dias.5,
                           BbQ600 = ensayo_12$X.Muertas.5,
                           tiempoBbQ750= ensayo_12$tiempo_dias.6,
                           BbQ750 = ensayo_12$X.Muertas.6)


# Añadir una columna de nido y combinar los data frames
df_ensayo_10$nido <- "RECS"
df_ensayo_11$nido <- "Bernal"
df_ensayo_12$nido <- "Gonnet"

# Combinar los data frames
df_combined <- bind_rows(df_ensayo_10,df_ensayo_11,df_ensayo_12)

# Organizar el data frame en formato largo
df_largo <- df_combined %>%
  pivot_longer(cols = c(tiempoControlSac, tiempoControlPGBQ, tiempoBbQ150, tiempoBbQ300, tiempoBbQ450, tiempoBbQ600, tiempoBbQ750),
               names_to = "tiempo_tratamiento",
               values_to = "tiempo") %>%
  pivot_longer(cols = c(controlSac, controlPGBQ, BbQ150, BbQ300, BbQ450, BbQ600, BbQ750),
               names_to = "tipo_tratamiento",
               values_to = "status") %>%
  mutate(status = as.numeric(status > 0)) # Define status binario: 1 para evento (muerte), 0 para censura


# Crear una nueva columna que clasifique los tratamientos y controles
df_largo <- df_largo %>%
  mutate(tipo_tratamiento = case_when(
    grepl("controlSac", tipo_tratamiento) ~ "ControlSac",
    grepl("controlPGBQ", tipo_tratamiento) ~ "ControlPGBQ",
    grepl("BbQ150", tipo_tratamiento) ~ "BbQ150",
    grepl("BbQ300", tipo_tratamiento) ~ "BbQ300",
    grepl("BbQ450", tipo_tratamiento) ~ "BbQ450",
    grepl("BbQ600", tipo_tratamiento) ~ "BbQ600",
    grepl("BbQ750", tipo_tratamiento) ~ "BbQ750",
    TRUE ~ "Desconocido"
  ))

# Crear una nueva columna que identifique cada combinación de nido y tratamiento/control
df_largo <- df_largo %>%
  mutate(grupo = paste(nido, tipo_tratamiento, sep = "_"))

# Cargar las librerías necesarias
library(dplyr)

#Prueba Kruskal-Wallis

df_largo_evento_1 <- df_largo %>%
  filter(status == 1)

# Realizar la prueba de Kruskal-Wallis
kruskal_resultado <- kruskal.test(tiempo ~ nido, data = df_largo_evento_1)

# Mostrar resultados
print(kruskal_resultado)

# Filtrar por cada tratamiento y realizar la prueba de Kruskal-Wallis
tratamientos <- c("BbQ150", "BbQ300", "BbQ450", "BbQ600", "BbQ750")

resultados <- lapply(tratamientos, function(trat) {
  df_filtrado <- df_largo %>% filter(tipo_tratamiento == trat)
  prueba <- kruskal.test(status ~ nido, data = df_filtrado)
  return(list(tratamiento = trat, resultado = prueba))
})

# Mostrar los resultados
for (res in resultados) {
  cat("\nTratamiento:", res$tratamiento, "\n")
  print(res$resultado)
}

# Realizar comparaciones de Wilcoxon entre pares de nidos para cada tratamiento

# Definir tratamientos y pares de nidos a comparar
tratamientos <- c("BbQ150", "BbQ300", "BbQ450", "BbQ600", "BbQ750")
pares_nidos <- list(
  c("Gonnet", "Bernal"),
  c("Bernal", "RECS"),
  c("RECS", "Gonnet")
)

# Lista para almacenar los resultados
resultados_kruskal <- list()

for (trat in tratamientos) {
  for (par in pares_nidos) {
    df_filtrado <- df_largo %>% filter(tipo_tratamiento == trat & nido %in% par)  # Filtrar por tratamiento y par de nidos

    if (nrow(df_filtrado) > 1) {  # Verificar que haya datos suficientes
      prueba_kruskal <- kruskal.test(status ~ nido, data = df_filtrado)  # Prueba Kruskal-Wallis

      # Guardar resultado en la lista
      resultados_kruskal[[paste(trat, paste(par, collapse = " vs. "), sep = " - ")]] <- prueba_kruskal
    }
  }
}

# Mostrar los resultados
for (nombre_comparacion in names(resultados_kruskal)) {
  cat("\nComparación:", nombre_comparacion, "\n")
  print(resultados_kruskal[[nombre_comparacion]])
}
