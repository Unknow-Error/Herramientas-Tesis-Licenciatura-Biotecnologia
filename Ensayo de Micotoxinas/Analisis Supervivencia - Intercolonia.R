#Análisis de Efectos de supervivencia Mixtos

#Efecto mixto - Cox
library(coxme)
library(dplyr)  # Para combinar los data frames de forma eficiente
library(tidyr)
library(survival)
library(bdsmatrix)

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

# Ajustar el modelo de Cox con efectos aleatorios
modelo_coxme <- coxme(Surv(tiempo, status) ~ tipo_tratamiento + (1 | nido), data = df_largo)
print(summary(modelo_coxme))

#Kaplan-Meier Estratificado

library(survival)
library(survminer)

surv_obj <- Surv(time = df_largo$tiempo, event = df_largo$status)
km_fit <- survfit(surv_obj ~ grupo, data = df_largo) #Si dividir en tratamiento y nido
km_fit <- survfit(surv_obj ~ nido, data = df_largo)
ggsurvplot(km_fit, data = df_largo,
           conf.int = TRUE,             # Intervalo de confianza
           pval = TRUE,                 # Mostrar el p-valor
           risk.table = TRUE,           # Mostrar la tabla de riesgo
           legend.title = "Nidos y Tratamientos",  # Título de la leyenda
           xlab = "Días de Supervivencia",  # Etiqueta del eje x
           ylab = "Probabilidad de Supervivencia")  # Etiqueta del eje y

# Prueba de log-rank
log_rank_test <- survdiff(surv_obj ~ nido, data = df_largo)
print(log_rank_test)

#Prueba Kruskal-Wallis
library(FSA)  # Para pruebas post hoc
df_largo_evento_1 <- df_largo %>%
  filter(status == 1)

# Realizar la prueba de Kruskal-Wallis
kruskal_resultado <- kruskal.test(tiempo ~ nido, data = df_largo_evento_1)

# Mostrar resultados
print(kruskal_resultado)

# Realizar prueba de Dunn
dunn_resultados <- dunnTest(tiempo ~ nido, data = df_largo_evento_1, method = "bh")

# Mostrar resultados de Dunn
print(dunn_resultados)


# Análisis de variancias de Milton Friedman

# Paso 1: Convertir cada data frame en formato largo y agregar el nido
df_ensayo_10_long <- df_ensayo_10 %>%
  pivot_longer(cols = starts_with("tiempo"), names_to = "tratamiento", values_to = "tiempo") %>%
  pivot_longer(cols = starts_with("control") | starts_with("BbQ"), names_to = "evento", values_to = "muertes") %>%
  filter(muertes > 0) %>%
  uncount(muertes) %>%
  mutate(nido = "RECS")

df_ensayo_11_long <- df_ensayo_11 %>%
  pivot_longer(cols = starts_with("tiempo"), names_to = "tratamiento", values_to = "tiempo") %>%
  pivot_longer(cols = starts_with("control") | starts_with("BbQ"), names_to = "evento", values_to = "muertes") %>%
  filter(muertes > 0) %>%
  uncount(muertes) %>%
  mutate(nido = "Bernal")

df_ensayo_12_long <- df_ensayo_12 %>%
  pivot_longer(cols = starts_with("tiempo"), names_to = "tratamiento", values_to = "tiempo") %>%
  pivot_longer(cols = starts_with("control") | starts_with("BbQ"), names_to = "evento", values_to = "muertes") %>%
  filter(muertes > 0) %>%
  uncount(muertes) %>%
  mutate(nido = "Gonnet")

# Paso 2: Combinar los data frames
df_combined <- bind_rows(df_ensayo_10_long, df_ensayo_11_long, df_ensayo_12_long)

# Crear una variable que combine nido y tratamiento
df_combined <- df_combined %>%
  mutate(tratamiento = case_when(
    grepl("ControlSac", tratamiento) ~ "ControlSac",
    grepl("ControlPGBQ", tratamiento) ~ "ControlPGBQ",
    grepl("BbQ150", tratamiento) ~ "BbQ150",
    grepl("BbQ300", tratamiento) ~ "BbQ300",
    grepl("BbQ450", tratamiento) ~ "BbQ450",
    grepl("BbQ600", tratamiento) ~ "BbQ600",
    grepl("BbQ750", tratamiento) ~ "BbQ750",
    TRUE ~ tratamiento
  )) # Sumar las muertes para cada combinación

# Paso 3: Asegurarnos de que cada fila tiene una observación por tratamiento y nido
df_friedman <- df_combined %>%
  group_by(nido, tratamiento) %>%
  summarise(muertes = n()) %>%  # Contar las muertes para cada combinación
  pivot_wider(names_from = tratamiento, values_from = muertes, values_fill = 0)


# Realizar la prueba de Friedman
friedman_test <- friedman.test(as.matrix(df_friedman[, -1]))  # Omitir la columna de 'nido'

# Mostrar el resultado
print(friedman_test)
