#Librerias
library(ggplot2)
library(dplyr)
library(tidyr)
library(extrafont)  # Para cargar fuentes adicionales en Linux
library(ggtext)  # Para etiquetas enriquecidas
library(tidyverse)
library(gridExtra)

#Carga de Datos
datos_crecimiento_fungico <- read.csv(file = "Ensayos de Micotoxinas de B. bassiana en L. humile - Crecimiento fúngico.csv", header = T, sep=",", skip=1)
df <- read.csv("Ensayos de Micotoxinas de B. bassiana en L. humile - Crecimiento fúngico.csv")

df_crecimiento_fungico <- data.frame(ensayo = datos_crecimiento_fungico$Ensayo,
                           tratamiento = datos_crecimiento_fungico$Tratamiento,
                           dia = datos_crecimiento_fungico$Dia,
                           sin_muerte = datos_crecimiento_fungico$Evento_Sin_muerte,
                           muerte_sin_hongo= datos_crecimiento_fungico$Evento_Muerte_sin_hongo,
                           muerte_con_Bb = datos_crecimiento_fungico$Evento_Muerte_con_Bb,
                           muerte_con_hongo = datos_crecimiento_fungico$Evento_Muerte_con_hongo,
                           muertes = datos_crecimiento_fungico$Muertes,
                           supervivientes = datos_crecimiento_fungico$Supervivientes)

# Contar el total de cada evento en todos los ensayos y tratamientos
df_totales <- df %>%
  select(Evento_Muerte_sin_hongo, Evento_Muerte_con_Bb, Evento_Muerte_con_hongo) %>%
  pivot_longer(cols = everything(), names_to = "Evento", values_to = "Ocurrencia") %>%
  filter(Ocurrencia == 1) %>%
  count(Evento)

# Etiquetas más legibles
df_totales$Evento <- recode(df_totales$Evento,
                            "Evento_Muerte_sin_hongo" = "Muerte sin crecimiento fúngico",
                            "Evento_Muerte_con_Bb" = "Crecimiento de B. bassiana",
                            "Evento_Muerte_con_hongo" = "Crecimiento de hongo observado")

# Graficar los totales en barras
ggplot(df_totales, aes(x = Evento, y = n, fill = Evento)) +
  geom_bar(stat = "identity") +
  labs(title = "Frecuencia total de muertes con o sin crecimiento fúngico",
       x = "Tipo de crecimiento fúngico",
       y = "Frecuencia total") +
  theme_bw() +
  theme(
    text = element_text(family = "DejaVu Sans"),  # Cambiar la fuente
    plot.title = element_text(size = 36, hjust = 0.5, face = "bold"),  # Tamaño del título
    axis.title.x = element_text(size = 28, face = "bold"),  # Tamaño del título del eje X
    axis.title.y = element_text(size = 28, face = "bold"),  # Tamaño del título del eje Y
    axis.text.x = element_text(size = 26, hjust = 0.5, face = "bold"),  # Tamaño del texto del eje X
    axis.text.y = element_text(size = 24),  # Tamaño del texto del eje Y
    legend.position = "none"  # Ocultar la leyenda ya que los nombres están en el eje X
  )


# Prueba de proporciones binomial (prueba unilateral: menor de lo esperado)
prop.test(7, 392, p = 293/392, alternative = "less")


#Modelo Bayesiano para predicir probaiblidad de que B. bassiana sea de micoflora natural
# Parámetros a priori
alpha <- 1  # Puedes cambiarlo según la información previa
beta <- 1

# Datos observados
X <- 7      # Casos de B. bassiana
n <- 392    # Total de cámaras húmedas

# Posterior Beta actualizada
alpha_post <- alpha + X
beta_post <- beta + (n - X)

# Muestreo de la distribución posterior
posterior_samples <- rbeta(10000, alpha_post, beta_post)

# Estimación de la media y del intervalo de credibilidad del 95%
mean(posterior_samples)  # Estimación bayesiana de la probabilidad
quantile(posterior_samples, c(0.025, 0.975))  # Intervalo de credibilidad

# Muestreo para otro hongo
X_otro <- 293
alpha_post_otro <- alpha + X_otro
beta_post_otro <- beta + (n - X_otro)
posterior_samples_otro <- rbeta(10000, alpha_post_otro, beta_post_otro)

# Calcular Odds Ratio
odds_ratio_samples <- posterior_samples / posterior_samples_otro

# Estimar intervalo de credibilidad para la odds ratio
quantile(odds_ratio_samples, c(0.025, 0.975))

# Datos observados
X_bb <- 7      # Casos de B. bassiana
X_otro <- 7584  # Casos de otro hongo
n <- 392       # Total de cámaras húmedas

# Definir diferentes parámetros a priori para Beta
priors <- list(
  "Beta(1,100)" = c(1, 100),
  "Beta(2,50)" = c(2, 50),
  "Beta(5,50)" = c(5, 50),
  "Beta(10,40)" = c(10, 40)
)

# Almacenar datos de densidad para graficar todas juntas
df_plot <- data.frame()
annotations <- data.frame()

for (prior in names(priors)) {
  alpha <- priors[[prior]][1]
  beta <- priors[[prior]][2]

  # Posterior para B. bassiana
  alpha_post_bb <- alpha + X_bb
  beta_post_bb <- beta + (n - X_bb)
  posterior_samples_bb <- rbeta(10000, alpha_post_bb, beta_post_bb)
  mean_bb <- mean(posterior_samples_bb)
  ci_bb <- quantile(posterior_samples_bb, c(0.025, 0.975))

  # Odds Ratio con respecto a otro hongo
  alpha_post_otro <- alpha + X_otro
  beta_post_otro <- beta + (n - X_otro)
  posterior_samples_otro <- rbeta(10000, alpha_post_otro, beta_post_otro)
  odds_ratio_samples <- posterior_samples_bb / posterior_samples_otro
  ci_or <- quantile(odds_ratio_samples, c(0.025, 0.975))

  # Obtener datos de densidad
  dens_data <- density(posterior_samples_bb)
  df_temp <- data.frame(Probabilidad = dens_data$x, Densidad = dens_data$y, Prior = prior)
  df_plot <- rbind(df_plot, df_temp)

  # Guardar datos para anotaciones
  annotations <- rbind(annotations, data.frame(
    Prior = prior,
    Media = mean_bb,
    IC_Low = ci_bb[1],
    IC_High = ci_bb[2],
    OR_Low = ci_or[1],
    OR_High = ci_or[2],
    Max_Y = max(dens_data$y)  # Para posicionar anotaciones en la parte superior
  ))
}

# Asignar colores
colors <- c("Beta(1,100)" = "blue", "Beta(2,50)" = "red", "Beta(5,50)" = "green", "Beta(10,40)" = "purple")
line_types <- c("Beta(1,100)" = "solid", "Beta(2,50)" = "dashed", "Beta(5,50)" = "dotted", "Beta(10,40)" = "dotdash")

# Crear gráfico superpuesto
p <- ggplot(df_plot, aes(x = Probabilidad, y = Densidad, color = Prior)) +
  geom_line(size = 1.6) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = line_types) +
  theme_minimal() +
  labs(title = "Distribuciones Beta Posteriores con Diferentes Priors",
       x = "Probabilidad de que B. bassiana esté en la micoflora", y = "Densidad") +
  theme(
    text = element_text(family = "DejaVu Sans"),  # Cambiar la fuente
    plot.title = element_text(size = 36, hjust = 0.5, face = "bold"),  # Tamaño del título
    axis.title.x = element_text(size = 28, face = "bold"),  # Tamaño del título del eje X
    axis.title.y = element_text(size = 28, face = "bold"),  # Tamaño del título del eje Y
    axis.text.x = element_text(size = 26, hjust = 0.5, face = "bold"),  # Tamaño del texto del eje X
    axis.text.y = element_text(size = 24),  # Tamaño del texto del eje Y
    legend.title = element_blank(),
    legend.text = element_text(size = 24, face = "bold")
  )

# Agregar anotaciones centradas respecto a la media de cada distribución
for (i in 1:nrow(annotations)) {
  p <- p +
    annotate("text", x = annotations$Media[i], y = annotations$Max_Y[i] * 1.05,
             label = sprintf("Probabilidad Media: %.4f\nIC 95%%: [%.4f, %.4f]\nOR 95%%: [%.4f, %.4f]",
                             annotations$Media[i], annotations$IC_Low[i], annotations$IC_High[i],
                             annotations$OR_Low[i], annotations$OR_High[i]),
             color = colors[annotations$Prior[i]], size = 7.25, hjust = 0.00000000001)
}

# Mostrar gráfico
print(p)

#Simulacion de obtener un cadaver de L. humile infectado con B.bassiana
set.seed(123)

# Parámetros
n_hormigas <- 2553  # Número total de hormigas analizadas = 7591, Solo de Bernal = 2553
observado_bb <- 6   # Casos reales de B. bassiana = 7, solo 6 en Bernal

# Definir diferentes probabilidades reales posibles de infección
p_vals <- c(0.0001, 0.0005, 0.001, 0.0025, 0.01)  # Probabilidades desde 0.01% hasta 1%
n_simulaciones <- 10000  # Número de simulaciones por cada p

# DataFrame para almacenar resultados
sim_results <- data.frame()

# Simular para cada valor de p
for (p in p_vals) {
  simulaciones <- rbinom(n_simulaciones, size = n_hormigas, prob = p)  # Simulación binomial
  df_temp <- data.frame(Probabilidad = as.factor(p), Aislamientos = simulaciones)
  sim_results <- rbind(sim_results, df_temp)
}

# Crear gráfico
ggplot(sim_results, aes(x = Aislamientos, fill = Probabilidad)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.6) +
  geom_vline(xintercept = observado_bb, linetype = "dashed", color = "black", size = 1.2) +
  labs(title = "Simulación de Aislamientos de B. bassiana",
       x = "Número de aislamientos simulados",
       y = "Frecuencia",
       fill = "Probabilidad\n asumida\n de infección") +
  theme_minimal() +
  theme(
    text = element_text(family = "DejaVu Sans"),  # Cambiar la fuente
    plot.title = element_text(size = 36, hjust = 0.5, face = "bold"),  # Tamaño del título
    axis.title.x = element_text(size = 28, face = "bold"),  # Tamaño del título del eje X
    axis.title.y = element_text(size = 28, face = "bold"),  # Tamaño del título del eje Y
    axis.text.x = element_text(size = 26, hjust = 0.5, face = "bold"),  # Tamaño del texto del eje X
    axis.text.y = element_text(size = 24),  # Tamaño del texto del eje Y
    legend.title = element_text(size = 24, face = "bold", hjust= 0.5),
    legend.text = element_text(size = 20)
  )

# Parámetros generales
n_hormigas <- 392
observado_bb <- 7  # Casos reales detectados
n_simulaciones <- 10000  # Número de simulaciones

# Definir distribuciones Beta
beta_params <- list(
  "Beta(1,100)" = c(1, 100),
  "Beta(2,50)"  = c(2, 50),
  "Beta(5,50)"  = c(5, 50),
  "Beta(10,40)" = c(10, 40)
)

# DataFrame para almacenar resultados
sim_results <- data.frame()

# Simular para cada distribución beta
for (dist in names(beta_params)) {
  alpha <- beta_params[[dist]][1]
  beta <- beta_params[[dist]][2]

  # Muestreo de probabilidad real desde Beta(alpha, beta)
  prob_reales <- rbeta(n_simulaciones, alpha, beta)

  # Simulación binomial usando esas probabilidades
  aislamientos <- rbinom(n_simulaciones, size = n_hormigas, prob = prob_reales)

  # Guardar resultados
  df_temp <- data.frame(Distribucion = dist, Aislamientos = aislamientos)
  sim_results <- rbind(sim_results, df_temp)
}

# Crear gráfico de distribuciones de aislamientos
ggplot(sim_results, aes(x = Aislamientos, fill = Distribucion)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  geom_vline(xintercept = observado_bb, linetype = "dashed", color = "black", size = 1.2) +
  labs(title = "Simulación de Aislamientos de B. bassiana con Probabilidades Beta",
       x = "Número de aislamientos simulados",
       y = "Frecuencia",
       fill = "Distribución Beta") +
  theme_minimal()

# Filtrar ensayos desde el 8 en adelante y considerar solo eventos con crecimiento fúngico
df_analisis_por_tratamiento <- df

# Transformar nombres en la columna Tratamiento
df_analisis_por_tratamiento$Tratamiento <- gsub("Control Sacarosa 25%m/V - Bernal", "Control Sacarosa 25%m/V", df_analisis_por_tratamiento$Tratamiento)
df_analisis_por_tratamiento$Tratamiento <- gsub("Control Sacarosa 25%m/V - Gonnet", "Control Sacarosa 25%m/V", df_analisis_por_tratamiento$Tratamiento)
df_analisis_por_tratamiento$Tratamiento <- gsub("Bb 100%", "Bb 150 μl", df_analisis_por_tratamiento$Tratamiento)
df_analisis_por_tratamiento$Tratamiento <- gsub("BbQ 100%", "BbQ 150 μl", df_analisis_por_tratamiento$Tratamiento)
df_analisis_por_tratamiento$Tratamiento <- gsub("BbQ 450 μL - Bernal", "BbQ 450 μl", df_analisis_por_tratamiento$Tratamiento)
df_analisis_por_tratamiento$Tratamiento <- gsub("BbQ 450 μL- Gonnet", "BbQ 450 μl", df_analisis_por_tratamiento$Tratamiento)
df_analisis_por_tratamiento$Tratamiento <- gsub("Control Sacarosa 25%m/V","Control Sacarosa", df_analisis_por_tratamiento$Tratamiento)

df_analisis_por_tratamiento_filtrado <- df_analisis_por_tratamiento %>%
  filter(Evento_Muerte_con_hongo == 1) %>%
  group_by(Tratamiento) %>%
  summarise(Apariciones_hongo = sum(Evento_Muerte_con_hongo)) %>%
  arrange(desc(Apariciones_hongo))  # Ordenar de mayor a menor

# Gráfico de barras
ggplot(df_analisis_por_tratamiento_filtrado, aes(x = reorder(Tratamiento, Apariciones_hongo), y = Apariciones_hongo, fill = Tratamiento)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8) +  # Colores con Viridis
  labs(title = "Cantidad de Apariciones de Hongos tras la Muerte por Tratamiento",
       x = "Tratamiento",
       y = "Número de Apariciones") +
  theme_minimal() +
  theme(
    text = element_text(family = "Liberation Sans"),  # Cambiar la fuente
    plot.title = element_text(size = 36, hjust = 0.5, face = "bold"),  # Tamaño del título
    axis.title.x = element_text(size = 28, face = "bold"),  # Tamaño del título del eje X
    axis.title.y = element_text(size = 28, face = "bold"),  # Tamaño del título del eje Y
    axis.text.x = element_text(size = 18, angle = 90, hjust = 0.5, face = "bold"),  # Tamaño del texto del eje X
    axis.text.y = element_text(size = 24),  # Tamaño del texto del eje Y
    legend.position = "none"  # Ocultar la leyenda ya que los nombres están en el eje X
  )

df_analisis_por_tratamiento_2 <- df_analisis_por_tratamiento

df_analisis_por_tratamiento_2$Tratamiento <- gsub("Control Sacarosa 25%m/V","Control Sacarosa", df_analisis_por_tratamiento_2$Tratamiento)
df_analisis_por_tratamiento_2$Tratamiento <- gsub("Bb .*", "Bb", df_analisis_por_tratamiento_2$Tratamiento)
df_analisis_por_tratamiento_2$Tratamiento <- gsub("BbQ.*", "BbQ", df_analisis_por_tratamiento_2$Tratamiento)

df_analisis_por_tratamiento_filtrado_2 <- df_analisis_por_tratamiento_2 %>%
  filter(Evento_Muerte_con_hongo == 1) %>%
  group_by(Tratamiento) %>%
  summarise(Apariciones_hongo = sum(Evento_Muerte_con_hongo)) %>%
  arrange(desc(Apariciones_hongo))  # Ordenar de mayor a menor

# Gráfico de barras
ggplot(df_analisis_por_tratamiento_filtrado_2, aes(x = reorder(Tratamiento, Apariciones_hongo), y = Apariciones_hongo, fill = Tratamiento)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8) +  # Colores con Viridis
  labs(title = "Cantidad de Apariciones de Hongos tras la Muerte por Tratamiento",
       x = "Tratamiento",
       y = "Número de Apariciones") +
  theme_minimal() +
  theme(
    text = element_text(family = "Liberation Sans"),  # Cambiar la fuente
    plot.title = element_text(size = 36, hjust = 0.5, face = "bold"),  # Tamaño del título
    axis.title.x = element_text(size = 28, face = "bold"),  # Tamaño del título del eje X
    axis.title.y = element_text(size = 28, face = "bold"),  # Tamaño del título del eje Y
    axis.text.x = element_text(size = 20, angle = 90, hjust = 0.5, face = "bold"),  # Tamaño del texto del eje X
    axis.text.y = element_text(size = 24),  # Tamaño del texto del eje Y
    legend.position = "none"  # Ocultar la leyenda ya que los nombres están en el eje X
  )

#Contar el número total de tratamientos realizados para cada tipo (Bb o BbQ)
conteo_tratamientos <- df_analisis_por_tratamiento_2 %>%
  group_by(Tratamiento) %>%
  summarise(Total_Tratamientos = n())

# Unir los datos de apariciones de hongo con el total de tratamientos realizados
df_final <- left_join(df_analisis_por_tratamiento_filtrado_2, conteo_tratamientos, by = "Tratamiento")

# Crear gráfico de barras
ggplot(df_final, aes(x = Tratamiento, y = Apariciones_hongo, fill = Tratamiento)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = paste0(Apariciones_hongo, " (", Total_Tratamientos, ")")),
            vjust = -0.5, size = 8, fontface = "bold") +  # Etiquetas con apariciones + total de tratamientos
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8) +
  labs(title = "Cantidad de Apariciones de Hongos tras la Muerte por Tratamiento",
       x = "Tratamiento",
       y = "Número de Apariciones",
       caption = "Entre paréntesis, el número total de tratamientos realizados") +
  theme_minimal() +
  theme(
    text = element_text(family = "Liberation Sans"),
    plot.title = element_text(size = 36, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 28, face = "bold"),
    axis.title.y = element_text(size = 28, face = "bold"),
    axis.text.x = element_text(size = 24, face = "bold"),
    axis.text.y = element_text(size = 24),
    legend.position = "none"
  )

library(broom)
library(tibble)
# Crear tabla de contingencia: Tratamiento vs. Evento de crecimiento fúngico
# Filtrar solo tratamientos "Bb" y "BbQ"
df_Bb_vs_BbQ <- df_analisis_por_tratamiento_2 %>%
  filter(Tratamiento %in% c("Bb", "BbQ"))

# Crear tabla de contingencia solo con Bb y BbQ
tabla_Bb_BbQ <- df_Bb_vs_BbQ %>%
  group_by(Tratamiento) %>%
  summarise(Apariciones_hongo = sum(Evento_Muerte_con_hongo),
            Total_Tratamientos = n()) %>%
  mutate(No_Apariciones = Total_Tratamientos - Apariciones_hongo) %>%
  select(Tratamiento, Apariciones_hongo, No_Apariciones) %>%  # Incluir Tratamiento
  column_to_rownames("Tratamiento") %>%
  as.matrix()

# 1️⃣ Prueba de proporciones para Bb vs BbQ
prop_result_Bb_BbQ <- prop.test(tabla_Bb_BbQ[,1], rowSums(tabla_Bb_BbQ))
print(prop_result_Bb_BbQ)

# 2️⃣ Prueba de chi-cuadrado para Bb vs BbQ
chisq_result_Bb_BbQ <- chisq.test(tabla_Bb_BbQ)
print(chisq_result_Bb_BbQ)


# 3️⃣ Regresión logística: Modelamos la probabilidad de crecimiento fúngico según el tratamiento
df_Bb_vs_BbQ <- df_Bb_vs_BbQ %>%
  mutate(Tratamiento_binario = ifelse(Tratamiento == "BbQ", 1, 0))

modelo_BbQ_vs_Bb <- glm(Evento_Muerte_con_hongo ~ Tratamiento_binario,
                         data = df_Bb_vs_BbQ,
                         family = binomial)
summary(modelo_BbQ_vs_Bb)

# Filtrar solo tratamientos "BbQ" y "Control Sacarosa"
df_BbQ_vs_Control <- df_analisis_por_tratamiento_2 %>%
  filter(Tratamiento %in% c("BbQ", "Control Sacarosa"))

# Crear tabla de contingencia para BbQ vs Control Sacarosa
tabla_BbQ_Control <- df_BbQ_vs_Control %>%
  group_by(Tratamiento) %>%
  summarise(Apariciones_hongo = sum(Evento_Muerte_con_hongo),
            Total_Tratamientos = n()) %>%
  mutate(No_Apariciones = Total_Tratamientos - Apariciones_hongo) %>%
  select(Tratamiento, Apariciones_hongo, No_Apariciones) %>%
  column_to_rownames("Tratamiento") %>%
  as.matrix()

# 1️⃣ Prueba de proporciones para BbQ vs Control Sacarosa
prop_result_BbQ_Control <- prop.test(tabla_BbQ_Control[,1], rowSums(tabla_BbQ_Control))
print(prop_result_BbQ_Control)

# 2️⃣ Prueba de chi-cuadrado para BbQ vs Control Sacarosa
chisq_result_BbQ_Control <- chisq.test(tabla_BbQ_Control)
print(chisq_result_BbQ_Control)


# 3️⃣ Regresión logística: Modelamos la probabilidad de crecimiento fúngico según el tratamiento
df_BbQ_vs_Control <- df_BbQ_vs_Control %>%
  mutate(Tratamiento_binario = ifelse(Tratamiento == "BbQ", 1, 0))

modelo_BbQ_vs_Control <- glm(Evento_Muerte_con_hongo ~ Tratamiento_binario,
                         data = df_BbQ_vs_Control,
                         family = binomial)
summary(modelo_BbQ_vs_Control)

# Filtrar solo tratamientos "Bb" y "Control Sacarosa"
df_Bb_vs_Control <- df_analisis_por_tratamiento_2 %>%
  filter(Tratamiento %in% c("Bb", "Control Sacarosa"))

# Crear tabla de contingencia para Bb vs Control Sacarosa
tabla_Bb_Control <- df_Bb_vs_Control %>%
  group_by(Tratamiento) %>%
  summarise(Apariciones_hongo = sum(Evento_Muerte_con_hongo),
            Total_Tratamientos = n()) %>%
  mutate(No_Apariciones = Total_Tratamientos - Apariciones_hongo) %>%
  select(Tratamiento, Apariciones_hongo, No_Apariciones) %>%
  column_to_rownames("Tratamiento") %>%
  as.matrix()

# 1️⃣ Prueba de proporciones para BbQ vs Control Sacarosa
prop_result_Bb_Control <- prop.test(tabla_Bb_Control[,1], rowSums(tabla_Bb_Control))
print(prop_result_Bb_Control)

# 2️⃣ Prueba de chi-cuadrado para BbQ vs Control Sacarosa
chisq_result_Bb_Control <- chisq.test(tabla_Bb_Control)
print(chisq_result_Bb_Control)

# 3️⃣ Regresión logística: Modelamos la probabilidad de crecimiento fúngico según el tratamiento
df_Bb_vs_Control <- df_Bb_vs_Control %>%
  mutate(Tratamiento_binario = ifelse(Tratamiento == "Bb", 1, 0))

modelo_Bb_vs_Control <- glm(Evento_Muerte_con_hongo ~ Tratamiento_binario,
                         data = df_Bb_vs_Control,
                         family = binomial)
summary(modelo_Bb_vs_Control)
