#Rutina para realizar modelado de curvas de Supervivencia por método de Kaplan-Meier y su respectivo gráfico

#Librerias
library(survival)
library(survMisc)
library(KMsurv)
library(survminer)
library(ggfortify)
library(flexsurv)
library(actuar)

library(survival)
library(ggplot2)
library(dplyr)

#Carga de Datos
ensayo_14 <- read.csv(file = "Ensayos de Micotoxinas de B. bassiana en L. humile - Ensayo #14 - Para analisis R.csv", header = T, sep=",", skip=1)


df_ensayo_14 <- data.frame(tiempoControlSac = ensayo_14$tiempo_dias,
                           controlSac = ensayo_14$X.Muertas,
                           tiempoControlPGB = ensayo_14$tiempo_dias.1,
                           controlPGB = ensayo_14$X.Muertas.1,
                           tiempoBbQ450= ensayo_14$tiempo_dias.2,
                           BbQ450 = ensayo_14$X.Muertas.2,
                           tiempoBb150= ensayo_14$tiempo_dias.3,
                           Bb150 = ensayo_14$X.Muertas.3,
                           tiempoBb300= ensayo_14$tiempo_dias.4,
                           Bb300 = ensayo_14$X.Muertas.4,
                           tiempoBb450= ensayo_14$tiempo_dias.5,
                           Bb450 = ensayo_14$X.Muertas.5)

library(ggpubr)

# Prueba de Shapiro-Wilk para cada grupo de tratamiento
shapiro_test_controlSac <- shapiro.test(df_ensayo_14$controlSac)
shapiro_test_controlPGB <- shapiro.test(df_ensayo_14$controlPGB)
shapiro_test_Bb150 <- shapiro.test(df_ensayo_14$Bb150)
shapiro_test_Bb300 <- shapiro.test(df_ensayo_14$Bb300)
shapiro_test_Bb450 <- shapiro.test(df_ensayo_14$Bb450)
shapiro_test_BbQ450 <- shapiro.test(df_ensayo_14$BbQ450)


# Imprimir los resultados
print(shapiro_test_controlSac)
print(shapiro_test_controlPGB)
print(shapiro_test_Bb150)
print(shapiro_test_Bb300)
print(shapiro_test_Bb450)
print(shapiro_test_BbQ450)

# Prueba de Levene

library(car)

levene_test <- leveneTest(
  c(controlSac = df_ensayo_14$controlSac,
    controlPGB = df_ensayo_14$controlPGB,
    BbQ450 = df_ensayo_14$BbQ450,
    Bb150 = df_ensayo_14$Bb150,
    Bb300 = df_ensayo_14$Bb300,
    Bb450 = df_ensayo_14$Bb450) ~ 
    factor(rep(c("controlSac", "controlPGB","Bb150", "Bb300", "Bb450", "BbQ450"),
               each = nrow(df_ensayo_14))),
  center = median
)

levene_test_2 <- leveneTest(
  c(Bb150 = df_ensayo_14$Bb150,
    Bb300 = df_ensayo_14$Bb300,
    Bb450 = df_ensayo_14$Bb450) ~ 
    factor(rep(c("Bb150", "Bb300", "Bb450"),
               each = nrow(df_ensayo_14))),
  center = median
)

levene_test_3 <- leveneTest(
  c(controlSac = df_ensayo_14$controlSac,
    controlPGB = df_ensayo_14$controlPGB) ~ 
    factor(rep(c("controlSac", "controlPGB"),
               each = nrow(df_ensayo_14))),
  center = median
)

levene_test_4 <- leveneTest(
  c(Bb450 = df_ensayo_14$Bb300,
    BbQ450 = df_ensayo_14$BbQ450) ~ 
    factor(rep(c("Bb450", "BbQ450"),
               each = nrow(df_ensayo_14))),
  center = median
)

print(levene_test)
print(levene_test_2)
print(levene_test_3)
print(levene_test_4)

# Definir los modelos Kaplan-Meier
survfits <- list(
  "Control Sac" = survfit(Surv(tiempoControlSac,controlSac) ~ 1, data = df_ensayo_14, type = "kaplan-meier", conf.type = "log-log"),
  "Control PGB" = survfit(Surv(tiempoControlPGB,controlPGB) ~ 1, data = df_ensayo_14, type = "kaplan-meier", conf.type = "log-log"),
  "Bb 150 µL" = survfit(Surv(tiempoBb150,Bb150) ~ 1, data = df_ensayo_14, type = "kaplan-meier", conf.type = "log-log"),
  "Bb 300 µL" = survfit(Surv(tiempoBb300,Bb300) ~ 1, data = df_ensayo_14, type = "kaplan-meier", conf.type = "log-log"),
  "Bb 450 µL" = survfit(Surv(tiempoBb450,Bb450) ~ 1, data = df_ensayo_14, type = "kaplan-meier", conf.type = "log-log"),
  "BbQ 450 µL" = survfit(Surv(tiempoBbQ450,BbQ450) ~ 1, data = df_ensayo_14, type = "kaplan-meier", conf.type = "log-log")
)


# Definir colores, tipos de línea y formas para los puntos
colores <- c("coral3","coral3", "coral3", "palegreen2", "darkblue", "darkblue")
linetypes <- c("solid", "dotted", "dashed", "solid", "solid", "dotted")
shapes <- c(15, 16, 17, 15, 15, 16)

# Convertir la lista de modelos a un solo data frame
df_surv <- bind_rows(lapply(names(survfits), function(name) {
  df <- as.data.frame(surv_summary(survfits[[name]]))
  df$strata <- name  # Agregar nombre del tratamiento
  return(df)
}))

# Calcular los tiempos medianos de supervivencia
median_times <- sapply(survfits, function(survfits) {
  if (!is.null(survfits$surv) && any(survfits$surv <= 0.5)) {
    return(surv_median(survfits)$median)
  } else {
    return(NA)  # No hay tiempo mediano definido
  }
})

# Crear un dataframe con los segmentos para los tiempos medianos
df_median <- data.frame(
  strata = names(survfits),
  x1 = median_times,
  x2 = median_times,
  y1 = rep(0, length(median_times)),
  y2 = rep(0.5, length(median_times))
) %>% na.omit()  # Eliminar valores NA

# Definir las anotaciones
annotaciones <- c("A", "B")
y_positions <- c(0.89, 0.25)

# Crear el gráfico Kaplan-Meier combinado
grafico <- ggplot(df_surv, aes(x = time, y = surv, color = strata, linetype = strata, shape = strata)) +
  geom_step(size = 3) +  # Línea de supervivencia
  geom_point(size = 9.5) + # Puntos en las curvas
  scale_color_manual(values = colores) +
  scale_linetype_manual(values = linetypes) +
  scale_shape_manual(values = shapes) +
  labs(
    title = "Curva de Supervivencia - Ensayo #14",
    x = "Tiempo (días)",
    y = "Probabilidad de Supervivencia = S(t)",
    color = "Tratamiento",
    linetype = "Tratamiento",
    shape = "Tratamiento") +
  scale_x_continuous(breaks = seq(0, max(df_surv$time), by = 1)) +
  theme_minimal() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +  # Línea horizontal en y = 0.5
  geom_segment(data = df_median, aes(x = 0, y = 0.5, xend = max(df_median$x1), yend = 0.5),
               linetype = "dashed", col = "darkblue", size = 0.75) + # Segmento horizontal
  geom_segment(data = df_median, aes(x = x1, y = y1, xend = x2, yend = y2),
               linetype = "dashed", col = "darkblue", size = 0.75) + # Segmentos verticales
  annotate("text", x = max(df_surv$time) * 0.95, y = y_positions,
           label = annotaciones, hjust = 0, size = 12, color = "black", face = "bold", family = "Helvetica") +
  theme(
    plot.title = element_text(size = 40, face = "bold", family = "Helvetica", hjust = 0.5, colour = "darkblue"),  # Tamaño y negrita del título
    axis.title.x = element_text(size = 30, face = "bold", family = "Helvetica"),  # Tamaño del título del eje X
    axis.title.y = element_text(size = 30, face = "bold", family = "Helvetica"),  # Tamaño del título del eje Y
    axis.text = element_text(size = 28, family = "Helvetica"),  # Tamaño de los valores en los ejes
    legend.text = element_text(size = 30, family = "Helvetica", face = "bold"),  # Tamaño y tipo de fuente en la leyenda
    legend.title = element_text(size = 32, face = "bold", family = "Helvetica", hjust = 0.5),  # Título de la leyenda
    legend.key.size = unit(2.5, "cm") # Ajustar tamaño de los símbolos
  )

print(grafico)

#Estadísticos de cada tratamiento

print(summary(ensayo_14_controlSac.km))
print(summary(ensayo_14_controlPGB.km))
print(summary(ensayo_14_BbQ450.km))
print(summary(ensayo_14_Bb150.km))
print(summary(ensayo_14_Bb300.km))
print(summary(ensayo_14_Bb450.km))

print(summary(ensayo_14_controlSac.km)$table["median"])
print(summary(ensayo_14_controlPGB.km)$table["median"])
print(summary(ensayo_14_BbQ450.km)$table["median"])
print(summary(ensayo_14_Bb300.km)$table["median"])
print(summary(ensayo_14_Bb450.km)$table["median"])


#Test de Cox-Matlan

df_controlSac <- data.frame(tiempo = ensayo_14$tiempo_dias,
                            muertas = ensayo_14$X.Muertas,
                            tratamiento = "ControlSac")

df_controlPGB <- data.frame(tiempo = ensayo_14$tiempo_dias.1,
                            muertas = ensayo_14$X.Muertas.1,
                            tratamiento = "ControlPGBQ")

df_BbQ450 <- data.frame(tiempo = ensayo_14$tiempo_dias.2,
                       muertas = ensayo_14$X.Muertas.2,
                       tratamiento = "BbQ450")

df_Bb150 <- data.frame(tiempo = ensayo_14$tiempo_dias.3,
                       muertas = ensayo_14$X.Muertas.3,
                       tratamiento = "Bb150")

df_Bb300 <- data.frame(tiempo = ensayo_14$tiempo_dias.4,
                       muertas = ensayo_14$X.Muertas.4,
                       tratamiento = "Bb300")

df_Bb450 <- data.frame(tiempo = ensayo_14$tiempo_dias.5,
                       muertas = ensayo_14$X.Muertas.5,
                       tratamiento = "Bb450")



df_combinado_todos <- rbind(df_controlPGB, df_controlSac, df_BbQ450, df_Bb150, df_Bb300, df_Bb450)
df_combinado_1 <- rbind(df_controlSac, df_controlPGB)
df_combinado_2 <- rbind(df_controlPGB, df_BbQ450)
df_combinado_3 <- rbind(df_controlPGB, df_Bb150)
df_combinado_4 <- rbind(df_Bb150, df_Bb300)
df_combinado_5 <- rbind(df_Bb300, df_Bb450)
df_combinado_6 <- rbind(df_Bb150, df_Bb300, df_BbQ450, df_Bb450)

print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_todos))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_1))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_2))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_3))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_4))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_5))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_6))
