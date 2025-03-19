#Rutina para realizar modelado de curvas de Supervivencia por método de Kaplan-Meier y su respectivo gráfico

#Librerias
library(survival)
library(survMisc)
library(KMsurv)
library(survminer)
library(ggfortify)
library(flexsurv)
library(actuar)
library(ggsurvfit)


#Carga de Datos
ensayo_6 <- read.csv(file = "Ensayos de Micotoxinas de B. bassiana en L. humile - Ensayo #6 - Para analisis R.csv", header = T, sep=",", skip=1)
#head(ensayo_6) #Para verificar

df_ensayo_6 <- data.frame(tiempoControl= ensayo_6$tiempo_dias,
                          control = ensayo_6$X.Muertas,
                          tiempoBbQ300= ensayo_6$tiempo_dias.1,
                          BbQ300 = ensayo_6$X.Muertas.1,
                          tiempoBbQ450= ensayo_6$tiempo_dias.2,
                          BbQ450 = ensayo_6$X.Muertas.2,
                          tiempoBbQ600= ensayo_6$tiempo_dias.3,
                          BbQ600 = ensayo_6$X.Muertas.3,
                          tiempoBbQ750= ensayo_6$tiempo_dias.4,
                          BbQ750 = ensayo_6$X.Muertas.4)


library(ggpubr)

# Prueba de Shapiro-Wilk para cada grupo de tratamiento
shapiro_test_controlSac <- shapiro.test(df_ensayo_6$control)
shapiro_test_BbQ300 <- shapiro.test(df_ensayo_6$BbQ300)
shapiro_test_BbQ450 <- shapiro.test(df_ensayo_6$BbQ450)
shapiro_test_BbQ600 <- shapiro.test(df_ensayo_6$BbQ600)
shapiro_test_BbQ750 <- shapiro.test(df_ensayo_6$BbQ750)

# Imprimir los resultados
print(shapiro_test_controlSac)
print(shapiro_test_BbQ300)
print(shapiro_test_BbQ450)
print(shapiro_test_BbQ600)
print(shapiro_test_BbQ750)

# Prueba de Levene

library(car)

levene_test <- leveneTest(
  c(controlSac = df_ensayo_6$control,
    BbQ300 = df_ensayo_6$BbQ300,
    BbQ450 = df_ensayo_6$BbQ450,
    BbQ600 = df_ensayo_6$BbQ600,
    BbQ750 = df_ensayo_6$BbQ750) ~ 
    factor(rep(c("controlSac", "BbQ300", "BbQ450", "BbQ600", "BbQ750"),
               each = nrow(df_ensayo_6))),
  # center = median
)

levene_test_2 <- leveneTest(
  c(BbQ300 = df_ensayo_6$BbQ300,
    BbQ450 = df_ensayo_6$BbQ450,
    BbQ600 = df_ensayo_6$BbQ600,
    BbQ750 = df_ensayo_6$BbQ750) ~ 
    factor(rep(c("BbQ300", "BbQ450", "BbQ600", "BbQ750"),
               each = nrow(df_ensayo_6))),
  center = median
)

print(levene_test)
print(levene_test_2)


# Definir los modelos Kaplan-Meier
survfits <- list(
  "Control Sac" = survfit(Surv(tiempoControl,control) ~ 1, data = df_ensayo_6, type = "kaplan-meier", conf.type = "log-log"),
  "BbQ 300 µL" = survfit(Surv(tiempoBbQ300,BbQ300) ~ 1, data = df_ensayo_6, type = "kaplan-meier", conf.type = "log-log"),
  "BbQ 450 µL" = survfit(Surv(tiempoBbQ450,BbQ450) ~ 1, data = df_ensayo_6, type = "kaplan-meier", conf.type = "log-log"),
  "BbQ 600 µL" = survfit(Surv(tiempoBbQ600,BbQ600) ~ 1, data = df_ensayo_6, type = "kaplan-meier", conf.type = "log-log"),
  "BbQ 750 µL" = survfit(Surv(tiempoBbQ750,BbQ750) ~ 1, data = df_ensayo_6, type = "kaplan-meier", conf.type = "log-log")
)


# Definir colores, tipos de línea y formas para los puntos
colores <- c("coral2", "coral2","coral2", "coral2", "darkblue")
linetypes <- c("solid", "dotted", "dashed", "longdash", "solid")
shapes <- c(15, 16, 17, 18, 15)

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
annotaciones <- c("A", "B", "C")
y_positions <- c(0.93, 0.425, 0.15)

# Crear el gráfico Kaplan-Meier combinado
grafico <- ggplot(df_surv, aes(x = time, y = surv, color = strata, linetype = strata, shape = strata)) +
  geom_step(size = 3) +  # Línea de supervivencia
  geom_point(size = 9.5) + # Puntos en las curvas
  scale_color_manual(values = colores) +
  scale_linetype_manual(values = linetypes) +
  scale_shape_manual(values = shapes) +
  labs(
    title = "Curva de Supervivencia - Ensayo #6",
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
  annotate("text", x = max(df_surv$time) * 0.75, y = y_positions,
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

print(summary(ensayo_6_controlSac.km))
print(summary(ensayo_6_BbQ300.km))
print(summary(ensayo_6_BbQ450.km))
print(summary(ensayo_6_BbQ600.km ))
print(summary(ensayo_6_BbQ750.km ))

print(summary(ensayo_6_controlSac.km)$table["median"])
print(summary(ensayo_6_BbQ300.km)$table["median"])
print(summary(ensayo_6_BbQ450.km)$table["median"])
print(summary(ensayo_6_BbQ600.km)$table["median"])
print(summary(ensayo_6_BbQ750.km)$table["median"])

#Test de Cox-Matlan

df_control <- data.frame(tiempo = ensayo_6$tiempo_dias,
                        muertas = ensayo_6$X.Muertas,
                        tratamiento = "Control")

df_Bb300 <- data.frame(tiempo = ensayo_6$tiempo_dias.1,
                       muertas = ensayo_6$X.Muertas.1,
                       tratamiento = "BbQ300")

df_Bb450 <- data.frame(tiempo = ensayo_6$tiempo_dias.2,
                       muertas = ensayo_6$X.Muertas.2,
                       tratamiento = "BbQ450")

df_Bb600 <- data.frame(tiempo = ensayo_6$tiempo_dias.3,
                       muertas = ensayo_6$X.Muertas.3,
                       tratamiento = "BbQ600")

df_Bb750 <- data.frame(tiempo = ensayo_6$tiempo_dias.4,
                       muertas = ensayo_6$X.Muertas.4,
                       tratamiento = "BbQ750")

df_combinado_todos <- rbind(df_control, df_Bb300, df_Bb450, df_Bb600, df_Bb750)
df_combinado_1 <- rbind(df_control, df_Bb300, df_Bb600, df_Bb750)
df_combinado_2 <- rbind(df_Bb300, df_Bb450)
df_combinado_3 <- rbind(df_Bb450, df_Bb600, df_Bb750)                        

print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_todos))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_1))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_2))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_3))

#Gráfica de Supervivencia
#png("Ensayo-6-Supervivencia.png",width = 5, height = 5, units = 'in', res = 600,pointsize=7.5)


#par(xpd=T, mar=par()$mar+c(0,0,0,10))
#plot(ensayo_6_controlSac.km, xlab="Tiempo (días)", ylab="Supervivencia = S(t)", main="Curva de Supervivencia - Ensayo #6", conf.int = F, col = "seagreen1", lwd=2, mark.time = TRUE, las = 1)
#lines(ensayo_6_BbQ300.km, conf.int = F, col = "darkolivegreen3" , lwd = 2, mark.time = TRUE)
#lines(ensayo_6_BbQ450.km, conf.int = F, col = "coral2" , lwd = 2, mark.time = TRUE)
#lines(ensayo_6_BbQ600.km, conf.int = F, col = "chocolate3" , lwd = 2, mark.time = TRUE)
#lines(ensayo_6_BbQ750.km, conf.int = F, col = "deeppink3" , lwd = 2, mark.time = TRUE)

#legend(7.5, 0.75, legend=c("Control  Sacarosa 25%m/V", "Bb 300 µL", "Bb 450 µL", "Bb 600 µL", "Bb 750 µL"), lty = 2, lwd = 3, col =c ("seagreen1", "darkolivegreen3", "coral2", "chocolate3", "deeppink2"), bty ="", cex=0.8)

#dev.off(3)
