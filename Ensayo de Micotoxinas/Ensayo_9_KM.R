#Rutina para realizar modelado de curvas de Supervivencia por método de Kaplan-Meier y su respectivo gráfico

#Librerias
library(survival)
library(survMisc)
library(KMsurv)
library(survminer)
library(ggfortify)
library(flexsurv)
library(actuar)

#Carga de Datos
ensayo_9 <- read.csv(file = "Ensayos de Micotoxinas de B. bassiana en L. humile - Ensayo #9 - Para analisis R.csv", header = T, sep=",", skip=1)
#head(ensayo_9) #Para verificar

df_ensayo_9 <- data.frame(tiempoControlSac = ensayo_9$tiempo_dias,
                          controlSac = ensayo_9$X.Muertas,
                          tiempoControlPGBQ = ensayo_9$tiempo_dias.1,
                          controlPGBQ = ensayo_9$X.Muertas.1,
                          tiempoBbQ150= ensayo_9$tiempo_dias.2,
                          BbQ150 = ensayo_9$X.Muertas.2,
                          tiempoBbQ300= ensayo_9$tiempo_dias.3,
                          BbQ300 = ensayo_9$X.Muertas.3,
                          tiempoBbQ450= ensayo_9$tiempo_dias.4,
                          BbQ450 = ensayo_9$X.Muertas.4,
                          tiempoBbQ600= ensayo_9$tiempo_dias.5,
                          BbQ600 = ensayo_9$X.Muertas.5,
                          tiempoBbQ750= ensayo_9$tiempo_dias.6,
                          BbQ750 = ensayo_9$X.Muertas.6)



library(ggpubr)

# Prueba de Shapiro-Wilk para cada grupo de tratamiento
shapiro_test_controlSac <- shapiro.test(df_ensayo_9$controlSac)
shapiro_test_controlPGBQ <- shapiro.test(df_ensayo_9$controlPGBQ)
shapiro_test_BbQ150 <- shapiro.test(df_ensayo_9$BbQ150)
shapiro_test_BbQ300 <- shapiro.test(df_ensayo_9$BbQ300)
shapiro_test_BbQ450 <- shapiro.test(df_ensayo_9$BbQ450)
shapiro_test_BbQ600 <- shapiro.test(df_ensayo_9$BbQ600)
shapiro_test_BbQ750 <- shapiro.test(df_ensayo_9$BbQ750)


# Imprimir los resultados
print(shapiro_test_controlSac)
print(shapiro_test_controlPGBQ)
print(shapiro_test_BbQ150)
print(shapiro_test_BbQ300)
print(shapiro_test_BbQ450)
print(shapiro_test_BbQ600)
print(shapiro_test_BbQ750)

# Prueba de Levene

library(car)

levene_test <- leveneTest(
  c(controlSac = df_ensayo_9$controlSac,
    controlPGBQ = df_ensayo_9$controlPGBQ,
    BbQ150 = df_ensayo_9$BbQ150,
    BbQ300 = df_ensayo_9$BbQ300,
    BbQ450 = df_ensayo_9$BbQ450,
    BbQ600 = df_ensayo_9$BbQ600,
    BbQ750 = df_ensayo_9$BbQ750) ~ 
    factor(rep(c("controlSac", "controlPGBQ","BbQ150", "BbQ300", "BbQ450", "BbQ600", "BbQ750"),
               each = nrow(df_ensayo_9))),
  center = median
)

levene_test_2 <- leveneTest(
  c(BbQ150 = df_ensayo_9$BbQ150,
    BbQ300 = df_ensayo_9$BbQ300,
    BbQ450 = df_ensayo_9$BbQ450,
    BbQ600 = df_ensayo_9$BbQ600,
    BbQ750 = df_ensayo_9$BbQ750) ~ 
    factor(rep(c("BbQ150", "BbQ300", "BbQ450", "BbQ600", "BbQ750"),
               each = nrow(df_ensayo_9))),
  center = median
)

levene_test_3 <- leveneTest(
  c(controlSac = df_ensayo_9$controlSac,
    controlPGBQ = df_ensayo_9$controlPGBQ) ~ 
    factor(rep(c("controlSac", "controlPGBQ"),
               each = nrow(df_ensayo_9))),
  center = median
)

print(levene_test)
print(levene_test_2)
print(levene_test_3)

#Modelo de Kaplan-Meier
ensayo_9_controlSac.km <- survfit(Surv(tiempoControlSac,controlSac) ~ 1, data = df_ensayo_9, type = "kaplan-meier")
ensayo_9_controlPGBQ.km  <- survfit(Surv(tiempoControlPGBQ,controlPGBQ) ~ 1, data = df_ensayo_9, type = "kaplan-meier")
ensayo_9_BbQ150.km <- survfit(Surv(tiempoBbQ150,BbQ150) ~ 1, data = df_ensayo_9, type = "kaplan-meier")
ensayo_9_BbQ300.km <- survfit(Surv(tiempoBbQ300,BbQ300) ~ 1, data = df_ensayo_9, type = "kaplan-meier")
ensayo_9_BbQ450.km <- survfit(Surv(tiempoBbQ450,BbQ450) ~ 1, data = df_ensayo_9, type = "kaplan-meier")
ensayo_9_BbQ600.km <- survfit(Surv(tiempoBbQ600,BbQ600) ~ 1, data = df_ensayo_9, type = "kaplan-meier")
ensayo_9_BbQ750.km <- survfit(Surv(tiempoBbQ750,BbQ750) ~ 1, data = df_ensayo_9, type = "kaplan-meier")

ensayo_9_survfits <- list(controlSac = ensayo_9_controlSac.km,
                          controlPGBQ = ensayo_9_controlPGBQ.km,
                          BbQ150 = ensayo_9_BbQ150.km,
                          BbQ300 = ensayo_9_BbQ300.km,
                          BbQ450 = ensayo_9_BbQ450.km,
                          BbQ600 = ensayo_9_BbQ600.km,
                          BbQ750 = ensayo_9_BbQ750.km)

# Definir los modelos Kaplan-Meier
survfits <- list(
  "Control Sac" = survfit(Surv(tiempoControlSac,controlSac) ~ 1, data = df_ensayo_9, type = "kaplan-meier", conf.type = "log-log"),
  "Control PGQ" = survfit(Surv(tiempoControlPGBQ,controlPGBQ) ~ 1, data = df_ensayo_9, type = "kaplan-meier", conf.type = "log-log"),
  "BbQ 150 µL" = survfit(Surv(tiempoBbQ150,BbQ150) ~ 1, data = df_ensayo_9, type = "kaplan-meier", conf.type = "log-log"),
  "BbQ 300 µL" = survfit(Surv(tiempoBbQ300,BbQ300) ~ 1, data = df_ensayo_9, type = "kaplan-meier", conf.type = "log-log"),
  "BbQ 450 µL" = survfit(Surv(tiempoBbQ450,BbQ450) ~ 1, data = df_ensayo_9, type = "kaplan-meier", conf.type = "log-log"),
  "BbQ 600 µL" = survfit(Surv(tiempoBbQ600,BbQ600) ~ 1, data = df_ensayo_9, type = "kaplan-meier", conf.type = "log-log"),
  "BbQ 750 µL" = survfit(Surv(tiempoBbQ750,BbQ750) ~ 1, data = df_ensayo_9, type = "kaplan-meier", conf.type = "log-log")
)


# Definir colores, tipos de línea y formas para los puntos
colores <- c("coral2", "coral2","coral2", "coral2", "coral2", "darkblue", "darkblue")
linetypes <- c("solid", "dotted", "dashed", "longdash", "twodash", "solid", "dotted")
shapes <- c(15, 16, 17, 18, 5, 15, 16)

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
annotaciones <- c("A", "B", "C", "D")
y_positions <- c(0.92, 0.74, 0.5, 0.25)

# Crear el gráfico Kaplan-Meier combinado
grafico <- ggplot(df_surv, aes(x = time, y = surv, color = strata, linetype = strata, shape = strata)) +
  geom_step(size = 3) +  # Línea de supervivencia
  geom_point(size = 9.5) + # Puntos en las curvas
  scale_color_manual(values = colores) +
  scale_linetype_manual(values = linetypes) +
  scale_shape_manual(values = shapes) +
  labs(
    title = "Curva de Supervivencia - Ensayo #9",
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

print(summary(ensayo_9_controlSac.km))
print(summary(ensayo_9_controlPGBQ.km))
print(summary(ensayo_9_BbQ300.km))
print(summary(ensayo_9_BbQ450.km))
print(summary(ensayo_9_BbQ600.km ))
print(summary(ensayo_9_BbQ750.km ))

print(summary(ensayo_9_controlSac.km)$table["median"])
print(summary(ensayo_9_controlPGBQ.km)$table["median"])
print(summary(ensayo_9_BbQ300.km)$table["median"])
print(summary(ensayo_9_BbQ450.km)$table["median"])
print(summary(ensayo_9_BbQ600.km)$table["median"])
print(summary(ensayo_9_BbQ750.km)$table["median"])

#Test de Cox-Matlan

df_controlSac <- data.frame(tiempo = ensayo_9$tiempo_dias,
                            muertas = ensayo_9$X.Muertas,
                            tratamiento = "ControlSac")

df_controlPGBQ <- data.frame(tiempo = ensayo_9$tiempo_dias.1,
                             muertas = ensayo_9$X.Muertas.1,
                             tratamiento = "ControlPGBQ")

df_Bb150 <- data.frame(tiempo = ensayo_9$tiempo_dias.2,
                       muertas = ensayo_9$X.Muertas.2,
                       tratamiento = "Bb150")

df_Bb300 <- data.frame(tiempo = ensayo_9$tiempo_dias.3,
                       muertas = ensayo_9$X.Muertas.3,
                       tratamiento = "Bb300")

df_Bb450 <- data.frame(tiempo = ensayo_9$tiempo_dias.4,
                       muertas = ensayo_9$X.Muertas.4,
                       tratamiento = "Bb450")

df_Bb600 <- data.frame(tiempo = ensayo_9$tiempo_dias.5,
                       muertas = ensayo_9$X.Muertas.5,
                       tratamiento = "Bb650")

df_Bb750 <- data.frame(tiempo = ensayo_9$tiempo_dias.6,
                       muertas = ensayo_9$X.Muertas.6,
                       tratamiento = "Bb750")

df_combinado_todos <- rbind(df_controlSac, df_controlPGBQ, df_Bb150, df_Bb300, df_Bb450, df_Bb600, df_Bb750)
df_combinado_1 <- rbind(df_controlSac, df_controlPGBQ)
df_combinado_2 <- rbind(df_Bb150, df_Bb300)
df_combinado_3 <- rbind(df_Bb300, df_Bb450)
df_combinado_4 <- rbind(df_Bb450, df_Bb600)
df_combinado_5 <- rbind(df_Bb600, df_Bb750)
df_combinado_6 <- rbind(df_controlSac, df_Bb150)
df_combinado_7 <- rbind(df_Bb150, df_Bb300, df_Bb750)

print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_todos))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_1))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_2))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_3))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_4))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_5))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_6))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_7))


#Gráfica de Supervivencia
#png("Ensayo-9-Supervivencia.png",width = 5, height = 5, units = 'in', res = 600,pointsize=7.5)
#par(xpd=T, mar=par()$mar+c(0,0,0,10))
#plot(ensayo_9_controlSac.km, xlab="Tiempo (días)", ylab="Supervivencia = S(t)", main="Curva de Supervivencia - Ensayo #9", conf.int = F, col = "seagreen1", lwd=2, mark.time = TRUE, las = 1)
#lines(ensayo_9_controlPGBQ.km, conf.int = F, col = "darkblue" , lwd = 2, mark.time = TRUE)
#lines(ensayo_9_BbQ150.km, conf.int = F, col = "violet" , lwd = 2, mark.time = TRUE)
#lines(ensayo_9_BbQ300.km, conf.int = F, col = "darkolivegreen3" , lwd = 2, mark.time = TRUE)
#lines(ensayo_9_BbQ450.km, conf.int = F, col = "coral2" , lwd = 2, mark.time = TRUE)
#lines(ensayo_9_BbQ600.km, conf.int = F, col = "chocolate3" , lwd = 2, mark.time = TRUE)
#lines(ensayo_9_BbQ750.km, conf.int = F, col = "deeppink3" , lwd = 2, mark.time = TRUE)

#legend(7.5, 0.75, legend=c("Control  Sacarosa 25%m/V", "Control  PGB+Q", "Bb 150 µL", "Bb 300 µL", "Bb 450 µL", "Bb 600 µL", "Bb 750 µL"), lty = 2, lwd = 3, col =c ("seagreen1","darkblue", "violet", "darkolivegreen3", "coral2", "chocolate3", "deeppink2"), bty ="", cex=0.8)
#dev.off(3)
