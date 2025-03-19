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
ensayo_7 <- read.csv(file = "Ensayos de Micotoxinas de B. bassiana en L. humile - Ensayo #7 - Para analisis R.csv", header = T, sep=",", skip=1)
#head(ensayo_7) #Para verificar

#Modelo de Kaplan-Meier
df_ensayo_7 <- data.frame(tiempocontrolSacPatri = ensayo_7$tiempo_dias,
                          controlSacPatri = ensayo_7$X.Muertas,
                          tiempocontrolSacDani = ensayo_7$tiempo_dias.1,
                          controlSacDani = ensayo_7$X.Muertas.1,
                          tiempoBbQ450Patri = ensayo_7$tiempo_dias.2,
                          BbQ450Patri = ensayo_7$X.Muertas.2,
                          tiempoBbQ450Dani = ensayo_7$tiempo_dias.3,
                          BbQ450Dani = ensayo_7$X.Muertas.3)

ensayo_7_controlSacPatri.km <- survfit(Surv(tiempocontrolSacPatri,controlSacPatri) ~ 1, data = df_ensayo_7, type = "kaplan-meier")
ensayo_7_controlSacDani.km  <- survfit(Surv(tiempocontrolSacDani,controlSacDani) ~ 1, data = df_ensayo_7, type = "kaplan-meier")
ensayo_7_BbQ450Patri.km <- survfit(Surv(tiempoBbQ450Patri,BbQ450Patri) ~ 1, data = df_ensayo_7, type = "kaplan-meier")
ensayo_7_BbQ450Dani.km <- survfit(Surv(tiempoBbQ450Dani,BbQ450Dani) ~ 1, data = df_ensayo_7, type = "kaplan-meier")

ensayo_7_survfits <- list(controlPatri = ensayo_7_controlSacPatri.km,
                          controlDani = ensayo_7_controlSacDani.km,
                          BbQ450Patri = ensayo_7_BbQ450Patri.km,
                          BbQ450Dani = ensayo_7_BbQ450Dani.km)

library(ggpubr)

# Prueba de Shapiro-Wilk para cada grupo de tratamiento
shapiro_test_controlSacPatri <- shapiro.test(df_ensayo_7$controlSacPatri)
shapiro_test_controlSacDani <- shapiro.test(df_ensayo_7$controlSacDani)
shapiro_test_BbQ450Patri <- shapiro.test(df_ensayo_7$BbQ450Patri)
shapiro_test_BbQ450Dani  <- shapiro.test(df_ensayo_7$BbQ450Dani)


# Imprimir los resultados
print(shapiro_test_controlSacPatri)
print(shapiro_test_controlSacDani)
print(shapiro_test_BbQ450Patri)
print(shapiro_test_BbQ450Dani)

# Prueba de Levene

library(car)

levene_test <- leveneTest(
  c(controlSacPatri = df_ensayo_7$controlSacPatri,
    controlSacDani = df_ensayo_7$controlSacDani,
    BbQ450Patri = df_ensayo_7$BbQ450Patri,
    BbQ450Dani = df_ensayo_7$BbQ450Dani) ~ 
    factor(rep(c("controlSacPatri", "controlSacDani", "BbQ450Patri", "BbQ450Dani"),
               each = nrow(df_ensayo_7))),
   center = median
)

levene_test_2 <- leveneTest(
  c(controlSacPatri = df_ensayo_7$controlSacPatri,
    controlSacDani = df_ensayo_7$controlSacDani) ~ 
    factor(rep(c("controlSacPatri", "controlSacDani"),
               each = nrow(df_ensayo_7))),
  center = median
)

levene_test_3 <- leveneTest(
  c(BbQ450Patri = df_ensayo_7$BbQ450Patri,
    BbQ450Dani = df_ensayo_7$BbQ450Dani) ~ 
    factor(rep(c("BbQ450Patri", "BbQ450Dani"),
               each = nrow(df_ensayo_7))),
  center = median
)

print(levene_test)
print(levene_test_2)
print(levene_test_3)


# Definir los modelos Kaplan-Meier
survfits <- list(
  "Control Sac - Bernal" = survfit(Surv(tiempocontrolSacPatri,controlSacPatri) ~ 1, data = df_ensayo_7, type = "kaplan-meier", conf.type = "log-log"),
  "Control Sac - Gonnet" = survfit(Surv(tiempocontrolSacDani,controlSacDani) ~ 1, data = df_ensayo_7, type = "kaplan-meier", conf.type = "log-log"),
  "BbQ 450 µL - Bernal" = survfit(Surv(tiempoBbQ450Patri,BbQ450Patri) ~ 1, data = df_ensayo_7, type = "kaplan-meier", conf.type = "log-log"),
  "BbQ 450 µL - Gonnet" = survfit(Surv(tiempoBbQ450Dani,BbQ450Dani) ~ 1, data = df_ensayo_7, type = "kaplan-meier", conf.type = "log-log")
)


# Definir colores, tipos de línea y formas para los puntos
colores <- c("palegreen2", "coral2", "darkblue", "darkblue")
linetypes <- c("solid", "dashed", "solid", "dotted")
shapes <- c(15, 16, 15, 16)

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
y_positions <- c(0.80, 0.42, 0.16)

# Crear el gráfico Kaplan-Meier combinado
grafico <- ggplot(df_surv, aes(x = time, y = surv, color = strata, linetype = strata, shape = strata)) +
  geom_step(size = 3) +  # Línea de supervivencia
  geom_point(size = 9.5) + # Puntos en las curvas
  scale_color_manual(values = colores) +
  scale_linetype_manual(values = linetypes) +
  scale_shape_manual(values = shapes) +
  labs(
    title = "Curva de Supervivencia - Ensayo #7",
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

print(summary(ensayo_7_controlSacPatri.km))
print(summary(ensayo_7_controlSacDani.km))
print(summary(ensayo_7_BbQ450Patri.km))
print(summary(ensayo_7_BbQ450Dani.km ))

print(summary(ensayo_7_controlSacPatri.km)$table["median"])
print(summary(ensayo_7_controlSacDani.km)$table["median"])
print(summary(ensayo_7_BbQ450Patri.km)$table["median"])
print(summary(ensayo_7_BbQ450Dani.km)$table["median"])


#Test de Cox-Matlan

df_controlPatri <- data.frame(tiempo = ensayo_7$tiempo_dias,
                         muertas = ensayo_7$X.Muertas,
                         tratamiento = "ControlPatri")

df_controlDani <- data.frame(tiempo = ensayo_6$tiempo_dias.1,
                       muertas = ensayo_6$X.Muertas.1,
                       tratamiento = "ControlDani")

df_BbQ450Patri <- data.frame(tiempo = ensayo_6$tiempo_dias.2,
                       muertas = ensayo_6$X.Muertas.2,
                       tratamiento = "Bb450Patri")

df_BbQ450Dani <- data.frame(tiempo = ensayo_6$tiempo_dias.3,
                       muertas = ensayo_6$X.Muertas.3,
                       tratamiento = "Bb450Dani")

df_combinado_todos <- rbind(df_controlPatri, df_controlDani, df_BbQ450Patri, df_BbQ450Dani)
df_combinado_1 <- rbind(df_controlPatri, df_controlDani)
df_combinado_2 <- rbind(df_BbQ450Patri, df_BbQ450Dani)
df_combinado_3 <- rbind(df_controlPatri, df_BbQ450Patri)
df_combinado_4 <- rbind(df_controlDani, df_BbQ450Dani)

print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_todos))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_1))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_2))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_3))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_4))


#Gráfica de Supervivencia
#png("Ensayo-7-Supervivencia.png",width = 5, height = 5, units = 'in', res = 600,pointsize=7.5)
#par(xpd=T, mar=par()$mar+c(0,0,0,10))
#plot(ensayo_7_controlSacPatri.km, xlab="Tiempo (días)", ylab="Supervivencia = S(t)", main="Curva de Supervivencia - Ensayo #7", conf.int = F, col = "seagreen1", lwd=2, mark.time = TRUE, las = 1)
#lines(ensayo_7_controlSacDani.km, conf.int = F, col = "cyan" , lwd = 2, mark.time = TRUE)
#lines(ensayo_7_BbQ450Patri.km, conf.int = F, col = "darkorchid3" , lwd = 2, mark.time = TRUE)
#lines(ensayo_7_BbQ450Dani.km, conf.int = F, col = "coral2" , lwd = 2, mark.time = TRUE)

#legend(8.5, 0.75, legend=c("Control Sac - Casa Patri", "Control Sac - Casa Dani", "Bb 450 µL - Casa Patri", "Bb 450 µL - Casa Dani"), lty = 2, lwd = 3, col =c ("seagreen1","cyan", "darkorchid3", "coral2"), bty ="", cex=0.8)
#dev.off(3)
