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
ensayo_2 <- read.csv(file = "Ensayos de Micotoxinas de B. bassiana en L. humile - Ensayo #2 - Para analisis R.csv", header = T, sep=",", skip=1)
#head(ensayo_2) #Para verificar

df_ensayo_2 <- data.frame(tiempoControlAgua = ensayo_2$tiempo_dias,
                          controlAgua = ensayo_2$X.Muertas,
                          tiempoControlSac = ensayo_2$tiempo_dias.1,
                          controlSac = ensayo_2$X.Muertas.1,
                          tiempoBb100= ensayo_2$tiempo_dias.2,
                          Bb100 = ensayo_2$X.Muertas.2,
                          tiempoBbQ100= ensayo_2$tiempo_dias.3,
                          BbQ100 = ensayo_2$X.Muertas.3)


library(ggpubr)

# Prueba de Shapiro-Wilk para cada grupo de tratamiento
shapiro_test_controlAgua <- shapiro.test(df_ensayo_2$controlAgua)
shapiro_test_controlSac <- shapiro.test(df_ensayo_2$controlSac)
shapiro_test_Bb100 <- shapiro.test(df_ensayo_2$Bb100)
shapiro_test_BbQ100 <- shapiro.test(df_ensayo_2$BbQ100)

# Imprimir los resultados
print(shapiro_test_controlAgua)
print(shapiro_test_controlSac)
print(shapiro_test_Bb100)
print(shapiro_test_BbQ100)


# Prueba de Levene

library(car)

levene_test <- leveneTest(
  c(controlSac = df_ensayo_2$controlAgua,
    Bb100 = df_ensayo_2$controlSac,
    Bb75 = df_ensayo_2$Bb100,
    BbQ100 = df_ensayo_2$BbQ100) ~ 
    factor(rep(c("controlAgua", "controlSac", "Bb100", "BbQ100"),
               each = nrow(df_ensayo_2))),
  center = median
)

print(levene_test)

# Definir los modelos Kaplan-Meier
survfits <- list(
  "Control Agua" = survfit(Surv(tiempoControlAgua,controlAgua) ~ 1, data = df_ensayo_2, type = "kaplan-meier"),
  "Control Sac" = survfit(Surv(tiempoControlSac,controlSac) ~ 1, data = df_ensayo_2, type = "kaplan-meier"),
  "Bb 150 μL" = survfit(Surv(tiempoBb100,Bb100) ~ 1, data = df_ensayo_2, type = "kaplan-meier"),
  "BbQ 150 μL" = survfit(Surv(tiempoBbQ100,BbQ100) ~ 1, data = df_ensayo_2, type = "kaplan-meier")
)


# Definir colores, tipos de línea y formas para los puntos
colores <- c("coral3","palegreen2", "darkblue", "darkblue")
linetypes <- c("solid", "solid", "solid", "dotted")
shapes <- c(15, 16, 17, 18)

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
y_positions <- c(0.96, 0.86)

# Crear el gráfico Kaplan-Meier combinado
grafico <- ggplot(df_surv, aes(x = time, y = surv, color = strata, linetype = strata, shape = strata)) +
  geom_step(size = 3) +  # Línea de supervivencia
  geom_point(size = 9.5) + # Puntos en las curvas
  scale_color_manual(values = colores) +
  scale_linetype_manual(values = linetypes) +
  scale_shape_manual(values = shapes) +
  labs(
    title = "Curva de Supervivencia - Ensayo #2",
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

print(summary(ensayo_2_controlAgua.km))
print(summary(ensayo_2_controlSac.km))
print(summary(ensayo_2_Bb100.km))
print(summary(ensayo_2_BbQ100.km))

print(summary(ensayo_2_controlAgua.km)$table["median"])
print(summary(ensayo_2_controlSac.km)$table["median"])
print(summary(ensayo_2_Bb100.km)$table["median"])
print(summary(ensayo_2_BbQ100.km)$table["median"])


#Gráfica de Supervivencia
#png("Ensayo-2-Supervivencia.png",width = 5, height = 5, units = 'in', res = 600,pointsize=7.5)
#par(xpd=T, mar=par()$mar+c(0,0,0,11))
#plot(ensayo_2_controlSac.km, xlab="Tiempo (días)", ylab="Supervivencia = S(t)", main="Curva de Supervivencia - Ensayo #2", conf.int = F, col = "seagreen1", lwd=2, mark.time = TRUE, las = 1)
#lines(ensayo_2_controlAgua.km, conf.int = F, col = "darkblue" , lwd = 2, mark.time = TRUE)
#lines(ensayo_2_Bb100.km, conf.int = F, col = "coral1" , lwd = 2, mark.time = TRUE)
#lines(ensayo_2_BbQ100.km, conf.int = F, col = "deeppink2" , lwd = 2, mark.time = TRUE)

#legend(11.5, 0.75, legend=c("Control  Sacarosa 25%", "Control Agua", "Bb 100%", "Bb+Q 100%"), lty = 2, lwd = 3, col =c ("seagreen1","darkblue", "coral1", "deeppink2"), bty ="", cex=0.8)
#dev.off(3)
