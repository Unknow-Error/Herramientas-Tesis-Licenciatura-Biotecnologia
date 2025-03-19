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
ensayo_4 <- read.csv(file = "Ensayos de Micotoxinas de B. bassiana en L. humile - Ensayo #4 - Para analisis R.csv", header = T, sep=",", skip=1)
#head(ensayo_4) #Para verificar

df_ensayo_4 <- data.frame(tiempoControlSac = ensayo_4$tiempo_dias,
                          controlSac = ensayo_4$X.Muertas,
                          tiempoBbQ1= ensayo_4$tiempo_dias.1,
                          BbQ1 = ensayo_4$X.Muertas.1,
                          tiempoBbQ2= ensayo_4$tiempo_dias.2,
                          BbQ2 = ensayo_4$X.Muertas.2)


library(ggpubr)

# Prueba de Shapiro-Wilk para cada grupo de tratamiento
shapiro_test_controlSac <- shapiro.test(df_ensayo_4$controlSac)
shapiro_test_BbQ1 <- shapiro.test(df_ensayo_4$BbQ1)
shapiro_test_BbQ2 <- shapiro.test(df_ensayo_4$BbQ2)

# Imprimir los resultados
print(shapiro_test_controlSac)
print(shapiro_test_BbQ1)
print(shapiro_test_BbQ2)


# Prueba de Levene

library(car)

levene_test <- leveneTest(
  c(controlSac = df_ensayo_4$controlSac,
    BbQ1 = df_ensayo_4$BbQ1,
    BbQ2 = df_ensayo_4$BbQ2) ~ 
    factor(rep(c("controlSac", "BbQ1", "BbQ2"),
               each = nrow(df_ensayo_4))),
  center = median
)

levene_test_2 <- leveneTest(
  c(BbQ1 = df_ensayo_4$BbQ1,
    BbQ2 = df_ensayo_4$BbQ2) ~ 
    factor(rep(c("BbQ1", "BbQ2"),
               each = nrow(df_ensayo_4))),
  center = median
)

print(levene_test)
print(levene_test_2)

# Definir los modelos Kaplan-Meier
survfits <- list(
  "Control Sac" = survfit(Surv(tiempoControlSac,controlSac) ~ 1, data = df_ensayo_4, type = "kaplan-meier"),
  "BbQ N°1" = survfit(Surv(tiempoBbQ1,BbQ1) ~ 1, data = df_ensayo_4, type = "kaplan-meier"),
  "BbQ N°2" = survfit(Surv(tiempoBbQ2,BbQ2) ~ 1, data = df_ensayo_4, type = "kaplan-meier")
)


# Definir colores, tipos de línea y formas para los puntos
colores <- c("palegreen2", "palegreen2", "darkblue")
linetypes <- c("solid", "dotted", "solid")
shapes <- c(15, 16, 15)

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
y_positions <- c(0.90, 0.10)

# Crear el gráfico Kaplan-Meier combinado
grafico <- ggplot(df_surv, aes(x = time, y = surv, color = strata, linetype = strata, shape = strata)) +
  geom_step(size = 3) +  # Línea de supervivencia
  geom_point(size = 9.5) + # Puntos en las curvas
  scale_color_manual(values = colores) +
  scale_linetype_manual(values = linetypes) +
  scale_shape_manual(values = shapes) +
  labs(
    title = "Curva de Supervivencia - Ensayo #4",
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

print(summary(ensayo_4_controlSac.km))
print(summary(ensayo_4_BbQ1.km))
print(summary(ensayo_4_BbQ2.km))


print(summary(ensayo_4_controlSac.km)$table["median"])
print(summary(ensayo_4_BbQ1.km)$table["median"])
print(summary(ensayo_4_BbQ2.km)$table["median"])


#Test de Cox-Matlan

df_control <- data.frame(tiempo = ensayo_4$tiempo_dias,
                         muertas = ensayo_4$X.Muertas,
                         tratamiento = "Control")

df_BbQ1 <- data.frame(tiempo = ensayo_4$tiempo_dias.1,
                       muertas = ensayo_4$X.Muertas.1,
                       tratamiento = "BbQ1")

df_BbQ2<- data.frame(tiempo = ensayo_4$tiempo_dias.2,
                       muertas = ensayo_4$X.Muertas.2,
                       tratamiento = "BbQ2")

df_combinado_todos <- rbind(df_control, df_BbQ1, df_BbQ2)
df_combinado_1 <- rbind(df_control, df_BbQ1)
df_combinado_2 <- rbind(df_control, df_BbQ2)
df_combinado_3 <- rbind(df_BbQ1,df_BbQ2)

                   

print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_todos, rho=0))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_1, rho=0))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_2, rho=0))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_3, rho=0))


#Gráfica de Supervivencia
#png("Ensayo-4-Supervivencia.png",width = 5, height = 5, units = 'in', res = 600,pointsize=7.5)
#par(xpd=T, mar=par()$mar+c(0,0,0,11))
#plot(ensayo_4_controlSac.km, xlab="Tiempo (días)", ylab="Supervivencia = S(t)", main="Curva de Supervivencia - Ensayo #4", conf.int = F, col = "seagreen1", lwd=2, mark.time = TRUE, las = 1)
#lines(ensayo_4_BbQ1.km, conf.int = F, col = "coral1" , lwd = 2, mark.time = TRUE)
#lines(ensayo_4_BbQ2.km, conf.int = F, col = "cadetblue4" , lwd = 2, mark.time = TRUE)

#legend(19, 0.75, legend=c("Control  Sacarosa 25%", "Bb+Q N°1", "Bb+Q N°2"), lty = 2, lwd = 3, col =c ("seagreen1","coral1", "deeppink2"), bty ="", cex=0.8)
#dev.off(3)
