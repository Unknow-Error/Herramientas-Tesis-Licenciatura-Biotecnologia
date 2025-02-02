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
ensayo_13 <- read.csv(file = "Ensayos de Micotoxinas de B. bassiana en L. humile - Ensayo #13 - Para analisis R.csv", header = T, sep=",", skip=1)
#head(ensayo_13) #Para verificar


df_ensayo_13 <- data.frame(tiempoControlSac = ensayo_13$tiempo_dias,
                           controlSac = ensayo_13$X.Muertas,
                           tiempoControlPGB = ensayo_13$tiempo_dias.1,
                           controlPGB = ensayo_13$X.Muertas.1,
                           tiempoBb150= ensayo_13$tiempo_dias.2,
                           Bb150 = ensayo_13$X.Muertas.2,
                           tiempoBb300= ensayo_13$tiempo_dias.3,
                           Bb300 = ensayo_13$X.Muertas.3,
                           tiempoBb450= ensayo_13$tiempo_dias.4,
                           Bb450 = ensayo_13$X.Muertas.4,
                           tiempoBb600= ensayo_13$tiempo_dias.5,
                           Bb600 = ensayo_13$X.Muertas.5,
                           tiempoBb750= ensayo_13$tiempo_dias.6,
                           Bb750 = ensayo_13$X.Muertas.6)

library(ggpubr)

# Prueba de Shapiro-Wilk para cada grupo de tratamiento
shapiro_test_controlSac <- shapiro.test(df_ensayo_13$controlSac)
shapiro_test_controlPGB <- shapiro.test(df_ensayo_13$controlPGB)
shapiro_test_Bb150 <- shapiro.test(df_ensayo_13$Bb150)
shapiro_test_Bb300 <- shapiro.test(df_ensayo_13$Bb300)
shapiro_test_Bb450 <- shapiro.test(df_ensayo_13$Bb450)
shapiro_test_Bb600 <- shapiro.test(df_ensayo_13$Bb600)
shapiro_test_Bb750 <- shapiro.test(df_ensayo_13$Bb750)


# Imprimir los resultados
print(shapiro_test_controlSac)
print(shapiro_test_controlPGB)
print(shapiro_test_Bb150)
print(shapiro_test_Bb300)
print(shapiro_test_Bb450)
print(shapiro_test_Bb600)
print(shapiro_test_Bb750)

# Prueba de Levene

library(car)

levene_test <- leveneTest(
  c(controlSac = df_ensayo_13$controlSac,
    controlPGB = df_ensayo_13$controlPGB,
    Bb150 = df_ensayo_13$Bb150,
    Bb300 = df_ensayo_13$Bb300,
    Bb450 = df_ensayo_13$Bb450,
    Bb600 = df_ensayo_13$Bb600,
    Bb750 = df_ensayo_13$Bb750) ~ 
    factor(rep(c("controlSac", "controlPGB","Bb150", "Bb300", "Bb450", "Bb600", "Bb750"),
               each = nrow(df_ensayo_13))),
  center = median
)

levene_test_2 <- leveneTest(
  c(Bb150 = df_ensayo_13$Bb150,
    Bb300 = df_ensayo_13$Bb300,
    Bb450 = df_ensayo_13$Bb450,
    Bb600 = df_ensayo_13$Bb600,
    Bb750 = df_ensayo_13$Bb750) ~ 
    factor(rep(c("Bb150", "Bb300", "Bb450", "Bb600", "Bb750"),
               each = nrow(df_ensayo_13))),
  center = median
)

levene_test_3 <- leveneTest(
  c(controlSac = df_ensayo_13$controlSac,
    controlPGB = df_ensayo_13$controlPGB) ~ 
    factor(rep(c("controlSac", "controlPGB"),
               each = nrow(df_ensayo_13))),
  center = median
)

print(levene_test)
print(levene_test_2)
print(levene_test_3)

#Modelo de Kaplan-Meier
ensayo_13_controlSac.km <- survfit(Surv(tiempoControlSac,controlSac) ~ 1, data = df_ensayo_13, type = "kaplan-meier")
ensayo_13_controlPGB.km  <- survfit(Surv(tiempoControlPGB,controlPGB) ~ 1, data = df_ensayo_13, type = "kaplan-meier")
ensayo_13_Bb150.km <- survfit(Surv(tiempoBb150,Bb150) ~ 1, data = df_ensayo_13, type = "kaplan-meier")
ensayo_13_Bb300.km <- survfit(Surv(tiempoBb300,Bb300) ~ 1, data = df_ensayo_13, type = "kaplan-meier")
ensayo_13_Bb450.km <- survfit(Surv(tiempoBb450,Bb450) ~ 1, data = df_ensayo_13, type = "kaplan-meier")
ensayo_13_Bb600.km <- survfit(Surv(tiempoBb600,Bb600) ~ 1, data = df_ensayo_13, type = "kaplan-meier")
ensayo_13_Bb750.km <- survfit(Surv(tiempoBb750,Bb750) ~ 1, data = df_ensayo_13, type = "kaplan-meier")

ensayo_13_survfits <- list(controlSac = ensayo_13_controlSac.km,
                           controlPGB = ensayo_13_controlPGB.km,
                           Bb150 = ensayo_13_Bb150.km,
                           Bb300 = ensayo_13_Bb300.km,
                           Bb450 = ensayo_13_Bb450.km,
                           Bb600 = ensayo_13_Bb600.km,
                           Bb750 = ensayo_13_Bb750.km)


#Gráfico
grafico_multiple <- 
  ggsurvplot(ensayo_13_survfits,
             data = df_ensayo_13,
             combine = TRUE,
             lwd = 1.5,
             legend.title = "",
             font.legend = c(19, "bold", "black"),
             font.main = c(20, "bold", "darkblue"),
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             font.tickslab = c(18, "plain", "black"),
             legend.labs = c("A: Control  Sacarosa 25%m/V", 
                             "B: Control  PGB", 
                             "C: Bb 150 µL", 
                             "D: Bb 300 µL", 
                             "E: Bb 450 µL", 
                             "F: Bb 600 µL", 
                             "G: Bb 750 µL"),
             conf.int = FALSE,
             conf.int.alpha = c(0.1),
             pval = FALSE,
             palette = c("seagreen1","darkblue", "violet", "darkolivegreen3", "coral2", "cyan", "deeppink2"),
             ggtheme=theme_bw()+theme(plot.title=element_text(hjust=0.5)))+
  labs(
    x = "Días",
    y = "Supervivencia = S(t)", 
    title = "Curva de Supervivencia - Ensayo #13")

# Añadir letras identificatorias en cada curva
annotaciones <- data.frame(x = c(11.1, 11.1, 8.625, 8.625, 8.625, 8.625, 8.625), # Ajusta las coordenadas de x e y
                           y = c(0.85, 0.4, 0.02, 0.07, 0.16, 0.25, 0.12),
                           label = c("A", "B", "C", "D", "E", "F", "G"))

grafico_multiple$plot <- grafico_multiple$plot + 
  geom_text(data = annotaciones, aes(x = x, y = y, label = label), size = 5, fontface = "bold", color = "black")



for (km in ensayo_13_survfits){
  surv_median <- as.vector(summary(km)$table["median"])
  df <- data.frame(x1 = surv_median, x2=  surv_median,
                   y1 = rep(0, length(surv_median)), y2 = rep(0.5, length(surv_median)))
  
  grafico_multiple$plot <- grafico_multiple$plot + 
    geom_segment(aes(x = 0, y = 0.5, xend = max(surv_median), yend = 0.5),
                 linetype = "dashed", col = "darkblue", size = 0.75)+ # horizontal segment
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df,
                 linetype = "dashed", col = "darkblue", size = 0.75) # vertical segments
}

print(grafico_multiple)

#Estadísticos de cada tratamiento

print(summary(ensayo_13_controlSac.km))
print(summary(ensayo_13_controlPGB.km))
print(summary(ensayo_13_Bb300.km))
print(summary(ensayo_13_Bb450.km))
print(summary(ensayo_13_Bb600.km ))
print(summary(ensayo_13_Bb750.km ))

print(summary(ensayo_13_controlSac.km)$table["median"])
print(summary(ensayo_13_controlPGB.km)$table["median"])
print(summary(ensayo_13_Bb300.km)$table["median"])
print(summary(ensayo_13_Bb450.km)$table["median"])
print(summary(ensayo_13_Bb600.km)$table["median"])
print(summary(ensayo_13_Bb750.km)$table["median"])

#Test de Cox-Matlan

df_controlSac <- data.frame(tiempo = ensayo_13$tiempo_dias,
                            muertas = ensayo_13$X.Muertas,
                            tratamiento = "ControlSac")

df_controlPGB <- data.frame(tiempo = ensayo_13$tiempo_dias.1,
                             muertas = ensayo_13$X.Muertas.1,
                             tratamiento = "ControlPGB")

df_Bb150 <- data.frame(tiempo = ensayo_13$tiempo_dias.2,
                       muertas = ensayo_13$X.Muertas.2,
                       tratamiento = "Bb150")

df_Bb300 <- data.frame(tiempo = ensayo_13$tiempo_dias.3,
                       muertas = ensayo_13$X.Muertas.3,
                       tratamiento = "Bb300")

df_Bb450 <- data.frame(tiempo = ensayo_13$tiempo_dias.4,
                       muertas = ensayo_13$X.Muertas.4,
                       tratamiento = "Bb450")

df_Bb600 <- data.frame(tiempo = ensayo_13$tiempo_dias.5,
                       muertas = ensayo_13$X.Muertas.5,
                       tratamiento = "Bb650")

df_Bb750 <- data.frame(tiempo = ensayo_13$tiempo_dias.6,
                       muertas = ensayo_13$X.Muertas.6,
                       tratamiento = "Bb750")

df_combinado_todos <- rbind(df_controlPGB, df_controlSac, df_Bb150, df_Bb300, df_Bb450, df_Bb600, df_Bb750)
df_combinado_1 <- rbind(df_controlSac, df_controlPGB)
df_combinado_2 <- rbind(df_Bb150, df_Bb300)
df_combinado_3 <- rbind(df_Bb300, df_Bb450)
df_combinado_4 <- rbind(df_Bb450, df_Bb600)
df_combinado_5 <- rbind(df_Bb600, df_Bb750)
df_combinado_6 <- rbind(df_controlPGB, df_Bb150)

print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_todos))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_1))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_2))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_3))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_4))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_5))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_6))