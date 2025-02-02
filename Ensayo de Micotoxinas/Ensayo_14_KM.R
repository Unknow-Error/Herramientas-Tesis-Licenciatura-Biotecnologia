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

#Modelo de Kaplan-Meier
ensayo_14_controlSac.km <- survfit(Surv(tiempoControlSac,controlSac) ~ 1, data = df_ensayo_14, type = "kaplan-meier")
ensayo_14_controlPGB.km  <- survfit(Surv(tiempoControlPGB,controlPGB) ~ 1, data = df_ensayo_14, type = "kaplan-meier")
ensayo_14_BbQ450.km <- survfit(Surv(tiempoBbQ450,BbQ450) ~ 1, data = df_ensayo_14, type = "kaplan-meier")
ensayo_14_Bb150.km <- survfit(Surv(tiempoBb150,Bb150) ~ 1, data = df_ensayo_14, type = "kaplan-meier")
ensayo_14_Bb300.km <- survfit(Surv(tiempoBb300,Bb300) ~ 1, data = df_ensayo_14, type = "kaplan-meier")
ensayo_14_Bb450.km <- survfit(Surv(tiempoBb450,Bb450) ~ 1, data = df_ensayo_14, type = "kaplan-meier")


ensayo_14_survfits <- list(controlSac = ensayo_14_controlSac.km,
                           controlPGB = ensayo_14_controlPGB.km,
                           BbQ450 = ensayo_14_BbQ450.km,
                           Bb150 = ensayo_14_Bb150.km,
                           Bb300 = ensayo_14_Bb300.km,
                           Bb450 = ensayo_14_Bb450.km)


#Gráfico
grafico_multiple <- 
  ggsurvplot(ensayo_14_survfits,
             data = df_ensayo_14,
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
                             "C: BbQ 450 µL", 
                             "D: Bb 150 µL", 
                             "E: Bb 300 µL", 
                             "F: Bb 450 µL"),
             conf.int = FALSE,
             conf.int.alpha = c(0.1),
             pval = FALSE,
             palette = c("seagreen1","darkblue", "violet", "darkolivegreen3", "coral2","deeppink2"),
             ggtheme=theme_bw()+theme(plot.title=element_text(hjust=0.5)))+
  labs(
    x = "Días",
    y = "Supervivencia = S(t)", 
    title = "Curva de Supervivencia - Ensayo #14")

# Añadir letras identificatorias en cada curva
annotaciones <- data.frame(x = c(11.1, 11.1, 11.1, 11.1, 11.1, 11.1), # Ajusta las coordenadas de x e y
                           y = c(0.8, 0.75, 0.1, 0.35, 0.2, 0.15),
                           label = c("A", "B", "C", "D", "E", "F"))

grafico_multiple$plot <- grafico_multiple$plot + 
  geom_text(data = annotaciones, aes(x = x, y = y, label = label), size = 5, fontface = "bold", color = "black")



for (km in ensayo_14_survfits){
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
df_combinado_6 <- rbind(df_BbQ450, df_Bb450)

print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_todos))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_1))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_2))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_3))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_4))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_5))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_6))