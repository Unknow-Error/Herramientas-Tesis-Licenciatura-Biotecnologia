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
ensayo_11 <- read.csv(file = "Ensayos de Micotoxinas de B. bassiana en L. humile - Ensayo #11 - Para analisis R.csv", header = T, sep=",", skip=1)
#head(ensayo_11) #Para verificar


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

library(ggpubr)

# Prueba de Shapiro-Wilk para cada grupo de tratamiento
shapiro_test_controlSac <- shapiro.test(df_ensayo_11$controlSac)
shapiro_test_controlPGBQ <- shapiro.test(df_ensayo_11$controlPGBQ)
shapiro_test_BbQ150 <- shapiro.test(df_ensayo_11$BbQ150)
shapiro_test_BbQ300 <- shapiro.test(df_ensayo_11$BbQ300)
shapiro_test_BbQ450 <- shapiro.test(df_ensayo_11$BbQ450)
shapiro_test_BbQ600 <- shapiro.test(df_ensayo_11$BbQ600)
shapiro_test_BbQ750 <- shapiro.test(df_ensayo_11$BbQ750)


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
  c(controlSac = df_ensayo_11$controlSac,
    controlPGBQ = df_ensayo_11$controlPGBQ,
    BbQ150 = df_ensayo_11$BbQ150,
    BbQ300 = df_ensayo_11$BbQ300,
    BbQ450 = df_ensayo_11$BbQ450,
    BbQ600 = df_ensayo_11$BbQ600,
    BbQ750 = df_ensayo_11$BbQ750) ~ 
    factor(rep(c("controlSac", "controlPGBQ","BbQ150", "BbQ300", "BbQ450", "BbQ600", "BbQ750"),
               each = nrow(df_ensayo_11))),
  center = median
)

levene_test_2 <- leveneTest(
  c(BbQ150 = df_ensayo_11$BbQ150,
    BbQ300 = df_ensayo_11$BbQ300,
    BbQ450 = df_ensayo_11$BbQ450,
    BbQ600 = df_ensayo_11$BbQ600,
    BbQ750 = df_ensayo_11$BbQ750) ~ 
    factor(rep(c("BbQ150", "BbQ300", "BbQ450", "BbQ600", "BbQ750"),
               each = nrow(df_ensayo_11))),
  center = median
)

levene_test_3 <- leveneTest(
  c(controlSac = df_ensayo_11$controlSac,
    controlPGBQ = df_ensayo_11$controlPGBQ) ~ 
    factor(rep(c("controlSac", "controlPGBQ"),
               each = nrow(df_ensayo_11))),
  center = median
)

print(levene_test)
print(levene_test_2)
print(levene_test_3)


#Modelo de Kaplan-Meier
ensayo_11_controlSac.km <- survfit(Surv(tiempoControlSac,controlSac) ~ 1, data = df_ensayo_11, type = "kaplan-meier")
ensayo_11_controlPGBQ.km  <- survfit(Surv(tiempoControlPGBQ,controlPGBQ) ~ 1, data = df_ensayo_11, type = "kaplan-meier")
ensayo_11_BbQ150.km <- survfit(Surv(tiempoBbQ150,BbQ150) ~ 1, data = df_ensayo_11, type = "kaplan-meier")
ensayo_11_BbQ300.km <- survfit(Surv(tiempoBbQ300,BbQ300) ~ 1, data = df_ensayo_11, type = "kaplan-meier")
ensayo_11_BbQ450.km <- survfit(Surv(tiempoBbQ450,BbQ450) ~ 1, data = df_ensayo_11, type = "kaplan-meier")
ensayo_11_BbQ600.km <- survfit(Surv(tiempoBbQ600,BbQ600) ~ 1, data = df_ensayo_11, type = "kaplan-meier")
ensayo_11_BbQ750.km <- survfit(Surv(tiempoBbQ750,BbQ750) ~ 1, data = df_ensayo_11, type = "kaplan-meier")

ensayo_11_survfits <- list(controlSac = ensayo_11_controlSac.km,
                           controlPGBQ = ensayo_11_controlPGBQ.km,
                           BbQ150 = ensayo_11_BbQ150.km,
                           BbQ300 = ensayo_11_BbQ300.km,
                           BbQ450 = ensayo_11_BbQ450.km,
                           BbQ600 = ensayo_11_BbQ600.km,
                           BbQ750 = ensayo_11_BbQ750.km)


#Gráfico
grafico_multiple <- 
  ggsurvplot(ensayo_11_survfits,
             data = df_ensayo_11,
             combine = TRUE,
             lwd = 1.5,
             legend.title = "",
             font.legend = c(19, "bold", "black"),
             font.main = c(20, "bold", "darkblue"),
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             font.tickslab = c(18, "plain", "black"),
             legend.labs = c("A: Control  Sacarosa 25%m/V", 
                             "B: Control  PGB+Q", 
                             "C: BbQ 150 µL", 
                             "D: BbQ 300 µL", 
                             "E: BbQ 450 µL", 
                             "F: BbQ 600 µL", 
                             "G: BbQ 750 µL"),
             conf.int = FALSE,
             conf.int.alpha = c(0.1),
             pval = FALSE,
             palette = c("seagreen1","darkblue", "violet", "darkolivegreen3", "coral2", "cyan", "deeppink2"),
             ggtheme=theme_bw()+theme(plot.title=element_text(hjust=0.5)))+
  labs(
    x = "Días",
    y = "Supervivencia = S(t)", 
    title = "Curva de Supervivencia - Ensayo #11")

# Añadir letras identificatorias en cada curva
annotaciones <- data.frame(x = c(7.1, 7.1, 7.2, 7.1, 7.1, 7.1, 7.2), # Ajusta las coordenadas de x e y
                           y = c(0.9, 0.8, 0.46, 0.46, 0.4, 0.5, 0.5),
                           label = c("A", "B", "C", "D", "E", "F", "G"))

grafico_multiple$plot <- grafico_multiple$plot + 
  geom_text(data = annotaciones, aes(x = x, y = y, label = label), size = 5, fontface = "bold", color = "black")


for (km in ensayo_11_survfits){
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

print(summary(ensayo_11_controlSac.km))
print(summary(ensayo_11_controlPGBQ.km))
print(summary(ensayo_11_BbQ300.km))
print(summary(ensayo_11_BbQ450.km))
print(summary(ensayo_11_BbQ600.km ))
print(summary(ensayo_11_BbQ750.km ))

print(summary(ensayo_11_controlSac.km)$table["median"])
print(summary(ensayo_11_controlPGBQ.km)$table["median"])
print(summary(ensayo_11_BbQ300.km)$table["median"])
print(summary(ensayo_11_BbQ450.km)$table["median"])
print(summary(ensayo_11_BbQ600.km)$table["median"])
print(summary(ensayo_11_BbQ750.km)$table["median"])

#Test de Cox-Matlan

df_controlSac <- data.frame(tiempo = ensayo_11$tiempo_dias,
                            muertas = ensayo_11$X.Muertas,
                            tratamiento = "ControlSac")

df_controlPGBQ <- data.frame(tiempo = ensayo_11$tiempo_dias.1,
                             muertas = ensayo_11$X.Muertas.1,
                             tratamiento = "ControlPGBQ")

df_Bb150 <- data.frame(tiempo = ensayo_11$tiempo_dias.2,
                       muertas = ensayo_11$X.Muertas.2,
                       tratamiento = "Bb150")

df_Bb300 <- data.frame(tiempo = ensayo_11$tiempo_dias.3,
                       muertas = ensayo_11$X.Muertas.3,
                       tratamiento = "Bb300")

df_Bb450 <- data.frame(tiempo = ensayo_11$tiempo_dias.4,
                       muertas = ensayo_11$X.Muertas.4,
                       tratamiento = "Bb450")

df_Bb600 <- data.frame(tiempo = ensayo_11$tiempo_dias.5,
                       muertas = ensayo_11$X.Muertas.5,
                       tratamiento = "Bb650")

df_Bb750 <- data.frame(tiempo = ensayo_11$tiempo_dias.6,
                       muertas = ensayo_11$X.Muertas.6,
                       tratamiento = "Bb750")

df_combinado_todos <- rbind(df_controlPGBQ, df_controlSac, df_Bb150, df_Bb300, df_Bb450, df_Bb600, df_Bb750)
df_combinado_1 <- rbind(df_controlSac, df_controlPGBQ)
df_combinado_2 <- rbind(df_Bb150, df_Bb300)
df_combinado_3 <- rbind(df_Bb300, df_Bb450)
df_combinado_4 <- rbind(df_Bb450, df_Bb600)
df_combinado_5 <- rbind(df_Bb600, df_Bb750)
df_combinado_6 <- rbind(df_controlPGBQ, df_Bb150)

print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_todos))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_1))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_2))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_3))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_4))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_5))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_6))