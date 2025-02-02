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


#Gráfico
grafico_multiple <- 
  ggsurvplot(ensayo_9_survfits,
             data = df_ensayo_9,
             combine = TRUE,
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
             lwd=1.5,
             pval = FALSE,
             risk.table = FALSE,
             palette = c("seagreen1","darkblue", "violet", "darkolivegreen3", "coral2", "cyan", "deeppink2"),
             ggtheme=theme_bw()+theme(plot.title=element_text(hjust=0.5)))+
  labs(
    x = "Días",
    y = "Supervivencia = S(t)", 
    title = "Curva de Supervivencia - Ensayo #9")

# Añadir letras identificatorias en cada curva
annotaciones <- data.frame(x = c(7.1, 7.1, 7.1, 7.1, 7.1, 7.1, 7.1), # Ajusta las coordenadas de x e y
                           y = c(0.85, 0.75, 0.4, 0.55, 0.20, 0.28, 0.60),
                           label = c("A", "B", "C", "D", "E", "F", "G"))

grafico_multiple$plot <- grafico_multiple$plot + 
  geom_text(data = annotaciones, aes(x = x, y = y, label = label), size = 5, fontface = "bold", color = "black")


for (km in ensayo_9_survfits){
  surv_median <- as.vector(summary(km)$table["median"])
  df <- data.frame(x1 = surv_median, x2=  surv_median,
                   y1 = rep(0, length(surv_median)), y2 = rep(0.5, length(surv_median)))
  
  grafico_multiple$plot <- grafico_multiple$plot + 
    geom_segment(aes(x = 0, y = 0.5, xend = max(surv_median), yend = 0.5),
                 linetype = "dashed", col = "darkblue", size = 0.75)+ # horizontal segment
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df,
                 linetype = "dashed", col = "darkblue", size = 0.75) # vertical segments
}

surv_median <- as.vector(summary(ensayo_9_survfits$controlPGBQ)$table["median"])
df <- data.frame(x1 = surv_median, x2=  surv_median,
                 y1 = rep(0, length(surv_median)), y2 = rep(0.5, length(surv_median)))

grafico_multiple$plot <- grafico_multiple$plot + 
  geom_segment(aes(x = 0, y = 0.5, xend = max(surv_median), yend = 0.5),
               linetype = "dashed", col = "darkblue", size = 0.75)+ # horizontal segment
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df,
               linetype = "dashed", col = "darkblue", size = 0.75) # vertical segments

print(grafico_multiple)

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

print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_todos))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_1))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_2))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_3))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_4))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_5))
print(survdiff(Surv(tiempo, muertas) ~ tratamiento, data = df_combinado_6))

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