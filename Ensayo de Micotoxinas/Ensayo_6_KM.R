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


#Modelo de Kaplan-Meier
ensayo_6_controlSac.km <- survfit(Surv(tiempoControl,control) ~ 1, data = df_ensayo_6, type = "kaplan-meier", conf.type = "log-log")
ensayo_6_BbQ300.km <- survfit(Surv(tiempoBbQ300,BbQ300) ~ 1, data = df_ensayo_6, type = "kaplan-meier", conf.type = "log-log")
ensayo_6_BbQ450.km <- survfit(Surv(tiempoBbQ450,BbQ450) ~ 1, data = df_ensayo_6, type = "kaplan-meier", conf.type = "log-log")
ensayo_6_BbQ600.km <- survfit(Surv(tiempoBbQ600,BbQ600) ~ 1, data = df_ensayo_6, type = "kaplan-meier", conf.type = "log-log")
ensayo_6_BbQ750.km <- survfit(Surv(tiempoBbQ750,BbQ750) ~ 1, data = df_ensayo_6, type = "kaplan-meier", conf.type = "log-log")

ensayo_6_survfits <- list(control = ensayo_6_controlSac.km,
                          BbQ300 = ensayo_6_BbQ300.km,
                          BbQ450 = ensayo_6_BbQ450.km,
                          BbQ600 = ensayo_6_BbQ600.km,
                          BbQ750 = ensayo_6_BbQ750.km)

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

#Gráfico
grafico_multiple <- 
  ggsurvplot(ensayo_6_survfits,
             data = df_ensayo_6 ,
             combine = TRUE,
             legend.title = "",
             font.legend = c(19, "bold", "black"),
             font.main = c(20, "bold", "darkblue"),
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             font.tickslab = c(18, "plain", "black"),
             legend.labs = c("A: Control Sacarosa 25%m/V",
                             "B: BbQ 300 µL",
                             "C: BbQ 450 µL",
                             "D: BbQ 600 µL",
                             "E: BbQ 750 µL"),
             conf.int = TRUE,
             conf.int.alpha = c(0.1),
             lwd = 2,
             pval = FALSE,
             risk.table = FALSE,
             palette = c("seagreen1", "coral1", "darkorchid3", "cadetblue4", "brown4"),
             ggtheme=theme_bw()+theme(plot.title=element_text(hjust=0.5)))+
  labs(
    x = "Días",
    y = "Supervivencia = S(t)", 
    title = "Curva de Supervivencia - Ensayo #6")

# Añadir letras identificatorias en cada curva
annotaciones <- data.frame(x = c(13, 12, 13.2, 11.5, 13.2), # Ajusta las coordenadas de x e y
                           y = c(1, 0.4, 0.0, 0.13, 0.05),
                           label = c("A", "B", "C", "D", "E"))

grafico_multiple$plot <- grafico_multiple$plot + 
  geom_text(data = annotaciones, aes(x = x, y = y, label = label), size = 5, fontface = "bold", color = "black")


for (km in ensayo_6_survfits){
  surv_median <- as.vector(summary(km)$table["median"])
  df <- data.frame(x1 = surv_median, x2=  surv_median,
                   y1 = rep(0, length(surv_median)), y2 = rep(0.5, length(surv_median)))
  
  grafico_multiple$plot <- grafico_multiple$plot + 
    geom_segment(aes(x = 0, y = 0.5, xend = max(surv_median), yend = 0.5),
                 linetype = "dashed", col = "darkblue", size = 0.75)+ # horizontal segment
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df,
                 linetype = "dashed", col = "darkblue", size = 0.75) # vertical segments
}

surv_median <- as.vector(summary(ensayo_6_survfits$BbQ300)$table["median"])
df <- data.frame(x1 = surv_median, x2=  surv_median,
                 y1 = rep(0, length(surv_median)), y2 = rep(0.5, length(surv_median)))

grafico_multiple$plot <- grafico_multiple$plot + 
  geom_segment(aes(x = 0, y = 0.5, xend = max(surv_median), yend = 0.5),
               linetype = "dashed", col = "darkblue", size = 0.75)+ # horizontal segment
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df,
               linetype = "dashed", col = "darkblue", size = 0.75) # vertical segments

print(grafico_multiple)

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
df_combinado_1 <- rbind(df_control, df_Bb300)
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