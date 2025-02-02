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
ensayo_1 <- read.csv(file = "Ensayos de Micotoxinas de B. bassiana en L. humile - Ensayo #1 - Para analisis R.csv", header = T, sep=",", skip=1)
#head(ensayo_1) #Para verificar

df_ensayo_1 <- data.frame(tiempoControl= ensayo_1$tiempo_dias,
                          control = ensayo_1$X.Muertas,
                          tiempoBb100= ensayo_1$tiempo_dias.1,
                          Bb100 = ensayo_1$X.Muertas.1,
                          tiempoBb75= ensayo_1$tiempo_dias.2,
                          Bb75 = ensayo_1$X.Muertas.2,
                          tiempoBb50= ensayo_1$tiempo_dias.3,
                          Bb50 = ensayo_1$X.Muertas.3,
                          tiempoBb25= ensayo_1$tiempo_dias.4,
                          Bb25 = ensayo_1$X.Muertas.4,
                          tiempoBbQ100= ensayo_1$tiempo_dias.5,
                          BbQ100 = ensayo_1$X.Muertas.5,
                          tiempoBbQ75= ensayo_1$tiempo_dias.6,
                          BbQ75 = ensayo_1$X.Muertas.6,
                          tiempoBbQ50= ensayo_1$tiempo_dias.7,
                          BbQ50 = ensayo_1$X.Muertas.7,
                          tiempoBbQ25= ensayo_1$tiempo_dias.8,
                          BbQ25 = ensayo_1$X.Muertas.8)

#Modelo de Kaplan-Meier
ensayo_1_controlSac.km <- survfit(Surv(tiempoControl,control) ~ 1, data = df_ensayo_1, type = "kaplan-meier")
ensayo_1_Bb100.km <- survfit(Surv(tiempoBb100,Bb100) ~ 1, data = df_ensayo_1, type = "kaplan-meier")
ensayo_1_Bb75.km <- survfit(Surv(tiempoBb75,Bb75) ~ 1, data = df_ensayo_1, type = "kaplan-meier")
ensayo_1_Bb50.km <- survfit(Surv(tiempoBb50,Bb50) ~ 1, data = df_ensayo_1, type = "kaplan-meier")
ensayo_1_Bb25.km <- survfit(Surv(tiempoBb25,Bb25) ~ 1, data = df_ensayo_1, type = "kaplan-meier")
ensayo_1_BbQ100.km <- survfit(Surv(tiempoBbQ100,BbQ100) ~ 1, data = df_ensayo_1, type = "kaplan-meier")
ensayo_1_BbQ75.km <- survfit(Surv(tiempoBbQ75,BbQ75) ~ 1, data = df_ensayo_1, type = "kaplan-meier")
ensayo_1_BbQ50.km <- survfit(Surv(tiempoBbQ50,BbQ50) ~ 1, data = df_ensayo_1, type = "kaplan-meier")
ensayo_1_BbQ25.km <- survfit(Surv(tiempoBbQ25,BbQ25) ~ 1, data = df_ensayo_1, type = "kaplan-meier")

ensayo_1_survfits <- list(control = ensayo_1_controlSac.km,
                          Bb100 = ensayo_1_Bb100.km,
                          Bb75 = ensayo_1_Bb75.km,
                          Bb50 = ensayo_1_Bb50.km,
                          Bb25 = ensayo_1_Bb25.km,
                          BbQ100 = ensayo_1_BbQ100.km,
                          BbQ75 = ensayo_1_BbQ75.km,
                          BbQ50 = ensayo_1_BbQ50.km,
                          BbQ25 = ensayo_1_BbQ25.km)

library(ggpubr)

# Prueba de Shapiro-Wilk para cada grupo de tratamiento
shapiro_test_control <- shapiro.test(df_ensayo_1$control)
shapiro_test_Bb100 <- shapiro.test(df_ensayo_1$Bb100)
shapiro_test_Bb75 <- shapiro.test(df_ensayo_1$Bb75)
shapiro_test_Bb50 <- shapiro.test(df_ensayo_1$Bb50)
shapiro_test_Bb25 <- shapiro.test(df_ensayo_1$Bb25)
shapiro_test_BbQ100 <- shapiro.test(df_ensayo_1$BbQ100)
shapiro_test_BbQ75 <- shapiro.test(df_ensayo_1$BbQ75)
shapiro_test_BbQ50 <- shapiro.test(df_ensayo_1$BbQ50)
shapiro_test_BbQ25 <- shapiro.test(df_ensayo_1$BbQ25)

# Imprimir los resultados
print(shapiro_test_control)
print(shapiro_test_Bb100)
print(shapiro_test_Bb75)
print(shapiro_test_Bb50)
print(shapiro_test_Bb25)
print(shapiro_test_BbQ100)
print(shapiro_test_BbQ75)
print(shapiro_test_BbQ50)
print(shapiro_test_BbQ25)

# Prueba de Levene

library(car)

levene_test <- leveneTest(
  c(controlSac = df_ensayo_1$control,
    Bb100 = df_ensayo_1$Bb100,
    Bb75 = df_ensayo_1$Bb75,
    Bb50 = df_ensayo_1$Bb50,
    Bb25 = df_ensayo_1$Bb25,
    BbQ100 = df_ensayo_1$BbQ100,
    BbQ75 = df_ensayo_1$BbQ75,
    BbQ50 = df_ensayo_1$BbQ50,
    BbQ25 = df_ensayo_1$BbQ25) ~ 
    factor(rep(c("controlSac", "Bb100", "Bb75", "Bb50", "Bb25", "BbQ100", "BbQ75", "BbQ50", "BbQ25"),
               each = nrow(df_ensayo_1))),
  center = median
  )

print(levene_test)

#Gráfico
grafico_multiple <- 
  ggsurvplot(ensayo_1_survfits,
             data = df_ensayo_1,
             combine = TRUE,
             font.legend = c(16, "bold", "black"),
             legend.labs = c("A: Control Sacarosa 25%", 
                             "B: Bb 100%", 
                             "C: Bb:Sacarosa 75:25", 
                             "D: Bb:Sacarosa 50:50", 
                             "E: Bb:Sacarosa 25:75", 
                             "F: Bb+Q 100%", 
                             "G: Bb+Q:Sacarosa 75:25", 
                             "H: Bb+Q:Sacarosa 50:50", 
                             "I: Bb+Q:Sacarosa 25:75"),
             legend.title = "",
             font.main = c(20, "bold", "darkblue"),
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             font.tickslab = c(18, "plain", "black"),
             conf.int = FALSE,
             conf.int.alpha = 0.1,
             lwd = 2.5,
             pval = FALSE,
             risk.table = FALSE,
             palette = c("seagreen1","darkblue", "violet", "darkolivegreen3", 
                         "green3", "chocolate3", "coral1", "deeppink2", "cyan"),
             ggtheme = theme_bw() + theme(plot.title = element_text(hjust = 0.5))) +
             labs(x = "Días", y = "Supervivencia = S(t)", title = "Curva de Supervivencia - Ensayo #1")

# Añadir letras identificatorias en cada curva
annotaciones <- data.frame(x = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5), # Ajusta las coordenadas de x e y
                           y = c(0.2, 0.25, 0.3, 0.5, 0.45, 0.12, 0.57, 0.63, 0.4),
                           label = c("A", "B", "C", "D", "E", "F", "G", "H", "I"))

grafico_multiple$plot <- grafico_multiple$plot + 
  geom_text(data = annotaciones, aes(x = x, y = y, label = label), size = 5, fontface = "bold", color = "black")

for (km in ensayo_1_survfits){
  surv_median <- as.vector(summary(km)$table["median"])
  df <- data.frame(x1 = surv_median, x2=  surv_median,
                   y1 = rep(0, length(surv_median)), y2 = rep(0.5, length(surv_median)))
  
  grafico_multiple$plot <- grafico_multiple$plot + 
    geom_segment(aes(x = 0, y = 0.5, xend = max(surv_median), yend = 0.5),
                 linetype = "dashed", col = "darkblue", size = 0.75)+ # horizontal segment
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df,
                 linetype = "dashed", col = "darkblue", size = 0.75) # vertical segments
}

surv_median <- as.vector(summary(ensayo_1_survfits$BbQ50)$table["median"])
df <- data.frame(x1 = surv_median, x2=  surv_median,
                 y1 = rep(0, length(surv_median)), y2 = rep(0.5, length(surv_median)))

grafico_multiple$plot <- grafico_multiple$plot + 
  geom_segment(aes(x = 0, y = 0.5, xend = max(surv_median), yend = 0.5),
               linetype = "dashed", col = "darkblue", size = 0.75)+ # horizontal segment
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df,
               linetype = "dashed", col = "darkblue", size = 0.75) # vertical segments

print(grafico_multiple)

#Estadísticos de cada tratamiento

print(summary(ensayo_1_controlSac.km))
print(summary(ensayo_1_Bb100.km))
print(summary(ensayo_1_Bb75.km))
print(summary(ensayo_1_Bb50.km))
print(summary(ensayo_1_Bb25.km))
print(summary(ensayo_1_BbQ100.km))
print(summary(ensayo_1_BbQ75.km))
print(summary(ensayo_1_BbQ50.km))
print(summary(ensayo_1_BbQ25.km))

print(summary(ensayo_1_controlSac.km)$table["median"])
print(summary(ensayo_1_Bb100.km)$table["median"])
print(summary(ensayo_1_Bb75.km)$table["median"])
print(summary(ensayo_1_Bb50.km)$table["median"])
print(summary(ensayo_1_Bb25.km)$table["median"])
print(summary(ensayo_1_BbQ100.km)$table["median"])
print(summary(ensayo_1_BbQ75.km)$table["median"])
print(summary(ensayo_1_BbQ50.km)$table["median"])
print(summary(ensayo_1_BbQ25.km)$table["median"])


#Gráfica de Supervivencia
#png("Ensayo-1-Supervivencia.png",width = 5, height = 5, units = 'in', res = 600,pointsize=7.5)
#par(xpd=T, mar=par()$mar+c(0,0,0,11))
#plot(ensayo_1_control.km, xlab="Tiempo (días)", ylab="Supervivencia = S(t)", main="Curva de Supervivencia - Ensayo #1", conf.int = F, col = "seagreen1", lwd=2, mark.time = TRUE, las = 1)
#lines(ensayo_1_Bb100.km, conf.int = F, col = "darkblue" , lwd = 2, mark.time = TRUE)
#lines(ensayo_1_Bb75.km, conf.int = F, col = "violet" , lwd = 2, mark.time = TRUE)
#lines(ensayo_1_Bb50.km, conf.int = F, col = "darkolivegreen3" , lwd = 2, mark.time = TRUE)
#lines(ensayo_1_Bb25.km, conf.int = F, col = "green3" , lwd = 2, mark.time = TRUE)
#lines(ensayo_1_BbQ100.km, conf.int = F, col = "chocolate3" , lwd = 2, mark.time = TRUE)
#lines(ensayo_1_BbQ75.km, conf.int = F, col = "coral1" , lwd = 2, mark.time = TRUE)
#lines(ensayo_1_BbQ50.km, conf.int = F, col = "deeppink2" , lwd = 2, mark.time = TRUE)
#lines(ensayo_1_BbQ25.km, conf.int = F, col = "cyan" , lwd = 2, mark.time = TRUE)
#legend(3.5, 0.75, legend=c("Control  Sacarosa 25%", "Bb 100%", "Bb:Sacarosa 75:20", "Bb:Sacarosa 50:50", "Bb:Sacarosa 25:70", "Bb+Q 100%", "Bb+Q:Sacarosa 75:25", "Bb+Q:Sacarosa 50:50", "Bb+Q:Sacarosa 25:75"), lty = 2, lwd = 3, col =c ("seagreen1","darkblue", "violet", "darkolivegreen3", "green3", "chocolate3", "coral1", "deeppink2", "cyan"), bty ="", cex=0.8)
#dev.off(3)