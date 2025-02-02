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


#Modelo de Kaplan-Meier
ensayo_2_controlAgua.km <- survfit(Surv(tiempoControlAgua,controlAgua) ~ 1, data = df_ensayo_2, type = "kaplan-meier")
ensayo_2_controlSac.km <- survfit(Surv(tiempoControlSac,controlSac) ~ 1, data = df_ensayo_2, type = "kaplan-meier")
ensayo_2_Bb100.km <- survfit(Surv(tiempoBb100,Bb100) ~ 1, data = df_ensayo_2, type = "kaplan-meier")
ensayo_2_BbQ100.km <- survfit(Surv(tiempoBbQ100,BbQ100) ~ 1, data = df_ensayo_2, type = "kaplan-meier")

ensayo_2_survfits <- list(controlAgua = ensayo_2_controlAgua.km,
                          controlSac = ensayo_2_controlSac.km,
                          Bb100 = ensayo_2_Bb100.km,
                          BbQ100 = ensayo_2_BbQ100.km)



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

#Gráfico
grafico_multiple <- 
  ggsurvplot(ensayo_2_survfits,
             data = df_ensayo_2,
             combine = TRUE,
             font.legend = c(18, "bold", "black"),
             legend.labs = c("A: Control  Sacarosa 25%",
                             "B: Control Agua",
                             "C: Bb 100%",
                             "D: Bb+Q 100%"),
             legend.title = "",
             font.main = c(20, "bold", "darkblue"),
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             font.tickslab = c(18, "plain", "black"),
             conf.int = FALSE,
             conf.int.alpha = c(0.1),
             lwd = 2,
             pval = FALSE,
             risk.table = FALSE,
             palette = c("seagreen1","darkblue", "coral1", "deeppink2"),
             ggtheme=theme_bw()+theme(plot.title=element_text(hjust=0.5)))+
  labs(
    x = "Días",
    y = "Supervivencia = S(t)", 
    title = "Curva de Supervivencia - Ensayo #2")


# Añadir letras identificatorias en cada curva
annotaciones <- data.frame(x = c(10, 10, 10, 10), # Ajusta las coordenadas de x e y
                           y = c(0.9, 0.95, 0.85, 0.8),
                           label = c("A", "B", "C", "D"))

grafico_multiple$plot <- grafico_multiple$plot + 
  geom_text(data = annotaciones, aes(x = x, y = y, label = label), size = 5, fontface = "bold", color = "black")

for (km in ensayo_2_survfits){
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