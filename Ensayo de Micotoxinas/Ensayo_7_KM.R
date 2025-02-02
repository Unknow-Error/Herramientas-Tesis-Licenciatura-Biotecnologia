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


#Gráfico
grafico_multiple <- 
  ggsurvplot(ensayo_7_survfits,
             data = df_ensayo_7 ,
             combine = TRUE,
             legend.title = "",
             font.legend = c(19, "bold", "black"),
             font.main = c(20, "bold", "darkblue"),
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             font.tickslab = c(18, "plain", "black"),
             legend.labs = c("A: Control Sacarosa 25%m/V - Bernal",
                             "B: Control Sacarosa 25%m/V - Gonnet",
                             "C: BbQ 450 uL Bernal",
                             "D: BbQ 450 uL Gonnet"),
             conf.int = TRUE,
             conf.int.alpha = c(0.15),
             lwd = 1.5,
             pval = FALSE,
             risk.table = FALSE,
             palette = c("seagreen1", "coral1", "cadetblue4", "darkorchid4"),
             ggtheme=theme_bw()+theme(plot.title=element_text(hjust=0.5)))+
  labs(
    x = "Días",
    y = "Supervivencia = S(t)", 
    title = "Curva de Supervivencia - Ensayo #7")

# Añadir letras identificatorias en cada curva
annotaciones <- data.frame(x = c(8, 8, 8.1, 8.1), # Ajusta las coordenadas de x e y
                           y = c(0.87, 0.8, 0.15, 0.0),
                           label = c("A", "B", "C", "D"))

grafico_multiple$plot <- grafico_multiple$plot + 
  geom_text(data = annotaciones, aes(x = x, y = y, label = label), size = 5, fontface = "bold", color = "black")


for (km in ensayo_7_survfits){
  surv_median <- as.vector(summary(km)$table["median"])
  df <- data.frame(x1 = surv_median, x2=  surv_median,
                   y1 = rep(0, length(surv_median)), y2 = rep(0.5, length(surv_median)))
  
  grafico_multiple$plot <- grafico_multiple$plot + 
    geom_segment(aes(x = 0, y = 0.5, xend = max(surv_median), yend = 0.5),
                 linetype = "dashed", col = "darkblue", size = 0.75)+ # horizontal segment
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df,
                 linetype = "dashed", col = "darkblue", size = 0.75) # vertical segments
}

surv_median <- as.vector(summary(ensayo_7_survfits$BbQ450Patri)$table["median"])
df <- data.frame(x1 = surv_median, x2=  surv_median,
                 y1 = rep(0, length(surv_median)), y2 = rep(0.5, length(surv_median)))

grafico_multiple$plot <- grafico_multiple$plot + 
  geom_segment(aes(x = 0, y = 0.5, xend = max(surv_median), yend = 0.5),
               linetype = "dashed", col = "darkblue", size = 0.75)+ # horizontal segment
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df,
               linetype = "dashed", col = "darkblue", size = 0.75) # vertical segments

print(grafico_multiple)

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
