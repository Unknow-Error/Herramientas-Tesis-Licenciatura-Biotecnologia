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

#Modelo de Kaplan-Meier
ensayo_4_controlSac.km <- survfit(Surv(tiempoControlSac,controlSac) ~ 1, data = df_ensayo_4, type = "kaplan-meier")
ensayo_4_BbQ1.km <- survfit(Surv(tiempoBbQ1,BbQ1) ~ 1, data = df_ensayo_4, type = "kaplan-meier")
ensayo_4_BbQ2.km <- survfit(Surv(tiempoBbQ2,BbQ2) ~ 1, data = df_ensayo_4, type = "kaplan-meier")

ensayo_4_survfits <- list(controlSac = ensayo_4_controlSac.km,
                          BbQ1 = ensayo_4_BbQ1.km,
                          BbQ2 = ensayo_4_BbQ2.km)



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

#Gráfico
grafico_multiple <- 
  ggsurvplot(ensayo_4_survfits,
             data = df_ensayo_4,
             combine = TRUE,
             legend.title = "",
             font.legend = c(19, "bold", "black"),
             font.main = c(20, "bold", "darkblue"),
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             font.tickslab = c(18, "plain", "black"),
             legend.labs = c("A: Control Sacarosa 25%m/V", "B: BbQ N°1", "C: BbQ N°2"),
             conf.int = TRUE,
             conf.int.alpha = c(0.2),
             lwd = 2.5,
             pval = FALSE,
             risk.table = FALSE,
             palette = c("seagreen1", "coral1", "cadetblue4"),
             ggtheme=theme_bw()+theme(plot.title=element_text(hjust=0.5)))+
  labs(
    x = "Días",
    y = "Supervivencia = S(t)", 
    title = "Curva de Supervivencia - Ensayo #4")

# Añadir letras identificatorias en cada curva
annotaciones <- data.frame(x = c(18, 18, 18), # Ajusta las coordenadas de x e y
                           y = c(0.8, 1, 0.0),
                           label = c("A", "B", "C"))

grafico_multiple$plot <- grafico_multiple$plot + 
  geom_text(data = annotaciones, aes(x = x, y = y, label = label), size = 5, fontface = "bold", color = "black")


for (km in ensayo_4_survfits){
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