# Install packages if needed
# install.packages(c("knitr", "dplyr", "survival", "ggplot2", "here", "tibble"))
library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)

# devtools::install_github("zabore/ezfun")
#ezfun::set_ccf_palette("contrast")

# install.packages(c("lubridate", "ggsurvfit", "gtsummary", "tidycmprsk"))
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)

# devtools::install_github("zabore/condsurv")
library(condsurv)

#Data de ejemplo

lung <- 
  lung %>% 
  mutate(
    status = recode(status, `1` = 0, `2` = 1)
  )

head(lung[, c("time", "status", "sex")])

#Calcular Tiempos supervivencia

date_ex <- 
  tibble(
    sx_date = c("2007-06-22", "2004-02-13", "2010-10-27"), 
    last_fup_date = c("2017-04-15", "2018-07-04", "2016-10-31")
  )

date_ex

date_ex <-
  date_ex %>% 
  mutate(
    sx_date = ymd(sx_date), 
    last_fup_date = ymd(last_fup_date)
  )

date_ex

#Crear Curvas de supervivencia

Surv(lung$time, lung$status)[1:10]
s1 <- survfit(Surv(time, status) ~ 1, data = lung)
str(s1)

#Curvas Kaplan-Meier

gg_default <-
  survfit2(Surv(lung$time, lung$status) ~ 1, data = lung) %>%
  ggsurvfit() +
  add_confidence_interval() +
  scale_ggsurvfit() +
  labs(title = "Default")
#plot(s1, xlab="Tiempo (días)", ylab="Supervivencia = S(t)", main="Curva de Supervivencia ", conf.int = F, col = "seagreen1", lwd=2, mark.time = TRUE, las = 1)

