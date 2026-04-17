rm(list = setdiff(ls(), "global_seed"))

require(knitr)
require(survival)
require(dplyr)
require(kableExtra)
require(flextable)
require(ggplot2)
require(ggsci)
require(ggthemes)
require(gridExtra)
require(paletteer)
require(Publish)
require(simsurv)
require(splines)
require(survminer)
require(kableExtra)
require(tidyr)
require(dplyr)
require(tibble)
require(car)

path_cox <- "Your path" # Path to the file where figures will be saved, needs to finish with a /

# Parameters ----

# seed
if (!exists("global_seed")) global_seed <- 1024

# parameters of simulation
N <- 1000

shape_1 <- 2
shape_2 <- 3 # Needs different shape otherwise we just create new HRs

median_surv_strate1 <- 50
median_surv_strate2 <- 20
scale1 <- log(2)/(median_surv_strate1^shape_1)
scale2 <- log(2)/(median_surv_strate2^shape_2)

beta_with.effect <- 1 # we
beta_no.effect <- 0 # ne

lambda_cens <- 0.023 # log(2)/30

xaxis <- seq(0, 80, 0.01)
times_print <- seq(0, max(xaxis), by = 10)

# Data generation ----

set.seed(global_seed)
U <- runif(N)

strate <- rbinom(N, 1, 0.5) + 1
X_with.effect <- rbinom(N, 1, 0.5)
X_no.effect <- rnorm(N, 0, 3)

LP_both <- beta_no.effect * X_no.effect + beta_with.effect * X_with.effect

#####################################'

S_true_fun_with.effect_strate1 <- function(t) {exp(-scale1 * exp(beta_with.effect) * t^shape_1)}
S_true_with.effect_strate1 <- S_true_fun_with.effect_strate1(xaxis)
S_true_fun_no.effect_strate1 <- function(t) {exp(-scale1 * 1 * t^shape_1)}
S_true_no.effect_strate1 <- S_true_fun_no.effect_strate1(xaxis)

S_true_fun_with.effect_strate2 <- function(t) {exp(-scale2 * exp(beta_with.effect) * t^shape_2)}
S_true_with.effect_strate2 <- S_true_fun_with.effect_strate2(xaxis)
S_true_fun_no.effect_strate2 <- function(t) {exp(-scale2 * 1 * t^shape_2)}
S_true_no.effect_strate2 <- S_true_fun_no.effect_strate2(xaxis)

#####################################'

surv_time_both <- numeric(N)
for(i in 1:N) {
  surv_time_both[i] <- ((-log(U[i])) / (ifelse(strate[i]==1, scale1, scale2) * exp(LP_both[i])))^(1/ifelse(strate[i]==1, shape_1, shape_2))
}

censor_time_both <- rexp(N, lambda_cens)

time_both <- pmin(surv_time_both, censor_time_both)
status_both <- as.numeric(surv_time_both <= censor_time_both)
data_both <- data.frame(time = time_both,
                        status = status_both,
                        strate = factor(strate),
                        X_no.effect = X_no.effect,
                        X_with.effect = X_with.effect)

cox_strate <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect + strata(strate), data=data_both)
cox_Nostrate <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect, data=data_both)
cox_triple <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect + strate, data=data_both)

# publish
invisible(capture.output(pub_cox_Nostrate <- publish(cox_Nostrate)))
invisible(capture.output(pub_cox_strate <- publish(cox_strate)))

# HR
HR_Nostrate <- exp(coef(cox_Nostrate))
HR_strate <- exp(coef(cox_strate))

# Comparison of results ----
# TABLE: 6
comparison <- data.frame(
  Model = c("NE", "WE", "NE", "WE"),
  HR_vrai = c(exp(beta_no.effect), exp(beta_with.effect), exp(beta_no.effect), exp(beta_with.effect)),
  HR_estime = c(HR_Nostrate[1], HR_Nostrate[2], HR_strate[1], HR_strate[2]),
  IC_inf = c(round(pub_cox_Nostrate$rawTable$Lower[1],4),
             round(pub_cox_Nostrate$rawTable$Lower[2],4),
             round(pub_cox_strate$rawTable$Lower[1],4),
             round(pub_cox_strate$rawTable$Lower[2],4)),
  IC_sup = c(round(pub_cox_Nostrate$rawTable$Upper[1],4),
             round(pub_cox_Nostrate$rawTable$Upper[2],4),
             round(pub_cox_strate$rawTable$Upper[1],4),
             round(pub_cox_strate$rawTable$Upper[2],4))
)


comparison <- comparison %>%
  kable(digits = 7, align = "lcc") %>%
  kable_styling(full_width = FALSE) %>%
  pack_rows("Basic model", 1, 2) %>%
  pack_rows("Stratified model", 3, 4)

# COX No strate
# NE = Q3, WE = 1
est_ns_NEQ3_WE1 <- survfit(cox_Nostrate, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=1))
df0 <- data.frame(time = est_ns_NEQ3_WE1$time, S = est_ns_NEQ3_WE1$surv, model = "No strate, NE=Q3 WE=1")
# NE = Q3, WE = 0
est_ns_NEQ3_WE0 <- survfit(cox_Nostrate, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=0))
df0bis <- data.frame(time = est_ns_NEQ3_WE0$time, S = est_ns_NEQ3_WE0$surv, model = "No strate, NE=Q3 WE=0")

# COX strate1
# NE = Q3, WE = 1
est_s1_NEQ3_WE1 <- survfit(cox_strate, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=1, strate=factor(1)))
df1 <- data.frame(time = est_s1_NEQ3_WE1$time, S = est_s1_NEQ3_WE1$surv, model = "Strate1, NE=Q3 WE=1")
# NE = Q3, WE = 0
est_s1_NEQ3_WE0 <- survfit(cox_strate, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=0, strate=factor(1)))
df1bis <- data.frame(time = est_s1_NEQ3_WE0$time, S = est_s1_NEQ3_WE0$surv, model = "Strate1, NE=Q3 WE=0")

# COX strate2
# NE = Q3, WE = 1
est_s2_NEQ3_WE1 <- survfit(cox_strate, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=1, strate=factor(2)))
df2 <- data.frame(time = est_s2_NEQ3_WE1$time, S = est_s2_NEQ3_WE1$surv, model = "Strate2, NE=Q3 WE=1")
# NE = Q3, WE = 0
est_s2_NEQ3_WE0 <- survfit(cox_strate, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=0, strate=factor(2)))
df2bis <- data.frame(time = est_s2_NEQ3_WE0$time, S = est_s2_NEQ3_WE0$surv, model = "Strate2, NE=Q3 WE=0")



# Cox triple hr
est_th_WE1_ST1 <- survfit(cox_triple, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), 
                                                           X_with.effect=1, 
                                                           strate=factor(1)))
df1_th <- data.frame(time = est_th_WE1_ST1$time, S = est_th_WE1_ST1$surv, model = "Strate1, NE=Q3 WE=1")
# NE = Q3, WE = 0
est_th_WE0_ST1 <- survfit(cox_triple, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), 
                                                            X_with.effect=0, 
                                                            strate=factor(1)))
df1bis_th <- data.frame(time = est_th_WE0_ST1$time, S = est_th_WE0_ST1$surv, model = "Strate1, NE=Q3 WE=0")

# COX strate2
# NE = Q3, WE = 1
est_th_WE1_ST2 <- survfit(cox_triple, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), 
                                                           X_with.effect=1, 
                                                           strate=factor(2)))
df2_th <- data.frame(time = est_th_WE1_ST2$time, S = est_th_WE1_ST2$surv, model = "Strate2, NE=Q3 WE=1")
# NE = Q3, WE = 0
est_th_WE0_ST2 <- survfit(cox_triple, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), 
                                                            X_with.effect=0, 
                                                            strate=factor(2)))
df2bis_th <- data.frame(time = est_th_WE0_ST2$time, S = est_th_WE0_ST2$surv, model = "Strate2, NE=Q3 WE=0")

# one per strata ----

data_strate1 <- data_both[which(data_both$strate==1),]
data_strate2 <- data_both[which(data_both$strate==2),]
cox_strate1 <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect, data=data_strate1)
cox_strate2 <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect, data=data_strate2)

# publish
invisible(capture.output(pub_cox_strate1 <- publish(cox_strate1)))
invisible(capture.output(pub_cox_strate2 <- publish(cox_strate2)))

# HR
HR_strate1 <- exp(coef(cox_strate1))
HR_strate2 <- exp(coef(cox_strate2))

# Comparison of results ----
# TABLE: 7
comparison_PerStrate <- data.frame(
  Modele = c("NE", "WE", "NE", "WE"),
  HR_vrai = c(exp(beta_no.effect), exp(beta_with.effect), exp(beta_no.effect), exp(beta_with.effect)),
  HR_estime = c(HR_strate1[1], HR_strate1[2], HR_strate2[1], HR_strate2[2]),
  IC_inf = c(round(pub_cox_strate1$rawTable$Lower[1],4),
             round(pub_cox_strate1$rawTable$Lower[2],4),
             round(pub_cox_strate2$rawTable$Lower[1],4),
             round(pub_cox_strate2$rawTable$Lower[2],4)),
  IC_sup = c(round(pub_cox_strate1$rawTable$Upper[1],4),
             round(pub_cox_strate1$rawTable$Upper[2],4),
             round(pub_cox_strate2$rawTable$Upper[1],4),
             round(pub_cox_strate2$rawTable$Upper[2],4))
)

# comparison_PerStrate <- comparison %>%
#   kable(digits = 7, align = "lcc") %>%
#   kable_styling(full_width = FALSE) %>%
#   pack_rows("Modèles strate 1", 1, 2) %>%
#   pack_rows("Modèles strate 2", 3, 4)

# COX No strate
# NE = Q3, WE = 1
est_s1_NEQ3_WE1_split <- survfit(cox_strate1, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=1))
df0_split <- data.frame(time = est_s1_NEQ3_WE1$time, S = est_s1_NEQ3_WE1$surv, model = "Strate 1, NE=Q3 WE=1")
# NE = Q3, WE = 0
est_s1_NEQ3_WE0_split <- survfit(cox_strate1, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=0))
df0bis_split <- data.frame(time = est_s1_NEQ3_WE0$time, S = est_s1_NEQ3_WE0$surv, model = "Strate 1, NE=Q3 WE=0")

# COX strate1
# NE = Q3, WE = 1
est_s2_NEQ3_WE1_split <- survfit(cox_strate2, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=1))
df1_split <- data.frame(time = est_s2_NEQ3_WE1$time, S = est_s2_NEQ3_WE1$surv, model = "Strate 2, NE=Q3 WE=1")
# NE = Q3, WE = 0
est_s2_NEQ3_WE0_split <- survfit(cox_strate2, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=0))
df1bis_split <- data.frame(time = est_s2_NEQ3_WE0$time, S = est_s2_NEQ3_WE0$surv, model = "Strate2, NE=Q3 WE=0")

# Figure 6 ----

my_colors <- paletteer_d("rcartocolor::Pastel")

fit <- survfit(cox_strate, type = "aalen")
fit2 <- survfit(cox_triple, type = "aalen")
at_risk <- matrix(c(summary(fit, times_print)$n.risk, 0, 0, 0, 0, 0, summary(fit2, times_print)$n.risk),
                  nrow = 3,
                  byrow = TRUE)

jpeg(paste(path_cox, "FIGURE6.jpeg", sep = ""), width = 6, height = 3.25, units = "in", res = 1200, pointsize = 4)
par(mar = c(9, 7, 5, 2), xaxs = "i", yaxs = "i", tcl = 0.2, cex.lab = 1.4, cex.axis = 1.2, mfrow = c(1, 2))

# a. ###############################################################'

plot(1, type = "n",
     col = NULL,
     ylim = c(0, 1),
     xlim = c(0, 80),
     sub = "",
     xlab = "Time", ylab = "Survival", axes = FALSE,
     panel.first = grid(nx = NULL, ny = NULL,
                        col = rgb(0.8, 0.8, 0.8, 0.7),
                        lty = "dotted", lwd = 1))

title(main = "a.", adj = 0, line = 0.5)

lines(xaxis, S_true_with.effect_strate1, lty = 1, lwd = 1, col = "black")
lines(xaxis, S_true_no.effect_strate1, lty = 3, lwd = 1, col = "black")
lines(xaxis, S_true_with.effect_strate2, lty = 2, lwd = 1, col = "black")
lines(xaxis, S_true_no.effect_strate2, lty = "aa", lwd = 1, col = "black")

polygon(rep(c(est_s2_NEQ3_WE1$time, rev(est_s2_NEQ3_WE1$time)), each = 2)[-1],
        rep(c(est_s2_NEQ3_WE1$upper, rev(est_s2_NEQ3_WE1$lower)), each = 2)[-length(c(est_s2_NEQ3_WE1$upper, rev(est_s2_NEQ3_WE1$lower)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(df2$time, df2$S,
      lwd = 1,
      lty = 2,
      col = my_colors[1],
      type = "s")

polygon(rep(c(est_s2_NEQ3_WE0$time, rev(est_s2_NEQ3_WE0$time)), each = 2)[-1],
        rep(c(est_s2_NEQ3_WE0$upper, rev(est_s2_NEQ3_WE0$lower)), each = 2)[-length(c(est_s2_NEQ3_WE0$upper, rev(est_s2_NEQ3_WE0$lower)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(df2bis$time, df2bis$S,
      lwd = 1,
      lty = "aa",
      col = my_colors[1],
      type = "s")

polygon(rep(c(est_s1_NEQ3_WE1$time, rev(est_s1_NEQ3_WE1$time)), each = 2)[-1],
        rep(c(est_s1_NEQ3_WE1$upper, rev(est_s1_NEQ3_WE1$lower)), each = 2)[-length(c(est_s1_NEQ3_WE1$upper, rev(est_s1_NEQ3_WE1$lower)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(est_s1_NEQ3_WE1$time, est_s1_NEQ3_WE1$surv,
      lwd = 1,
      lty = 1,
      col = my_colors[1],
      type = "s")

polygon(rep(c(est_s1_NEQ3_WE0$time, rev(est_s1_NEQ3_WE0$time)), each = 2)[-1],
        rep(c(est_s1_NEQ3_WE0$upper, rev(est_s1_NEQ3_WE0$lower)), each = 2)[-length(c(est_s1_NEQ3_WE0$upper, rev(est_s1_NEQ3_WE0$lower)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(est_s1_NEQ3_WE0$time, est_s1_NEQ3_WE0$surv,
      lwd = 1,
      lty = 3,
      col = my_colors[1],
      type = "s")



lines(df2_th$time, df2_th$S,
      lwd = 1,
      lty = 2,
      col = my_colors[2],
      type = "s")


lines(df2bis_th$time, df2bis_th$S,
      lwd = 1,
      lty = "aa",
      col = my_colors[2],
      type = "s")


lines(df1_th$time, df1_th$S,
      lwd = 1,
      lty = 1,
      col = my_colors[2],
      type = "s")


lines(df1bis_th$time, df1bis_th$S,
      lwd = 1,
      lty = 3,
      col = my_colors[2],
      type = "s")


axis(side = 1, at = times_print, lwd = 0.7)
axis(side = 2, at = seq(0, 1, 0.1), lwd = 0.7)

mtext(side = 1, line = 4, at = times_print[1], text = expression(italic("Number at risk")), adj = 0, cex = 1.4)

mtext(side = 1, line = 5, at = times_print[1] - diff(times_print)[1] * 0.5,
      text = "Strata Z=1", adj = 1, cex = 1.2)
mtext(side = 1, line = 6, at = times_print[1] - diff(times_print)[1] * 0.5,
      text = "Strata Z=2", adj = 1, cex = 1.2)
mtext(side = 1, line = 7, at = times_print[1] - diff(times_print)[1] * 0.5,
      text = "Naive", adj = 1, cex = 1.2)
for (j in 1:nrow(at_risk))
{
  for(i in seq_along(times_print)){
    mtext(side = 1, line = 5+j-1, at = times_print[i],
          text = at_risk[j, i], cex = 1.2)
  }
}

# b. ###############################################################'

plot(1, type = "n",
     col = NULL,
     ylim = c(0, 1),
     xlim = c(0, 80),
     sub = "",
     xlab = "Time", ylab = "Survival", axes = FALSE,
     panel.first = grid(nx = NULL, ny = NULL,
                        col = rgb(0.8, 0.8, 0.8, 0.7),
                        lty = "dotted", lwd = 1))

title(main = "b.", adj = 0, line = 0.5)

lines(xaxis, S_true_with.effect_strate1, lty = 1, lwd = 1, col = "black")
lines(xaxis, S_true_no.effect_strate1, lty = 3, lwd = 1, col = "black")
lines(xaxis, S_true_with.effect_strate2, lty = 2, lwd = 1, col = "black")
lines(xaxis, S_true_no.effect_strate2, lty = "aa", lwd = 1, col = "black")

polygon(rep(c(est_s1_NEQ3_WE1_split$time, rev(est_s1_NEQ3_WE1_split$time)), each = 2)[-1],
        rep(c(est_s1_NEQ3_WE1_split$upper, rev(est_s1_NEQ3_WE1_split$lower)), each = 2)[-length(c(est_s1_NEQ3_WE1_split$upper, rev(est_s1_NEQ3_WE1_split$lower)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(est_s1_NEQ3_WE1_split$time, est_s1_NEQ3_WE1_split$surv,
      lwd = 1,
      lty = 1,
      col = my_colors[1],
      type = "s")

polygon(rep(c(est_s1_NEQ3_WE0_split$time, rev(est_s1_NEQ3_WE0_split$time)), each = 2)[-1],
        rep(c(est_s1_NEQ3_WE0_split$upper, rev(est_s1_NEQ3_WE0_split$lower)), each = 2)[-length(c(est_s1_NEQ3_WE0_split$upper, rev(est_s1_NEQ3_WE0_split$lower)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(est_s1_NEQ3_WE0_split$time, est_s1_NEQ3_WE0_split$surv,
      lwd = 1,
      lty = 3,
      col = my_colors[1],
      type = "s")

polygon(rep(c(est_s2_NEQ3_WE0_split$time, rev(est_s2_NEQ3_WE0_split$time)), each = 2)[-1],
        rep(c(est_s2_NEQ3_WE0_split$upper, rev(est_s2_NEQ3_WE0_split$lower)), each = 2)[-length(c(est_s2_NEQ3_WE0_split$upper, rev(est_s2_NEQ3_WE0_split$lower)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(est_s2_NEQ3_WE0_split$time, est_s2_NEQ3_WE0_split$surv,
      lwd = 1,
      lty = "aa",
      col = my_colors[1],
      type = "s")

polygon(rep(c(est_s2_NEQ3_WE1_split$time, rev(est_s2_NEQ3_WE1_split$time)), each = 2)[-1],
        rep(c(est_s2_NEQ3_WE1_split$upper, rev(est_s2_NEQ3_WE1_split$lower)), each = 2)[-length(c(est_s2_NEQ3_WE1_split$upper, rev(est_s2_NEQ3_WE1_split$lower)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(est_s2_NEQ3_WE1_split$time, est_s2_NEQ3_WE1_split$surv,
      lwd = 1,
      lty = 2,
      col = my_colors[1],
      type = "s")

axis(side = 1, at = times_print, lwd = 0.7)
axis(side = 2, at = seq(0, 1, 0.1), lwd = 0.7)

mtext(side = 1, line = 4, at = times_print[1], text = expression(italic("Number at risk")), adj = 0, cex = 1.4)

mtext(side = 1, line = 5, at = times_print[1] - diff(times_print)[1] * 0.5,
      text = "Z=1", adj = 1, cex = 1.2)
mtext(side = 1, line = 6, at = times_print[1] - diff(times_print)[1] * 0.5,
      text = "Z=2", adj = 1, cex = 1.2)
for (j in 1:2)
{
  for(i in seq_along(times_print)){
    mtext(side = 1, line = 5+j-1, at = times_print[i],
          text = at_risk[j, i], cex = 1.2)
  }
}

par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
lg <- legend("top", inset = c(0, -0.01), xpd = TRUE,
             legend = c("Adjusted estimator", "Naive estimator",
                        "Truth", "", 
                        "X=0/Z=1", "X=0/Z=2",
                        "X=1/Z=1", "X=1/Z=2"),
             col = c(my_colors[1], my_colors[2],
                     "black", "white", 
                     "black", "black",
                     "black", "black"),
             lty = c(1, 1,
                     1, 1,
                     3, 0,   # use 0 (blank) instead of "aa"
                     1, 2), 
             lwd = 1,
             bty = "n",
             ncol = 4,
             cex = 1.2)

# now draw a custom line where the "aa" entry should go
# (the 6th legend item in your vector)
segments(lg$rect$left + 0.48, lg$text$y[6],
         lg$rect$left + 0.58, lg$text$y[6],
         col = "black", lwd = 1, lty = "aa")  # choose your style

dev.off()

par(mfrow = c(1, 1))

# monte carlo ----

HR <- rep(NA, 1000)
CI_low <- rep(NA, 1000)
CI_up <- rep(NA, 1000)

for (j in 1:10000)
{
  
  N <- 1000
  
  shape_1 <- 2
  shape_2 <- 3
  
  median_surv_strate1 <- 50
  median_surv_strate2 <- 20
  scale1 <- log(2)/(median_surv_strate1^shape_1)
  scale2 <- log(2)/(median_surv_strate2^shape_2)
  
  beta_with.effect <- 1
  beta_no.effect <- 0
  
  lambda_cens <- 0.023 # log(2)/30
  
  xaxis <- seq(0, 80, 0.01)
  times_print <- seq(0, max(xaxis), by = 10)
  
  # Data generation ----
  
  set.seed(35138+j)
  U <- runif(N)
  
  strate <- rbinom(N, 1, 0.5) + 1
  X_with.effect <- rbinom(N, 1, 0.5)
  X_no.effect <- rnorm(N, 0, 3)
  
  LP_both <- beta_no.effect * X_no.effect + beta_with.effect * X_with.effect
  
  #####################################'
  
  surv_time_both <- numeric(N)
  for(i in 1:N) {
    surv_time_both[i] <- ((-log(U[i])) / (ifelse(strate[i]==1, scale1, scale2) * exp(LP_both[i])))^(1/ifelse(strate[i]==1, shape_1, shape_2))
  }
  
  censor_time_both <- rexp(N, lambda_cens)
  
  time_both <- pmin(surv_time_both, censor_time_both)
  status_both <- as.numeric(surv_time_both <= censor_time_both)
  data_both <- data.frame(time = time_both,
                          status = status_both,
                          strate = factor(strate),
                          X_no.effect = X_no.effect,
                          X_with.effect = X_with.effect)
  
  cox_triple <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect + strate, data=data_both)
  
  HR[j] <- exp(cox_triple$coefficients[2])
  CI_up[j] <- publish(cox_triple)$rawTable$Upper[2]
  CI_low[j] <- publish(cox_triple)$rawTable$Lower[2]
}

table(CI_up >= exp(1) & CI_low <= exp(1))

summary(HR)
