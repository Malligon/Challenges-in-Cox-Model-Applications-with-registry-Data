
rm(list = ls())

# library
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
library(car)

path_cox <- "Your path" # Path to the file where figures will be saved, needs to finish with a /

# Parameters ----

# Define a global seed
if (!exists("global_seed")) global_seed <- 1024

# parameters of simulation
N <- 1000
shape <- 2
median_surv <- 50
scale <- log(2)/(median_surv^shape)

beta_with.effect <- 1
beta_no.effect <- 0

lambda_cens <- 0.023 # log(2)/30

xaxis <- seq(0, 100, 0.01)
times_print <- seq(0, max(xaxis), by = 10)

# Data generation ----

set.seed(global_seed)
U <- runif(N)

X_with.effect <- rbinom(N, 1, 0.5) # WE
X_no.effect <- rnorm(N, 0, 3) # NE

# Linear predictor
LP_both <- beta_no.effect * X_no.effect + beta_with.effect * X_with.effect

# Times 
surv_time_both <- ((-log(U)) / (scale * exp(LP_both)))^(1/shape)

censor_time_both <- rexp(N, lambda_cens)

time_both <- pmin(surv_time_both, censor_time_both)
status_both <- as.numeric(surv_time_both <= censor_time_both)
data_both <- data.frame(time = time_both,
                        status = status_both,
                        X_no.effect = X_no.effect,
                        X_with.effect = X_with.effect)

cox_both <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect, data = data_both)
publish(cox_both)

# Cox models ----

# no covariate
cox_noCov <- coxph(Surv(time, status) ~ 1, data = data_both)

# survival calculation
est_noCov <- survfit(cox_noCov)

df1 <- data.frame(time = est_noCov$time, S = est_noCov$surv, model = "Model without covariate")


# Cox with just the no effect covariate
coX_no.effect <- coxph(Surv(time, status) ~ X_no.effect, data = data_both)

# NE = Q3
est_no.effect_Q3 <- survfit(coX_no.effect, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75))))
df2 <- data.frame(time = est_no.effect_Q3$time, S = est_no.effect_Q3$surv, model = "Model without effect, NE=Q3")

# NE = Q1
est_no.effect_Q1 <- survfit(coX_no.effect, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.25))))
df2bis <- data.frame(time = est_no.effect_Q1$time, S = est_no.effect_Q1$surv, model = "Model without effect, NE=Q1")


# Cox with just the covariate with effect
coX_with.effect <- coxph(Surv(time, status) ~ X_with.effect, data = data_both)

# WE = 1
est_with.effect_1 <- survfit(coX_with.effect, newdata = data.frame(X_with.effect=1))
df3 <- data.frame(time = est_with.effect_1$time, S = est_with.effect_1$surv, model = "Model with effect, WE=1")

# WE = 0
est_with.effect_0 <- survfit(coX_with.effect, newdata = data.frame(X_with.effect=0))
df3bis <- data.frame(time = est_with.effect_0$time, S = est_with.effect_0$surv, model = "Model with effect, WE=0")


# Cox with both covariates
cox_both <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect, data = data_both)

# NE = Q3, WE = 1
est_both_NEQ3_WE1 <- survfit(cox_both, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=1))
df4 <- data.frame(time = est_both_NEQ3_WE1$time, S = est_both_NEQ3_WE1$surv, model = "Model with both covariates, NE=Q3 WE=1")

# NE = Q3, WE = 0
est_both_NEQ3_WE0 <- survfit(cox_both, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=0))
df4bis <- data.frame(time = est_both_NEQ3_WE0$time, S = est_both_NEQ3_WE0$surv, model = "Model with both covariates, NE=Q3 WE=0")

# NE = Q1, WE = 1
est_both_NEQ1_WE1 <- survfit(cox_both, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.25)), X_with.effect=1))
df4ter <- data.frame(time = est_both_NEQ1_WE1$time, S = est_both_NEQ1_WE1$surv, model = "Model with both covariates, NE=Q1 WE=1")

# NE = Q1, WE = 0
est_both_NEQ1_WE0 <- survfit(cox_both, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.25)), X_with.effect=0))
df4quat <- data.frame(time = est_both_NEQ1_WE0$time, S = est_both_NEQ1_WE0$surv, model = "Model with both covariates, NE=Q1 WE=0")

# HR
HR_no.effect <- exp(coef(coX_no.effect))
HR_with.effect <- exp(coef(coX_with.effect))
HR_both <- exp(coef(cox_both))

# CI
s <- summary(coX_no.effect)$coef
beta_j <- coX_no.effect$coefficients["X_no.effect"]
se_j <- s["X_no.effect","se(coef)"]

hr_j <- exp(beta_j)
ic_j <- exp(beta_j + c(-1,1)*qnorm(0.975)*se_j)
round(c(HR=hr_j, IC_inf=ic_j[1], IC_sup=ic_j[2]), 2)

# publish
invisible(capture.output(pub_no.effect <- publish(coX_no.effect)))
invisible(capture.output(pub_with.effect <- publish(coX_with.effect)))
invisible(capture.output(pub_both <- publish(cox_both)))

# comparison of results ----

comparison <- data.frame(
  Modele = c("No effect", "With effect", "Both - No effect", "Both - With effect"),
  HR_vrai = c(exp(beta_no.effect), exp(beta_with.effect), exp(beta_no.effect), exp(beta_with.effect)),
  HR_estime = c(HR_no.effect, HR_with.effect, HR_both[1], HR_both[2]),
  IC_inf = c(round(pub_no.effect$rawTable$Lower,4),
             round(pub_with.effect$rawTable$Lower,4),
             round(pub_both$rawTable$Lower[1],4),
             round(pub_both$rawTable$Lower[2],4)),
  IC_sup = c(round(pub_no.effect$rawTable$Upper,4),
             round(pub_with.effect$rawTable$Upper,4),
             round(pub_both$rawTable$Upper[1],4),
             round(pub_both$rawTable$Upper[2],4))
)

# TABLE: 1
comparison <- comparison %>%
  kable(digits = 7, align = "lcc") %>%
  kable_styling(full_width = FALSE) %>%
  pack_rows("Models with one covariate", 1, 2) %>%
  pack_rows("Model with both covariates", 3, 4)


# Prep figures ----

S_true_fun_with.effect <- function(t) {
  exp(-scale * exp(beta_with.effect) * t^shape)
  # exp(-exp(beta_with.effect) * (t/scale)^shape)
}
S_true_with.effect <- S_true_fun_with.effect(xaxis)

S_true_fun_no.effect <- function(t) {
  exp(-scale * 1 * t^shape)
  # exp(-(t/scale)^shape)
}
S_true_no.effect <- S_true_fun_no.effect(xaxis)

df_true_with.effect <- data.frame(time = xaxis, S = S_true_with.effect, model = "Truth with WE = 1")
df_true_no.effect <- data.frame(time = xaxis, S = S_true_no.effect, model = "Truth with WE = 0")

# FIGURE 1 ----

my_colors <- paletteer_d("rcartocolor::Pastel")

jpeg(paste(path_cox, "FIGURE1.jpeg", sep = ""), width = 6, height = 3.25, units = "in", res = 1200, pointsize = 4)
par(mar = c(7, 5, 5, 2), xaxs = "i", yaxs = "i", tcl = 0.2, cex.lab = 1.4, cex.axis = 1.2
    ,mfrow = c(2, 2)
)

# a. ###############################################################'

fit <- survfit(cox_noCov, type = "aalen")
at_risk <- summary(fit, times_print)$n.risk

plot(1,
     col = NULL,
     ylim = c(0, 1),
     xlim = c(0, 100),
     sub = "",
     xlab = "Time", ylab = "Survival", axes = FALSE,
     panel.first = grid(nx = NULL, ny = NULL,
                        col = rgb(0.8, 0.8, 0.8, 0.7),
                        lty = "dotted", lwd = 1))
title(main = "a.", adj = 0, line = 0.5)

lines(xaxis, df_true_no.effect$S, lty = 1, lwd = 1, col = "black")
lines(xaxis, df_true_with.effect$S, lty = 3, lwd = 1, col = "black")

polygon(rep(c(fit$time, rev(fit$time)), each = 2)[-1],
        rep(c(fit$upper, rev(fit$lower)), each = 2)[-length(c(fit$upper, rev(fit$lower)))*2],
        col = adjustcolor(my_colors[7], alpha.f = 0.2),
        border = NA)
lines(df1$time, df1$S,
      lwd = 1,
      lty = 1,
      col = my_colors[7],
      type = "s")

axis(side = 1, at = times_print, lwd = 0.7)
axis(side = 2, at = seq(0, 1, 0.1), lwd = 0.7)

mtext(side = 1, line = 4, at = times_print[1], text = expression(italic("Number at risk")), adj = 0, cex = 1)

mtext(side = 1, line = 5, at = times_print[1] - diff(times_print)[1] * 0.5,
      text = rownames(at_risk)[1], adj = 1, cex = 1)
for(i in seq_along(times_print)){
  mtext(side = 1, line = 5, at = times_print[i],
        text = at_risk[i], cex = 1)
}

legend("top", inset = c(0, -0.17), xpd = TRUE,
       legend = c("Without covariates"),
       col = c(my_colors[7]),
       lty = c(1), 
       lwd = 1,
       bty = "n",
       ncol = 1,
       cex = 1.2)

# b. ###############################################################'

fit <- survfit(coX_no.effect, type = "aalen")
at_risk <- summary(fit, times_print)$n.risk


plot(1,
     col = NULL,
     ylim = c(0, 1),
     xlim = c(0, 100),
     sub = "",
     xlab = "Time", ylab = "Survival", axes = FALSE,
     panel.first = grid(nx = NULL, ny = NULL,
                        col = rgb(0.8, 0.8, 0.8, 0.7),
                        lty = "dotted", lwd = 1))
title(main = "b.", adj = 0, line = 0.5)

lines(xaxis, df_true_no.effect$S, lty = 1, lwd = 1, col = "black")
lines(xaxis, df_true_with.effect$S, lty = 3, lwd = 1, col = "black")

polygon(rep(c(est_no.effect_Q3$time, rev(est_no.effect_Q3$time)), each = 2)[-1],
        rep(c(est_no.effect_Q3$upper, rev(est_no.effect_Q3$lower)), each = 2)[-length(c(est_no.effect_Q3$upper, rev(est_no.effect_Q3$lower)))*2],
        col = adjustcolor(my_colors[8], alpha.f = 0.2),
        border = NA)
lines(df2$time, df2$S,
      lwd = 1,
      lty = 1,
      col = my_colors[8],
      type = "s")

polygon(rep(c(est_no.effect_Q1$time, rev(est_no.effect_Q1$time)), each = 2)[-1],
        rep(c(est_no.effect_Q1$upper, rev(est_no.effect_Q1$lower)), each = 2)[-length(c(est_no.effect_Q1$upper, rev(est_no.effect_Q1$lower)))*2],
        col = adjustcolor(my_colors[9], alpha.f = 0.2),
        border = NA, lty=3)
lines(df2bis$time, df2bis$S,
      lwd = 1,
      lty = 3,
      col = my_colors[9],
      type = "s")

axis(side = 1, at = times_print, lwd = 0.7)
axis(side = 2, at = seq(0, 1, 0.1), lwd = 0.7)

mtext(side = 1, line = 4, at = times_print[1], text = expression(italic("Number at risk")), adj = 0, cex = 1)

mtext(side = 1, line = 5, at = times_print[1] - diff(times_print)[1] * 0.5,
      text = rownames(at_risk)[1], adj = 1, cex = 1)
for(i in seq_along(times_print)){
  mtext(side = 1, line = 5, at = times_print[i],
        text = at_risk[i], cex = 1)
}

legend("top", inset = c(0, -0.17), xpd = TRUE,
       legend = c(expression(paste(tilde(X), "=Q1")), 
                  expression(paste(tilde(X), "=Q3"))),
       col = c(my_colors[8], my_colors[9]),
       lty = c(1, 3), 
       lwd = 1,
       bty = "n",
       ncol = 2,
       cex = 1.2)

# c. ###############################################################'

fit <- survfit(coX_with.effect, type = "aalen")
at_risk <- summary(fit, times_print)$n.risk

plot(1,
     col = NULL,
     ylim = c(0, 1),
     xlim = c(0, 100),
     sub = "",
     xlab = "Time", ylab = "Survival", axes = FALSE,
     panel.first = grid(nx = NULL, ny = NULL,
                        col = rgb(0.8, 0.8, 0.8, 0.7),
                        lty = "dotted", lwd = 1))
title(main = "c.", adj = 0, line = 0.5)

lines(xaxis, df_true_no.effect$S, lty = 1, lwd = 1, col = "black")
lines(xaxis, df_true_with.effect$S, lty = 3, lwd = 1, col = "black")

polygon(rep(c(est_with.effect_1$time, rev(est_with.effect_1$time)), each = 2)[-1],
        rep(c(est_with.effect_1$upper, rev(est_with.effect_1$lower)), each = 2)[-length(c(est_with.effect_1$upper, rev(est_with.effect_1$lower)))*2],
        col = adjustcolor(my_colors[2], alpha.f = 0.2),
        border = NA)
lines(df3$time, df3$S,
      lwd = 1,
      lty = 1,
      col = my_colors[2],
      type = "s")

polygon(rep(c(est_with.effect_0$time, rev(est_with.effect_0$time)), each = 2)[-1],
        rep(c(est_with.effect_0$upper, rev(est_with.effect_0$lower)), each = 2)[-length(c(est_with.effect_0$upper, rev(est_with.effect_0$lower)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(df3bis$time, df3bis$S,
      lwd = 1,
      lty = 1,
      col = my_colors[1],
      type = "s")

axis(side = 1, at = times_print, lwd = 0.7)
axis(side = 2, at = seq(0, 1, 0.1), lwd = 0.7)

mtext(side = 1, line = 4, at = times_print[1], text = expression(italic("Number at risk")), adj = 0, cex = 1)

mtext(side = 1, line = 5, at = times_print[1] - diff(times_print)[1] * 0.5,
      text = rownames(at_risk)[1], adj = 1, cex = 1)
for(i in seq_along(times_print)){
  mtext(side = 1, line = 5, at = times_print[i],
        text = at_risk[i], cex = 1)
}

legend("top", inset = c(0, -0.17), xpd = TRUE,
       legend = c("X=0", "X=1"),
       col = c(my_colors[1], my_colors[2]),
       lty = c(1, 1), 
       lwd = 1,
       bty = "n",
       ncol = 2,
       cex = 1.2)

# d. ###############################################################'

fit <- survfit(cox_both, type = "aalen")
at_risk <- summary(fit, times_print)$n.risk

plot(1,
     col = NULL,
     ylim = c(0, 1),
     xlim = c(0, 100),
     sub = "",
     xlab = "Time", ylab = "Survival", axes = FALSE,
     panel.first = grid(nx = NULL, ny = NULL,
                        col = rgb(0.8, 0.8, 0.8, 0.7),
                        lty = "dotted", lwd = 1))
title(main = "d.", adj = 0, line = 0.5)

lines(xaxis, df_true_no.effect$S, lty = 1, lwd = 1, col = "black")
lines(xaxis, df_true_with.effect$S, lty = 3, lwd = 1, col = "black")

polygon(rep(c(est_both_NEQ3_WE1$time, rev(est_both_NEQ3_WE1$time)), each = 2)[-1],
        rep(c(est_both_NEQ3_WE1$upper, rev(est_both_NEQ3_WE1$lower)), each = 2)[-length(c(est_both_NEQ3_WE1$upper, rev(est_both_NEQ3_WE1$lower)))*2],
        col = adjustcolor(my_colors[6], alpha.f = 0.2),
        border = NA)
lines(df4$time, df4$S,
      lwd = 1,
      lty = 1,
      col = my_colors[6],
      type = "s")

polygon(rep(c(est_both_NEQ3_WE0$time, rev(est_both_NEQ3_WE0$time)), each = 2)[-1],
        rep(c(est_both_NEQ3_WE0$upper, rev(est_both_NEQ3_WE0$lower)), each = 2)[-length(c(est_both_NEQ3_WE0$upper, rev(est_both_NEQ3_WE0$lower)))*2],
        col = adjustcolor(my_colors[5], alpha.f = 0.2),
        border = NA)
lines(df4bis$time, df4bis$S,
      lwd = 1,
      lty = 1,
      col = my_colors[5],
      type = "s")


polygon(rep(c(est_both_NEQ1_WE1$time, rev(est_both_NEQ1_WE1$time)), each = 2)[-1],
        rep(c(est_both_NEQ1_WE1$upper, rev(est_both_NEQ1_WE1$lower)), each = 2)[-length(c(est_both_NEQ1_WE1$upper, rev(est_both_NEQ1_WE1$lower)))*2],
        col = adjustcolor(my_colors[3], alpha.f = 0.2),
        border = NA,
        lty    = 3)
lines(df4ter$time, df4ter$S,
      lwd = 1,
      lty = 3,
      col = my_colors[3],
      type = "s")

polygon(rep(c(est_both_NEQ1_WE0$time, rev(est_both_NEQ1_WE0$time)), each = 2)[-1],
        rep(c(est_both_NEQ1_WE0$upper, rev(est_both_NEQ1_WE0$lower)), each = 2)[-length(c(est_both_NEQ1_WE0$upper, rev(est_both_NEQ1_WE0$lower)))*2],
        col = adjustcolor(my_colors[4], alpha.f = 0.2),
        border = NA,
        lty    = 3)
lines(df4quat$time, df4quat$S,
      lwd = 1,
      lty = 3,
      col = my_colors[4],
      type = "s")

axis(side = 1, at = times_print, lwd = 0.7)
axis(side = 2, at = seq(0, 1, 0.1), lwd = 0.7)

mtext(side = 1, line = 4, at = times_print[1], text = expression(italic("Number at risk")), adj = 0, cex = 1)

mtext(side = 1, line = 5, at = times_print[1] - diff(times_print)[1] * 0.5,
      text = rownames(at_risk)[1], adj = 1, cex = 1)
for(i in seq_along(times_print)){
  mtext(side = 1, line = 5, at = times_print[i],
        text = at_risk[i], cex = 1)
}

legend("top", inset = c(0, -0.17), xpd = TRUE,
       legend = c(expression(paste("X=0/", tilde(X), "=Q1")),
                  expression(paste("X=1/", tilde(X), "=Q1")),
                  expression(paste("X=0/", tilde(X), "=Q3")),
                  expression(paste("X=1/", tilde(X), "=Q3"))),
       col = c(my_colors[4], my_colors[3], my_colors[5], my_colors[6]),
       lty = c(3, 3, 1, 1), 
       lwd = 1,
       bty = "n",
       ncol = 2,
       cex = 1.2)

par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("top", inset = c(0, -0.01), xpd = TRUE,
       legend = c("Truth for X=0", "Truth for X=1"),
       col = c("black", "black"),
       lty = c(3, 1), 
       lwd = 1,
       bty = "n",
       ncol = 2,
       cex = 1.2)

dev.off()

par(mfrow = c(1, 1))
