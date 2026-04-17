
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
# require(prodlim)

path_cox <- "Your path" # Path to the file where figures will be saved, needs to finish with a /

# Constantes ----

if (!exists("global_seed")) global_seed <- 1024

N = 1000

shape = 2

median_surv1 <- 50
scale1 <- log(2)/(median_surv1^shape)

median_surv2 <- 40
scale2 <- log(2)/(median_surv2^shape)

lambda_cens <- 0.023

beta_no.effect1 = 0
beta_with.effect1 = 1

beta_no.effect2 = 0
beta_with.effect2 = 0

xaxis <- seq(0, 100, 0.01)
times_print <- seq(0, max(xaxis), by = 10)

# Data generation ----

set.seed(global_seed)

X_with.effect <- rbinom(N, 1, 0.5) # WE
X_no.effect <- rnorm(N, 0, 10) # NE

LP1 <- beta_no.effect1 * X_no.effect + beta_with.effect1 * X_with.effect
LP2 <- beta_no.effect2 * X_no.effect + beta_with.effect2 * X_with.effect

# Simulation of times for outcome 1 ----

U <- runif(N)
surv_time1 <- ((-log(U)) / (scale1 * exp(LP1)))^(1/shape)

# Simulation of times for outcome 1 ----
# Covariates will have an effect on this outcome

U <- runif(N)
surv_time2 <- ((-log(U)) / (scale2 * exp(LP2)))^(1/shape)

# Observed time ----

T <- pmin(surv_time1, surv_time2)
C <- rexp(N, lambda_cens)

status <- as.numeric(surv_time2 < surv_time1) + 1 
status[C<T] <- 0 # 0 for censoring, 1 for event, 2 for competing risk

time_obs <- pmin(C, T)

# summary(time_obs)
# table(status)

DB <- data.frame(Tobs = time_obs, 
                 status = status,
                 X_with.effect = X_with.effect,
                 X_no.effect = X_no.effect)

# Boostrap ----

N_boot <- 1000

DB$id <- 1:1000

DB_boot <- list()

for (i in 1:N_boot)
{
  DB_boot[[i]] <- DB[sample(1:1000, N, replace = TRUE), ]
}

# True CIF ----

LP1_sE0_WE0 <- beta_no.effect1 * 0 + beta_with.effect1 * 0
LP1_sE0_WE1 <- beta_no.effect1 * 0 + beta_with.effect1 * 1
LP1_sE1_WE0 <- beta_no.effect1 * 1 + beta_with.effect1 * 0
LP1_sE1_WE1 <- beta_no.effect1 * 1 + beta_with.effect1 * 1

LP2_sE0_WE0 <- beta_no.effect2 * 0 + beta_with.effect2 * 0
LP2_sE0_WE1 <- beta_no.effect2 * 0 + beta_with.effect2 * 1
LP2_sE1_WE0 <- beta_no.effect2 * 1 + beta_with.effect2 * 0
LP2_sE1_WE1 <- beta_no.effect2 * 1 + beta_with.effect2 * 1

true_surv1_sE0_WE0 <- (1-exp(-(scale1*exp(LP1_sE0_WE0) + scale2*exp(LP2_sE0_WE0))*xaxis^shape))*scale1*exp(LP1_sE0_WE0)/(scale1*exp(LP1_sE0_WE0) + scale2*exp(LP2_sE0_WE0))
true_surv2_sE0_WE0 <- (1-exp(-(scale1*exp(LP1_sE0_WE0) + scale2*exp(LP2_sE0_WE0))*xaxis^shape))*scale2*exp(LP2_sE0_WE0)/(scale1*exp(LP1_sE0_WE0) + scale2*exp(LP2_sE0_WE0))

true_surv1_sE0_WE1 <- (1-exp(-(scale1*exp(LP1_sE0_WE1) + scale2*exp(LP2_sE0_WE1))*xaxis^shape))*scale1*exp(LP1_sE0_WE1)/(scale1*exp(LP1_sE0_WE1) + scale2*exp(LP2_sE0_WE1))
true_surv2_sE0_WE1 <- (1-exp(-(scale1*exp(LP1_sE0_WE1) + scale2*exp(LP2_sE0_WE1))*xaxis^shape))*scale2*exp(LP2_sE0_WE1)/(scale1*exp(LP1_sE0_WE1) + scale2*exp(LP2_sE0_WE1))

true_surv1_sE1_WE0 <- (1-exp(-(scale1*exp(LP1_sE1_WE0) + scale2*exp(LP2_sE1_WE0))*xaxis^shape))*scale1*exp(LP1_sE1_WE0)/(scale1*exp(LP1_sE1_WE0) + scale2*exp(LP2_sE1_WE0))
true_surv2_sE1_WE0 <- (1-exp(-(scale1*exp(LP1_sE1_WE0) + scale2*exp(LP2_sE1_WE0))*xaxis^shape))*scale2*exp(LP2_sE1_WE0)/(scale1*exp(LP1_sE1_WE0) + scale2*exp(LP2_sE1_WE0))

true_surv1_sE1_WE1 <- (1-exp(-(scale1*exp(LP1_sE1_WE1) + scale2*exp(LP2_sE1_WE1))*xaxis^shape))*scale1*exp(LP1_sE1_WE1)/(scale1*exp(LP1_sE1_WE1) + scale2*exp(LP2_sE1_WE1))
true_surv2_sE1_WE1 <- (1-exp(-(scale1*exp(LP1_sE1_WE1) + scale2*exp(LP2_sE1_WE1))*xaxis^shape))*scale2*exp(LP2_sE1_WE1)/(scale1*exp(LP1_sE1_WE1) + scale2*exp(LP2_sE1_WE1))

# Cox models ----

cox1 <- coxph(Surv(Tobs, status==1) ~ X_no.effect + X_with.effect, data = DB)
round(exp(cox1$coefficients), 2)
publish(cox1)

cox2 <- coxph(Surv(Tobs, status==2) ~ X_no.effect + X_with.effect, data = DB)
round(exp(cox2$coefficients), 2)
publish(cox2)

cox_cr <- coxph(Surv(Tobs, as.factor(status)) ~ X_no.effect + X_with.effect, data = DB, id = 1:nrow(DB))

# TABLE: 3
cox_cr_ft <- as.data.frame(round(summary(cox_cr)$coefficients[, c(2, 6)], 3))
cox_cr_ft[cox_cr_ft[, 2] == 0, 2] <- "<0.001"
cox_cr_ft[, 1] <- round(cox_cr_ft[, 1], 2)
cox_cr_ft[, 3] <- c("X_se, outcome 1",
                    "X_ae, outcome 1",
                    "X_se, outcome 2",
                    "X_ae, outcome 2")

cox_cr_ft[, 4] <- NA

for (i in 1:nrow(summary(cox_cr)$coefficients))
{
  cox_cr_ft[i, 4] <- paste("[", 
                           round(exp(summary(cox_cr)$coefficients[i, 1] - qnorm(0.975)*summary(cox_cr)$coefficients[i, 3]), 2),
                           ";",
                           round(exp(summary(cox_cr)$coefficients[i, 1] + qnorm(0.975)*summary(cox_cr)$coefficients[i, 3]), 2),
                           "]",
                           sep = "")
}



cox_cr_ft <- cox_cr_ft[, c(3, 1, 4, 2)]
names(cox_cr_ft) <- c("Variable", "HazardRatio", "CI.95", "p-value")

ft <- autofit(align(flextable(cox_cr_ft), align = "left"))

# Cox bootstrap ----

fit <- survfit(cox_cr, newdata = c(X_no.effect=0,
                                   X_with.effect=0))

DB_boot_surv1_0 <- as.data.frame(matrix(nrow = N_boot, ncol = length(fit$time)))
DB_boot_surv2_0 <- as.data.frame(matrix(nrow = N_boot, ncol = length(fit$time)))

DB_boot_surv1_1 <- as.data.frame(matrix(nrow = N_boot, ncol = length(fit$time)))
DB_boot_surv2_1 <- as.data.frame(matrix(nrow = N_boot, ncol = length(fit$time)))

for (i in 1:N_boot)
{
  cox_cr_boot <- coxph(Surv(Tobs, as.factor(status)) ~ X_no.effect + X_with.effect, data = DB_boot[[i]], id = 1:nrow(DB))
  fit_boot <- survfit(cox_cr_boot, newdata = c(X_no.effect=as.numeric(quantile(X_no.effect, 0.75)),
                                          X_with.effect=0))
  
  fun1 <- stepfun(fit_boot$time, c(0, as.data.frame(fit_boot$pstate)[, 2]))
  fun2 <- stepfun(fit_boot$time, c(0, as.data.frame(fit_boot$pstate)[, 3]))
  
  DB_boot_surv1_0[i, ] <- fun1(fit$time)
  DB_boot_surv2_0[i, ] <- fun2(fit$time)
  
  fit_boot <- survfit(cox_cr_boot, newdata = c(X_no.effect=as.numeric(quantile(X_no.effect, 0.75)),
                                               X_with.effect=1))
  
  fun1 <- stepfun(fit_boot$time, c(0, as.data.frame(fit_boot$pstate)[, 2]))
  fun2 <- stepfun(fit_boot$time, c(0, as.data.frame(fit_boot$pstate)[, 3]))
  
  DB_boot_surv1_1[i, ] <- fun1(fit$time)
  DB_boot_surv2_1[i, ] <- fun2(fit$time)
}

CI_low1_0 <- apply(DB_boot_surv1_0, 2, quantile, 0.025)
CI_up1_0 <- apply(DB_boot_surv1_0, 2, quantile, 0.975)
CI_low2_0 <- apply(DB_boot_surv2_0, 2, quantile, 0.025)
CI_up2_0 <- apply(DB_boot_surv2_0, 2, quantile, 0.975)

CI_low1_1 <- apply(DB_boot_surv1_1, 2, quantile, 0.025)
CI_up1_1 <- apply(DB_boot_surv1_1, 2, quantile, 0.975)
CI_low2_1 <- apply(DB_boot_surv2_1, 2, quantile, 0.025)
CI_up2_1 <- apply(DB_boot_surv2_1, 2, quantile, 0.975)

# FIGURE3 ----

my_colors <- paletteer_d("rcartocolor::Pastel")

at_risk <- summary(survfit(cox1, type = "aalen"), times_print)$n.risk

jpeg(paste(path_cox, "FIGURE3.jpeg", sep = ""), width = 6, height = 3.25, units = "in", res = 1200, pointsize = 4)
par(mar = c(7, 5, 5, 2), xaxs = "i", yaxs = "i", tcl = 0.2, cex.lab = 1.4, cex.axis = 1.2, mfrow = c(1, 2))

# a. ###############################################################'

plot(1, type = "n",
     col = NULL,
     ylim = c(0, 1),
     xlim = c(0, 80),
     sub = "",
     xlab = "Time", ylab = "CIF", axes = FALSE,
     panel.first = grid(nx = NULL, ny = NULL,
                        col = rgb(0.8, 0.8, 0.8, 0.7),
                        lty = "dotted", lwd = 1))

title(main = "a.", adj = 0, line = 0.5)

lines(xaxis, true_surv1_sE0_WE0, lty = 1, lwd = 1, col = "black")
lines(xaxis, true_surv2_sE0_WE0, lty = 3, lwd = 1, col = "black")

lines(survfit(cox1, newdata = c(X_no.effect = as.numeric(quantile(X_no.effect, 0.75)),
                                X_with.effect = 0))$time,
      1-survfit(cox1, newdata = c(X_no.effect = as.numeric(quantile(X_no.effect, 0.75)),
                                  X_with.effect = 0))$surv,
      lwd = 1,
      lty = 1,
      col = my_colors[2],
      type = "s")
lines(survfit(cox2, newdata = c(X_no.effect = as.numeric(quantile(X_no.effect, 0.75)),
                                X_with.effect = 0))$time,
      1-survfit(cox2, newdata = c(X_no.effect = as.numeric(quantile(X_no.effect, 0.75)),
                                  X_with.effect = 0))$surv,
      lwd = 1,
      lty = 3,
      col = my_colors[2],
      type = "s")

polygon(rep(c(survfit(cox_cr, newdata = c(X_no.effect=as.numeric(quantile(X_no.effect, 0.75)),
                                          X_with.effect=0))$time, 
              rev(survfit(cox_cr, newdata = c(X_no.effect=as.numeric(quantile(X_no.effect, 0.75)),
                                              X_with.effect=0))$time)), each = 2)[-1],
        rep(c(CI_up1_0, rev(CI_low1_0)), each = 2)[-length(c(CI_up1_0, rev(CI_low1_0)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
polygon(rep(c(survfit(cox_cr, newdata = c(X_no.effect=as.numeric(quantile(X_no.effect, 0.75)),
                                          X_with.effect=0))$time, 
              rev(survfit(cox_cr, newdata = c(X_no.effect=as.numeric(quantile(X_no.effect, 0.75)),
                                              X_with.effect=0))$time)), each = 2)[-1],
        rep(c(CI_up2_0, rev(CI_low2_0)), each = 2)[-length(c(CI_up2_0, rev(CI_low2_0)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(survfit(cox_cr, newdata = c(X_no.effect=as.numeric(quantile(X_no.effect, 0.75)),
                                  X_with.effect=0)),
      lwd = 1,
      lty = c(1, 3),
      col = my_colors[1],
      type = "s")

axis(side = 1, at = times_print, lwd = 0.7)
axis(side = 2, at = seq(0, 1, 0.1), lwd = 0.7)

mtext(side = 1, line = 4, at = times_print[1], text = expression(italic("Number at risk")), adj = 0, cex = 1.4)

mtext(side = 1, line = 5, at = times_print[1] - diff(times_print)[1] * 0.5,
      text = rownames(at_risk)[1], adj = 1, cex = 1.2)
for(i in seq_along(times_print)){
  mtext(side = 1, line = 5, at = times_print[i],
        text = at_risk[i], cex = 1.2)
}

# b. ###############################################################'

plot(1, type = "n",
     col = NULL,
     ylim = c(0, 1),
     xlim = c(0, 80),
     sub = "",
     xlab = "Time", ylab = "CIF", axes = FALSE,
     panel.first = grid(nx = NULL, ny = NULL,
                        col = rgb(0.8, 0.8, 0.8, 0.7),
                        lty = "dotted", lwd = 1))

title(main = "b.", adj = 0, line = 0.5)

lines(xaxis, true_surv1_sE0_WE1, lty = 1, lwd = 1, col = "black")
lines(xaxis, true_surv2_sE0_WE1, lty = 3, lwd = 1, col = "black")

lines(survfit(cox1, newdata = c(X_no.effect = as.numeric(quantile(X_no.effect, 0.75)),
                                X_with.effect = 1))$time,
      1-survfit(cox1, newdata = c(X_no.effect = as.numeric(quantile(X_no.effect, 0.75)),
                                  X_with.effect = 1))$surv,
      lwd = 1,
      lty = 1,
      col = my_colors[2],
      type = "s")
lines(survfit(cox2, newdata = c(X_no.effect = as.numeric(quantile(X_no.effect, 0.75)),
                                X_with.effect = 1))$time,
      1-survfit(cox2, newdata = c(X_no.effect = as.numeric(quantile(X_no.effect, 0.75)),
                                  X_with.effect = 1))$surv,
      lwd = 1,
      lty = 3,
      col = my_colors[2],
      type = "s")

polygon(rep(c(survfit(cox_cr, newdata = c(X_no.effect=as.numeric(quantile(X_no.effect, 0.75)),
                                          X_with.effect=1))$time, 
              rev(survfit(cox_cr, newdata = c(X_no.effect=as.numeric(quantile(X_no.effect, 0.75)),
                                              X_with.effect=1))$time)), each = 2)[-1],
        rep(c(CI_up1_1, rev(CI_low1_1)), each = 2)[-length(c(CI_up1_1, rev(CI_low1_1)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
polygon(rep(c(survfit(cox_cr, newdata = c(X_no.effect=as.numeric(quantile(X_no.effect, 0.75)),
                                          X_with.effect=1))$time, 
              rev(survfit(cox_cr, newdata = c(X_no.effect=as.numeric(quantile(X_no.effect, 0.75)),
                                              X_with.effect=1))$time)), each = 2)[-1],
        rep(c(CI_up2_1, rev(CI_low2_1)), each = 2)[-length(c(CI_up2_1, rev(CI_low2_1)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)

lines(survfit(cox_cr, newdata = c(X_no.effect=as.numeric(quantile(X_no.effect, 0.75)),
                                  X_with.effect=1)),
      lwd = 1,
      lty = c(1, 3),
      col = my_colors[1],
      type = "s")

axis(side = 1, at = times_print, lwd = 0.7)
axis(side = 2, at = seq(0, 1, 0.1), lwd = 0.7)

mtext(side = 1, line = 4, at = times_print[1], text = expression(italic("Number at risk")), adj = 0, cex = 1.4)

mtext(side = 1, line = 5, at = times_print[1] - diff(times_print)[1] * 0.5,
      text = rownames(at_risk)[1], adj = 1, cex = 1.2)
for(i in seq_along(times_print)){
  mtext(side = 1, line = 5, at = times_print[i],
        text = at_risk[i], cex = 1.2)
}



par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("top", inset = c(0, -0.01), xpd = TRUE,
       legend = c("Adjusted estimator", "Naive estimator",
                  "Truth", "", 
                  "Outcome 1", "Outcome 2"),
       col = c(my_colors[1], my_colors[2],
               "black", "white", 
               "black", "black"),
       lty = c(1,1,
               1,1,
               1,3), 
       lwd = 1,
       bty = "n",
       ncol = 4,
       cex = 1.2)

dev.off()

par(mfrow = c(1, 1))