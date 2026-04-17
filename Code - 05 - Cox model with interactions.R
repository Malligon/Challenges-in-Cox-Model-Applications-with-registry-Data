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
shape <- 2
median_surv <- 50
scale <- log(2)/(median_surv^shape)

beta_avec.effet <- 1
beta_sans.effet <- 0

lambda_cens <- 0.023 # log(2)/30

beta_int = log(4)

xaxis <- seq(0, 100, 0.01)
times_print <- seq(0, max(xaxis), by = 10)

# Data generation ----

set.seed(global_seed)
U <- runif(N)

X_avec.effet <- rbinom(N, 1, 0.5)
X_sans.effet <- rnorm(N, 0, 3)

LP_int <- beta_sans.effet * X_sans.effet + beta_avec.effet * X_avec.effet + beta_int*(X_avec.effet*X_sans.effet)

# Simulation ----

surv_time_int <- ((-log(U)) / (scale * exp(LP_int)))^(1/shape)
censor_time_both <- rexp(N, lambda_cens)

time_both <- pmin(surv_time_int, censor_time_both)
status_both <- as.numeric(surv_time_int <= censor_time_both)
data_both <- data.frame(time = time_both,
                        status = status_both,
                        X_sans.effet = X_sans.effet,
                        X_avec.effet = X_avec.effet)

# Cox 1 ----

cox_both <- coxph(Surv(time, status) ~ X_sans.effet + X_avec.effet, data = data_both)
est_NEQ3_WE1 <- survfit(cox_both, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.75)), X_avec.effet=1))
df4 <- data.frame(time = est_NEQ3_WE1$time, S = est_NEQ3_WE1$surv, model = "NE=Q3 WE=1")


est_NEQ3_WE0 <- survfit(cox_both, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.75)), X_avec.effet=0))
df4bis <- data.frame(time = est_NEQ3_WE0$time, S = est_NEQ3_WE0$surv, model = "NE=Q3 WE=0")


est_NEQ1_WE1 <- survfit(cox_both, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.25)), X_avec.effet=1))
df4ter <- data.frame(time = est_NEQ1_WE1$time, S = est_NEQ1_WE1$surv, model = "NE=Q1 WE=1")


est_NEQ1_WE0 <- survfit(cox_both, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.25)), X_avec.effet=0))
df4quat <- data.frame(time = est_NEQ1_WE0$time, S = est_NEQ1_WE0$surv, model = "NE=Q1 WE=0")

# Cox 2 ----

cox_both_int  <- coxph(Surv(time, status) ~ X_sans.effet * as.factor(X_avec.effet), data = data_both)
est_NEQ3_WE1_int <- survfit(cox_both_int, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.75)), X_avec.effet=1))
df5 <- data.frame(time = est_NEQ3_WE1_int$time, S = est_NEQ3_WE1_int$surv, model = "Modèle int, NE=Q3 WE=1")


est_NEQ3_WE0_int <- survfit(cox_both_int, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.75)), X_avec.effet=0))
df5bis <- data.frame(time = est_NEQ3_WE0_int$time, S = est_NEQ3_WE0_int$surv, model = "Modèle int, NE=Q3 WE=0")


est_NEQ1_WE1_int <- survfit(cox_both_int, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.25)), X_avec.effet=1))
df5ter <- data.frame(time = est_NEQ1_WE1_int$time, S = est_NEQ1_WE1_int$surv, model = "Modèle int, NE=Q1 WE=1")


est_NEQ1_WE0_int <- survfit(cox_both_int, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.25)), X_avec.effet=0))
df5quat <- data.frame(time = est_NEQ1_WE0_int$time, S = est_NEQ1_WE0_int$surv, model = "Modèle int, NE=Q1 WE=0")

# Comparaison ----

HR_both <- exp(coef(cox_both))
HR_both_int <- exp(coef(cox_both_int))
names(HR_both) <- c("Xne", 
                    "Xwe")
names(HR_both_int) <- c("Xne, Xwe=0", 
                        "Xwe, Xne=Q1", 
                        "Xwe=1 & Xne=Q3")

invisible(capture.output(tab_both <- publish(cox_both)))
invisible(capture.output(tab_both_int <- publish(cox_both_int)))

comparaison <- data.frame(
  Modele    = c("Both - NE", "Both - WE", "Both - NE (int)", "Both - WE (int)", "Both - NE/WE (int)"),
  HR_vrai   = c(exp(beta_sans.effet), exp(beta_avec.effet), exp(beta_sans.effet), exp(beta_avec.effet), exp(beta_int)),
  HR_estime = c(HR_both[1], HR_both[2], HR_both_int[1], HR_both_int[2], HR_both_int[3]),
  IC_inf    = c(tab_both$rawTable$Lower[1],
                tab_both$rawTable$Lower[2],
                tab_both_int$rawTable$Lower[3],
                tab_both_int$rawTable$Lower[1],
                summary(cox_both_int)$conf.int[3,3]),
  IC_sup    = c(tab_both$rawTable$Upper[1],
                tab_both$rawTable$Upper[2],
                tab_both_int$rawTable$Upper[3],
                tab_both_int$rawTable$Upper[1],
                summary(cox_both_int)$conf.int[3,4])
)

comparaison <- comparaison %>%
  kable(digits = 7, align = "lcc") %>%
  kable_styling(full_width = FALSE) %>%
  pack_rows("Basic model", 1, 2) %>%
  pack_rows("Model with interaction", 3, 4)

# Truth ----

S_true_fun <- function(t, beta = 0) { exp(-scale * exp(beta) * t^shape) }

df_true_00 <- data.frame(time = xaxis, S = S_true_fun(xaxis, quantile(X_sans.effet, 0.25) * beta_sans.effet), model = "NE=Q1/WE=0")
df_true_10 <- data.frame(time = xaxis, S = S_true_fun(xaxis, quantile(X_sans.effet, 0.75) * beta_sans.effet), model = "NE=Q3/WE=0")
df_true_01 <- data.frame(time = xaxis, S = S_true_fun(xaxis, beta_sans.effet + beta_avec.effet + quantile(X_sans.effet, 0.25) * beta_int), model = "NE=Q1/WE=1")
df_true_11_int <- data.frame(time = xaxis, S = S_true_fun(xaxis, beta_sans.effet + beta_avec.effet + quantile(X_sans.effet, 0.75) * beta_int), model = "NE=Q3/WE=1")

# Prep publish ----

publish_int <- publish(cox_both_int)$regressionTable

publish_int$Variable <- c("X_ne, X_we=0",
                          "X_ne, X_we=1")

# Expected values : 1 and exp(beta_sans.effet * 1 + beta_avec.effet * 1 + beta_int*(1*1))/exp(beta_sans.effet * 0 + beta_avec.effet * 1 + beta_int*(1*0))

# TABLE: 4
ft <- autofit(flextable(publish(cox_both)$regressionTable))
# TABLE: 5
ft <- autofit(flextable(publish_int))

# calculation for the tables
beta <- coef(cox_both_int)
covm <- vcov(cox_both_int)

# Linear combination for Xne
lc <- c(1, 0, 1)

HR <- exp(sum(beta * lc))
NE <- sqrt(t(lc) %*% covm %*% lc)
CI <- exp(sum(beta * lc) + c(-1,1)*1.96*NE)

# Linear combination for Xwe
lc <- c(0, 1, 1)

HR <- exp(sum(beta * lc))
NE <- sqrt(t(lc) %*% covm %*% lc)
CI <- exp(sum(beta * lc) + c(-1,1)*1.96*NE)



# FIGURE5 ----

my_colors <- paletteer_d("rcartocolor::Pastel")

fit <- survfit(cox_both, type = "aalen")
at_risk <- summary(fit, times_print)$n.risk

jpeg(paste(path_cox, "FIGURE5.jpeg", sep = ""), width = 6, height = 3.25, units = "in", res = 1200, pointsize = 4)
par(mar = c(7, 7, 5, 2), xaxs = "i", yaxs = "i", tcl = 0.2, cex.lab = 1.4, cex.axis = 1.2)

# a. ###############################################################'

plot(1,
     col = NULL,
     ylim = c(0, 1),
     xlim = c(0, 100),
     sub = "",
     xlab = "Time", ylab = "Survival", axes = FALSE,
     panel.first = grid(nx = NULL, ny = NULL,
                        col = rgb(0.8, 0.8, 0.8, 0.7),
                        lty = "dotted", lwd = 1))

lines(xaxis, df_true_10$S, lty = 1, lwd = 1.5, col = "black")
lines(xaxis, df_true_01$S, lty = 3, lwd = 1.5, col = "black")
lines(xaxis, df_true_11_int$S, lty = 5, lwd = 1.5, col = "black")

polygon(rep(c(est_NEQ3_WE1_int$time[which(!is.na(est_NEQ3_WE1_int$upper))], rev(est_NEQ3_WE1_int$time[which(!is.na(est_NEQ3_WE1_int$upper))])), each = 2)[-1],
        rep(c(est_NEQ3_WE1_int$upper[which(!is.na(est_NEQ3_WE1_int$upper))], rev(est_NEQ3_WE1_int$lower[which(!is.na(est_NEQ3_WE1_int$upper))])), each = 2)[-length(c(est_NEQ3_WE1_int$upper[which(!is.na(est_NEQ3_WE1_int$upper))], rev(est_NEQ3_WE1_int$lower[which(!is.na(est_NEQ3_WE1_int$upper))])))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(df5$time, df5$S,
      lwd = 1,
      lty = 5,
      col = my_colors[1],
      type = "s")

polygon(rep(c(est_NEQ1_WE1_int$time, rev(est_NEQ1_WE1_int$time)), each = 2)[-1],
        rep(c(est_NEQ1_WE1_int$upper, rev(est_NEQ1_WE1_int$lower)), each = 2)[-length(c(est_NEQ1_WE1_int$upper, rev(est_NEQ1_WE1_int$lower)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(df5ter$time, df5ter$S,
      lwd = 1,
      lty = 3,
      col = my_colors[1],
      type = "s")

polygon(rep(c(est_NEQ3_WE0_int$time, rev(est_NEQ3_WE0_int$time)), each = 2)[-1],
        rep(c(est_NEQ3_WE0_int$upper, rev(est_NEQ3_WE0_int$lower)), each = 2)[-length(c(est_NEQ3_WE0_int$upper, rev(est_NEQ3_WE0_int$lower)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(df5bis$time, df5bis$S,
      lwd = 1,
      lty = 1,
      col = my_colors[1],
      type = "s")

lines(df4$time, df4$S,
      lwd = 1,
      lty = 5,
      col = my_colors[2],
      type = "s")

lines(df4bis$time, df4bis$S,
      lwd = 1,
      lty = 1,
      col = my_colors[2],
      type = "s")

lines(df4ter$time, df4ter$S,
      lwd = 1.5,
      lty = 3,
      col = my_colors[2],
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

legend("top", inset = c(0, -0.12), xpd = TRUE,
       legend = c("Adjusted estimator", "Naive estimator",
                  "Truth", "", 
                  expression(paste("X=0/", tilde(X), "=Q3")), expression(paste("X=1/", tilde(X), "=Q1")),
                  expression(paste("X=1/", tilde(X), "=Q3"))),
       col = c(my_colors[1], my_colors[2],
               "black", "white", 
               "black", "black",
               "black"),
       lty = c(1,1,
               1,1,
               1,3,
               5), 
       lwd = 1,
       bty = "n",
       ncol = 4,
       cex = 1.2)

dev.off()

# Cox 1 ----
cox_both_low <- coxph(Surv(time, status) ~ X_sans.effet + X_avec.effet, data = data_both)
est_NEQ3_WE1 <- survfit(cox_both_low, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.75)), X_avec.effet=1))
df4 <- data.frame(time = est_NEQ3_WE1$time, S = est_NEQ3_WE1$surv, model = "NE=Q3 WE=1")
est_NEQ3_WE0 <- survfit(cox_both_low, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.75)), X_avec.effet=0))
df4bis <- data.frame(time = est_NEQ3_WE0$time, S = est_NEQ3_WE0$surv, model = "NE=Q3 WE=0")
est_NEQ1_WE1 <- survfit(cox_both_low, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.25)), X_avec.effet=1))
df4ter <- data.frame(time = est_NEQ1_WE1$time, S = est_NEQ1_WE1$surv, model = "NE=Q1 WE=1")
est_NEQ1_WE0 <- survfit(cox_both_low, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.25)), X_avec.effet=0))
df4quat <- data.frame(time = est_NEQ1_WE0$time, S = est_NEQ1_WE0$surv, model = "NE=Q1 WE=0")
# Cox 2 ----
cox_both_int_low  <- coxph(Surv(time, status) ~ X_sans.effet * as.factor(X_avec.effet), data = data_both)
est_NEQ3_WE1_int <- survfit(cox_both_int_low, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.75)), X_avec.effet=1))
df5 <- data.frame(time = est_NEQ3_WE1_int$time, S = est_NEQ3_WE1_int$surv, model = "Modèle int, NE=Q3 WE=1")
est_NEQ3_WE0_int <- survfit(cox_both_int_low, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.75)), X_avec.effet=0))
df5bis <- data.frame(time = est_NEQ3_WE0_int$time, S = est_NEQ3_WE0_int$surv, model = "Modèle int, NE=Q3 WE=0")
est_NEQ1_WE1_int <- survfit(cox_both_int_low, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.25)), X_avec.effet=1))
df5ter <- data.frame(time = est_NEQ1_WE1_int$time, S = est_NEQ1_WE1_int$surv, model = "Modèle int, NE=Q1 WE=1")
est_NEQ1_WE0_int <- survfit(cox_both_int_low, newdata = data.frame(X_sans.effet=as.numeric(quantile(data_both$X_sans.effet, 0.25)), X_avec.effet=0))
df5quat <- data.frame(time = est_NEQ1_WE0_int$time, S = est_NEQ1_WE0_int$surv, model = "Modèle int, NE=Q1 WE=0")
# Comparaison ----
HR_both <- exp(coef(cox_both_low))
HR_both_int <- exp(coef(cox_both_int))
names(HR_both) <- c("Xne", 
                    "Xwe")
names(HR_both_int) <- c("Xne, Xwe=0", 
                        "Xwe, Xne=Q1", 
                        "Xwe=1 & Xne=Q3")

invisible(capture.output(tab_both <- publish(cox_both_low)))
invisible(capture.output(tab_both_int <- publish(cox_both_int_low)))

comparaison2 <- data.frame(
  Modele    = c("Both - NE", "Both - WE", "Both - NE (int)", "Both - WE (int)", "Both - NE/WE (int)"),
  HR_vrai   = c(exp(beta_sans.effet), exp(beta_avec.effet), exp(beta_sans.effet), exp(beta_avec.effet), exp(beta_int)),
  HR_estime = c(HR_both[1], HR_both[2], HR_both_int[1], HR_both_int[2], HR_both_int[3]),
  IC_inf    = c(tab_both$rawTable$Lower[1],
                tab_both$rawTable$Lower[2],
                tab_both_int$rawTable$Lower[3],
                tab_both_int$rawTable$Lower[1],
                summary(cox_both_int_low)$conf.int[3,3]),
  IC_sup    = c(tab_both$rawTable$Upper[1],
                tab_both$rawTable$Upper[2],
                tab_both_int$rawTable$Upper[3],
                tab_both_int$rawTable$Upper[1],
                summary(cox_both_int_low)$conf.int[3,4])
)

# Truth ----
df_true_00 <- data.frame(time = xaxis, S = S_true_fun(xaxis, quantile(X_sans.effet, 0.25) * beta_sans.effet), model = "NE=Q1/WE=0")
df_true_10 <- data.frame(time = xaxis, S = S_true_fun(xaxis, quantile(X_sans.effet, 0.75) * beta_sans.effet), model = "NE=Q3/WE=0")
df_true_01 <- data.frame(time = xaxis, S = S_true_fun(xaxis, beta_sans.effet + beta_avec.effet + quantile(X_sans.effet, 0.25) * beta_avec.effet), model = "NE=Q1/WE=1")
df_true_11_int <- data.frame(time = xaxis, S = S_true_fun(xaxis, beta_sans.effet + beta_avec.effet + quantile(X_sans.effet, 0.75) * beta_avec.effet), model = "NE=Q3/WE=1")
# Prep publish ----
publish_int2 <- publish(cox_both_int_low)$regressionTable
publish_int2$Variable <- c("X_ne, X_we=0",
                          "X_ne, X_we=1")
# Expected values : 1 and exp(beta_sans.effet * 1 + beta_avec.effet * 1 + beta_int*(1*1))/exp(beta_sans.effet * 0 + beta_avec.effet * 1 + beta_int*(1*0))