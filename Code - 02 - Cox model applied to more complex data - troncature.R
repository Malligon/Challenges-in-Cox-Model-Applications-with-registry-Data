
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

beta_with.effect <- 1
beta_no.effect <- 0

lambda_cens <- 0.023 # log(2)/30
lambda_trunc <- 0.03 # 1/30

xaxis <- seq(0, 100, 0.01)
times_print <- seq(0, max(xaxis), by = 10)

# Data generation ----

set.seed(global_seed)
U <- runif(N)

X_with.effect <- rbinom(N, 1, 0.5) # WE
X_no.effect <- rnorm(N, 0, 3) # NE

LP_both <- beta_no.effect * X_no.effect + beta_with.effect * X_with.effect

#####################################

S_true_fun_with.effect <- function(t) {exp(-scale * exp(beta_with.effect) * t^shape)}
S_true_with.effect <- S_true_fun_with.effect(xaxis)
S_true_fun_no.effect <- function(t) {exp(-scale * 1 * t^shape)}
S_true_no.effect <- S_true_fun_no.effect(xaxis)

#####################################

surv_time_both <- ((-log(U)) / (scale * exp(LP_both)))^(1/shape)

censor_time_both <- rexp(N, lambda_cens)
trunc_time <- rexp(N, lambda_trunc)

time_both <- pmin(surv_time_both, censor_time_both)
status_both <- as.numeric(surv_time_both <= censor_time_both)

# Truncation
time_both2 <- time_both[time_both >= trunc_time]
status_both2 <- status_both[time_both >= trunc_time]
trunc_time2 <- trunc_time[time_both >= trunc_time]
X_no.effect2 <- X_no.effect[time_both >= trunc_time]
X_with.effect2 <- X_with.effect[time_both >= trunc_time]
data_both <- data.frame(time = time_both2,
                        t_trunc = trunc_time2,
                        status = status_both2,
                        X_no.effect = X_no.effect2,
                        X_with.effect = X_with.effect2)

#####################################

# adjusted cox model
cox_both <- coxph(Surv(t_trunc, time, status) ~ X_no.effect + X_with.effect, data = data_both)

# NE = Q3, WE = 1
est_both_NEQ3_WE1 <- survfit(cox_both, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=1))
df4 <- data.frame(time = est_both_NEQ3_WE1$time, S = est_both_NEQ3_WE1$surv, model = "Adjusted model, NE=Q3 WE=1")

# NE = Q3, WE = 0
est_both_NEQ3_WE0 <- survfit(cox_both, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=0))
df4bis <- data.frame(time = est_both_NEQ3_WE0$time, S = est_both_NEQ3_WE0$surv, model = "Adjusted model, NE=Q3 WE=0")

# NE = Q1, WE = 1
est_both_NEQ1_WE1 <- survfit(cox_both, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.25)), X_with.effect=1))
df4ter <- data.frame(time = est_both_NEQ1_WE1$time, S = est_both_NEQ1_WE1$surv, model = "Adjusted model, NE=Q1 WE=1")

# NE = Q1, WE = 0
est_both_NEQ1_WE0 <- survfit(cox_both, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.25)), X_with.effect=0))
df4quat <- data.frame(time = est_both_NEQ1_WE0$time, S = est_both_NEQ1_WE0$surv, model = "Adjusted model, NE=Q1 WE=0")


#####################################


# Biased Cox model
cox_both_biased <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect, data = data_both)

# NE = Q3, WE = 1
est_both_NEQ3_WE1_biased <- survfit(cox_both_biased, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=1))
df4_biased <- data.frame(time = est_both_NEQ3_WE1_biased$time, S = est_both_NEQ3_WE1_biased$surv, model = "Biased model, NE=Q3 WE=1")

# NE = Q3, WE = 0
est_both_NEQ3_WE0_biased <- survfit(cox_both_biased, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.75)), X_with.effect=0))
df4bis_biased <- data.frame(time = est_both_NEQ3_WE0_biased$time, S = est_both_NEQ3_WE0_biased$surv, model = "Biased model, NE=Q3 WE=0")

# NE = Q1, WE = 1
est_both_NEQ1_WE1_biased <- survfit(cox_both_biased, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.25)), X_with.effect=1))
df4ter_biased <- data.frame(time = est_both_NEQ1_WE1_biased$time, S = est_both_NEQ1_WE1_biased$surv, model = "Biased model, NE=Q1 WE=1")

# NE = Q1, WE = 0
est_both_NEQ1_WE0_biased <- survfit(cox_both_biased, newdata = data.frame(X_no.effect=as.numeric(quantile(data_both$X_no.effect, 0.25)), X_with.effect=0))
df4quat_biased <- data.frame(time = est_both_NEQ1_WE0_biased$time, S = est_both_NEQ1_WE0_biased$surv, model = "Biased model, NE=Q1 WE=0")


HR_both <- exp(coef(cox_both))
HR_both_biased <- exp(coef(cox_both_biased))

invisible(capture.output(tab_both <- publish(cox_both)))
invisible(capture.output(tab_both_biased <- publish(cox_both_biased)))

comparison <- data.frame(
  Modele    = c("Both (adjusted)- NE", "Both (adjusted)- WE", "Both (biased)- NE", "Both (biased)- WE"),
  HR_vrai   = c(exp(beta_no.effect),
                exp(beta_with.effect),
                exp(beta_no.effect),
                exp(beta_with.effect)),
  HR_estime = c(HR_both[1], HR_both[2],
                HR_both_biased[1], HR_both_biased[2]),
  IC_inf    = c(tab_both$rawTable$Lower[1],
                tab_both$rawTable$Lower[2],
                tab_both_biased$rawTable$Lower[1],
                tab_both_biased$rawTable$Lower[2]),
  IC_sup    = c(tab_both$rawTable$Upper[1],
                tab_both$rawTable$Upper[2],
                tab_both_biased$rawTable$Upper[1],
                tab_both_biased$rawTable$Upper[2])
)

# TABLE: 2
comparison <- comparison %>%
  kable(digits = 7, align = "lcc") %>%
  kable_styling(full_width = FALSE) %>%
  pack_rows("Biased model", 1, 2) %>%
  pack_rows("Adjusted model", 3, 4)



fit <- survfit(cox_both, type = "aalen", entry = TRUE)
fit_biased <- survfit(cox_both_biased, type = "aalen", entry = TRUE)

# FIGURE 2 ----

my_colors <- paletteer_d("rcartocolor::Pastel")

at_risk_biased <- summary(fit_biased, times_print)$n.risk
at_risk <- summary(fit, times_print)$n.risk

jpeg(paste(path_cox, "FIGURE2.jpeg", sep = ""), width = 6, height = 3.25, units = "in", res = 1200, pointsize = 4)
par(mar = c(7, 7, 5, 2), xaxs = "i", yaxs = "i", tcl = 0.2, cex.lab = 1.4, cex.axis = 1.2
    # ,mfrow = c(1, 2)
)

plot(1,
     col = NULL,
     ylim = c(0, 1),
     xlim = c(0, 100),
     sub = "",
     xlab = "Time", ylab = "Survival", axes = FALSE,
     panel.first = grid(nx = NULL, ny = NULL,
                        col = rgb(0.8, 0.8, 0.8, 0.7),
                        lty = "dotted", lwd = 1))

lines(xaxis, S_true_no.effect, lty = 1, lwd = 1, col = "black")
lines(xaxis, S_true_with.effect, lty = 3, lwd = 1, col = "black")

polygon(rep(c(est_both_NEQ3_WE1$time, rev(est_both_NEQ3_WE1$time)), each = 2)[-1],
        rep(c(est_both_NEQ3_WE1$upper, rev(est_both_NEQ3_WE1$lower)), each = 2)[-length(c(est_both_NEQ3_WE1$upper, rev(est_both_NEQ3_WE1$lower)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(df4$time, df4$S,
      lwd = 1,
      lty = 3,
      col = my_colors[1],
      type = "s")

polygon(rep(c(est_both_NEQ3_WE0$time, rev(est_both_NEQ3_WE0$time)), each = 2)[-1],
        rep(c(est_both_NEQ3_WE0$upper, rev(est_both_NEQ3_WE0$lower)), each = 2)[-length(c(est_both_NEQ3_WE0$upper, rev(est_both_NEQ3_WE0$lower)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(df4bis$time, df4bis$S,
      lwd = 1,
      lty = 1,
      col = my_colors[1],
      type = "s")

lines(df4_biased$time, df4_biased$S,
      lwd = 1,
      lty = 3,
      col = my_colors[2],
      type = "s")

lines(df4bis_biased$time, df4bis_biased$S,
      lwd = 1,
      lty = 1,
      col = my_colors[2],
      type = "s")

axis(side = 1, at = times_print, lwd = 0.7)
axis(side = 2, at = seq(0, 1, 0.1), lwd = 0.7)

mtext(side = 1, line = 4, at = times_print[1], text = expression(italic("Number at risk")), adj = 0, cex = 1.4)

mtext(side = 1, line = 5, at = -7,
      text = "Adj. est.", adj = 0, cex = 1.2)
mtext(side = 1, line = 6, at = -7,
      text = "Naive est.", adj = 0, cex = 1.2)
for(i in seq_along(times_print)){
  mtext(side = 1, line = 5, at = times_print[i],
        text = at_risk[i], cex = 1.2)
  mtext(side = 1, line = 6, at = times_print[i],
        text = at_risk_biased[i], cex = 1.2)
}

legend("top", inset = c(0, -0.12), xpd = TRUE,
       legend = c("Adjusted estimator", "Naive estimator",
                  "Truth", "", 
                  "X=0", "X=1"),
       col = c(my_colors[1], my_colors[2],
               "black", "white", 
               "black", "black"),
       lty = c(1,1,
               1,1,
               1,3), 
       lwd = 1,
       bty = "n",
       ncol = 3,
       cex = 1.2)

dev.off()
