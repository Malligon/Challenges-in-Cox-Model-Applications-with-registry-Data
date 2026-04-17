
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
require(prodlim)

path_cox <- "Your path" # Path to the file where figures will be saved, needs to finish with a /

# Constantes ----

if (!exists("global_seed")) global_seed <- 1024

set.seed(global_seed+20)

N <- 1000
shape <- 2
median_surv <- 50
scale <- log(2)/(median_surv^shape)

beta0 <- 0
delta <- 0.05

time_function <- function(t) t # the HR will be multiplied by exp(delta) for each unit of increase in time_function.

beta0_with.effect <- beta0 # WE
beta0_no.effect <- 0 # NE

lambda_cens <- 0.023

xaxis <- seq(0, 100, 0.01)
times_print <- seq(0, max(xaxis), by = 10)

# Function for PH model ----

h0 <- function(t) {
  scale * shape * t^(shape - 1)
}

S0_fun <- function(t) {
  exp(- scale * t^shape)
}

h <- function(t) {
  h0(t) * exp(beta0 + delta * time_function(t))
}

H1_fun <- function(t) {
  integrate(function(u) h(u), lower=0, upper=t)$value
}

S1_fun <- function(t) {
  exp(-H1_fun(t))
}

# Truth ----

S0_values <- S0_fun(xaxis)
S1_values <- sapply(xaxis, S1_fun)

df_true_sans <- data.frame(time = xaxis, S = S0_values)
df_true_avec <- data.frame(time = xaxis, S = S1_values)

# Simulation ----

set.seed(global_seed)

X <- data.frame(id=1:N,
                X_with.effect = rbinom(N, 1, 0.5),
                X_no.effect = rnorm(N, 0, 3))

sim <- simsurv(
  lambdas = scale,
  gammas = shape,
  x = X,
  tde = c(X_with.effect = delta), # The increase for each unit of time gained in time_function
  tdefunction = time_function, # function of time
  beta = c(X_no.effect = beta0_no.effect, X_with.effect = beta0_with.effect)
)
sim_data <- merge(sim, X, by = "id")

cens_time <- rexp(N, lambda_cens)
sim_data$status <- as.numeric(sim_data$eventtime <= cens_time)
sim_data$time <- pmin(cens_time, sim_data$eventtime)

# Building of the most optimal Cox model ----

# We propose 5 different time_function to test everytime. Linear, logarithmic, exponential and splines (natural and B)
time_function1 <- function(x, t, ...) x * t
time_function2 <- function(x, t, ...) x * log(t + 1) # We put +1 to manage t=0
time_function3 <- function(x, t, ...) x * exp(t) # We put +1 to manage t=0
time_function4 <- function(x, t, ...) x * ns(t, df = 3)
time_function5 <- function(x, t, ...) x * bs(t, df = 3)

fit1 <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect + tt(X_with.effect),
              data = sim_data, tt = time_function1)
fit2 <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect + tt(X_with.effect),
              data = sim_data, tt = time_function2)
fit3 <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect + tt(X_with.effect),
              data = sim_data, tt = time_function3)
fit4 <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect + tt(X_with.effect),
              data = sim_data, tt = time_function4)
fit5 <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect + tt(X_with.effect),
              data = sim_data, tt = time_function5)

# fit1 ; fit2 ; fit3 ; fit4 ; fit5
# TABLE: 14
# Here we see that one likelihood is far from other: the exponential time function.
table_loglikeAIC <- data.frame("Fonction"= c("Linéaire",
                                             "Logarithmique",
                                             "Exponentielle",
                                             "Natural splines",
                                             "B-splines"),
                               "Likelihood ratio test"=c(round(summary(fit1)$logtest[["test"]], 2),
                                                         round(summary(fit2)$logtest[["test"]], 2),
                                                         round(summary(fit3)$logtest[["test"]], 2),
                                                         round(summary(fit4)$logtest[["test"]], 2),
                                                         round(summary(fit5)$logtest[["test"]], 2)),
                               AIC=AIC(fit1, fit2, fit3, fit4, fit5)$AIC)

# And so here we see this AIC is higher than others. The best is the first one, but everything is pretty close.
# In practice, if splines are very close to a more comprehensive time_function, aka linear, logarithmic or exponential, we advise to use the comprehensive function instead.

# Now that we found our time_function, we can structure the dataset to draw curves

events_times <- sort(unique(sim_data$time[sim_data$status==1]))
sim_data_split <- survSplit(Surv(time, status) ~ id + X_with.effect + X_no.effect, sim_data, cut=events_times)
sim_data_split$X_with.effect_timeD <- time_function1(sim_data_split$X_with.effect, sim_data_split$tstart)

# Cox model
cox_model <- coxph(Surv(tstart, time, status) ~ X_no.effect + X_with.effect + X_with.effect_timeD, 
                   data=sim_data_split, id=id)

invisible(capture.output(publish_cox <- publish(cox_model)$regressionTable))

publish_cox$Variable <- c("X_ne",
                          "X_we",
                          "X_we_timeD")
# survfit cannot be applied on cox models using tt that is why we need to do this step.
# Coefficient of X_with.effect_timeD should be very close to delta since we are using linear function. X_no.effect and X_with.effect should be equal at there respective beta at baseline, here 0.

# fit1 ; cox_model
# Both Cox models are very similar

# Definition of two patients
newdata1 <- data.frame(
  tstart=c(0, events_times[1:length(events_times)-1]),
  time=events_times,
  X_with.effect = rep(0, length(events_times)),
  X_no.effect = rep(0, length(events_times)))

newdata1$status <- 0
newdata1$X_with.effect_timeD <- time_function1(newdata1$X_with.effect, newdata1$tstart) # will be 0 everywhere


newdata2 <- newdata1
newdata2$X_with.effect <- 1
newdata2$X_with.effect_timeD <- time_function1(newdata2$X_with.effect, newdata2$tstart)


newdata1$id <- "1"
newdata2$id <- "2"

# We use conf.type = "log-log" because for small survival, the confidence interval glitches due to asymptotic method.
fit_surv0.lin <- survfit(cox_model, newdata = newdata1, id = newdata1$id, se.fit = TRUE, conf.type = "log-log")
fit_surv1.lin <- survfit(cox_model, newdata = newdata2, id = newdata1$id, se.fit = TRUE, conf.type = "log-log")


cox_model_notime <- coxph(Surv(tstart, time, status) ~ X_no.effect + X_with.effect, 
                          data=sim_data_split, id=id)
fit_surv0.nt <- survfit(cox_model_notime, newdata = newdata1, id = newdata1$id, se.fit = TRUE, conf.type = "log-log")
fit_surv1.nt <- survfit(cox_model_notime, newdata = newdata2, id = newdata1$id, se.fit = TRUE, conf.type = "log-log")

# FIGURE8 ----

my_colors <- paletteer_d("rcartocolor::Pastel")

fit <- survfit(cox_model, type = "aalen")
at_risk <- summary(fit, times_print)$n.risk

jpeg(paste(path_cox, "FIGURE8.jpeg", sep = ""), width = 6, height = 3.25, units = "in", res = 1200, pointsize = 5)
par(mar = c(7, 7, 5, 2), xaxs = "i", yaxs = "i", tcl = 0.2, cex.lab = 1.4, cex.axis = 1.2)

plot(1, col = NULL, lwd = 2,
     xlim = c(0, 100), ylim = c(0, 1),
     xlab = "Time", ylab = "Survival", axes = FALSE,
     # mark.time = FALSE, conf.int = FALSE,
     panel.first = grid(nx = NULL, ny = NULL,
                        col = rgb(0.8, 0.8, 0.8, 0.7),
                        lty = "dotted", lwd = 1))

lines(df_true_sans$time, df_true_sans$S, lwd = 1, lty = 1, col = "black")
lines(df_true_avec$time, df_true_avec$S, lwd = 1, lty = 3, col = "black")

polygon(rep(c(fit_surv0.lin$time, rev(fit_surv0.lin$time)), each = 2)[-1],
        rep(c(fit_surv0.lin$upper, rev(fit_surv0.lin$lower)), each = 2)[-length(c(fit_surv0.lin$upper, rev(fit_surv0.lin$lower)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(fit_surv0.lin$time, fit_surv0.lin$surv, type = "s",
      lwd = 1,
      lty = 1,
      col = my_colors[1])

polygon(rep(c(fit_surv1.lin$time[!is.na(fit_surv1.lin$upper)], rev(fit_surv1.lin$time[!is.na(fit_surv1.lin$upper)])), each = 2)[-1],
        rep(c(fit_surv1.lin$upper[!is.na(fit_surv1.lin$upper)], rev(fit_surv1.lin$lower[!is.na(fit_surv1.lin$upper)])), each = 2)[-length(c(fit_surv1.lin$upper[!is.na(fit_surv1.lin$upper)], rev(fit_surv1.lin$lower[!is.na(fit_surv1.lin$upper)])))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(fit_surv1.lin$time, fit_surv1.lin$surv, type = "s",
      lwd = 1,
      lty = 3,
      col = my_colors[1])


lines(fit_surv0.nt$time, fit_surv0.nt$surv, type = "s",
      lwd = 1,
      lty = 1,
      col = my_colors[2])
lines(fit_surv1.nt$time, fit_surv1.nt$surv, type = "s",
      lwd = 1,
      lty = 3,
      col = my_colors[2])

axis(side = 1, at = times_print, lwd = 0.7)
axis(side = 2, at = seq(0, 1, 0.1), lwd = 0.7)

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

mtext(side = 1, line = 4, at = times_print[1], text = expression(italic("Number at risk")), adj = 0, cex = 1.4)

mtext(side = 1, line = 5, at = times_print[1] - diff(times_print)[1] * 0.5,
      text = rownames(at_risk)[1], adj = 1, cex = 1.2)
for(i in seq_along(times_print)){
  mtext(side = 1, line = 5, at = times_print[i],
        text = at_risk[i], cex = 1.2)
}

dev.off()

# Schoenfeld residuals ----


test.ph_notime <- cox.zph(cox_model_notime)

test.ph <- cox.zph(cox_model)

# Spline time function ----

# When linear, logarithmic or exponential is the right time_function, results are pretty straight forward. The risk will increased linearly, logarithmically or exponentially over time.

# When spline function are really better fits than usual time functions, then the interpretation is more difficult. Indeed, we do not know the form of the splines. Imagine if the following model was the best one, here is the method we would apply:

time_function4 <- function(x, t, ...) x * ns(t, df = 3)
fit4 <- coxph(Surv(time, status) ~ X_no.effect + X_with.effect + tt(X_with.effect),
              data = sim_data, tt = time_function4)
# Here we can see there i a significant time effect on the first spline.

# To draw curves, we advise to use the exact same methodology, wihile keeping the models with the splines:

events_times <- sort(unique(sim_data$time[sim_data$status==1]))
sim_data_split <- survSplit(Surv(time, status) ~ id + X_with.effect + X_no.effect, sim_data, cut=events_times)
sim_data_split$X_with.effect_timeD <- sim_data_split$X_with.effect * time_function4(sim_data_split$X_with.effect, sim_data_split$tstart)

# Cox model
cox_model <- coxph(Surv(tstart, time, status) ~ X_no.effect + X_with.effect + X_with.effect_timeD, 
                   data=sim_data_split, id=id)



# Definition of two patients
newdata1 <- data.frame(
  tstart=c(0, events_times[1:length(events_times)-1]),
  time=events_times,
  X_with.effect = rep(0, length(events_times)),
  X_no.effect = rep(0, length(events_times)))

newdata1$status <- 0
newdata1$X_with.effect_timeD <- time_function4(newdata1$X_with.effect, newdata1$tstart) # will be 0 everywhere


newdata2 <- newdata1
newdata2$X_with.effect <- 1
newdata2$X_with.effect_timeD <- time_function4(newdata2$X_with.effect, newdata2$tstart)


newdata1$id <- "1"
newdata2$id <- "2"

# We use conf.type = "log-log" and the summary because for small survival, the confidence interval glitches due to asymptotic method.
fit_surv0.ns <- survfit(cox_model, newdata = newdata1, id = newdata1$id, se.fit = TRUE, conf.type = "log-log")
fit_surv1.ns <- survfit(cox_model, newdata = newdata2, id = newdata1$id, se.fit = TRUE, conf.type = "log-log")

my_colors <- pal_jco("default", alpha = 1)(6)
my_colors <- paletteer_d("ggthemes::Classic_10_Light")
my_colors <- paletteer_d("rcartocolor::Pastel")

fit <- survfit(cox_model, type = "aalen")
at_risk <- summary(fit, times_print)$n.risk

par(mar = c(7, 5, 5, 2), xaxs = "i", yaxs = "i", tcl = 0.2, cex.lab = 1.4, cex.axis = 1.2)

plot(fit_surv0.ns$time, fit_surv0.ns$surv, col = NULL, lwd = 2,
     xlim = c(0, 60), ylim = c(0, 1),
     xlab = "Time", ylab = "Survival", axes = FALSE,
     # mark.time = FALSE, conf.int = FALSE,
     panel.first = grid(nx = NULL, ny = NULL,
                        col = rgb(0.8, 0.8, 0.8, 0.7),
                        lty = "dotted", lwd = 1))
title("Survie : Théorie vs Modèle de Cox", line = 3)

lines(df_true_sans$time, df_true_sans$S, lwd = 1.5, lty = 1, col = "black")
lines(df_true_avec$time, df_true_avec$S, lwd = 1.5, lty = 3, col = "black")


polygon(c(fit_surv0.ns$time, rev(fit_surv0.ns$time)),
        c(fit_surv0.ns$upper, rev(fit_surv0.ns$lower)),
        col = adjustcolor(my_colors[1], alpha.f = 0.1),
        border = my_colors[1])
lines(fit_surv0.ns$time, fit_surv0.ns$surv, type = "s",
      lwd = 1.5,
      lty = 1,
      col = my_colors[1])


polygon(rep(c(fit_surv1.ns$time[!is.na(fit_surv1.ns$upper)], rev(fit_surv1.ns$time[!is.na(fit_surv1.ns$upper)])), each = 2)[-1],
        rep(c(fit_surv1.ns$upper[!is.na(fit_surv1.ns$upper)], rev(fit_surv1.ns$lower[!is.na(fit_surv1.ns$upper)])), each = 2)[-length(c(fit_surv1.ns$upper[!is.na(fit_surv1.ns$upper)], rev(fit_surv1.ns$lower[!is.na(fit_surv1.ns$upper)]))) * 2],
        col = adjustcolor(my_colors[2], alpha.f = 0.1),
        border = my_colors[2])
lines(fit_surv1.ns$time, fit_surv1.ns$surv, type = "s",
      lwd = 1.5,
      lty = 1,
      col = my_colors[2])

legend("top",
       inset = c(0, -0.08),
       xpd = TRUE,
       legend = c("Cox: Xwe = 0", "Cox: Xwe = 1", "Truth Xwe = 0", "Truth Xwe = 1"),
       col = c(my_colors[1:2], "black", "black"),
       lty = c(1, 1, 1, 3), 
       lwd = 1.5,
       bty = "n",
       ncol = 2)

axis(side = 1, at = times_print, lwd = 0.7)
axis(side = 2, at = seq(0, 1, 0.1), lwd = 0.7)

mtext(side = 1, line = 4, at = times_print[1], text = expression(italic("Number at risk")), adj = 0, cex = 1.4)

mtext(side = 1, line = 5, at = times_print[1] - diff(times_print)[1] * 0.5,
      text = rownames(at_risk)[1], adj = 1, cex = 1.4)
for(i in seq_along(times_print)){
  mtext(side = 1, line = 5, at = times_print[i],
        text = at_risk[i], cex = 1.4)
}

# Here we can see that NS is less accurate than the linear function, which is expected seeing the difference between likelihood and AIC.

# Now, in order to present comprehensive HR, we advise praticians to do a Piecewise model. This model will force the HR to be constant on interval of time. It will allow us to show an increasing/decreasing effect over time. Note this method should be used to make HR comprehensive only, one should not draw curves based on this method.

sim_data_split <- survSplit(
  data = sim_data_split,
  cut = c(10, 20, 30, 40),
  start = "tstart",
  end = "time",
  event = "status",
  id = "split_id"
)

sim_data_split$time_interval <- cut(sim_data_split$tstart,
                                    breaks = c(0, 10, 20, 30, 40, Inf),
                                    labels = c("0-10", 
                                               "10-20",
                                               "20-30",
                                               "30-40",
                                               "40+"))

cox_PW.ns <- coxph(Surv(tstart, time, status) ~ X_no.effect + X_with.effect * time_interval, data=sim_data_split, id=id)

invisible(capture.output(publish_cox_PW <- publish(cox_PW.ns)$regressionTable))

publish_cox_PW$Variable <- c("X_ne",
                             "X_we: (0-10)",
                             "X_we: (10-20)", 
                             "X_we: (20-30)",
                             "X_we: (30-40)", 
                             "X_we: (40+)")

# publish(cox_PW.ns)
# CI cannot be calculated using publish package.We need to calculate them manually.

# Values should be relatively close to exp(mean(interval) * delta) because we used a linear time function for simulation.
exp(5*delta)
exp(15*delta)
exp(25*delta)
exp(35*delta)
exp(45*delta)

intervals <- levels(sim_data_split$time_interval)

# Extract baseline term
beta_base <- cox_PW.ns$coefficients["X_with.effect"]
var_base <- vcov(cox_PW.ns)["X_with.effect", "X_with.effect"]


# Initialize result table
hr_df <- tibble(
  interval = intervals,
  logHR = NA_real_,
  se = NA_real_
)


for (i in 1:length(intervals)) {
  interval <- intervals[i]
  
  if (interval == "0-10") {
    logHR <- beta_base
    se <- sqrt(var_base)
  } else {
    interaction_term <- paste0("X_with.effect:time_interval", interval)
    
    # Check if the interaction exists in model
    if (interaction_term %in% rownames(vcov(cox_PW.ns))) {
      beta_inter <- coef(cox_PW.ns)[interaction_term]
      var_inter <- vcov(cox_PW.ns)[interaction_term, interaction_term]
      covar <- vcov(cox_PW.ns)["X_with.effect", interaction_term]
      
      logHR <- beta_base + beta_inter
      se <- sqrt(var_base + var_inter + 2 * covar)
    } else {
      logHR <- NA
      se <- NA
    }
  }
  
  hr_df$logHR[i] <- logHR
  hr_df$se[i] <- se
}

hr_df <- hr_df %>%
  mutate(
    HR = exp(logHR),
    CI_lower = exp(logHR - 1.96 * se),
    CI_upper = exp(logHR + 1.96 * se)
  )

# Add expected HRs to hr_df
hr_df <- hr_df %>%
  mutate(
    midpoint = c(5, 15, 25, 35, 50),
    expected_HR = exp(midpoint * delta)
  )
# pval
beta <- coef(cox_PW.ns)
V    <- vcov(cox_PW.ns)
main_name <- "X_with.effect"

L <- beta[main_name]
SE <- sqrt(V[main_name, main_name])
z <- L / SE
pval <- 2 * (1 - pnorm(abs(z)))



main_name <- "X_with.effect"
int_name  <- "X_with.effect:time_interval40+"   # use the exact name printed above

lc_vec <- rep(0, length(beta))
names(lc_vec) <- names(beta)
lc_vec[main_name] <- 1
lc_vec[int_name]  <- 1

L <- sum(beta * lc_vec, na.rm = TRUE)                  # log(HR) for X at 10-20
SE <- sqrt(t(lc_vec) %*% V %*% lc_vec)
z <- L / SE
pval <- 2 * (1 - pnorm(abs(z)))

# TABLE: 16
hr_df$interval <- factor(hr_df$interval,
                         levels = hr_df$interval)

# FIGURE9 ----
# GRAPH: 12
plot_df <- hr_df %>%
  select(interval, HR, CI_lower, CI_upper, expected_HR) %>%
  pivot_longer(cols = c(HR, expected_HR), names_to = "Type", values_to = "Value") %>%
  mutate(
    Type = dplyr::recode(Type, 
                  HR = "Estimated HR", 
                  expected_HR = "Expected HR")
  )
# Plot with legend
forestplot <- ggplot() +
  geom_errorbar(
    data = hr_df,
    aes(x = interval, ymin = CI_lower, ymax = CI_upper),
    width = 0.2
  ) +
  geom_point(
    data = plot_df,
    aes(x = interval, y = Value, color = Type, shape = Type),
    size = 4
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgray", size = 1) +
  scale_y_log10() +
  scale_color_manual(values = c("Estimated HR" = "black", "Expected HR" = "red")) +
  scale_shape_manual(values = c("Estimated HR" = 16, "Expected HR" = 18)) +
  labs(
    title = "Time-varying effect of X",
    x = "Time interval",
    y = "Hazard Ratio (log scale)",
    color = "Type",
    shape = "Type"
  ) +
  theme_minimal() + 
  coord_flip()

ggsave(paste(path_cox, "FIGURE9.jpeg", sep = ""), forestplot,
       dpi = 600,
       width = 30,
       height = 20,
       units = "cm")


# Validation Schoenfeld piecewise ----

test.ph_PW <- cox.zph(cox_PW.ns)
