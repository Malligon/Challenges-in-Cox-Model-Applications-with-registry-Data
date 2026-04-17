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
require(data.table)
require(parallel)

path_cox <- "Your path" # Path to the file where figures will be saved, needs to finish with a /

# Parameters ----

# seed
if (!exists("global_seed")) global_seed <- 1024

set.seed(global_seed)

N <- 1000

shape <- 2
scale <- 30
shape_event <- 2.5
scale_event <- 20

lambda_cens <- 0.023 # log(2)/30
xaxis <- seq(0, 100, 0.01)
times_print <- seq(0, max(xaxis), by = 10)

randef = rep(1, N) # when simulating a random effect, to cancel, randef = 1
# randef = rlnorm(N, 0, 0.7)
# randef = rep(1.7, N)
# randef = runif(N, 1, 2)
beta_we = 1
beta_ne = 0




# Simulation ----
set.seed(global_seed)

C <- rexp(N, lambda_cens)
T <- rweibull(N, shape, scale)

status <- T <= C
Tobs <- pmin(T, C)

we <- rbinom(N, 1, 0.5)
ne <- rnorm(N, 0, 3)

LP = beta_we*we + beta_ne*ne

U <- runif(N)
E1 <- scale_event*(-log(U)/(randef*exp(LP)))^(1/shape_event) # note that -log(runif(N)) = rexp(N, 1)
j = 1
nbevent <- rep(NA, N)
recdata <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(recdata) <- c("id", "start", "stop", "event", "terminal", "we", "ne")
for (i in 1:N)
{
  stop = E1[i]
  start = 0
  while(stop < Tobs[i])
  {
    recdata[j, ] <- c(i, start, stop, 1, 0, we[i], ne[i])
    start = stop
    stop = scale_event*((stop/scale_event)^(shape_event)-log(runif(1))/(randef[i]*exp(LP[i])))^(1/shape_event)
    # (stop/scale_event)^(shape_event) corresponds to the classical (t/scale)^shape where t = stop, meaning (stop/scale_event)^(shape_event) is equal to the cumulative hazard at stop aka last recurrent event.
    # Also, log(runif(1))/(randef[i]*exp(LP[i])) corresponds to the increment which will define the new cumulative hazard at the new stop of Ei.
    # because H(stop + t) - H(stop) = -log(U)/(randef*LP) donc H(stop + t) = H(stop) - log(U)/(randef*LP)
    # then, the inversion method is t = scale * H^(1/shape)
    j = j+1
  }
  nbevent[i] <- j-1
  recdata[j,] <- c(i, start, Tobs[i], 0, status[i], we[i], ne[i])
  j = j+1
}

# recdata <- fread("recdata.csv")


# Truth ----

h0 <- function(t) {
  (shape_event/scale_event) * (t/scale_event)^(shape_event-1)
}

h1 <- function(t, randef) {
  h0(t) * randef * exp(beta_we)
}

S_term <- function(t) {
  exp(-(t/scale)^shape)
}


integral_h1 <- function(u, randef) {
  integrate(f = function(t) h1(t, randef) * S_term(t), lower = 0, upper = u)$value
}

cl <- makeCluster(detectCores() - 2) # Use all available cores minus 2
clusterExport(cl, c("N", "h0", "h1", "S_term", "integral_h1", "xaxis", "randef", "we", "beta_we", "shape", "scale", "shape_event", "scale_event"))

N_event0 <- parSapply(cl, xaxis, function(u) {
  mean(sapply(which(we == 0), function(i) {
    integrate(f = function(t) h0(t) * randef[i] * S_term(t), lower = 0, upper = u)$value
  }))
})


N_event1 <- parSapply(cl, xaxis, function(u) {
  mean(sapply(which(we == 1), function(i) {
    integral_h1(u, randef[i])
  }))
})


stopCluster(cl)

# Recurrent event: estimation ----

Eventfit <- coxph(Surv(start, stop, event) ~ ne + we + cluster(id), data = recdata)
expecEvent <- survfit(Eventfit, type = "aalen")

Eventfit2 <- coxph(Surv(start, stop, event) ~ ne + we, data = recdata)
expecEvent2 <- survfit(Eventfit2, type = "aalen")

Terminalfit <- coxph(Surv(start, stop, terminal) ~ 1, data = recdata, cluster = id)
survTerm <- survfit(Terminalfit)

invisible(capture.output(pub_Eventfit <- publish(Eventfit)))
pub_Eventfit$regressionTable$Variable <- c("X_ne", "X_we")

invisible(capture.output(pub_Eventfit2 <- publish(Eventfit2)))
pub_Eventfit2$regressionTable$Variable <- c("X_ne", "X_we")

First_event <- coxph(Surv(start, stop, event) ~ ne + we, data = recdata[which(recdata$start == 0), ])

# Individuals ----

surv <- stepfun(survTerm$time, c(1, survTerm$surv))

bh <- basehaz(Eventfit, centered = FALSE)

H_we0q1 <- bh$hazard * exp(quantile(ne, 0.25) * coef(Eventfit)["ne"])
H_we0q3 <- bh$hazard * exp(quantile(ne, 0.75) * coef(Eventfit)["ne"])
H_we1q1 <- bh$hazard * exp(coef(Eventfit)["we"] + quantile(ne, 0.25) * coef(Eventfit)["ne"])
H_we1q3 <- bh$hazard * exp(coef(Eventfit)["we"] + quantile(ne, 0.75) * coef(Eventfit)["ne"])

dH_we0q1 <- c(H_we0q1[1], diff(H_we0q1))
dH_we0q3 <- c(H_we0q3[1], diff(H_we0q3))
dH_we1q1 <- c(H_we1q1[1], diff(H_we1q1))
dH_we1q3 <- c(H_we1q3[1], diff(H_we1q3))

MNRE_0q1 <- cumsum(surv(bh$time) * dH_we0q1)
MNRE_0q3 <- cumsum(surv(bh$time) * dH_we0q3)
MNRE_1q1 <- cumsum(surv(bh$time) * dH_we1q1)
MNRE_1q3 <- cumsum(surv(bh$time) * dH_we1q3)

compute_MNRE_CI <- function(Eventfit, survTerm, ne_quantile, we_value) {
  # Extract baseline hazard and model info
  bh <- basehaz(Eventfit, centered = FALSE)
  beta <- coef(Eventfit)
  vcov_mat <- vcov(Eventfit)

  # Define linear predictor
  lp <- beta["we"] * we_value + beta["ne"] * ne_quantile
  exp_lp <- exp(lp)

  # Gradient for delta method
  grad <- exp_lp * c(we = we_value, ne = ne_quantile)

  # Estimated cumulative hazard and its increments
  H_hat <- bh$hazard * exp_lp
  dH_hat <- c(H_hat[1], diff(H_hat))

  # Standard error of cumulative hazard at each time (delta method)
  ne_H <- sapply(1:nrow(bh), function(i) {
    h0 <- bh$hazard[i]
    sqrt(t(grad) %*% vcov_mat[c("we", "ne"), c("we", "ne")] %*% grad) * h0
  })
  d_ne_H <- c(ne_H[1], diff(ne_H))

  # Survival function from survTerm
  surv_fun <- stepfun(survTerm$time, c(1, survTerm$surv))
  S_t <- surv_fun(bh$time)

  # Compute mean number of recurrent events
  MNRE <- cumsum(S_t * dH_hat)

  # Compute standard error of MNRE (delta method)
  ne_dMNRE <- S_t * d_ne_H
  ne_MNRE <- sqrt(cumsum(ne_dMNRE^2))

  # Return as data.frame
  return(data.frame(
    time = bh$time,
    MNRE = MNRE,
    we = we_value,
    ne_quantile = ne_quantile
  ))
}

# For we = 0, ne at Q1
MNRE_0q1 <- compute_MNRE_CI(Eventfit, survTerm, quantile(ne, 0.25), 0)
MNRE_0q3 <- compute_MNRE_CI(Eventfit, survTerm, quantile(ne, 0.75), 0)
MNRE_1q1 <- compute_MNRE_CI(Eventfit, survTerm, quantile(ne, 0.25), 1)
MNRE_1q3 <- compute_MNRE_CI(Eventfit, survTerm, quantile(ne, 0.75), 1)

# CI bootstrap ----

recdatabis <- setDT(recdata)



N <- 1000

N_boot <- 1000



bootstrap_list_0q1 <- vector("list", N_boot)

bootstrap_list_0q3 <- vector("list", N_boot)

bootstrap_list_1q1 <- vector("list", N_boot)

bootstrap_list_1q3 <- vector("list", N_boot)



Event_time <- unique(recdata$stop)[order(unique(recdata$stop))]



ids <- unique(recdata$id)



for (i in seq_len(N_boot)) {

  set.seed(2353843 + i)

  boot_ids <- sample(ids, N, replace = TRUE)

  recdata_boot <- recdatabis[id %in% boot_ids]



  Eventfit_boot <- coxph(Surv(start, stop, event) ~ ne + we + cluster(id), data = recdata_boot)

  expecEvent_boot <- survfit(Eventfit_boot, type = "aalen")



  Terminalfit_boot <- coxph(Surv(start, stop, terminal) ~ 1, data = recdata_boot, cluster = id)

  survTerm_boot <- survfit(Terminalfit_boot)



  surv_boot <- stepfun(survTerm_boot$time, c(1, survTerm_boot$surv))



  bh_boot <- basehaz(Eventfit_boot, centered = FALSE)



  H_we0q1_boot <- bh_boot$hazard * exp(quantile(ne, 0.25) * coef(Eventfit_boot)["ne"])

  H_we0q3_boot <- bh_boot$hazard * exp(quantile(ne, 0.75) * coef(Eventfit_boot)["ne"])

  H_we1q1_boot <- bh_boot$hazard * exp(coef(Eventfit_boot)["we"] + quantile(ne, 0.25) * coef(Eventfit_boot)["ne"])

  H_we1q3_boot <- bh_boot$hazard * exp(coef(Eventfit_boot)["we"] + quantile(ne, 0.75) * coef(Eventfit_boot)["ne"])



  dH_we0q1_boot <- c(H_we0q1_boot[1], diff(H_we0q1_boot))

  dH_we0q3_boot <- c(H_we0q3_boot[1], diff(H_we0q3_boot))

  dH_we1q1_boot <- c(H_we1q1_boot[1], diff(H_we1q1_boot))

  dH_we1q3_boot <- c(H_we1q3_boot[1], diff(H_we1q3_boot))



  MNRE_0q1_boot <- cumsum(surv(bh_boot$time) * dH_we0q1_boot)

  MNRE_0q3_boot <- cumsum(surv(bh_boot$time) * dH_we0q3_boot)

  MNRE_1q1_boot <- cumsum(surv(bh_boot$time) * dH_we1q1_boot)

  MNRE_1q3_boot <- cumsum(surv(bh_boot$time) * dH_we1q3_boot)



  MNRE_0q1_boot <- compute_MNRE_CI(Eventfit_boot, survTerm_boot, quantile(ne, 0.25), 0)

  MNRE_0q3_boot <- compute_MNRE_CI(Eventfit_boot, survTerm_boot, quantile(ne, 0.75), 0)

  MNRE_1q1_boot <- compute_MNRE_CI(Eventfit_boot, survTerm_boot, quantile(ne, 0.25), 1)

  MNRE_1q3_boot <- compute_MNRE_CI(Eventfit_boot, survTerm_boot, quantile(ne, 0.75), 1)



  MNRE_0q1_boot <- stepfun(MNRE_0q1_boot$time, c(0, MNRE_0q1_boot$MNRE))

  MNRE_0q3_boot <- stepfun(MNRE_0q3_boot$time, c(0, MNRE_0q3_boot$MNRE))

  MNRE_1q1_boot <- stepfun(MNRE_1q1_boot$time, c(0, MNRE_1q1_boot$MNRE))

  MNRE_1q3_boot <- stepfun(MNRE_1q3_boot$time, c(0, MNRE_1q3_boot$MNRE))



  bootstrap_list_0q1[[i]] <- MNRE_0q1_boot(Event_time)

  bootstrap_list_0q3[[i]] <- MNRE_0q3_boot(Event_time)

  bootstrap_list_1q1[[i]] <- MNRE_1q1_boot(Event_time)

  bootstrap_list_1q3[[i]] <- MNRE_1q3_boot(Event_time)



  print(i)

}

bootstrap_0q1 <- as.data.frame(do.call(rbind, bootstrap_list_0q1))
bootstrap_0q3 <- as.data.frame(do.call(rbind, bootstrap_list_0q3))
bootstrap_1q1 <- as.data.frame(do.call(rbind, bootstrap_list_1q1))
bootstrap_1q3 <- as.data.frame(do.call(rbind, bootstrap_list_1q3))

CI_low_0 <- apply(bootstrap_0q3, 2, quantile, 0.025)
CI_up_0 <- apply(bootstrap_0q3, 2, quantile, 0.975)

CI_low_1 <- apply(bootstrap_1q3, 2, quantile, 0.025)
CI_up_1 <- apply(bootstrap_1q3, 2, quantile, 0.975)

Event_time <- unique(recdata$stop)[order(unique(recdata$stop))]

# FIGURE4 ----

at_risk <- summary(survTerm, times_print)$n.risk

my_colors <- paletteer_d("rcartocolor::Pastel")

jpeg(paste(path_cox, "FIGURE4.jpeg", sep = ""), width = 6, height = 3.25, units = "in", res = 1200, pointsize = 4)
par(mar = c(7, 7, 5, 2), xaxs = "i", yaxs = "i", tcl = 0.2, cex.lab = 1.4, cex.axis = 1.2
    # ,mfrow = c(1, 2)
)

plot(1,
     col = NULL,
     ylim = c(0, 10),
     xlim = c(0, 70),
     sub = "",
     xlab = "Time", ylab = "Number of events", axes = FALSE,
     panel.first = grid(nx = NULL, ny = NULL,
                        col = rgb(0.8, 0.8, 0.8, 0.7),
                        lty = "dotted", lwd = 1))

lines(xaxis, N_event0, lty = 1, lwd = 1, col = "black")
lines(xaxis, N_event1, lty = 3, lwd = 1, col = "black")

axis(side = 1, at = times_print, lwd = 0.7)
axis(side = 2, at = seq(0, 15, 1), lwd = 0.7)


polygon(rep(c(Event_time, rev(Event_time)), each = 2)[-1],
        rep(c(CI_up_0, rev(CI_low_0)), each = 2)[-length(c(CI_up_0, rev(CI_low_0)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(bh$time, MNRE_0q3$MNRE,
      lwd = 1,
      lty = 1,
      col = my_colors[1],
      type = "s")

polygon(rep(c(Event_time, rev(Event_time)), each = 2)[-1],
        rep(c(CI_up_1, rev(CI_low_1)), each = 2)[-length(c(CI_up_1, rev(CI_low_1)))*2],
        col = adjustcolor(my_colors[1], alpha.f = 0.2),
        border = NA)
lines(bh$time, MNRE_1q3$MNRE,
      lwd = 1,
      lty = 3,
      col = my_colors[1],
      type = "s")

mtext(side = 1, line = 4, at = times_print[1], text = expression(italic("Number at risk")), adj = 0, cex = 1.4)

mtext(side = 1, line = 5, at = -7,
      text = "", adj = 0, cex = 1.2)
for(i in seq_along(times_print)){
  mtext(side = 1, line = 5, at = times_print[i],
        text = at_risk[i], cex = 1.2)
}

legend("top", inset = c(0, -0.12), xpd = TRUE,
       legend = c("MNRE estimator", "Truth",
                  "X=0", "X=1"),
       col = c(my_colors[1], "black",
               "black", "black"),
       lty = c(1,1,
               1,3), 
       lwd = 1,
       bty = "n",
       ncol = 3,
       cex = 1.2)

dev.off()