
# ============================================================
# PART 2 — FITTING STUDY
# Fit a simple 1-compartment model to TMDD-generated data and
# show it is structurally misspecified (systematic residuals).
# ============================================================

library(nlmixr2)

# --- Generate synthetic "observed" data from the FULL model ---
set.seed(123)   # reproducible noise

# Realistic blood-draw times (days)
samp_times <- c(0.25, 0.5, 1, 2, 4, 7, 10, 14, 21, 28)

# Pull the full-model free-drug value at those times
# (uses sim_full from Part 1; times that match the 0.1-day grid are kept)
truth <- sim_full[sim_full$time %in% samp_times, c("time", "L")]

# Add 15% proportional measurement noise
truth$DV <- truth$L * (1 + rnorm(nrow(truth), 0, 0.15))

# Plot: smooth true curve + noisy "observed" dots
plot(sim_full$time, sim_full$L, type = "l", log = "y", lwd = 2,
     xlab = "Time (days)", ylab = "Free drug L (nM, log)",
     main = "Synthetic observed data (dots) vs truth (line)")
points(truth$time, truth$DV, pch = 19, col = "red", cex = 1.3)

# --- Format data for nlmixr2 (NONMEM-style) ---
# ID as integer; dose row uses DV = NA; column order ID,TIME,AMT,DV,EVID,CMT
obs <- data.frame(
  ID = 1L, TIME = truth$time, AMT = 0, DV = truth$DV, EVID = 0, CMT = 1
)
dose <- data.frame(
  ID = 1L, TIME = 0, AMT = 100, DV = NA_real_, EVID = 1, CMT = 1
)
fit_data <- rbind(dose, obs)
fit_data <- fit_data[order(fit_data$TIME, -fit_data$EVID), ]

# --- 1-compartment model, fixed effects only (single subject) ---
# NOTE: SAEM fails with a single subject ("uninformed etas"); with one
# profile, between-subject variability is not estimable. We therefore
# drop the random effect and estimate fixed effects with FOCEi.
# linCmt() also mis-behaved on this single-subject data, so the model is
# written as an explicit ODE instead.
one_cmt_single <- function() {
  ini({
    tCL <- log(15)     # clearance guess (L/day)
    tV  <- log(50)     # volume guess (L)
    prop.err <- 0.15
  })
  model({
    CL <- exp(tCL)
    V  <- exp(tV)
    ke <- CL / V
    d/dt(centr) = -ke * centr      # IV bolus into central compartment
    cp = centr / V                 # concentration = amount / volume
    cp ~ prop(prop.err)
  })
}

fit3 <- nlmixr2(one_cmt_single, fit_data, est = "focei")
print(fit3)

# --- Overlay: observed data (red dots) vs the 1-cmt fit (blue line) ---
pred <- fit3[, c("TIME", "IPRED")]

plot(truth$time, truth$DV, pch = 19, col = "red", cex = 1.3, log = "y",
     xlab = "Time (days)", ylab = "Free drug (nM, log)",
     main = "1-compartment fit to TMDD data: misses early phase")
lines(pred$TIME, pred$IPRED, col = "blue", lwd = 2)
legend("topright", c("TMDD data (truth)", "1-cmt model fit"),
       col = c("red", "blue"), pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2))
