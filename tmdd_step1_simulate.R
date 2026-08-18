library(rxode2)

# --- FULL TMDD MODEL (Mager & Jusko structure) ---
full_tmdd <- rxode2({
  # L  = free drug (central), Lt = drug in tissue
  # R  = free target, RC = drug-target complex
  d/dt(L)  = -(kel + kon*R)*L + koff*RC - kpt*L + ktp*Lt
  d/dt(Lt) =  kpt*L - ktp*Lt
  d/dt(R)  =  ksyn - kdeg*R - kon*L*R + koff*RC
  d/dt(RC) =  kon*L*R - (koff + kint)*RC
})
params <- c(
  kel  = 0.15,   # drug elimination (1/day)
  kpt  = 0.30,   # central -> tissue
  ktp  = 0.20,   # tissue -> central
  kon  = 0.50,   # binding association (1/(nM*day))  <- the fast one
  koff = 0.10,   # dissociation (1/day)
  kint = 0.10,   # complex internalization (1/day)
  ksyn = 0.55,   # target synthesis (nM/day)
  kdeg = 0.11    # target degradation (1/day) -> R0 = ksyn/kdeg = 5 nM
)

# Baseline: system at steady state before drug. R starts at R0 = ksyn/kdeg.
init <- c(L = 100, Lt = 0, R = 5, RC = 0)   # L=100 nM IV bolus at t=0
ev <- et(seq(0, 30, by = 0.1))   # 30 days, fine grid
sim_full <- rxSolve(full_tmdd, params, ev, inits = init)

plot(sim_full$time, sim_full$L, type = "l", log = "y",
     xlab = "Time (days)", ylab = "Free drug L (nM, log scale)",
     main = "Full TMDD: free drug over time")


# --- QSS APPROXIMATION ---
# Tracks TOTAL drug (Ltot) and TOTAL target (Rtot) instead of free species.
# Free drug L is recovered algebraically each step via the binding balance.
qss_tmdd <- rxode2({
  # Kss is the steady-state binding constant: (koff + kint)/kon
  Kss = (koff + kint) / kon
  
  # Solve the quadratic for free drug L given total drug & total target
  Lfree = 0.5 * ((Ltot - Rtot - Kss) +
                   sqrt((Ltot - Rtot - Kss)^2 + 4*Kss*Ltot))
  
  RC   = Rtot * Lfree / (Kss + Lfree)   # complex, from the balance
  Lt2c = Ltot - Lfree                   # (not used in ODE; just bookkeeping)
  
  d/dt(Ltot) = -kel*Lfree - kint*RC - kpt*Lfree + ktp*Ltc
  d/dt(Ltc)  =  kpt*Lfree - ktp*Ltc
  d/dt(Rtot) =  ksyn - kdeg*Rtot - (kint - kdeg)*RC
})
# --- MICHAELIS-MENTEN APPROXIMATION ---
# No target at all. Linear elimination + a saturable Vmax/Km term.
mm_tmdd <- rxode2({
  d/dt(L)  = -kel*L - (Vmax*L)/(Km + L) - kpt*L + ktp*Lt
  d/dt(Lt) =  kpt*L - ktp*Lt
})
# QSS starts with total drug = dose, total target = R0
init_qss <- c(Ltot = 100, Ltc = 0, Rtot = 5)
sim_qss  <- rxSolve(qss_tmdd, params, ev, inits = init_qss)

# MM needs Vmax and Km. Derive them from the binding params:
params_mm <- c(params, Vmax = 0.10*5, Km = (0.10+0.10)/0.50)  # kint*R0 , (koff+kint)/kon
init_mm   <- c(L = 100, Lt = 0)
sim_mm    <- rxSolve(mm_tmdd, params_mm, ev, inits = init_mm)

# Plot full (black), QSS (red dashed), MM (blue dotted)
plot(sim_full$time, sim_full$L, type="l", log="y", lwd=2,
     xlab="Time (days)", ylab="Free drug L (nM, log scale)",
     main="Full vs QSS vs MM")
lines(sim_qss$time, sim_qss$Lfree, col="red",  lty=2, lwd=2)
lines(sim_mm$time,  sim_mm$L,      col="blue", lty=3, lwd=2)
legend("topright", c("Full TMDD","QSS","MM"),
       col=c("black","red","blue"), lty=c(1,2,3), lwd=2)


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
# (uses sim_full from Part 1; times matching the 0.1-day grid are kept)
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
