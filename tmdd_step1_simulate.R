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