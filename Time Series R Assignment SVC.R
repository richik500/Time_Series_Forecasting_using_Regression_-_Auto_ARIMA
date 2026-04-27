###############################################################################
# TIME SERIES ANALYSIS ASSIGNMENT
# Dataset: USAccDeaths (Monthly accidental deaths in the USA, 1973-1978)
###############################################################################

#Install packages if not already installed

library(tseries)   # ADF, KPSS tests
library(forecast)  # auto.arima, Acf, Pacf, decompose helpers
library(ggplot2)   # optional pretty plots
library(urca)      # ur.kpss for KPSS test
###############################################################################
# LOAD DATA
###############################################################################

data("USAccDeaths")          # Built-in R dataset
ts_data <- USAccDeaths       # Monthly accidental deaths: Jan 1973 – Dec 1978

# Quick look at the series
print(ts_data)
cat("\nClass:", class(ts_data),
    "\nFrequency:", frequency(ts_data),
    "\nStart:", start(ts_data),
    "\nEnd:", end(ts_data), "\n")

#Base plot
plot(ts_data,
     main = "Monthly Accidental Deaths in the USA (1973–1978)",
     ylab = "Deaths", xlab = "Year", col = "steelblue", lwd = 2)
# INTERPRETATION: Clear seasonal peaks every summer (higher deaths in warm
# months). No obvious long-run upward/downward trend, but variance seems
# roughly constant — mild heteroscedasticity possible.


###############################################################################
# SECTION A: DEEP STATIONARITY ANALYSIS
###############################################################################
#
# A1. Decompose non-stationarity into trend, seasonal, and variance components
#
# Classical additive decomposition (appropriate when seasonal amplitude is
# roughly constant — which is the case here visually)
decomp_add <- decompose(ts_data, type = "additive")
plot(decomp_add)
# Components extracted:
#   Trend     : Slow downward drift from ~8000 to ~7500 deaths/month
#   Seasonal  : Peak in July/August, trough in February — amplitude ~2000
#   Random    : Residuals; irregular, mean ≈ 0, but some clustering visible

# Check whether multiplicative decomposition fits better
decomp_mul <- decompose(ts_data, type = "multiplicative")
plot(decomp_mul)


# Compare random (residual) variance
var_add <- var(decomp_add$random, na.rm = TRUE)
var_mul <- var(decomp_mul$random, na.rm = TRUE)
cat("Residual variance — Additive:", round(var_add, 2),
    "  Multiplicative:", round(var_mul, 6), "\n")
# INTERPRETATION: Additive model is preferred here because the seasonal
# amplitude does NOT grow proportionally with the level.
# Trend component reveals a gentle decline in accidental deaths over the
# 6-year window. The seasonal component shows strong within-year periodicity
# (period = 12). The remainder (random) component captures irregular shocks.

# Seasonal subseries plot — visualise each month's average level
monthplot(ts_data,
          main = "Seasonal Sub-Series Plot (USAccDeaths)",
          ylab = "Deaths", xlab = "Month")
# INTERPRETATION: Confirms July/August are systematically highest;
# February is lowest. The horizontal lines inside each panel show each
# month's mean, revealing systematic seasonal structure.

# Box-plot by month — another view of seasonal variation
month_factor <- cycle(ts_data)
boxplot(ts_data ~ month_factor,
        names = month.abb,
        main  = "Monthly Box-Plots of Accidental Deaths",
        xlab  = "Month", ylab  = "Deaths",
        col   = "lightblue")
# INTERPRETATION: Median deaths in June–August consistently higher.
# Within-month variability is moderate but similar across months (supporting
# additive rather than multiplicative seasonality).


###############################################################################
# A2. Justify whether first differencing and seasonal differencing are required
###############################################################################

# ADF test on the original series (H0: unit root / non-stationary)
adf_orig <- adf.test(ts_data, alternative = "stationary")
print(adf_orig)
# INTERPRETATION: If p-value > 0.05 → fail to reject H0 → unit root present
# (non-stationary in mean). ADF on raw USAccDeaths typically gives p ≈ 0.01
# due to seasonal structure inflating the test — so we check KPSS as well.

# KPSS test on original series (H0: stationary)
kpss_orig <- ur.kpss(as.numeric(ts_data), type = "mu", lags = "short")
summary(kpss_orig)
# INTERPRETATION: If test-stat > critical value → reject H0 → non-stationary.

# Apply FIRST DIFFERENCING (removes linear trend / unit root)
ts_diff1 <- diff(ts_data, differences = 1)
plot(ts_diff1,
     main = "After First Differencing (d=1)",
     ylab = "Δ Deaths", col = "darkgreen", lwd = 1.5)
# The trend is removed; seasonal oscillation is still very visible.

adf_diff1 <- adf.test(ts_diff1, alternative = "stationary")
print(adf_diff1)

# Apply SEASONAL DIFFERENCING (lag=12, removes annual seasonality)
ts_sdiff <- diff(ts_data, lag = 12, differences = 1)
plot(ts_sdiff,
     main = "After Seasonal Differencing (D=1, s=12)",
     ylab = "Seasonal Δ Deaths", col = "darkorange", lwd = 1.5)

adf_sdiff <- adf.test(ts_sdiff, alternative = "stationary")
print(adf_sdiff)

# Apply BOTH differences: first + seasonal (as in Airline Model)
ts_both <- diff(diff(ts_data, lag = 12), differences = 1)
plot(ts_both,
     main = "After Both Differences (d=1, D=1)",
     ylab = "Differenced Deaths", col = "purple", lwd = 1.5)
abline(h = 0, lty = 2, col = "red")

adf_both <- adf.test(ts_both, alternative = "stationary")
print(adf_both)
kpss_both <- ur.kpss(as.numeric(ts_both), type = "mu", lags = "short")
summary(kpss_both)
kpss.test(ts_both)
# THEORETICAL JUSTIFICATION:
# The raw series has:
#   (i)  A slow-moving trend → first differencing (d=1) removes the unit root.
#   (ii) Strong periodic seasonality at lag 12 → seasonal differencing (D=1)
#        removes the seasonal unit root.
# Together, SARIMA(p,1,q)(P,1,Q)[12] is the standard framework.
# After both differences the series should be stationary (ADF: p<0.05,
# KPSS: stat < critical value).


###############################################################################
# A3. Over-differencing: concept and ACF impact
###############################################################################

# Apply an EXTRA first difference (d=2) to show over-differencing
ts_overdiff <- diff(ts_both, differences = 1)  # d=2, D=1 total

par(mfrow = c(2, 1))
acf(ts_both,    lag.max = 40, main = "ACF: d=1, D=1 (correct differencing)")
acf(ts_overdiff,lag.max = 40, main = "ACF: d=2, D=1 (over-differenced)")
par(mfrow = c(1, 1))

#ACF & PACF Analysis
par(mfrow=c(1,2))
acf(ts_both, main="ACF Plot")
pacf(ts_both, main="PACF Plot")
#Explanation
#ACF helps identify MA terms
#PACF helps identify AR terms

#KEY OBSERVATIONS:
#Lag 2 spike → MA(2)
#Both ACF and PACF show gradual decomposition pattern 
#Interaction exists. ARIMA(1,1,2)

# INTERPRETATION / THEORY:
# Over-differencing occurs when more differences are applied than necessary.
# Consequences on ACF:
#  • The ACF of an over-differenced series shows a large NEGATIVE spike at lag 1
#    (ρ₁ ≈ -0.5 in the extreme). This is because differencing once too often
#    induces a non-invertible MA root at θ = 1 (unit root in the MA polynomial).
#  • The ACF decays extremely quickly but with oscillating signs.
#  • The spectral density near frequency 0 is zero — the series has been
#    "over-smoothed" and important low-frequency information is destroyed.
# RULE: Apply the minimum number of differences needed to achieve stationarity.
# If ACF at lag 1 < -0.5 after differencing, suspect over-differencing.


###############################################################################
# A4. Why ADF and KPSS tests can conflict
###############################################################################

# Run both on the original series
cat("\n=== ADF on raw series ===\n")
print(adf.test(ts_data, alternative = "stationary"))

cat("\n=== KPSS on raw series ===\n")
kpss_raw <- ur.kpss(as.numeric(ts_data), type = "mu", lags = "short")
summary(kpss_raw)

# INTERPRETATION:
# ADF (Augmented Dickey-Fuller):
#   H₀ = unit root (non-stationary) ; H₁ = stationary.
#   Reject H₀ if p < 0.05 → evidence FOR stationarity.
#
# KPSS (Kwiatkowski-Phillips-Schmidt-Shin):
#   H₀ = stationary ; H₁ = unit root.
#   Reject H₀ if test-stat > critical value → evidence AGAINST stationarity.
#
# Possible conflict scenarios:
# (a) ADF: fail to reject H₀  +  KPSS: reject H₀
#     → Both agree: non-stationary. (Most common, no conflict.)
# (b) ADF: reject H₀  +  KPSS: fail to reject H₀
#     → Both agree: stationary.
# (c) ADF: reject H₀  +  KPSS: reject H₀
#     → CONFLICT. The series is "trend-stationary" (stationary around a
#       deterministic trend). ADF rejects because the deterministic trend
#       absorbs the power of the test, while KPSS detects the non-constant mean.
# (d) ADF: fail to reject H₀  +  KPSS: fail to reject H₀
#     → CONFLICT. Low power — small sample or the series is near-integrated.
#
# option B is what we interpret from the console output.Both agree: stationary.
# For USAccDeaths: the seasonal structure can confuse ADF (seasonal unit root
# looks like a regular unit root). KPSS is sensitive to the shifting seasonal
# mean. Resolution: use OCSB or CH seasonal unit-root tests in addition.


###############################################################################
# SECTION B: ACF–PACF AND STRUCTURAL IDENTIFICATION
###############################################################################

# Plot ACF and PACF of the doubly-differenced series
par(mfrow = c(2, 1))
Acf(ts_both,  lag.max = 40, main = "ACF of Doubly-Differenced USAccDeaths")
Pacf(ts_both, lag.max = 40, main = "PACF of Doubly-Differenced USAccDeaths")
par(mfrow = c(1, 1))


###############################################################################
# B1. Expected ACF and PACF patterns for MA(1) × SMA(1)
###############################################################################

# MATHEMATICAL DERIVATION:
# The SARIMA(0,1,1)(0,1,1)[12] model has MA component: (1 + θB)(1 + ΘB¹²)
# where θ = non-seasonal MA parameter, Θ = seasonal MA parameter.
#
# After expanding: εₜ + θεₜ₋₁ + Θεₜ₋₁₂ + θΘεₜ₋₁₃
#
# For a pure MA(q) process, ACF cuts off after lag q and PACF decays
# exponentially. For the product MA(1)×SMA(1):
#
# Non-zero ACF autocorrelations (theoretical):
#   ρ(1)  = θ(1 + Θ²) / [(1+θ²)(1+Θ²)]  → spike at lag 1
#   ρ(11) = -θΘ / [(1+θ²)(1+Θ²)]         → small spike at lag 11
#   ρ(12) = Θ(1 + θ²) / [(1+θ²)(1+Θ²)]  → spike at lag 12
#   ρ(13) = θΘ / [(1+θ²)(1+Θ²)]          → spike at lag 13
#   ρ(k)  = 0  for all other k > 0
#
# So ACF has exactly 4 non-zero spikes: lags 1, 11, 12, 13.
# PACF shows exponential decay at seasonal and non-seasonal lags (infinite AR).

cat("\nB1: Theoretical non-zero ACF lags for MA(1)×SMA(1):\n")
cat("Lags with non-zero ACF: 1, 11, 12, 13\n")
cat("ACF at all other lags: 0 (theoretical cutoff)\n")
cat("PACF: exponential / geometric decay — no sharp cutoff\n")


# ─────────────────────────────────────────────────────────────────────────────
# B2. Why lag 13 appears in the ACF
# ─────────────────────────────────────────────────────────────────────────────

# MATHEMATICAL PROOF:
# The moving average polynomial is:
#   Θ(B) = (1 + θB)(1 + ΘB¹²)
#         = 1 + θB + ΘB¹² + θΘB¹³
#
# The autocovariance of a MA process γ(k) = σ² × Σ θⱼθⱼ₊ₖ
# where the θ's are the MA coefficients.
# Coefficients here: θ₀=1, θ₁=θ, θ₁₂=Θ, θ₁₃=θΘ, all others = 0.
#
# γ(13) = σ² × (θ₀×θ₁₃ + θ₁×θ₁₂+θ₁₂×θ₁ + θ₁₃×θ₀)
#       = σ² × (1×θΘ + θ×Θ)   [only cross-products at distance 13]
#   Wait — simpler: γ(k) = σ² Σⱼ θⱼ θⱼ₊ₖ
#   At k=13: j=0 → θ₀θ₁₃ = 1×θΘ = θΘ
#   No other j gives both θⱼ≠0 and θⱼ₊₁₃≠0.
#   Therefore γ(13) = σ²θΘ ≠ 0  provided θ≠0 and Θ≠0.
#
# Conclusion: The INTERACTION between the non-seasonal MA(1) shock at lag 1
# and the seasonal MA shock at lag 12 propagates to create a cross-correlation
# at lag 1+12 = 13. This is the hallmark of multiplicative SARIMA models.

cat("\nB2: Lag 13 in ACF arises because:\n")
cat("  MA polynomial = (1+θB)(1+ΘB^12) = 1 + θB + ΘB^12 + θΘB^13\n")
cat("  The θΘ coefficient at B^13 gives γ(13) = σ²θΘ ≠ 0\n")
cat("  This is the product of the non-seasonal and seasonal MA parameters.\n")

# Visualise: compute ACF and highlight lags 1, 12, 13
acf_vals <- acf(ts_both, lag.max = 40, plot = FALSE)
plot(acf_vals, main = "ACF — Highlighting key lags (1, 12, 13)",
     col = "grey50")
abline(v = c(1, 12, 13)/12, col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("Key lags: 1, 12, 13"), col = "red",
       lty = 2, lwd = 2)


# ─────────────────────────────────────────────────────────────────────────────
# B3. Prove MA has finite cutoff; AR has exponential decay
# ─────────────────────────────────────────────────────────────────────────────

# MATHEMATICAL PROOF:
#
# --- MA(q): Finite ACF cutoff ---
# An MA(q) process: Xₜ = εₜ + θ₁εₜ₋₁ + ... + θqεₜ₋q
# Autocovariance: γ(k) = σ² Σⱼ₌₀^{q-k} θⱼθⱼ₊ₖ   for k ≤ q
#                       = 0                          for k > q
# Because when k > q there are NO overlapping indices j and j+k both
# in {0,1,...,q}. Hence ACF ρ(k) = 0 for k > q. ■ (FINITE CUTOFF)
#
# --- AR(p): Exponential ACF decay ---
# An AR(p): Xₜ = φ₁Xₜ₋₁ + ... + φpXₜ₋p + εₜ
# The Yule-Walker equations give:
#   ρ(k) = φ₁ρ(k-1) + φ₂ρ(k-2) + ... + φpρ(k-p)   for k ≥ 1
# This is a p-th order linear recurrence. Its solution is:
#   ρ(k) = Σᵢ Aᵢ rᵢᵏ
# where rᵢ are roots of the characteristic polynomial (= reciprocals of
# AR polynomial roots). For a STATIONARY AR(p) all |rᵢ| < 1, so:
#   |ρ(k)| → 0 exponentially as k → ∞. ■ (EXPONENTIAL DECAY)
#
# Implication for model identification:
#   ACF cuts off after q  → MA(q)
#   PACF cuts off after p → AR(p)
#   Both decay slowly    → ARMA(p,q)

# Simulate and plot to verify
set.seed(42)
ma1_sim <- arima.sim(model = list(ma = c(0.7)), n = 300)
ar1_sim <- arima.sim(model = list(ar = c(0.7)), n = 300)

par(mfrow = c(2, 2))
acf(ma1_sim,  main = "MA(1): ACF — cuts off at lag 1")
acf(ar1_sim,  main = "AR(1): ACF — exponential decay")
pacf(ma1_sim, main = "MA(1): PACF — exponential decay")
pacf(ar1_sim, main = "AR(1): PACF — cuts off at lag 1")
par(mfrow = c(1, 1))

