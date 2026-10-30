# =============================================================================
# 04_estimation.R - Section D: design-based estimation from the realised sample
#
# D1  survey design object: strata, no clusters, weights w_hi = N_h / r_h and
#     the finite population correction (stratum sizes N_h)
# D2-D3 estimate of T with variance, SE, 95% CI, CV and design effect, also with the
#       separate ratio estimator that uses the register's colony counts
# D5  naive (unweighted) and design-weights-only estimates against the truth
# D6  SRSWOR benchmark sample of the same issued size, same response process
# D7  theoretical SE of the stratified estimator for the realised r_h
# =============================================================================

source("R/00_setup.R")

pop    <- read.csv(PATH_POP)
smp    <- read.csv(PATH_SAMPLE)
T_true <- sum(pop$colonies_lost)
T0     <- read.csv(file.path(DIR_TAB, "plan_summary.csv"))$T0
H      <- strata_def$stratum
resp   <- subset(smp, responded == 1)

# ---- D1: design object -----------------------------------------------------------
des <- svydesign(ids = ~1, strata = ~stratum, fpc = ~N_h, weights = ~weight, data = resp)

# ---- D2-D3: estimate of the total ------------------------------------------------
est  <- svytotal(~colonies_lost, des, deff = TRUE)
ci   <- confint(est, level = CONF)                    # normal-theory 95% CI
st <- data.frame(
  method = "Stratified, adjusted weights", estimate = coef(est)[[1]], variance = vcov(est)[1, 1],
  se = SE(est)[[1]], lcl = ci[1, 1], ucl = ci[1, 2], cv = SE(est)[[1]] / coef(est)[[1]],
  deff = deff(est)[[1]], n_resp = nrow(resp))

# formula check: T_hat = sum N_h ybar_h, V = sum N_h^2 (1 - r_h/N_h) s_h^2 / r_h
by_h <- do.call(rbind, lapply(H, function(h) { d <- resp[resp$stratum == h, ]
  data.frame(stratum = h, N_h = d$N_h[1], r_h = nrow(d), ybar = mean(d$colonies_lost), s2 = var(d$colonies_lost)) }))
T_form <- sum(by_h$N_h * by_h$ybar)
V_form <- sum(by_h$N_h^2 * (1 - by_h$r_h / by_h$N_h) * by_h$s2 / by_h$r_h)
stopifnot(abs(T_form - st$estimate) < 1e-6, abs(V_form - st$variance) < 1e-3)

# ---- Separate ratio estimator with the register's colony counts M_h ---------------
# T_R = sum_h M_h ybar_h / mbar_h; variance as for T_hat with S_h^2 replaced by the
# variance of the residuals y - R_h m, times (M_h / Xhat_h)^2 as in svyratio()
M_h <- tapply(pop$colonies_oct, factor(pop$stratum, levels = H), sum)
rs  <- svyratio(~colonies_lost, ~colonies_oct, des, separate = TRUE)
pr  <- predict(rs, total = M_h[names(rs$ratios)])               # svyratio orders strata itself
T_r <- as.numeric(pr$total); se_r <- as.numeric(pr$se)
cf_r <- sapply(H, function(h) { d <- resp[resp$stratum == h, ]; N <- d$N_h[1]; r <- nrow(d)
  R <- mean(d$colonies_lost) / mean(d$colonies_oct); e <- d$colonies_lost - R * d$colonies_oct
  c(T = R * M_h[[h]], V = (M_h[[h]] / (N * mean(d$colonies_oct)))^2 * N^2 * (1 - r / N) * var(e) / r) })
stopifnot(abs(sum(cf_r["T", ]) - T_r) < 1e-6, abs(sum(cf_r["V", ]) - se_r^2) < 1e-3)   # formula check
V_srs_ref <- st$variance / st$deff                               # SRSWOR variance behind svytotal's deff
ra <- data.frame(
  method = "Stratified, ratio estimator", estimate = T_r, variance = se_r^2, se = se_r,
  lcl = T_r - Z * se_r, ucl = T_r + Z * se_r, cv = se_r / T_r, deff = se_r^2 / V_srs_ref, n_resp = nrow(resp))

# ---- D5: naive and design-weights-only estimates ---------------------------------
naive_T  <- N * mean(resp$colonies_lost)                         # ignores the design
design_T <- sum(resp$design_weight * resp$colonies_lost)         # ignores non-response

# ---- D6: SRSWOR benchmark of the same issued size ---------------------------------
set.seed(SEED_BENCH)
n_issued <- nrow(smp)
bench <- pop[sample(N, n_issued), ]
bench$responded <- rbinom(n_issued, 1, RESP_PLAN)
bench$incl_prob <- n_issued / N
bench$weight    <- ifelse(bench$responded == 1, N / sum(bench$responded), NA)
bench$colonies_lost[bench$responded == 0] <- NA
bench$migratory[bench$responded == 0]     <- NA
bench_out <- bench[order(bench$bk_id), c("bk_id", "stratum", "colonies_oct", "incl_prob", "responded", "weight",
                                          "migratory", "colonies_lost")]
write.csv(bench_out, "data/sample_srswor_benchmark.csv", row.names = FALSE)

bresp <- subset(bench, responded == 1)
bresp$N <- N
bdes <- svydesign(ids = ~1, fpc = ~N, weights = ~weight, data = bresp)
best <- svytotal(~colonies_lost, bdes, deff = TRUE)
bci  <- confint(best, level = CONF)
sr <- data.frame(
  method = "SRSWOR benchmark", estimate = coef(best)[[1]], variance = vcov(best)[1, 1],
  se = SE(best)[[1]], lcl = bci[1, 1], ucl = bci[1, 2], cv = SE(best)[[1]] / coef(best)[[1]],
  deff = deff(best)[[1]], n_resp = nrow(bresp))

# ---- D7: theoretical SE for the realised respondents (true S_h, unknown in practice) ---
S2_true <- tapply(pop$colonies_lost, factor(pop$stratum, levels = H), var)
r_h     <- setNames(by_h$r_h, by_h$stratum)[H]
Nh      <- table(factor(pop$stratum, levels = H))
se_theory <- sqrt(sum(Nh^2 * (1 - r_h / Nh) * S2_true / r_h))
S2_all  <- var(pop$colonies_lost); r_b <- nrow(bresp)
se_theory_srs <- sqrt(N^2 * (1 - r_b / N) * S2_all / r_b)

# ---- D8: respondents the plan would have needed with the true S_h (an evaluation;
#      the planner cannot know these). Same rule as Section B: commercial stratum
#      taken in full with floor(R N_h) respondents, Neyman allocation for the rest.
alloc  <- read.csv(file.path(DIR_TAB, "plan_allocation.csv"))
V0     <- (REL_MOE * T0 / Z)^2
r_com  <- floor(RESP_PLAN * Nh[["COM"]])
V_com  <- Nh[["COM"]]^2 * (1 - r_com / Nh[["COM"]]) * S2_true[["COM"]] / r_com
rest   <- c("HOB", "SID"); NS <- Nh[rest] * sqrt(S2_true[rest])
n_rest <- sum(NS)^2 / (V0 - V_com + sum(Nh[rest] * S2_true[rest]))
n_true_plan <- sum(ceiling(n_rest * NS / sum(NS))) + r_com

# ---- Save ---------------------------------------------------------------------------------
estimates <- rbind(st, ra, sr)
estimates$halfwidth     <- Z * estimates$se
estimates$rel_half_T0   <- estimates$halfwidth / T0
estimates$error_pct     <- 100 * (estimates$estimate / T_true - 1)
estimates$covers_true   <- estimates$lcl <= T_true & T_true <= estimates$ucl
write.csv(estimates, file.path(DIR_TAB, "estimates.csv"), row.names = FALSE)

compare <- data.frame(
  method   = c("Naive (unweighted)", "Design weights only", "Adjusted weights (final)", "True total"),
  estimate = c(naive_T, design_T, st$estimate, T_true))
compare$error_pct <- 100 * (compare$estimate / T_true - 1)
write.csv(compare, file.path(DIR_TAB, "estimates_compare.csv"), row.names = FALSE)

by_h$weight <- by_h$N_h / by_h$r_h
by_h$T_hat  <- by_h$N_h * by_h$ybar
by_h$T_R    <- as.numeric(cf_r["T", by_h$stratum])         # ratio estimate by stratum, M_h ybar_h / mbar_h
by_h$T_true <- as.numeric(tapply(pop$colonies_lost, factor(pop$stratum, levels = H), sum))
write.csv(by_h, file.path(DIR_TAB, "estimates_by_stratum.csv"), row.names = FALSE)
write.csv(data.frame(T_true = T_true, T0 = T0, se_theory = se_theory, se_theory_srs = se_theory_srs,
                     se_ratio = sr$se / st$se, var_ratio = st$variance / sr$variance,
                     M = sum(pop$colonies_oct), n_true_plan = n_true_plan,
                     n_plan = sum(alloc$n_resp), V_com_share_true = V_com / V0),
          file.path(DIR_TAB, "estimates_extra.csv"), row.names = FALSE)

cat(sprintf("True T = %d | baseline T0 = %.0f\n", T_true, T0))
print(estimates, digits = 5)
print(compare, digits = 6)
cat(sprintf("Theoretical SE (true S_h, realised r_h): stratified %.0f, SRSWOR %.0f\n", se_theory, se_theory_srs))
