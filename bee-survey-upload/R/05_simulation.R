# =============================================================================
# 05_simulation.R - Section E (and D7): Monte Carlo evaluation, R = 10,000
#
# Each repetition repeats the complete procedure on the fixed population:
#   stratified design: SRSWOR of n_h* = 1,565 / 567 / 320 in each stratum,
#                      response with probability RESP_TRUE[h], estimator
#                      T_hat = sum N_h ybar_h and variance sum N_h^2 (1 - r_h/N_h) s_h^2 / r_h,
#                      and from the same sample the separate ratio estimator
#                      T_R = sum M_h ybar_h / mbar_h with the svyratio() variance
#   SRSWOR benchmark:  SRSWOR of the same issued size, the same response process,
#                      T_hat = N ybar and variance N^2 (1 - r/N) s^2 / r
# These are the estimators svytotal() and svyratio() give in 04_estimation.R.
# For each repetition the estimate, estimated SE and 95% CI limits are kept.
# D8 context: the stratified procedure is also repeated under a plan built on the
# true S_h, to show how often the 5% bound is met when the plan is right, and with
# larger hobby and side-line samples, to find the size that meets it in 80% or 90%.
# Limitations: the procedure is also repeated with response that depends on the losses.
# =============================================================================

source("R/00_setup.R")
set.seed(SEED_SIM)

R_SIM  <- 10000
pop    <- read.csv(PATH_POP)
alloc  <- read.csv(file.path(DIR_TAB, "plan_allocation.csv"))
T0     <- read.csv(file.path(DIR_TAB, "plan_summary.csv"))$T0
H      <- strata_def$stratum
T_true <- sum(pop$colonies_lost)
y_h    <- split(pop$colonies_lost, factor(pop$stratum, levels = H))
m_h    <- split(pop$colonies_oct,  factor(pop$stratum, levels = H))   # register colonies
N_h    <- sapply(y_h, length)
M_h    <- sapply(m_h, sum)
n_iss  <- setNames(alloc$issued, alloc$stratum)[H]
n_srs  <- sum(n_iss)

one_stratified <- function(sizes = n_iss) {
  parts <- lapply(H, function(h) {
    i <- sample.int(N_h[[h]], sizes[[h]])                   # SRSWOR within stratum (same draws as sample())
    i <- i[rbinom(length(i), 1, RESP_TRUE[[h]]) == 1]       # respondents
    y <- y_h[[h]][i]; m <- m_h[[h]][i]; N <- N_h[[h]]; r <- length(y)
    R <- mean(y) / mean(m); e <- y - R * m                  # separate ratio estimator, as svyratio()
    c(N = N, r = r, ybar = mean(y), s2 = var(y), T_r = R * M_h[[h]],
      V_r = (M_h[[h]] / (N * mean(m)))^2 * N^2 * (1 - r / N) * var(e) / r) })
  p <- do.call(rbind, parts)
  est <- sum(p[, "N"] * p[, "ybar"])
  se  <- sqrt(sum(p[, "N"]^2 * (1 - p[, "r"] / p[, "N"]) * p[, "s2"] / p[, "r"]))
  c(estimate = est, se = se, n_resp = sum(p[, "r"]), estimate_r = sum(p[, "T_r"]), se_r = sqrt(sum(p[, "V_r"])))
}
one_srswor <- function() {
  y <- sample(pop$colonies_lost, n_srs)
  y <- y[rbinom(n_srs, 1, RESP_PLAN) == 1]
  r <- length(y)
  c(estimate = N * mean(y), se = sqrt(N^2 * (1 - r / N) * var(y) / r), n_resp = r)
}

sim <- do.call(rbind, lapply(seq_len(R_SIM), function(b) {
  s <- one_stratified()                                     # both estimators from the same sample
  rbind(data.frame(rep = b, design = "Stratified (Neyman)", estimate = s[["estimate"]], se = s[["se"]], n_resp = s[["n_resp"]]),
        data.frame(rep = b, design = "Stratified, ratio",   estimate = s[["estimate_r"]], se = s[["se_r"]], n_resp = s[["n_resp"]]),
        data.frame(rep = b, design = "SRSWOR", t(one_srswor()))) }))
sim$lcl <- sim$estimate - Z * sim$se
sim$ucl <- sim$estimate + Z * sim$se
write.csv(sim, file.path(DIR_TAB, "sim_reps.csv"), row.names = FALSE)

# ---- E3-E5: performance measures ---------------------------------------------------
perf <- do.call(rbind, lapply(split(sim, sim$design), function(d) {
  cover <- mean(d$lcl <= T_true & T_true <= d$ucl)
  data.frame(design = d$design[1], reps = nrow(d), mean_estimate = mean(d$estimate),
             bias = mean(d$estimate) - T_true, rel_bias = mean(d$estimate) / T_true - 1,
             emp_se = sd(d$estimate), avg_se = mean(d$se), se_ratio = mean(d$se) / sd(d$estimate),
             rmse = sqrt(mean((d$estimate - T_true)^2)), coverage = cover,
             mc_se_coverage = sqrt(cover * (1 - cover) / nrow(d)),
             avg_width = mean(d$ucl - d$lcl), meets_5pct = mean(Z * d$se <= REL_MOE * T0),
             avg_resp = mean(d$n_resp), mc_se_bias = sd(d$estimate) / sqrt(nrow(d)),
             emp_var = var(d$estimate), miss_low = sum(d$ucl < T_true), miss_high = sum(d$lcl > T_true),
             cor_est_se = cor(d$estimate, d$se))
}))
perf <- perf[match(c("Stratified (Neyman)", "Stratified, ratio", "SRSWOR"), perf$design), ]
perf$emp_deff <- perf$emp_se^2 / perf$emp_se[perf$design == "SRSWOR"]^2
write.csv(perf, file.path(DIR_TAB, "sim_summary.csv"), row.names = FALSE)

# ---- D8 context: a plan built on the true S_h ---------------------------------------
# Same rule as Section B (commercial stratum taken in full with floor(R N_h) respondents,
# Neyman allocation for the rest) but with the true S_h, which the planner cannot know.
# V(T_hat) is then close to V0, yet the estimated SE still varies from sample to sample.
V0     <- (REL_MOE * T0 / Z)^2
S2_h   <- sapply(y_h, var)
rest   <- c("HOB", "SID")
r_com  <- floor(RESP_PLAN * N_h[["COM"]])
V_com  <- N_h[["COM"]]^2 * (1 - r_com / N_h[["COM"]]) * S2_h[["COM"]] / r_com
NS     <- N_h[rest] * sqrt(S2_h[rest])
n_rest <- sum(NS)^2 / (V0 - V_com + sum(N_h[rest] * S2_h[rest]))
n_true <- c(ceiling(n_rest * NS / sum(NS)), COM = r_com)
stopifnot(sum(n_true) == read.csv(file.path(DIR_TAB, "estimates_extra.csv"))$n_true_plan)   # as in 04
n_iss_true <- c(ceiling(n_true[rest] / RESP_PLAN), COM = N_h[["COM"]])
true_reps  <- as.data.frame(t(replicate(R_SIM, one_stratified(n_iss_true))))
meets      <- mean(Z * true_reps$se <= REL_MOE * T0)
write.csv(data.frame(reps = R_SIM, n_resp_plan = sum(n_true), n_issued = sum(n_iss_true),
                     emp_se = sd(true_reps$estimate), target_se = sqrt(V0), meets_5pct = meets,
                     mc_se_meets = sqrt(meets * (1 - meets) / R_SIM)),
          file.path(DIR_TAB, "sim_true_plan.csv"), row.names = FALSE)

# ---- Figure: sampling distributions of T_hat ----------------------------------------
simf <- subset(sim, design %in% c("Stratified, ratio", "SRSWOR"))  # primary estimator against SRSWOR; the expansion estimator is in the table
simf$design <- factor(simf$design, levels = c("Stratified, ratio", "SRSWOR"))
cols_sim <- c("Stratified, ratio" = "#0d3b78", "SRSWOR" = "#b9b8b2")
fig4 <- ggplot(mapping = aes(estimate / 1000, fill = design)) +          # SRSWOR drawn first, behind
  geom_histogram(data = subset(simf, design == "SRSWOR"), binwidth = 0.5, boundary = 0,
                 colour = "white", linewidth = 0.2) +
  geom_histogram(data = subset(simf, design == "Stratified, ratio"), binwidth = 0.5, boundary = 0,
                 colour = "white", linewidth = 0.2, alpha = 0.92) +
  geom_vline(xintercept = T_true / 1000, colour = "#eb6834", linewidth = 0.8) +
  annotate("text", x = T_true / 1000, y = Inf, label = "True T", vjust = 1.3, hjust = -0.1, size = 3.2, colour = "#eb6834") +
  scale_fill_manual(values = cols_sim, name = NULL,
                    labels = c("Stratified, ratio" = "Stratified, ratio estimator", "SRSWOR" = "SRSWOR")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Estimated total colonies lost (thousands)", y = "Repetitions") +
  theme_report(base_size = 11) +
  theme(legend.position = "top", legend.justification = "left", panel.grid.major.x = element_blank())
save_fig(fig4, "fig4_sim_distribution", width = 6.5, height = 2.6)

# ---- D8 context: issued size for which the bound is met in 80% or 90% of samples -----
# Hobby and side-line issued sizes scaled up together (capped at N_h), all commercial
# beekeepers issued as before; 2,000 repetitions at each size. Run last, so the
# results above keep their random numbers.
assure <- do.call(rbind, lapply(seq(1, 3, by = 0.1), function(k) {
  sz <- pmin(round(n_iss * c(HOB = k, SID = k, COM = 1)), N_h)
  r  <- as.data.frame(t(replicate(2000, one_stratified(sz))))
  data.frame(k = k, issued = sum(sz), meets_expansion = mean(Z * r$se <= REL_MOE * T0),
             meets_ratio = mean(Z * r$se_r <= REL_MOE * T0)) }))
write.csv(assure, file.path(DIR_TAB, "sim_assurance.csv"), row.names = FALSE)

# ---- Sensitivity to non-response that depends on the losses (MNAR) -----------------------
# logit P(respond) = a + b * y_i / m_i: with b = -/+ log(2) / 0.5, a beekeeper who lost half
# their colonies has half / twice the odds of responding of one who lost none [G20 report both
# tendencies]. a is set within each stratum (overall for SRSWOR) so that the mean response
# probability stays at the planned rate; the weighting still assumes MAR. 2,000 repetitions each.
lr_h  <- lapply(H, function(h) y_h[[h]] / m_h[[h]]); names(lr_h) <- H
lr    <- pop$colonies_lost / pop$colonies_oct
calib <- function(b, lr_x, target) { a <- uniroot(function(a) mean(plogis(a + b * lr_x)) - target, c(-10, 10))$root
  plogis(a + b * lr_x) }
mnar <- do.call(rbind, lapply(c(-1, 1) * log(2) / 0.5, function(b) {
  pr_h <- lapply(H, function(h) calib(b, lr_h[[h]], RESP_TRUE[[h]])); names(pr_h) <- H
  pr   <- calib(b, lr, RESP_PLAN)
  r <- t(replicate(2000, {
    s <- colSums(do.call(rbind, lapply(H, function(h) {
      i <- sample.int(N_h[[h]], n_iss[[h]]); i <- i[rbinom(length(i), 1, pr_h[[h]][i]) == 1]
      y <- y_h[[h]][i]; m <- m_h[[h]][i]; N <- N_h[[h]]; r <- length(y); R <- mean(y) / mean(m); e <- y - R * m
      c(T = N * mean(y), V = N^2 * (1 - r / N) * var(y) / r,
        T_r = R * M_h[[h]], V_r = (M_h[[h]] / (N * mean(m)))^2 * N^2 * (1 - r / N) * var(e) / r) })))
    i <- sample.int(N, n_srs); i <- i[rbinom(n_srs, 1, pr[i]) == 1]; y <- pop$colonies_lost[i]; r <- length(y)
    c(s, T_s = N * mean(y), V_s = N^2 * (1 - r / N) * var(y) / r) }))
  cover <- function(t, v) mean(abs(r[, t] - T_true) <= Z * sqrt(r[, v]))
  data.frame(b = b, odds_ratio_half_lost = exp(b * 0.5),
             rel_bias_ratio = mean(r[, "T_r"]) / T_true - 1, cover_ratio = cover("T_r", "V_r"),
             rel_bias_exp = mean(r[, "T"]) / T_true - 1,     cover_exp = cover("T", "V"),
             rel_bias_srs = mean(r[, "T_s"]) / T_true - 1,   cover_srs = cover("T_s", "V_s")) }))
write.csv(mnar, file.path(DIR_TAB, "sim_mnar.csv"), row.names = FALSE)

cat(sprintf("R = %d repetitions | true T = %d\n", R_SIM, T_true))
print(perf, digits = 4)
