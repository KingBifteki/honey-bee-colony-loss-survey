# =============================================================================
# 02_sample_size.R - Section B: sample size, allocation and sensitivity
#
# The planner knows the register (N_h and each beekeeper's colonies m_i) and
# the COLOSS loss rates by size class [G23], but not the colonies lost y_i.
# Prior quantities come from the beta-binomial model of A5:
#   E(y_i) = p_h m_i,   Var(y_i) = m_i p_h (1 - p_h) + m_i (m_i - 1) sigma^2
# with sigma = LOSS_SD, the assumed between-beekeeper SD of the loss risk.
# =============================================================================

source("R/00_setup.R")

pop <- read.csv(PATH_POP)
H   <- strata_def$stratum
Nh  <- setNames(strata_def$N_h, H)

# ---- Planner's prior: baseline total and stratum SDs --------------------------
p_i  <- PRIOR_LOSS_RATE[pop$stratum]
m_i  <- pop$colonies_oct
Ey   <- p_i * m_i
Vy   <- m_i * p_i * (1 - p_i) + m_i * (m_i - 1) * LOSS_SD^2
prior <- data.frame(
  stratum = H, N_h = as.numeric(Nh),
  M_h     = sapply(H, function(h) sum(m_i[pop$stratum == h])),
  p_h     = as.numeric(PRIOR_LOSS_RATE[H]),
  T0_h    = sapply(H, function(h) sum(Ey[pop$stratum == h])),
  S_h     = sapply(H, function(h) { k <- pop$stratum == h; sqrt(var(Ey[k]) + mean(Vy[k])) }))
T0      <- sum(prior$T0_h)                    # baseline value of T
S_all   <- sqrt(var(Ey) + mean(Vy))           # prior overall SD (for SRSWOR)
S_prior <- setNames(prior$S_h, H)

# ---- Neyman allocation with take-all strata, a minimum and non-response ------
# n_h are respondents; issued n_h* = n_h / R. A take-all stratum issues all N_h
# beekeepers, so with response rate R it keeps the variance
# N_h^2 (1 - r_h/N_h) S_h^2 / r_h with r_h = floor(R N_h) whole respondents,
# which the other strata must make up for.
neyman_plan <- function(rel = REL_MOE, theta = 1, R = 1, S = S_prior) {
  V0 <- (rel * theta * T0 / Z)^2
  take_all <- character(0)
  repeat {
    rest   <- setdiff(H, take_all)
    r_ta   <- floor(R * Nh[take_all])
    V_ta   <- sum(Nh[take_all]^2 * (1 - r_ta / Nh[take_all]) * S[take_all]^2 / r_ta)
    budget <- V0 - V_ta
    if (budget <= 0 || !length(rest)) return(NULL)            # target unattainable
    NS <- Nh[rest] * S[rest]
    n_rest <- sum(NS)^2 / (budget + sum(Nh[rest] * S[rest]^2))
    n_h <- n_rest * NS / sum(NS)
    over <- n_h / (R * Nh[rest])
    if (all(over <= 1)) break
    take_all <- c(take_all, rest[which.max(over)])
  }
  n_h <- pmax(ceiling(n_h), MIN_NH)                           # whole respondents
  n   <- setNames(numeric(length(H)), H)
  n[rest] <- n_h; n[take_all] <- floor(R * Nh[take_all])       # whole respondents
  issued <- setNames(ceiling(n / R), H); issued[take_all] <- Nh[take_all]
  list(n = n, issued = issued, take_all = take_all, V0 = V0, V_ta = V_ta,
       n_total = sum(n), issued_total = sum(issued))
}

srs_plan <- function(rel = REL_MOE, theta = 1, R = 1) {
  V0 <- (rel * theta * T0 / Z)^2
  n0 <- N^2 * S_all^2 / V0
  n  <- ceiling(n0 / (1 + n0 / N))
  list(n0 = n0, n = n, issued = ceiling(n / R), feasible = n / R <= N)
}

# ---- B3-B5: plans at the 5% target --------------------------------------------
st_full <- neyman_plan(R = 1)
st_resp <- neyman_plan(R = RESP_PLAN)
srs     <- srs_plan(R = 1)
srs_r   <- srs_plan(R = RESP_PLAN)
V0      <- st_full$V0

# first pass of Neyman allocation over all strata (full response), before the
# take-all rule: shows how many commercial respondents it would ask for
n_first     <- sum(Nh * S_prior)^2 / (V0 + sum(Nh * S_prior^2))
n_com_first <- n_first * Nh[["COM"]] * S_prior[["COM"]] / sum(Nh * S_prior)

# fpc: Neyman without fpc or take-all (n0) versus the full calculation
NS_all  <- Nh * S_prior
n0_st   <- sum(NS_all)^2 / V0

# proportional allocation for comparison (n_h proportional to N_h, full response)
prop_plan <- function(rel = REL_MOE) {
  V0 <- (rel * T0 / Z)^2
  n  <- N * sum(Nh * S_prior^2) / (V0 + sum(Nh * S_prior^2))
  setNames(pmax(ceiling(n * Nh / N), MIN_NH), H)
}
n_prop <- prop_plan()

# finer strata on colony counts (A6): under the planner's model, the best single split of
# the hobby class on m_i, with the same rules (Neyman, take-all, minimum, full response)
plan_n_groups <- function(g) {
  gs <- unique(g); Ng <- sapply(gs, function(k) sum(g == k))
  Sg <- sapply(gs, function(k) sqrt(var(Ey[g == k]) + mean(Vy[g == k]))); ta <- integer(0)
  repeat {
    rest <- setdiff(seq_along(gs), ta); NS <- Ng[rest] * Sg[rest]
    nh <- sum(NS)^2 / (V0 + sum(Ng[rest] * Sg[rest]^2)) * NS / sum(NS)
    if (all(nh <= Ng[rest])) break
    ta <- c(ta, rest[which.max(nh / Ng[rest])])
  }
  sum(pmax(ceiling(nh), MIN_NH)) + sum(Ng[ta])
}
split_n <- sapply(1:49, function(k) { g <- pop$stratum; g[g == "HOB" & m_i <= k] <- "HOB_small"; plan_n_groups(g) })
split_hob   <- which.min(split_n)                      # hobby operations with up to split_hob colonies
n_split_hob <- min(split_n)
stopifnot(plan_n_groups(pop$stratum) == st_full$n_total)   # same rules as neyman_plan()

# the ratio estimator planned under the same model (D8): within a stratum E(y_i) = p_h m_i is
# proportional to m_i, so the prior residual SD about R_h m_i is sqrt(mean Var(y_i))
Se_prior     <- sapply(H, function(h) sqrt(mean(Vy[pop$stratum == h])))
n_ratio_full <- neyman_plan(R = 1, S = Se_prior)$n_total     # B3 minimum for the primary estimator
n_ratio_resp <- neyman_plan(R = RESP_PLAN, S = Se_prior)$n_total
# safety margin (C1): the largest factor by which the prior residual SDs could be too small
# with the planned respondents still meeting the bound for T_R
infl_grid <- seq(1, 1.5, by = 0.01)
n_ratio_k <- sapply(infl_grid, function(k) {
  pl <- neyman_plan(R = RESP_PLAN, S = Se_prior * k); if (is.null(pl)) Inf else pl$n_total })
ratio_margin_sd <- max(infl_grid[n_ratio_k <= st_resp$n_total]) - 1

# design effect at the planned number of respondents (full response)
V_st  <- function(n) sum(Nh^2 * (1 - n / Nh) * S_prior^2 / n)
V_srs <- function(n) N^2 * (1 - n / N) * S_all^2 / n
deff_plan <- V_st(st_full$n) / V_srs(sum(st_full$n))

# minimum response rate in the take-all stratum for the target to be reachable:
# its variance N^2 (1 - r/N) S^2 / r stays within V0 when r >= N^2 S^2 / (V0 + N S^2),
# counted in whole respondents as in neyman_plan()
r_min <- ceiling(Nh[["COM"]]^2 * S_prior[["COM"]]^2 / (V0 + Nh[["COM"]] * S_prior[["COM"]]^2))
R_min <- r_min / Nh[["COM"]]
# the same threshold for the ratio estimator, with the residual SD about the colony counts
R_min_ratio <- ceiling(Nh[["COM"]]^2 * Se_prior[["COM"]]^2 / (V0 + Nh[["COM"]] * Se_prior[["COM"]]^2)) / Nh[["COM"]]

# smallest margin reachable at the planned response: with every beekeeper issued,
# only the non-response variance sum N_h^2 (1 - r_h/N_h) S_h^2 / r_h is left
r_all        <- floor(RESP_PLAN * Nh)
rel_min_resp <- Z * sqrt(sum(Nh^2 * (1 - r_all / Nh) * S_prior^2 / r_all)) / T0

# B10: prior SDs that are too small (all S_h raised together, planned response)
prior_infl <- do.call(rbind, lapply(c(1.2, 1.3), function(k) {
  pl <- neyman_plan(R = RESP_PLAN, S = S_prior * k); s <- S_prior[["COM"]] * k
  data.frame(inflation = k, n_resp = pl$n_total, issued = pl$issued_total,
             R_min_com = ceiling(Nh[["COM"]]^2 * s^2 / (V0 + Nh[["COM"]] * s^2)) / Nh[["COM"]]) }))

allocation <- data.frame(
  stratum        = H,
  N_h            = as.numeric(Nh),
  S_h_prior      = as.numeric(S_prior),
  Se_h_prior     = as.numeric(Se_prior[H]),                    # residual SD about p_h m_i (ratio estimator)
  NS_share       = as.numeric(NS_all / sum(NS_all)),
  n_full         = as.numeric(st_full$n),
  n_resp         = as.numeric(st_resp$n),
  issued         = as.numeric(st_resp$issued),
  frac_issued    = as.numeric(st_resp$issued / Nh),
  take_all       = H %in% st_resp$take_all)

plan_summary <- data.frame(
  T0 = T0, rel_moe = REL_MOE, halfwidth = REL_MOE * T0, V0 = V0, z = Z,
  S_all_prior = S_all, resp_plan = RESP_PLAN,
  n0_strat = n0_st, n_strat_full = st_full$n_total,
  n_strat_resp = st_resp$n_total, issued_strat = st_resp$issued_total,
  V_takeall_resp = st_resp$V_ta, share_V0_takeall = st_resp$V_ta / V0,
  n0_srs = srs$n0, n_srs_full = srs$n, issued_srs = srs_r$issued, n_prop_full = sum(n_prop), n_com_first = n_com_first,
  split_hob = split_hob, n_split_hob = n_split_hob, n_ratio_full = n_ratio_full, n_ratio_resp = n_ratio_resp, ratio_margin_sd = ratio_margin_sd,
  deff_plan = deff_plan, n_ratio = srs$n / st_full$n_total, R_min_com = R_min, R_min_com_ratio = R_min_ratio,
  rel_min_resp = rel_min_resp)

write.csv(prior,        file.path(DIR_TAB, "plan_prior.csv"),      row.names = FALSE)
write.csv(prior_infl,   file.path(DIR_TAB, "plan_prior_inflation.csv"), row.names = FALSE)
write.csv(allocation,   file.path(DIR_TAB, "plan_allocation.csv"), row.names = FALSE)
write.csv(plan_summary, file.path(DIR_TAB, "plan_summary.csv"),    row.names = FALSE)

# ---- B6: required respondents against the relative margin of error ----------
curve_st <- function(rel, theta = 1, R = 1) sapply(rel, function(e) {
  pl <- neyman_plan(rel = e, theta = theta, R = R); if (is.null(pl)) NA else pl$n_total })
curve_srs <- function(rel) sapply(rel, function(e) srs_plan(rel = e)$n)

LAB_RESP <- sprintf("Stratified, %d%% response", round(100 * RESP_PLAN))   # e.g. "Stratified, 77% response"
b6 <- rbind(
  data.frame(rel = REL_GRID, n = curve_st(REL_GRID),                  design = "Stratified, full response"),
  data.frame(rel = REL_GRID, n = curve_st(REL_GRID, R = RESP_PLAN),   design = LAB_RESP),
  data.frame(rel = REL_GRID, n = curve_srs(REL_GRID),                 design = "SRSWOR, full response"))
b6$design <- factor(b6$design, levels = c("SRSWOR, full response", LAB_RESP, "Stratified, full response"))
rel_min_resp <- plan_summary$rel_min_resp

# ---- B7: sensitivity to the assumed value of T --------------------------------
LAB_R77 <- sprintf("%d%% response", round(100 * RESP_PLAN))
b7 <- do.call(rbind, lapply(THETA, function(th) rbind(
  data.frame(rel = REL_GRID, n = curve_st(REL_GRID, theta = th),                theta = th, response = "Full response"),
  data.frame(rel = REL_GRID, n = curve_st(REL_GRID, theta = th, R = RESP_PLAN), theta = th, response = LAB_R77))))
b7$label <- factor(paste0(100 * b7$theta, "% of baseline"), levels = paste0(100 * rev(THETA), "% of baseline"))
sens_at_target <- data.frame(theta = THETA, n = sapply(THETA, function(th) neyman_plan(theta = th)$n_total),
                             take_all = sapply(THETA, function(th) paste(neyman_plan(theta = th)$take_all, collapse = " ")),
                             rel_min_resp = rel_min_resp / THETA)   # smallest margin reachable at the planned response
write.csv(b6, file.path(DIR_TAB, "plan_curve_margin.csv"), row.names = FALSE)
write.csv(b7, file.path(DIR_TAB, "plan_curve_sensitivity.csv"), row.names = FALSE)
write.csv(sens_at_target, file.path(DIR_TAB, "plan_sensitivity_at_5pct.csv"), row.names = FALSE)

# ---- Figures (sized to sit side by side in the report) ---------------------------
pct_axis <- scale_x_continuous(labels = function(x) paste0(100 * x, "%"), breaks = seq(0.02, 0.14, 0.04),
                               expand = expansion(mult = c(0.02, 0.02)))
n_axis   <- scale_y_continuous(labels = scales::comma, limits = c(0, N), breaks = seq(0, N, 2500),
                               expand = expansion(mult = c(0, 0.02)))
ref_lines <- list(
  geom_hline(yintercept = N, colour = INK_2, linewidth = 0.3),
  annotate("text", x = 0.15, y = N, label = "N = 10,000", hjust = 1, vjust = 1.6, size = 3, colour = INK_2),
  geom_vline(xintercept = REL_MOE, linetype = "dashed", colour = INK_2, linewidth = 0.4))
legend_top <- theme(legend.position = "top", legend.justification = "left", legend.title = element_text(size = 9),
                    legend.text = element_text(size = 9), legend.key.width = unit(1.4, "lines"),
                    legend.margin = margin(0, 0, -4, 0))

cols_b6 <- setNames(c("#0d3b78", "#2a78d6", INK_2), c("Stratified, full response", LAB_RESP, "SRSWOR, full response"))
lty_b6  <- setNames(c("solid", "dashed", "solid"),   c("Stratified, full response", LAB_RESP, "SRSWOR, full response"))
b6$design <- factor(b6$design, levels = names(cols_b6))
pts_b6 <- data.frame(rel = REL_MOE, n = c(st_full$n_total, st_resp$n_total, srs$n),
                     design = factor(names(cols_b6), levels = names(cols_b6)),
                     hj = c(1.2, 1.25, -0.25), vj = c(1.5, -0.4, -0.6))
fig2 <- ggplot(b6, aes(rel, n, colour = design, linetype = design)) +
  ref_lines +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(data = pts_b6, size = 2, show.legend = FALSE) +
  geom_text(data = pts_b6, aes(label = scales::comma(round(n)), hjust = hj, vjust = vj), size = 3.2,
            colour = INK, show.legend = FALSE) +
  annotate("curve", x = 0.0575, y = 3700, xend = 0.0405, yend = 4700, curvature = 0.25, colour = "#2a78d6",
           linewidth = 0.4, arrow = arrow(length = unit(0.08, "inches"))) +
  annotate("text", x = 0.0585, y = 3700, hjust = 0, size = 2.9, colour = "#2a78d6", lineheight = 0.9,
           label = sprintf("%d%% response:\nunattainable\nbelow %.1f%%", round(100 * RESP_PLAN), 100 * rel_min_resp)) +
  scale_colour_manual(values = cols_b6, name = NULL) +
  scale_linetype_manual(values = lty_b6, name = NULL) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  pct_axis + n_axis +
  labs(x = "Relative margin of error (% of baseline)", y = "Respondents needed") +
  theme_report(base_size = 11) + legend_top
save_fig(fig2, "fig2_n_vs_margin", width = 4.3, height = 3.2)

blues <- c("100%" = "#0d3b78", "75%" = "#2a78d6", "50%" = "#5f9fe6", "25%" = "#93c0f0")
b7$share <- factor(paste0(100 * b7$theta, "%"), levels = names(blues))
pts_b7 <- data.frame(rel = REL_MOE, n = sens_at_target$n,
                     share = factor(paste0(100 * THETA, "%"), levels = names(blues)),
                     response = factor("Full response", levels = c("Full response", LAB_R77)))
b7$response <- factor(b7$response, levels = c("Full response", LAB_R77))
fig3 <- ggplot(b7, aes(rel, n, colour = share, linetype = response)) +
  ref_lines +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_point(data = pts_b7, size = 2) +
  geom_label(data = pts_b7, aes(label = scales::comma(round(n)), hjust = ifelse(share == "100%", 1.2, -0.25),
                                vjust = ifelse(share == "100%", 1.5, 0.5)), size = 3.2, colour = INK,
             fill = "white", border.colour = NA, label.padding = unit(0.08, "lines"),
             show.legend = FALSE) +                       # boxes keep labels clear of the dashed curves
  scale_colour_manual(values = blues, name = "Assumed T, as % of baseline") +
  scale_linetype_manual(values = setNames(c("solid", "dashed"), c("Full response", LAB_R77)), name = NULL) +
  guides(colour = guide_legend(title.position = "top", nrow = 1, order = 1),
         linetype = guide_legend(nrow = 1, order = 2, override.aes = list(colour = INK_2),
                                 theme = theme(legend.key.width = unit(2.4, "lines")))) +
  pct_axis + n_axis +
  labs(x = "Relative margin of error (% of assumed T)", y = "Respondents needed") +
  theme_report(base_size = 11) + legend_top +
  theme(legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(0, "pt"))
save_fig(fig3, "fig3_sensitivity", width = 4.3, height = 3.2)

# ---- Console summary ---------------------------------------------------------------
cat(sprintf("Baseline T0 = %.0f, half-width target = %.0f colonies, V0 = %.0f\n", T0, REL_MOE * T0, V0))
cat("Prior S_h:", round(S_prior, 2), "| prior overall S:", round(S_all, 2), "\n")
cat(sprintf("Stratified: n0 (no fpc) = %.0f | full response n = %.0f | %d%% response: respondents %.0f, issued %d\n",
            n0_st, st_full$n_total, round(100 * RESP_PLAN), st_resp$n_total, st_resp$issued_total))
cat(sprintf("SRSWOR: n0 = %.0f, n = %.0f, issued at %d%% = %.0f | deff = %.3f | ratio = %.1f\n",
            srs$n0, srs$n, round(100 * RESP_PLAN), srs_r$issued, deff_plan, srs$n / st_full$n_total))
cat(sprintf("Take-all variance at %d%% response = %.0f (%.0f%% of V0); minimum COM response = %.3f; min margin at %d%% = %.2f%%\n",
            round(100 * RESP_PLAN), st_resp$V_ta, 100 * st_resp$V_ta / V0, R_min, round(100 * RESP_PLAN), 100 * rel_min_resp))
print(allocation, digits = 4)
print(sens_at_target)
