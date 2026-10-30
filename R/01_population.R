# =============================================================================
# 01_population.R - Section A5 / C2: generate the finite population (N = 10,000)
#
# Population: 10,000 registered beekeepers with at least one colony going into
# winter on 1 October 2025. Outcome: colonies lost by 31 March 2026 as dead or
# empty colonies or through unsolvable queen problems (natural disasters
# excluded), the loss measure Gray et al. (2023) use for operation size.
#
# Model, for beekeeper i in size class h:
#   colonies m_i  : discrete distribution on the class range whose mean equals
#                   the COLOSS class mean [G23] (shape: ASSUMPTION)
#   migratory     : Bernoulli(0.183); 81.7% of COLOSS beekeepers did not migrate [G23]
#   loss risk p_i ~ Beta with mean = COLOSS loss rate of the class [G23], SD 0.15
#   colonies lost y_i ~ Binomial(m_i, p_i)
# =============================================================================

source("R/00_setup.R")
set.seed(SEED_POP)

# Geometric-shaped probabilities on lo..hi with a given mean (solved exactly)
shaped_probs <- function(lo, hi, target) {
  k <- lo:hi
  f <- function(q) sum(k * q^(k - lo)) / sum(q^(k - lo)) - target
  q <- uniroot(f, c(1e-9, 1 - 1e-9))$root
  p <- q^(k - lo)
  list(k = k, p = p / sum(p))
}

draw_colonies <- function(s, target, n) {
  if (s$stratum == "COM") {
    # long right tail: 151 + lognormal with mean (target - 151)
    sdlog <- 1.0
    pmin(s$lower + round(rlnorm(n, log(target - s$lower) - sdlog^2 / 2, sdlog)), s$upper)
  } else {
    d <- shaped_probs(s$lower, s$upper, target)
    sample(d$k, n, replace = TRUE, prob = d$p)
  }
}

beta_draw <- function(n, mu, sd) {
  k <- mu * (1 - mu) / sd^2 - 1
  rbeta(n, mu * k, (1 - mu) * k)
}

simulate_stratum <- function(h) {
  s     <- strata_def[strata_def$stratum == h, ]
  cl    <- coloss[coloss$stratum == h, ]
  n     <- s$N_h
  colonies  <- draw_colonies(s, cl$mean_colonies, n)
  migratory <- rbinom(n, 1, P_MIGRATORY)
  risk      <- beta_draw(n, cl$loss_rate, LOSS_SD)
  lost      <- rbinom(n, colonies, risk)
  data.frame(stratum = h, operation = s$operation, migratory,
             colonies_oct = colonies, colonies_lost = lost,
             loss_pct = round(100 * lost / colonies, 1))
}

population <- do.call(rbind, lapply(strata_def$stratum, simulate_stratum))

# Register order: shuffle, then number beekeepers BK00001..BK10000
population <- population[sample(nrow(population)), ]
population <- cbind(bk_id = sprintf("BK%05d", seq_len(N)), population)
rownames(population) <- NULL
stopifnot(nrow(population) == N, !anyNA(population),
          all(population$colonies_lost <= population$colonies_oct))
write.csv(population, PATH_POP, row.names = FALSE)

# ---- True finite-population values ------------------------------------------
M_TOTAL <- sum(population$colonies_oct)       # colonies on the register (known)
T_TRUE  <- sum(population$colonies_lost)      # the parameter: total colonies lost

strata_pop <- population |>
  group_by(stratum) |>
  summarise(N_h           = n(),
            mean_colonies = mean(colonies_oct),
            M_h           = sum(colonies_oct),
            pct_migratory = 100 * mean(migratory),
            Ybar_h        = mean(colonies_lost),
            S_h           = sd(colonies_lost),
            S2_h          = var(colonies_lost),
            T_h           = sum(colonies_lost),
            avg_loss_pct  = mean(loss_pct),
            .groups = "drop") |>
  mutate(total_loss_pct   = 100 * T_h / M_h,
         share_beekeepers = 100 * N_h / N,
         share_colonies   = 100 * M_h / M_TOTAL,
         share_lost       = 100 * T_h / T_TRUE) |>
  left_join(strata_def[, c("stratum", "operation", "W_h")], by = "stratum") |>
  left_join(transmute(coloss, stratum, coloss_pct_beekeepers = pct_beekeepers,
                      coloss_mean_colonies = mean_colonies,
                      coloss_pct_colonies = pct_colonies,
                      coloss_loss_pct = 100 * loss_rate), by = "stratum")
strata_pop <- strata_pop[match(strata_def$stratum, strata_pop$stratum), ]

ss_within <- sum((strata_pop$N_h - 1) * strata_pop$S2_h)
ss_total  <- (N - 1) * var(population$colonies_lost)
overall_pop <- data.frame(
  N = N, M = M_TOTAL, T = T_TRUE, Ybar = T_TRUE / N,
  S = sd(population$colonies_lost), S2 = var(population$colonies_lost),
  CV = sd(population$colonies_lost) / (T_TRUE / N),
  total_loss_pct = 100 * T_TRUE / M_TOTAL,
  avg_loss_pct   = mean(population$loss_pct),
  pct_migratory  = 100 * mean(population$migratory),
  eta2 = 1 - ss_within / ss_total)

write.csv(strata_pop,  file.path(DIR_TAB, "strata_population.csv"), row.names = FALSE)
write.csv(overall_pop, file.path(DIR_TAB, "overall_population.csv"), row.names = FALSE)

# ---- Figure 1: who holds the colonies, and who loses them -------------------
fig_df <- rbind(
  data.frame(measure = "Beekeepers",            stratum = strata_pop$stratum, pct = strata_pop$share_beekeepers),
  data.frame(measure = "Colonies going into winter", stratum = strata_pop$stratum, pct = strata_pop$share_colonies),
  data.frame(measure = "Colonies lost",         stratum = strata_pop$stratum, pct = strata_pop$share_lost))
fig_df$measure <- factor(fig_df$measure, levels = c("Beekeepers", "Colonies going into winter", "Colonies lost"))
fig_df$stratum <- factor(fig_df$stratum, levels = rev(strata_def$stratum))
lab_map <- setNames(gsub("(\\d)-(\\d)", "\\1\u2013\\2", strata_def$operation), strata_def$stratum)

# label positions: hobby stacks at the bottom, commercial at the top
fig_df <- fig_df[order(fig_df$measure, -as.integer(fig_df$stratum)), ]
fig_df$ymid  <- ave(fig_df$pct, fig_df$measure, FUN = cumsum) - fig_df$pct / 2
fig_df$label <- sprintf("%.0f%%", fig_df$pct)
fig_df$small <- fig_df$pct < 5           # too thin for a label inside the bar
BAR_W <- 0.62

fig1 <- ggplot(fig_df, aes(measure, pct, fill = stratum)) +
  geom_col(width = BAR_W, colour = "white", linewidth = 0.6) +
  geom_text(data = subset(fig_df, !small), aes(y = ymid, label = label),
            size = 3.2, colour = INK, fontface = "bold") +
  geom_text(data = subset(fig_df, small), aes(x = as.integer(measure) + BAR_W / 2 + 0.04, y = ymid, label = label),
            size = 3.2, colour = INK, fontface = "bold", hjust = 0) +
  scale_fill_manual(values = STRATUM_COLS, labels = lab_map, name = NULL,
                    breaks = rev(strata_def$stratum)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Share of the population total") +
  theme_report() +
  theme(panel.grid.major.x = element_blank(), legend.position = "right")
save_fig(fig1, "fig1_population_shares", width = 6.5, height = 2.6)

cat(sprintf("Population written to %s\n", PATH_POP))
cat(sprintf("Colonies M = %s | colonies lost T = %s | total loss %.1f%% | average loss per beekeeper %.1f%% | migratory %.1f%% | eta^2 = %.3f\n",
            format(M_TOTAL, big.mark = ","), format(T_TRUE, big.mark = ","),
            overall_pop$total_loss_pct, overall_pop$avg_loss_pct, overall_pop$pct_migratory,
            overall_pop$eta2))
print(as.data.frame(strata_pop[, c("stratum", "N_h", "coloss_pct_beekeepers", "mean_colonies",
                                   "coloss_mean_colonies", "share_colonies", "coloss_pct_colonies",
                                   "total_loss_pct", "coloss_loss_pct", "S_h", "T_h", "share_lost")]),
      digits = 3)
