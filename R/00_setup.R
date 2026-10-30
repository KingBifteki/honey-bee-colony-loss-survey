# =============================================================================
# 00_setup.R - shared settings, sources of every input, and helpers
# National Colony Loss Survey 2026 (STAT4102 Mini Project)
#
# Every other script sources this file. Run all scripts from the project root
# (open BeeSurvey.Rproj in RStudio, or run `Rscript run_all.R`).
#
# Sources used for the simulation inputs:
#   [G23] Gray et al. (2023). Honey bee colony loss rates in 37 countries using
#         the COLOSS survey for winter 2019-2020: the combined effects of
#         operation size, migration and queen replacement. Journal of
#         Apicultural Research 62(2): 204-210.
#   [Z13] van der Zee et al. (2013). Standard survey methods for estimating
#         colony losses and explanatory risk factors in Apis mellifera.
#         Journal of Apicultural Research 52(4): 1-36.
# Anything not taken from a source is marked ASSUMPTION.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(survey)
  library(sampling)
})

# ---- Random seeds (one per stage so each stage is reproducible on its own) --
RNGkind("Mersenne-Twister", "Inversion", "Rejection")   # R's defaults since 3.6, fixed so any setup gives the same draws
SEED_POP    <- 4102   # Section C2: generating the finite population
SEED_SAMPLE <- 2026   # Section C3: drawing the realised sample (+ fieldwork)
SEED_BENCH  <- 2027   # Section D6: drawing the SRSWOR benchmark sample
SEED_SIM    <- 1021   # Section E:  Monte Carlo simulation

# ---- Survey constants --------------------------------------------------------
N         <- 10000    # finite population size (fixed by the project brief)
CONF      <- 0.95
Z         <- qnorm(1 - (1 - CONF) / 2)
REL_MOE   <- 0.05     # CI half-width must be <= 5% of the baseline value
RESP_PLAN <- 0.77     # planned response rate: highest rate in the surveys reviewed by [Z13]
                      # (Scotland 2006, with personal contact; Box 5). Applying it here is an ASSUMPTION
MIN_NH    <- 30       # minimum respondents per stratum

# ---- COLOSS winter 2019/20 results by operation size [G23] -------------------
# 30,491 beekeepers: 90.2%, 6.6% and 3.2% had up to 50, 51-150 and over 150
# colonies, managing 296,763, 177,501 and 362,817 colonies. Loss rates from
# dead/empty colonies and unsolvable queen problems (natural disasters
# excluded): 16.6%, 14.6% and 15.2%. 81.7% of beekeepers did not migrate.
COLOSS_N_BEEKEEPERS <- 30491
coloss <- data.frame(
  stratum        = c("HOB", "SID", "COM"),
  operation      = c("Hobby (up to 50 colonies)", "Side-line (51-150 colonies)",
                     "Commercial (more than 150 colonies)"),
  pct_beekeepers = c(90.2, 6.6, 3.2),
  colonies       = c(296763, 177501, 362817),
  loss_rate      = c(0.166, 0.146, 0.152),
  stringsAsFactors = FALSE
)
coloss$mean_colonies <- coloss$colonies / (coloss$pct_beekeepers / 100 * COLOSS_N_BEEKEEPERS)
coloss$pct_colonies  <- 100 * coloss$colonies / sum(coloss$colonies)
P_MIGRATORY <- 1 - 0.817            # 81.7% did not migrate [G23]; the rest includes 163 "Don't know"

# ---- Strata for the simulated population (size classes as in [G23]) -------
strata_def <- data.frame(
  stratum   = coloss$stratum,
  operation = coloss$operation,
  N_h       = round(N * coloss$pct_beekeepers / 100),
  lower     = c(1, 51, 151),
  upper     = c(50, 150, 5000),     # upper cap for commercial: ASSUMPTION
  stringsAsFactors = FALSE
)
strata_def$N_h[1] <- N - sum(strata_def$N_h[-1])
stopifnot(sum(strata_def$N_h) == N)
strata_def$W_h <- strata_def$N_h / N

# ---- Loss model --------------------------------------------------------------
# Each beekeeper's loss risk has mean equal to the [G23] loss rate of their
# size class. Colony loss data vary more than a binomial model allows, so a
# beta-binomial model is used ([G23], [Z13]); the between-beekeeper SD of the
# risk is an ASSUMPTION.
LOSS_SD <- 0.15

# ---- Prior information for planning [G23] -------------------------------------
# The planner uses the winter 2019/20 loss rates by size class together with
# the colony counts on the register.
PRIOR_LOSS_RATE <- setNames(coloss$loss_rate, coloss$stratum)

# Response propensities used to SIMULATE fieldwork: the planned rate in every
# stratum, with response random within strata (ASSUMPTION).
RESP_TRUE <- c(HOB = RESP_PLAN, SID = RESP_PLAN, COM = RESP_PLAN)

# ---- Sample-size plots (Section B) ---------------------------------------------
REL_GRID <- seq(0.02, 0.15, by = 0.0025)   # relative margins of error for B6/B7
THETA    <- c(0.25, 0.50, 0.75, 1.00)      # assumed T as a share of the baseline (B7)

# ---- File paths --------------------------------------------------------------
PATH_POP    <- "data/population.csv"
PATH_SAMPLE <- "data/sample.csv"
DIR_FIG     <- "output/figures"
DIR_TAB     <- "output/tables"
for (d in c(dirname(PATH_POP), DIR_FIG, DIR_TAB)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Plot styling ------------------------------------------------------------
STRATUM_COLS <- c(HOB = "#2a78d6", SID = "#eb6834", COM = "#1baf7a")
INK   <- "#0b0b0b"
INK_2 <- "#52514e"
GRID  <- "#e1e0d9"

theme_report <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      text             = element_text(colour = INK),
      axis.text        = element_text(colour = INK_2),
      panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold", hjust = 0, colour = INK),
      plot.title.position = "plot",
      legend.position  = "top"
    )
}

save_fig <- function(plot, name, width = 6.5, height = 4) {
  ggsave(file.path(DIR_FIG, paste0(name, ".png")), plot, device = ragg::agg_png,
         width = width, height = height, dpi = 300, bg = "white")
}
