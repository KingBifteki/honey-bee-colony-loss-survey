# =============================================================================
# 03_sample.R - Section C: draw one stratified sample and simulate fieldwork
#
# 1. Sort the register by stratum and draw an independent SRSWOR of the issued
#    size n_h* in each stratum with sampling::strata() (take-all: n_h* = N_h).
# 2. Inclusion probability pi_hi = n_h*/N_h and design weight d_hi = 1/pi_hi.
# 3. Simulate response: each issued beekeeper responds with probability
#    RESP_TRUE[h], independently (random within strata, Section B4).
# 4. Non-response-adjusted weight w_hi = N_h / r_h for respondents.
# 5. Save the realised sample (issued beekeepers, with a response flag).
# =============================================================================

source("R/00_setup.R")
set.seed(SEED_SAMPLE)

pop   <- read.csv(PATH_POP)
alloc <- read.csv(file.path(DIR_TAB, "plan_allocation.csv"))
H     <- strata_def$stratum
n_iss <- setNames(alloc$issued, alloc$stratum)[H]
N_h   <- table(factor(pop$stratum, levels = H))

# ---- C5: every beekeeper on the register has a known, positive pi -------------
pi_h <- n_iss / as.numeric(N_h)
stopifnot(all(pi_h > 0), all(pi_h <= 1))

# ---- C3: stratified SRSWOR with sampling::strata() -----------------------------
frame <- pop[order(match(pop$stratum, H), pop$bk_id), ]          # sorted by stratum
s <- sampling::strata(frame, stratanames = "stratum", size = as.numeric(n_iss), method = "srswor")
smp <- sampling::getdata(frame, s)
stopifnot(nrow(smp) == sum(n_iss), !anyDuplicated(smp$bk_id),
          all(table(factor(smp$stratum, levels = H)) == n_iss))

smp$N_h           <- as.numeric(N_h[smp$stratum])
smp$n_issued      <- as.numeric(n_iss[smp$stratum])
smp$incl_prob     <- smp$Prob                                   # from sampling::strata()
stopifnot(isTRUE(all.equal(smp$incl_prob, smp$n_issued / smp$N_h)))
smp$design_weight <- 1 / smp$incl_prob

# ---- Simulated fieldwork ----------------------------------------------------------
smp$responded <- rbinom(nrow(smp), 1, RESP_TRUE[smp$stratum])
r_h <- tapply(smp$responded, factor(smp$stratum, levels = H), sum)
smp$n_resp_h <- as.numeric(r_h[smp$stratum])
smp$weight   <- ifelse(smp$responded == 1, smp$N_h / smp$n_resp_h, NA)   # w_hi = N_h / r_h

# survey answers are only observed for respondents
smp$colonies_lost <- ifelse(smp$responded == 1, smp$colonies_lost, NA)
smp$migratory     <- ifelse(smp$responded == 1, smp$migratory, NA)

sample_out <- smp[, c("bk_id", "stratum", "operation", "colonies_oct", "N_h", "n_issued",
                      "incl_prob", "design_weight", "responded", "n_resp_h", "weight",
                      "migratory", "colonies_lost")]
sample_out <- sample_out[order(match(sample_out$stratum, H), sample_out$bk_id), ]
rownames(sample_out) <- NULL

# ---- Checks -------------------------------------------------------------------------
stopifnot(abs(sum(sample_out$design_weight) - N) < 1e-6,                       # sum d = N
          abs(sum(sample_out$weight, na.rm = TRUE) - N) < 1e-6,                # sum w = N
          all(!is.na(sample_out$colonies_lost[sample_out$responded == 1])))
write.csv(sample_out, PATH_SAMPLE, row.names = FALSE)

# ---- Summary for the report (C1, C4, C5) ---------------------------------------------------
sample_summary <- data.frame(
  stratum       = H,
  N_h           = as.numeric(N_h),
  planned_resp  = setNames(alloc$n_resp, alloc$stratum)[H],
  issued        = as.numeric(n_iss),
  incl_prob     = as.numeric(pi_h),
  design_weight = 1 / as.numeric(pi_h),
  respondents   = as.numeric(r_h),
  resp_rate     = as.numeric(r_h / n_iss),
  weight        = as.numeric(N_h / r_h))
write.csv(sample_summary, file.path(DIR_TAB, "sample_summary.csv"), row.names = FALSE)

cat(sprintf("Sample written to %s: %d issued, %d respondents (%.1f%%)\n", PATH_SAMPLE,
            nrow(sample_out), sum(sample_out$responded), 100 * mean(sample_out$responded)))
cat(sprintf("Sum of inclusion probabilities over the register = %.0f (= issued)\n", sum(pi_h * N_h)))
print(sample_summary, digits = 4)
