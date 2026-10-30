# =============================================================================
# run_all.R - reproduce every result, table, figure and the report
# National Colony Loss Survey 2026 (STAT4102 Mini Project)
#
# Run from the project root:   Rscript run_all.R
# (or open BeeSurvey.Rproj in RStudio and source this file)
#
# Each script runs in its own R process, as when run one by one, and sets its
# own seed (R/00_setup.R), so the outputs are identical on every run:
#   R/01_population.R   finite population            -> data/population.csv
#   R/02_sample_size.R  sample size and allocation   -> output/tables/plan_*.csv, figures 2-3
#   R/03_sample.R       stratified sample, fieldwork -> data/sample.csv
#   R/04_estimation.R   design-based estimation      -> output/tables/estimates*.csv
#   R/05_simulation.R   Monte Carlo study (about 1.5 minutes) -> output/tables/sim_*.csv, figure 4
# Then the report is rendered with Quarto           -> report/report.pdf
# =============================================================================

needed  <- c("dplyr", "ggplot2", "survey", "sampling", "ragg", "scales", "knitr", "kableExtra")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these packages first: install.packages(c(",
                          paste0('"', missing, '"', collapse = ", "), "))")
if (!file.exists("R/00_setup.R")) stop("Run this script from the project root (the folder holding R/ and report/).")

rscript <- file.path(R.home("bin"), "Rscript")
for (s in c("01_population", "02_sample_size", "03_sample", "04_estimation", "05_simulation")) {
  message("Running R/", s, ".R ...")
  if (system2(rscript, file.path("R", paste0(s, ".R"))) != 0) stop("R/", s, ".R failed")
}

# ---- Render the report (needs Quarto with a LaTeX installation) -----------------
quarto <- Sys.which("quarto")
bundled <- "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto"   # RStudio's copy on macOS
if (!nzchar(quarto) && file.exists(bundled)) quarto <- bundled
if (nzchar(quarto)) {
  if (system2(quarto, c("render", "report/report.qmd", "--to", "pdf")) != 0) stop("Rendering report/report.qmd failed")
  message("Done: report/report.pdf")
} else {
  message("Results reproduced. Quarto was not found: open report/report.qmd in RStudio and click Render.")
}
