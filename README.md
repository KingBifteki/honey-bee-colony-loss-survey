# Estimating Winter Honey Bee Colony Losses

Design and simulation of a stratified sample survey of registered beekeepers.
STAT4102 Applied Sampling mini project, University of the Witwatersrand, 2026 (Jun Xi Yip, 2554427).

The project plans a probability survey to estimate the total number of colonies lost over winter 2025/26
by the 10,000 registered beekeepers of a hypothetical country that takes part in the COLOSS monitoring.
The population is simulated and calibrated to the COLOSS survey of winter 2019/20 (Gray et al., 2023).
The design is stratified random sampling by operation size (hobby, side-line, commercial) with Neyman
allocation and a take-all commercial stratum, a planned 77% response rate and non-response weighting within
strata. The total is estimated with a separate ratio estimator that uses the colony counts on the register,
and the design is compared with simple random sampling without replacement (SRSWOR) in one realised sample
and over 10,000 simulated repetitions.

## Files

| Path | Contents |
|---|---|
| `report/report.qmd` | Quarto source of the report; every result in the text is read from the data and output files |
| `report/report.pdf` | The rendered report |
| `report/ai_declaration.pdf` | The signed generative-AI declaration; `report.qmd` adds it as Appendix D when the file is present and renders without it otherwise |
| `run_all.R` | Runs the whole analysis and renders the report |
| `R/00_setup.R` | Packages, seeds, constants and the COLOSS calibration values (sourced by every script) |
| `R/01_population.R` | Generates the finite population of 10,000 beekeepers |
| `R/02_sample_size.R` | Section B: baseline, sample size, allocation, sensitivity analysis, Figures 2 and 3 |
| `R/03_sample.R` | Section C: draws the stratified sample and simulates fieldwork |
| `R/04_estimation.R` | Section D: design-based estimation with the `survey` package and the SRSWOR benchmark |
| `R/05_simulation.R` | Section E: Monte Carlo study (10,000 repetitions) and Figure 4 |
| `data/population.csv` | The complete finite population (10,000 rows) |
| `data/sample.csv` | The realised stratified sample (2,452 issued beekeepers) |
| `data/sample_srswor_benchmark.csv` | The SRSWOR benchmark sample (2,452 issued beekeepers) |
| `output/tables/` | Every table the report uses, as CSV files |
| `output/figures/` | Figures 1 to 4 |
| `BeeSurvey.Rproj` | RStudio project file |

The course brief is not included, as it is the lecturer's document.

## How to run

From the project folder:

```bash
Rscript run_all.R
```

or open `BeeSurvey.Rproj` in RStudio and source `run_all.R`. Each script runs in its own R session, and the
whole run takes about two minutes. The report is rendered with Quarto, which needs a LaTeX installation
(for example TinyTeX); if Quarto is not found, `run_all.R` reproduces all results and asks you to click
Render on `report/report.qmd` in RStudio.

Requirements: R 4.5 with the packages `dplyr`, `ggplot2`, `survey`, `sampling`, `ragg`, `scales`, `knitr`
and `kableExtra` (`run_all.R` lists any that are missing). The results were produced with R 4.5.2,
survey 4.5, sampling 2.11, ggplot2 4.0.2 and Quarto 1.8.25.

The random seeds are set in `R/00_setup.R`: 4102 for the population, 2026 for the sample and fieldwork,
2027 for the SRSWOR benchmark and 1021 for the simulation. A fresh run reproduces every data and output
file exactly.

## Data

`data/population.csv`, one row per beekeeper:

| Column | Meaning |
|---|---|
| `bk_id` | Beekeeper identifier |
| `stratum` | `HOB`, `SID` or `COM` |
| `operation` | Operation-size class: hobby (up to 50 colonies), side-line (51 to 150) or commercial (more than 150) |
| `migratory` | 1 if the beekeeper migrates colonies, 0 otherwise |
| `colonies_oct` | Colonies going into winter on 1 October 2025 ($m_i$) |
| `colonies_lost` | Colonies lost by 31 March 2026 ($y_i$) |
| `loss_pct` | Colonies lost as a percentage of colonies going into winter |

`data/sample.csv`, one row per issued beekeeper:

| Column | Meaning |
|---|---|
| `bk_id`, `stratum`, `operation` | As in the population file |
| `colonies_oct` | Colonies on the register ($m_i$), an auxiliary variable |
| `N_h` | Number of beekeepers in the stratum |
| `n_issued` | Number of beekeepers issued in the stratum ($n_h^*$) |
| `incl_prob` | Inclusion probability $\pi_{hi}=n_h^*/N_h$ |
| `design_weight` | Design weight $d_{hi}=1/\pi_{hi}$ |
| `responded` | 1 if the beekeeper responded, 0 otherwise |
| `n_resp_h` | Number of respondents in the stratum ($r_h$) |
| `weight` | Non-response-adjusted weight $w_{hi}=N_h/r_h$ |
| `migratory`, `colonies_lost` | Survey answers |

`weight`, `migratory` and `colonies_lost` are blank for the 531 non-respondents, whose answers were not
obtained. `data/sample_srswor_benchmark.csv` holds the SRSWOR benchmark with the columns `bk_id`, `stratum`,
`colonies_oct`, `incl_prob`, `responded`, `weight`, `migratory` and `colonies_lost`, and a single weight
$N/r$ for its 1,903 respondents.

All data are simulated; no real beekeepers are involved.
