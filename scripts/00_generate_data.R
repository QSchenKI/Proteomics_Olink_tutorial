###############################################################################
#  00_generate_data.R
#  Generate a simulated SCAPIS-like dataset (n = 1000)
#  -----------------------------------------------------------------------
#  Variables: demographics, anthropometry, blood pressure, lipids, glucose,
#             inflammatory markers, renal, lifestyle, medications, pulmonary,
#             coronary artery calcium score (CACS), 25 Olink proteins,
#             survival / CVD events, and engineered missing data.
###############################################################################

library(MASS)      # mvrnorm()
set.seed(42)

n <- 1000

# ── 1. Demographics ─────────────────────────────────────────────────────────
id       <- sprintf("SCAPIS_%04d", seq_len(n))
sex      <- sample(c("Male", "Female"), n, replace = TRUE, prob = c(0.49, 0.51))
is_male  <- as.integer(sex == "Male")
age      <- round(runif(n, 50, 64), 1)

# ── 2. Lifestyle ─────────────────────────────────────────────────────────────
smoking <- sample(c("Never", "Former", "Current"), n, replace = TRUE,
                  prob = c(0.45, 0.35, 0.20))
physical_activity <- sample(c("Low", "Moderate", "High"), n, replace = TRUE,
                            prob = c(0.25, 0.50, 0.25))

# ── 3. Anthropometry ────────────────────────────────────────────────────────
bmi     <- round(rnorm(n, mean = 27, sd = 4.3), 1)
bmi     <- pmax(bmi, 16)                        # floor at 16
waist_cm <- round(ifelse(is_male,
                         rnorm(n, 98, 12),
                         rnorm(n, 88, 13)), 0)
waist_cm <- pmax(waist_cm, 55)

# ── 4. Blood pressure ──────────────────────────────────────────────────────
sbp <- round(rnorm(n, 130, 18), 0)
sbp <- pmax(sbp, 85)
dbp <- round(rnorm(n, 80, 10), 0)
dbp <- pmax(dbp, 50)

# ── 5. Lipids ───────────────────────────────────────────────────────────────
total_chol    <- round(rnorm(n, 5.5, 1.0), 2)
total_chol    <- pmax(total_chol, 2.5)
hdl           <- round(ifelse(is_male,
                              rnorm(n, 1.3, 0.35),
                              rnorm(n, 1.6, 0.40)), 2)
hdl           <- pmax(hdl, 0.5)
triglycerides <- round(exp(rnorm(n, log(1.4), 0.5)), 2)
triglycerides <- pmax(triglycerides, 0.3)
# Friedewald-derived LDL
ldl           <- round(total_chol - hdl - triglycerides / 2.2, 2)
ldl           <- pmax(ldl, 0.5)

# ── 6. Glucose ──────────────────────────────────────────────────────────────
fasting_glucose <- round(rnorm(n, 5.6, 0.9), 2)
fasting_glucose <- pmax(fasting_glucose, 3.0)
hba1c           <- round(rnorm(n, 38, 6), 0)
hba1c           <- pmax(hba1c, 20)

# ── 7. Inflammatory ────────────────────────────────────────────────────────
hs_crp <- round(exp(rnorm(n, log(1.5), 0.9)), 2)
hs_crp <- pmax(hs_crp, 0.1)

# ── 8. Renal ────────────────────────────────────────────────────────────────
creatinine <- round(ifelse(is_male,
                           rnorm(n, 88, 14),
                           rnorm(n, 70, 12)), 0)
creatinine <- pmax(creatinine, 35)

# ── 9. Medications ──────────────────────────────────────────────────────────
antihypertensive <- rbinom(n, 1, 0.30)
statin           <- rbinom(n, 1, 0.20)
diabetes_med     <- rbinom(n, 1, 0.08)

# ── 10. Pulmonary ──────────────────────────────────────────────────────────
fev1_pct <- round(rnorm(n, 95, 14), 1)
fev1_pct <- pmax(fev1_pct, 30)
fvc_pct  <- round(rnorm(n, 98, 13), 1)
fvc_pct  <- pmax(fvc_pct, 35)

# ── 11. CACS via latent propensity ─────────────────────────────────────────
smoke_num <- ifelse(smoking == "Current", 2,
             ifelse(smoking == "Former", 1, 0))

# Latent score  (higher → more calcium)
# Intercept calibrated so ~42% have CACS > 0 (58% = 0)
latent <- -0.90 +
  0.06  * (age - 57) +
  0.50  * is_male +
  0.30  * smoke_num +
  0.010 * (sbp - 130) +
  0.15  * (total_chol - 5.5) +
  0.20  * (fasting_glucose - 5.6) +
  rnorm(n, 0, 1)

# Probability of CACS > 0
prob_pos <- plogis(latent)
cacs_positive <- rbinom(n, 1, prob_pos)

# Among those with CACS > 0, sample a continuous score
cacs <- rep(0, n)
pos_idx <- which(cacs_positive == 1)
# Log-normal among positives, with heavier tail for >400 group
cacs[pos_idx] <- round(exp(rnorm(length(pos_idx), log(80), 1.3)))
cacs <- pmax(cacs, 0)

# Categorical CACS
cacs_cat <- cut(cacs,
                breaks = c(-Inf, 0, 100, 400, Inf),
                labels = c("0", "1-100", "101-400", ">400"))
cacs_binary <- as.integer(cacs > 0)

# ── 12. Proteomics (25 proteins, Olink NPX units) ─────────────────────────
#  10 CACS-associated + 15 null
protein_names_assoc <- c("NT_proBNP", "GDF15", "IL6", "IL18", "MPO",
                         "MMP9", "MMP12", "SPP1", "LGALS3", "PTX3")
protein_names_null  <- c("RETN", "LEP", "REN", "PCSK9", "FABP4",
                         "F3", "THBD", "SELP", "DCN", "CHIT1",
                         "SORT1", "SELE", "CSTB", "GRN", "CDH5")
all_proteins <- c(protein_names_assoc, protein_names_null)
n_prot <- length(all_proteins)

# Base correlation matrix: block correlations among associated proteins
Sigma <- diag(n_prot)
# Mild positive correlation among the 10 associated proteins
for (i in 1:10) {
  for (j in 1:10) {
    if (i != j) Sigma[i, j] <- 0.20
  }
}
# Small correlations among some null proteins (biological realism)
for (i in 11:15) {
  for (j in 11:15) {
    if (i != j) Sigma[i, j] <- 0.10
  }
}

# Simulate base protein levels (mean ~5 NPX, SD ~1.2)
base_prot <- mvrnorm(n, mu = rep(5, n_prot), Sigma = Sigma * 1.44)

# Add CACS-dependent signal to the 10 associated proteins
# Effect sizes (NPX shift per unit of cacs_binary)
effect_sizes <- c(0.45, 0.55, 0.40, 0.35, 0.30,
                  0.38, 0.32, 0.42, 0.28, 0.36)
for (k in 1:10) {
  base_prot[, k] <- base_prot[, k] +
    effect_sizes[k] * cacs_binary +
    0.02 * (age - 57) +        # mild age confounding
    0.10 * is_male +            # mild sex confounding
    rnorm(n, 0, 0.25)           # extra noise
}

# Add some clinical confounding to a few null proteins (realistic nuisance)
base_prot[, 11] <- base_prot[, 11] + 0.3 * (bmi - 27) / 4.3  # RETN ~ BMI
base_prot[, 12] <- base_prot[, 12] + 0.4 * (bmi - 27) / 4.3  # LEP  ~ BMI

prot_df <- as.data.frame(base_prot)
colnames(prot_df) <- all_proteins
prot_df <- round(prot_df, 3)

# ── 13. Survival data ─────────────────────────────────────────────────────
# Baseline hazard ~ Weibull, ~8% event rate over ~5 yr follow-up
max_followup <- 6
lambda0 <- 0.007   # baseline scale (calibrated for ~8% event rate)
shape   <- 1.2

# Linear predictor for hazard
lp_surv <- 0.04 * (age - 57) +
  0.30 * is_male +
  0.50 * (cacs_cat == "1-100") +
  1.00 * (cacs_cat == "101-400") +
  1.60 * (cacs_cat == ">400") +
  0.20 * smoke_num +
  0.15 * ((prot_df$GDF15 - 5) / 1.2)   # GDF15 effect

# Weibull event times
U <- runif(n)
event_time <- (-log(U) / (lambda0 * exp(lp_surv)))^(1 / shape)

# Censoring: uniform(3, max_followup) to get median ~4.5 yr follow-up
censor_time <- runif(n, 3, max_followup)

time_to_event <- pmin(event_time, censor_time)
time_to_event <- round(pmin(time_to_event, max_followup), 2)
cvd_event     <- as.integer(event_time <= censor_time & event_time <= max_followup)

# Event type (among events): MI ~60%, Stroke ~40%
event_type <- rep(NA_character_, n)
n_events <- sum(cvd_event == 1)
if (n_events > 0) {
  event_type[cvd_event == 1] <- sample(c("MI", "Stroke"), n_events,
                                       replace = TRUE, prob = c(0.60, 0.40))
}

# ── 14. Assemble data frame ───────────────────────────────────────────────
scapis <- data.frame(
  id                = id,
  age               = age,
  sex               = sex,
  bmi               = bmi,
  waist_cm          = waist_cm,
  sbp               = sbp,
  dbp               = dbp,
  total_chol        = total_chol,
  hdl               = hdl,
  ldl               = ldl,
  triglycerides     = triglycerides,
  fasting_glucose   = fasting_glucose,
  hba1c             = hba1c,
  hs_crp            = hs_crp,
  creatinine        = creatinine,
  smoking           = smoking,
  physical_activity = physical_activity,
  antihypertensive  = antihypertensive,
  statin            = statin,
  diabetes_med      = diabetes_med,
  fev1_pct          = fev1_pct,
  fvc_pct           = fvc_pct,
  cacs              = cacs,
  cacs_binary       = cacs_binary,
  cacs_cat          = as.character(cacs_cat),
  prot_df,
  time_to_event     = time_to_event,
  cvd_event         = cvd_event,
  event_type        = event_type,
  stringsAsFactors  = FALSE
)

# ── 15. Introduce missing data ─────────────────────────────────────────────
introduce_mcar <- function(x, prop) {
  miss <- sample(seq_along(x), size = round(prop * length(x)))
  x[miss] <- NA
  x
}

introduce_mar <- function(x, indicator, prob_if_1, prob_if_0) {
  # Higher missingness when indicator == 1
  p <- ifelse(indicator == 1, prob_if_1, prob_if_0)
  miss <- rbinom(length(x), 1, p) == 1
  x[miss] <- NA
  x
}

# MCAR
scapis$bmi             <- introduce_mcar(scapis$bmi, 0.03)
scapis$triglycerides   <- introduce_mcar(scapis$triglycerides, 0.05)
scapis$fasting_glucose <- introduce_mcar(scapis$fasting_glucose, 0.04)

# MAR
scapis$hdl   <- introduce_mar(scapis$hdl,
                              indicator = as.integer(scapis$sex == "Male"),
                              prob_if_1 = 0.10, prob_if_0 = 0.04)
scapis$hs_crp <- introduce_mar(scapis$hs_crp,
                               indicator = as.integer(scapis$smoking == "Current"),
                               prob_if_1 = 0.14, prob_if_0 = 0.06)
scapis$hba1c <- introduce_mar(scapis$hba1c,
                              indicator = as.integer(age > 60),
                              prob_if_1 = 0.10, prob_if_0 = 0.04)

# Protein missingness (3-4 columns, ~5-8%)
scapis$GDF15  <- introduce_mcar(scapis$GDF15, 0.06)
scapis$IL6    <- introduce_mcar(scapis$IL6, 0.05)
scapis$MPO    <- introduce_mcar(scapis$MPO, 0.07)
scapis$FABP4  <- introduce_mcar(scapis$FABP4, 0.08)

# ── 16. Write data ─────────────────────────────────────────────────────────
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
write.csv(scapis, "data/raw/scapis_simulated.csv", row.names = FALSE)
saveRDS(scapis, "data/raw/scapis_simulated.rds")

# ── 17. Quick summary ─────────────────────────────────────────────────────
cat("Dataset written: data/raw/scapis_simulated.csv\n")
cat("Dataset written: data/raw/scapis_simulated.rds\n")
cat("Dimensions:", nrow(scapis), "x", ncol(scapis), "\n")
cat("CACS distribution:\n")
print(table(scapis$cacs_cat, useNA = "ifany"))
cat("\nCACS binary (0 vs >0):\n")
print(table(scapis$cacs_binary))
cat("\nCVD event rate:", mean(scapis$cvd_event), "\n")
cat("Median follow-up:", median(scapis$time_to_event), "years\n")
cat("\nMissing values per column:\n")
miss <- colSums(is.na(scapis))
print(miss[miss > 0])
