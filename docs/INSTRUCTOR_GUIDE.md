# Instructor Guide

Teaching notes for the SCAPIS Proteomics & Cardiovascular Risk Exercise.

## Dataset Design Rationale

The simulated dataset is deliberately constructed to provide realistic teaching scenarios:

1. **Known ground truth**: 10 of the 25 proteins are truly associated with CACS, and 15 are null. This allows students to evaluate the performance of their analyses against a known answer.

2. **Controlled confounding**: Age and sex confound the protein–CACS associations. RETN and LEP are confounded by BMI (they will show up as marginally associated with CACS through BMI, but are not directly associated). This teaches students about confounding and adjustment.

3. **Realistic missing data**: A mix of MCAR (BMI, triglycerides, fasting glucose, proteins) and MAR (HDL depends on sex, hs-CRP depends on smoking, HbA1c depends on age) patterns lets students explore missingness mechanisms and compare complete-case vs. multiple imputation approaches.

4. **Survival outcomes**: The Weibull-generated event times with ~8% event rate provide enough events for meaningful Cox regression while reflecting realistic cardiovascular event rates.

## Which Proteins Are Truly Associated

### CACS-associated (10 proteins — should be significant)

| Protein | Effect size (NPX) | Expected detectability |
|---------|-------------------|----------------------|
| GDF15 | 0.55 | Strongest signal — should always be detected |
| NT_proBNP | 0.45 | Strong signal |
| SPP1 | 0.42 | Strong signal |
| IL6 | 0.40 | Strong (but 5% missing data may slightly reduce power) |
| MMP9 | 0.38 | Moderate-strong |
| PTX3 | 0.36 | Moderate |
| IL18 | 0.35 | Moderate |
| MMP12 | 0.32 | Moderate |
| MPO | 0.30 | Moderate (7% missing data) |
| LGALS3 | 0.28 | Weakest — may not survive Bonferroni correction |

### Null proteins (15 — should NOT be significant)

- Most will have p > 0.05, but some may reach nominal significance by chance.
- **RETN** and **LEP** may show marginal associations due to BMI confounding (BMI is weakly associated with CACS through shared risk factors). This is a good teaching point about confounding vs. true associations.

## Expected Results by Section

### Section 3 — Hypothesis Testing
- **Q12**: BMI t-test by CACS: may or may not be significant (no direct simulated association).
- **Q14**: Chi-squared smoking vs CACS: should be significant (smoking drives CACS in the latent model).
- **Q15**: ANOVA GDF15 by CACS categories: should be highly significant with a dose-response pattern.

### Section 5 — Figures
- **Q28**: Correlation heatmap should show a clear block of correlated CACS-associated proteins (upper-left 10x10 block).
- **Q29**: Volcano plot should clearly separate the 10 associated proteins (right side, high -log10 p) from the 15 null proteins (center, low -log10 p).

### Section 6 — Regression
- **Q31 vs Q32**: The cacs_binary coefficient for GDF15 should decrease slightly after adjusting for age and sex (confounders), demonstrating the concept of confounding.
- **Q33**: Logistic regression — age, sex, smoking, SBP, total cholesterol, and fasting glucose should all be significant predictors of CACS binary (they are in the data-generating model).
- **Q37**: The quadratic age term is unlikely to be significant (age effect is simulated as linear).

### Section 7 — Multiple Testing
- **Q38**: Most of the 10 associated proteins should have p < 0.05. A few null proteins may also reach nominal significance.
- **Q40**: Bonferroni correction (threshold = 0.002) should retain ~8–10 of the associated proteins. LGALS3 may be lost due to its small effect size.
- **Q41**: FDR correction should retain all 10 associated proteins while rejecting most null proteins.

### Section 8 — Missing Data
- **Q45**: HDL missingness should be significantly associated with sex (MAR mechanism by design).
- **Q47**: Multiple imputation estimates should be very similar to complete-case estimates for this dataset, since the missingness rates are modest (~3–8%).

### Section 9 — Survival Analysis
- **Q50**: KM curves by CACS category should show clear separation with log-rank p < 0.001. The >400 group should have the worst survival.
- **Q51**: Cox HRs for CACS categories: expect ~1.5 for 1-100, ~2.5 for 101-400, ~5 for >400 (approximate — depends on simulation seed).
- **Q53**: GDF15 should be a significant predictor in the Cox model (HR ~1.15–1.25 per NPX unit). AIC comparison should favor the model with GDF15.

## Timing Suggestions

| Section | Estimated time | Notes |
|---------|---------------|-------|
| 1. Data Loading & Exploration | 15–20 min | Straightforward warm-up |
| 2. Statistical Description | 20–25 min | Good for reviewing summary stats |
| 3. Hypothesis Testing | 25–30 min | Core statistical concepts |
| 4. Table 1 | 25–30 min | Hands-on with gtsummary/table1/htmlTable |
| 5. Figures | 30–40 min | ggplot2 practice; volcano plot is the highlight |
| 6. Regression | 30–40 min | Key epidemiological modeling |
| 7. Multiple Testing | 20–25 min | Critical for omics analyses |
| 8. Missing Data | 25–30 min | mice imputation takes time to run |
| 9. Survival Analysis | 30–35 min | Kaplan-Meier and Cox models |
| **Total** | **~3.5–4.5 hours** | Can be split across multiple sessions |

## Suggested Session Splits

**Option A — Two sessions (2 hours each)**:
- Session 1: Sections 1–5 (data handling, description, testing, tables, plots)
- Session 2: Sections 6–9 (regression, multiple testing, missing data, survival)

**Option B — Three sessions (~90 min each)**:
- Session 1: Sections 1–3
- Session 2: Sections 4–6
- Session 3: Sections 7–9

## Common Student Questions

**Q: Why do some null proteins appear significant?**
A: With 15 null tests at alpha = 0.05, you expect ~0.75 false positives on average. This is the whole point of Section 7 on multiple testing.

**Q: Why does adjusting for covariates change the protein association?**
A: Age and sex are confounders (they affect both the protein level and CACS). Adjustment removes the confounded portion, revealing the direct effect.

**Q: Why are the imputed results so similar to complete-case?**
A: The missingness rates are low (3-8%), so the bias from complete-case analysis is small. With higher missingness, the differences would be more pronounced.

## Dependencies

This exercise requires R >= 4.0 and the following packages:
tidyverse, gtsummary, survival, survminer, mice, naniar, corrplot, broom, ggrepel, htmlTable, table1, rmarkdown, MASS.
