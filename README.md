# SCAPIS Proteomics & Cardiovascular Risk Exercise

A biostatistics tutorial using a simulated SCAPIS-like dataset (n = 1,000) with clinical variables, coronary artery calcium scores (CACS), 25 plasma proteins (Olink NPX), and survival data. Students work through 53 progressively structured questions covering descriptive statistics, hypothesis testing, regression, multiple testing correction, missing data handling, and survival analysis.

## Project Structure

```
proteiomics_tutorial/
├── README.md                          # This file
├── SCAPIS_Proteomics_Exercise.Rproj   # RStudio project file
├── .gitignore                         # Version control exclusions
├── .Rprofile                          # renv activation
├── renv.lock                          # Package lockfile
│
├── data/
│   ├── raw/                           # Simulated data (generated, not tracked)
│   └── processed/                     # Cleaned/analysis-ready data
│
├── scripts/
│   ├── 00_generate_data.R             # Simulate the SCAPIS-like dataset
│   └── run_all.R                      # Master pipeline
│
├── docs/
│   ├── data_dictionary.md             # Variable definitions (53 columns)
│   └── INSTRUCTOR_GUIDE.md            # Teaching notes & answer key
│
├── output/
│   ├── figures/                       # Generated figures
│   └── tables/                        # Generated tables
│
└── scapis_exercise.Rmd                # Student exercise (53 questions)
```

## Getting Started

1. **Open the project** in RStudio by double-clicking `SCAPIS_Proteomics_Exercise.Rproj`.

2. **Restore packages** (first time only):
   ```r
   renv::restore()
   ```

3. **Generate the simulated dataset**:
   ```r
   source("scripts/00_generate_data.R")
   ```

4. **Open `scapis_exercise.Rmd`** and work through the questions.

Alternatively, run the full pipeline:
```r
source("scripts/run_all.R")
```

## Dataset Overview

The simulated dataset (`data/raw/scapis_simulated.csv`) contains 1,000 observations and 53 variables:

- **Identifiers**: participant ID
- **Demographics**: age, sex
- **Anthropometry**: BMI, waist circumference
- **Blood pressure**: systolic (SBP), diastolic (DBP)
- **Lipids**: total cholesterol, HDL, LDL, triglycerides
- **Glucose metabolism**: fasting glucose, HbA1c
- **Inflammation**: hs-CRP
- **Renal**: creatinine
- **Lifestyle**: smoking status, physical activity
- **Medications**: antihypertensive, statin, diabetes medication
- **Pulmonary**: FEV1%, FVC%
- **CACS**: continuous score, binary (0 vs >0), categorical (4 groups)
- **Proteins**: 25 Olink proteins (10 CACS-associated, 15 null)
- **Survival**: time-to-event, CVD event indicator, event type

See `docs/data_dictionary.md` for full variable definitions.

## Exercise Topics

| Section | Topic | Questions |
|---------|-------|-----------|
| 1 | Data Loading & Exploration | Q1–Q5 |
| 2 | Statistical Description | Q6–Q11 |
| 3 | Hypothesis Testing | Q12–Q17 |
| 4 | Table 1 (gtsummary, table1, htmlTable) | Q18–Q23 |
| 5 | Figures with ggplot2 | Q24–Q30 |
| 6 | Regression Analysis | Q31–Q37 |
| 7 | Multiple Testing Correction | Q38–Q42 |
| 8 | Missing Data Handling | Q43–Q47 |
| 9 | Survival / Cox Regression | Q48–Q53 |

## Requirements

R >= 4.0 with the following packages: tidyverse, gtsummary, survival, survminer, mice, naniar, corrplot, broom, ggrepel, htmlTable, table1, rmarkdown, MASS.
