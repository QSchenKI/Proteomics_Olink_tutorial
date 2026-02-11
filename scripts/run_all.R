###############################################################################
#  run_all.R
#  Master pipeline: generate data → render exercise
###############################################################################

cat("=== SCAPIS Proteomics Exercise Pipeline ===\n\n")

# ── 1. Generate simulated data ──────────────────────────────────────────────
cat("Step 1: Generating simulated data...\n")
source("scripts/00_generate_data.R")
cat("\n")

# ── 2. Render exercise to HTML ──────────────────────────────────────────────
cat("Step 2: Rendering scapis_exercise.Rmd to HTML...\n")
rmarkdown::render("scapis_exercise.Rmd", quiet = TRUE)
cat("  Rendered: scapis_exercise.html\n\n")

# ── 3. Verify outputs ──────────────────────────────────────────────────────
cat("=== Output verification ===\n")
expected <- c(
  "data/raw/scapis_simulated.csv",
  "data/raw/scapis_simulated.rds",
  "scapis_exercise.html"
)
for (f in expected) {
  status <- if (file.exists(f)) "OK" else "MISSING"
  cat(sprintf("  [%s] %s\n", status, f))
}
cat("\nPipeline complete.\n")
