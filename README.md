# Cardiovascular Disease (CVD) Risk Factor Cohort Datasets

This repository contains the synthetic cohort datasets used for the comparative information-theoretic and quantum-inspired analysis of cardiovascular disease profiles in Indian and global populations.

## Datasets

The datasets consist of $N=8,000$ synthetic patient profiles each:

1. **`reference-india-cohort.csv`**: Modeled cohort representing high cardiometabolic risk-factor architectures, parameterized using prevalence statistics from the ICMR-INDIAB and PURE studies.
2. **`reference-global-cohort.csv`**: Modeled comparator cohort representing global risk-factor architectures, parameterized using prevalence statistics from NHANES, UK Biobank, and the WHO Global Health Observatory.

## Data Structure

Each dataset contains the following 6 discretized categorical risk-factor variables:

* **`age`**: 20-39, 40-59, 60-79, or >=80.
* **`sex`**: Female or Male.
* **`smoking`**: Non-Smoker or Current Smoker.
* **`diabetes`**: Normoglycemic or Diabetic.
* **`bp`**: Blood pressure status (Normal, Elevated, Stage 1 HTN, Stage 2 HTN).
* **`cholesterol`**: Cholesterol status (Desirable, Borderline High, High, Very High).
