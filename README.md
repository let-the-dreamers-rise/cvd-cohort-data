# Cardiovascular Disease (CVD) Risk Factor Cohort Datasets

This repository contains the synthetic cohort datasets used for the comparative information-theoretic and quantum-inspired analysis of cardiovascular disease profiles in Indian and global populations.

## Datasets

The datasets consist of $N=10,000$ synthetic patient profiles each:

1. **`reference-india-cohort.csv`**: Modeled cohort representing high cardiometabolic risk-factor architectures, parameterized using prevalence statistics from the ICMR-INDIAB and PURE studies.
2. **`reference-global-cohort.csv`**: Modeled comparator cohort representing global risk-factor architectures, parameterized using prevalence statistics from NHANES, UK Biobank, and the WHO Global Health Observatory.

## Data Structure

Each dataset contains the following 8 risk-factor variables:

* **`age`**: 20-39, 40-59, 60-79, or >=80.
* **`sex`**: Female or Male.
* **`smoking`**: Non-Smoker or Current Smoker.
* **`diabetes`**: Normoglycemic or Diabetic.
* **`bp`**: Blood pressure status (Normal, Elevated, Stage 1 HTN, Stage 2 HTN).
* **`cholesterol`**: Cholesterol status (Desirable, Borderline High, High, Very High).
* **`bmi`**: Body mass index in kg/m2, continuous, to one decimal place.
* **`obesity`**: Obese if `bmi` is 30.0 or above, otherwise Non-Obese.

## Notes on BMI and Obesity

The obesity threshold of 30.0 kg/m2 is the WHO cut-off, applied the same way to both cohorts so that comparisons are not driven by differing category definitions. No BMI value falls exactly on 30.0, so the "30 or above" and "greater than 30" rules give the same labels. The lower Asian-Indian cut-offs (23 and 25 kg/m2) can be applied directly to the continuous `bmi` column if needed.

BMI is generated conditionally on the six other risk factors rather than independently, so that it carries joint structure for the entropy, mutual information, and fidelity analyses. Effect directions follow standard epidemiology: BMI peaks in midlife and declines in the oldest age bin, current smokers average about 1 kg/m2 lower, and diabetes, higher blood pressure stage, and worse lipid status are each associated with higher BMI. The diabetes-BMI gradient is weaker in the Indian cohort to reflect the South Asian pattern of high cardiometabolic risk at comparatively lower body mass index.

The Indian cohort has a mean BMI of 24.6 with 10.0% obesity; the global cohort has a mean BMI of 27.6 with 27.0% obesity. These targets follow the sources already cited for each cohort.

Including the binary `obesity` variable extends the categorical state space from 512 to 1024 states. The continuous `bmi` column is not part of the categorical encoding.

## Reproducibility

The BMI and obesity columns are regenerated deterministically by running the following from the repository root:

```bash
python3 scripts/generate_bmi_obesity.py
```

The script uses fixed seeds and the Python standard library only. It appends to each row without modifying the six original columns, and can be re-run safely. Model parameters and source citations are documented inside the script.
