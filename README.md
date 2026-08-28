# Cardiovascular Disease (CVD) Risk Factor Cohort Datasets

This repository contains the synthetic cohort datasets used for the comparative information-theoretic and quantum-inspired analysis of cardiovascular disease profiles in Indian and global populations.

## Datasets

The datasets consist of $N=10,000$ synthetic patient profiles each:

1. **`reference-india-cohort.csv`**: Modeled cohort representing high cardiometabolic risk-factor architectures, parameterized using prevalence statistics from the ICMR-INDIAB and PURE studies.
2. **`reference-global-cohort.csv`**: Modeled comparator cohort representing global risk-factor architectures, parameterized using prevalence statistics from NHANES, UK Biobank, and the WHO Global Health Observatory.

## Data Structure

Each dataset contains the following 8 risk-factor variables (6 discretized categorical variables, plus one continuous adiposity measure and its derived binary category):

* **`age`**: 20-39, 40-59, 60-79, or >=80.
* **`sex`**: Female or Male.
* **`smoking`**: Non-Smoker or Current Smoker.
* **`diabetes`**: Normoglycemic or Diabetic.
* **`bp`**: Blood pressure status (Normal, Elevated, Stage 1 HTN, Stage 2 HTN).
* **`cholesterol`**: Cholesterol status (Desirable, Borderline High, High, Very High).
* **`bmi`**: Body mass index in kg/m$^2$, continuous, reported to one decimal place.
* **`obesity`**: Obese if `bmi` $\geq 30.0$ kg/m$^2$, otherwise Non-Obese.

## Adiposity Variables (`bmi` and `obesity`)

### Obesity threshold

Obesity uses the WHO class-I threshold of 30.0 kg/m$^2$, applied **identically to both cohorts**. This follows the study's stated discretization policy that thresholds are held common across cohorts so that comparisons are not driven by differing bin definitions. No BMI value in either file falls exactly on 30.0, so the `>= 30` and `> 30` rules produce identical labels.

Note that the Asian-Indian / WHO Asia-Pacific cut-offs for South Asian populations are lower (overweight $\geq 23$, generalized obesity $\geq 25$ kg/m$^2$). Applying a single WHO threshold here is deliberate: it keeps the two cohorts comparable, and it is what makes the low-BMI / high-metabolic-risk contrast visible in the data. Analyses using the Asian-Indian criterion can be derived directly from the continuous `bmi` column.

### Statistical construction

`bmi` is **not** sampled independently. It is drawn conditionally on the six risk factors already present in each profile, so that the column carries joint structure for the entropy, systemic mutual information (SMI), and fidelity analyses. An independently sampled BMI column would share approximately 0 bits of mutual information with every other factor and would be uninformative in that framework.

For profile $i$, the conditional mean is additive in the risk factors,

$$\mu_i = \text{base} + b_{\text{age}} + b_{\text{sex}} + b_{\text{smoking}} + b_{\text{diabetes}} + b_{\text{bp}} + b_{\text{chol}}$$

and the realized value is lognormal about that mean,

$$\mathrm{BMI}_i = \mu_i \exp\left(\sigma z_i - \sigma^2/2\right), \qquad z_i \sim \mathcal{N}(0,1)$$

which reproduces the right skew of empirical adult BMI distributions while preserving $\mathbb{E}[\mathrm{BMI}_i] = \mu_i$. Values are clipped to [15.0, 55.0] and rounded to one decimal place. `base` and $\sigma$ are solved numerically per cohort so the final rounded data reproduce the target mean BMI and target obesity prevalence exactly.

Effect directions follow established epidemiology: BMI peaks in midlife and falls in the oldest bin; current smokers average roughly 1 kg/m$^2$ lower; and diabetes, higher blood-pressure stage, and worse lipid status are each associated with higher BMI. The diabetes-BMI gradient is deliberately **weaker in the Indian cohort** (approximately 1.6 vs 3.1 kg/m$^2$) to reflect the "South Asian Paradox" — high cardiometabolic risk at comparatively lower body mass index.

### Cohort targets and sources

| | Indian cohort | Global cohort |
|---|---|---|
| Mean BMI (kg/m$^2$) | 24.6 | 27.6 |
| Obesity ($\geq 30$) | 10.0% | 27.0% |
| Overweight+ ($\geq 25$) | 43.2% | 69.9% |
| Underweight ($< 18.5$) | 4.4% | 0.7% |
| SD | 4.04 | 4.56 |

* **Indian cohort**: NCD-RisC and NFHS-5 report a mean adult BMI near 22.5 kg/m$^2$ and WHO-threshold obesity near 4-6%; ICMR-INDIAB reports generalized obesity (Asian-Indian cut-off, $\geq 25$) at 28.6%. Because this cohort is a deliberately high-cardiometabolic-risk scenario (diabetes ~28%, roughly 2.5x the ICMR-INDIAB national figure), the targets sit above these national averages while remaining well below the comparator.
* **Global cohort**: a blend of the high-resource sources used for that cohort — UK Biobank (mean BMI ~27.4, obesity ~24%) and NHANES (mean BMI ~29.8, obesity ~42%) — tempered toward the WHO Global Health Observatory's lower worldwide figure.

### Resulting information-theoretic quantities

System entropy $S(\rho_X)$ of the new categorical variable, in bits:

| Risk factor | India | Global |
|---|---|---|
| Obesity | 0.47 | 0.84 |

Systemic mutual information with obesity, in bits:

| Pair | India | Global |
|---|---|---|
| SMI(Obesity : Age) | 0.017 | 0.029 |
| SMI(Obesity : Sex) | 0.003 | 0.000 |
| SMI(Obesity : Smoking) | 0.002 | 0.000 |
| SMI(Obesity : Diabetes) | 0.073 | 0.085 |
| SMI(Obesity : Hypertension) | 0.075 | 0.073 |
| SMI(Obesity : Dyslipidemia) | 0.046 | 0.013 |

The two burdens move in **opposite directions across the cohorts**: the Indian cohort has 10.0% obesity with 28.0% diabetes, against 27.0% obesity with 16.2% diabetes globally — roughly 0.37x the obesity prevalence but 1.7x the diabetes prevalence. Diabetes prevalence among non-obese profiles is 22.9% in the Indian cohort versus 8.1% in the global comparator, and the mean BMI of diabetic profiles is 27.9 versus 31.8 kg/m$^2$. The obesity-diabetes SMI is correspondingly lower in the Indian cohort, consistent with metabolic risk that is less tightly coupled to generalized adiposity.

### Effect on the categorical state space

The six original variables span $4 \times 2 \times 2 \times 2 \times 4 \times 4 = 512$ basis states. Including the binary `obesity` variable extends this to $4 \times 2 \times 2 \times 2 \times 4 \times 4 \times 2 = 1024$, so the population state matrix $\rho$ becomes $1024 \times 1024$ if obesity is included in the tensor-product encoding. The continuous `bmi` column does not enter the categorical encoding directly; it is provided for descriptive statistics, for sensitivity analysis under alternative thresholds, and for any future continuous-variable extension.

## Reproducibility

The adiposity columns are regenerated deterministically (fixed seeds, Python standard library only, no third-party dependencies) by running from the repository root:

```bash
python3 scripts/generate_bmi_obesity.py
```

The script appends to each row without modifying the six original columns. All model parameters, cohort targets, and source citations are documented inline in `scripts/generate_bmi_obesity.py`.
