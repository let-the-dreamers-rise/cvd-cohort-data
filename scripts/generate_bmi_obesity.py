#!/usr/bin/env python3
"""Append `bmi` and `obesity` columns to the CVD risk-factor cohort datasets.

The two cohorts in this repository are modeled (synthetic) populations whose
categorical risk factors were parameterized from published epidemiological
sources.  This script extends them with the two adiposity variables using the
same modeling philosophy: BMI is *not* drawn independently, but conditionally
on the six risk factors already present in each row, so that the resulting
column carries the joint structure the information-theoretic analysis measures
(entropy, systemic mutual information, fidelity).  An independently sampled
BMI column would have ~0 bits of mutual information with every other factor
and would therefore be uninformative for that analysis.

Model
-----
For profile i, the conditional mean BMI is additive in the risk factors,

    mu_i = base + b_age[age_i] + b_sex[sex_i] + b_smoking[smoking_i]
                + b_diabetes[diabetes_i] + b_bp[bp_i] + b_chol[cholesterol_i]

and the realized value is drawn from a lognormal with that mean,

    BMI_i = mu_i * exp(sigma * z_i - sigma^2 / 2),   z_i ~ N(0, 1)

which reproduces the right skew of empirical adult BMI distributions while
preserving E[BMI_i] = mu_i.  Values are clipped to a physiologically plausible
range and rounded to one decimal place.  `base` and `sigma` are solved
numerically so that each cohort reproduces its target mean BMI and its target
obesity prevalence exactly on the final, rounded data.

Obesity is the WHO class-I threshold applied identically to both cohorts:
BMI >= 30.0 kg/m^2 -> "Obese", otherwise "Non-Obese".  Rounded values that
land exactly on 30.0 are nudged off the boundary so that the `>= 30` and
`> 30` rules give identical labels.

Parameterization
----------------
Cohort-level targets follow the sources already cited for each cohort:

* Indian cohort - NCD-RisC and NFHS-5 report a mean adult BMI near 22.5 kg/m^2
  and WHO-threshold (>=30) obesity near 4-6%; ICMR-INDIAB reports generalized
  obesity (Asian-Indian cut-off, >=25) at 28.6%.  Because this cohort is a
  deliberately high-cardiometabolic-risk scenario (diabetes ~28%, roughly
  2.5x the ICMR-INDIAB national figure), the targets are set above the
  national averages but still well below the comparator.
* Global cohort - a blend of the high-resource sources used for that cohort:
  UK Biobank (mean BMI ~27.4, obesity ~24%) and NHANES (mean BMI ~29.8,
  obesity ~42%), tempered toward the WHO Global Health Observatory's lower
  worldwide figure.

Within-cohort effect sizes follow well-established directions: BMI peaks in
midlife and falls in the oldest bin; current smokers average roughly 1 kg/m^2
lower; and diabetes, higher blood-pressure stage, and worse lipid status are
each associated with higher BMI.  The diabetes-BMI gradient is deliberately
weaker in the Indian cohort (~1.6 vs ~3.1 kg/m^2) to reflect the "South Asian
Paradox" the paper describes: high cardiometabolic risk at comparatively lower
body mass index.

Run from the repository root:

    python3 scripts/generate_bmi_obesity.py
"""

import csv
import math
import os
import random

OBESITY_THRESHOLD = 30.0  # kg/m^2, WHO class-I obesity, applied to both cohorts
BMI_MIN, BMI_MAX = 15.0, 55.0

# Risk-factor effects on mean BMI (kg/m^2), relative to the fitted cohort base.
COHORTS = {
    "reference-india-cohort.csv": {
        "seed": 20240517,
        "target_mean_bmi": 24.6,
        "target_obesity_rate": 0.10,
        "effects": {
            "age": {"20-39": -1.00, "40-59": 0.70, "60-79": 0.30, "80+": -1.30},
            "sex": {"Female": 0.35, "Male": -0.35},
            "smoking": {"Non-Smoker": 0.25, "Current Smoker": -0.75},
            # Narrow diabetes gradient: the South Asian Paradox.
            "diabetes": {"Normoglycemic": -0.45, "Diabetic": 1.15},
            "bp": {"Normal": -1.15, "Elevated": 0.10,
                   "Stage 1 HTN": 0.95, "Stage 2 HTN": 1.85},
            "cholesterol": {"Desirable": -0.85, "Borderline High": 0.05,
                            "High": 0.60, "Very High": 1.25},
        },
    },
    "reference-global-cohort.csv": {
        "seed": 20240518,
        "target_mean_bmi": 27.6,
        "target_obesity_rate": 0.27,
        "effects": {
            "age": {"20-39": -1.20, "40-59": 0.80, "60-79": 0.40, "80+": -1.50},
            "sex": {"Female": 0.25, "Male": -0.25},
            "smoking": {"Non-Smoker": 0.30, "Current Smoker": -1.00},
            # Wider diabetes gradient: obesity-driven metabolic risk.
            "diabetes": {"Normoglycemic": -0.50, "Diabetic": 2.60},
            "bp": {"Normal": -1.40, "Elevated": 0.15,
                   "Stage 1 HTN": 1.15, "Stage 2 HTN": 2.20},
            "cholesterol": {"Desirable": -1.00, "Borderline High": 0.05,
                            "High": 0.70, "Very High": 1.45},
        },
    },
}


def offsets(rows, effects):
    """Per-row sum of risk-factor effects on mean BMI."""
    return [sum(effects[col][row[col]] for col in effects) for row in rows]


def realize(offs, noise, base, sigma):
    """Final rounded BMI values for a candidate (base, sigma)."""
    shift = sigma * sigma / 2.0
    out = []
    for off, z in zip(offs, noise):
        raw = (base + off) * math.exp(sigma * z - shift)
        raw = min(max(raw, BMI_MIN), BMI_MAX)
        bmi = round(raw, 1)
        # Keep the threshold unambiguous: never land exactly on 30.0.
        if bmi == OBESITY_THRESHOLD:
            bmi = 30.1 if raw >= OBESITY_THRESHOLD else 29.9
        out.append(bmi)
    return out


def solve_base(offs, noise, sigma, target_mean):
    """Bisect on `base` so the realized mean BMI matches the target."""
    lo, hi = 5.0, 60.0
    for _ in range(80):
        mid = (lo + hi) / 2.0
        vals = realize(offs, noise, mid, sigma)
        if sum(vals) / len(vals) < target_mean:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def solve(offs, noise, target_mean, target_rate):
    """Solve for (base, sigma) hitting both the mean and obesity-rate targets.

    At a fixed mean below the obesity threshold, the realized share above the
    threshold increases monotonically with sigma, so a bisection on sigma with
    an inner bisection on base converges on both targets at once.
    """
    lo, hi = 0.01, 1.0
    base = solve_base(offs, noise, hi, target_mean)
    for _ in range(60):
        mid = (lo + hi) / 2.0
        base = solve_base(offs, noise, mid, target_mean)
        vals = realize(offs, noise, base, mid)
        rate = sum(v >= OBESITY_THRESHOLD for v in vals) / len(vals)
        if rate < target_rate:
            lo = mid
        else:
            hi = mid
    sigma = (lo + hi) / 2.0
    return solve_base(offs, noise, sigma, target_mean), sigma


NEW_COLUMNS = ("bmi", "obesity")


def strip_existing(header, raw_lines):
    """Drop previously generated columns so re-runs are idempotent."""
    names = header.split(",")
    keep = {i for i, name in enumerate(names) if name not in NEW_COLUMNS}
    if len(keep) == len(names):
        return header, raw_lines

    def trim(line):
        return ",".join(f for i, f in enumerate(line.split(",")) if i in keep)

    return trim(header), [trim(line) for line in raw_lines]


def process(path, spec):
    with open(path, newline="") as fh:
        header = fh.readline().rstrip("\n")
        # Keep the original six columns byte-for-byte: only append to each line.
        raw_lines = [line.rstrip("\n") for line in fh if line.strip()]

    header, raw_lines = strip_existing(header, raw_lines)
    rows = list(csv.DictReader([header] + raw_lines))

    offs = offsets(rows, spec["effects"])
    rng = random.Random(spec["seed"])
    noise = [rng.gauss(0.0, 1.0) for _ in rows]

    base, sigma = solve(offs, noise,
                        spec["target_mean_bmi"], spec["target_obesity_rate"])
    values = realize(offs, noise, base, sigma)

    with open(path, "w", newline="") as fh:
        fh.write(header + ",bmi,obesity\n")
        for line, bmi in zip(raw_lines, values):
            label = "Obese" if bmi >= OBESITY_THRESHOLD else "Non-Obese"
            fh.write(f"{line},{bmi:.1f},{label}\n")

    n = len(values)
    obese = sum(v >= OBESITY_THRESHOLD for v in values)
    print(f"{path}: n={n} base={base:.3f} sigma={sigma:.4f} "
          f"mean={sum(values)/n:.2f} obese={obese/n:.2%} "
          f"range={min(values):.1f}-{max(values):.1f}")


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for name, spec in COHORTS.items():
        process(os.path.join(root, name), spec)


if __name__ == "__main__":
    main()
