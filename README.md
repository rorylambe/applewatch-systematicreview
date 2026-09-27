# The Accuracy of Apple Watch Measurements: A Living Systematic Review and Meta-Analysis

Analysis code accompanying:

> Lambe, R. et al. The accuracy of Apple Watch measurements: a living systematic review and meta-analysis. *npj Digital Medicine* **9**, 63 (2026). [https://doi.org/10.1038/s41746-025-02238-1](https://doi.org/10.1038/s41746-025-02238-1)

[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41746--025--02238--1-blue)](https://doi.org/10.1038/s41746-025-02238-1)
[![PROSPERO](https://img.shields.io/badge/PROSPERO-CRD42023481841-green)](https://www.crd.york.ac.uk/PROSPERO/view/CRD42023481841)
[![OSF](https://img.shields.io/badge/OSF-osf.io%2Fv5d3k-lightgrey)](https://osf.io/v5d3k)

## About the review

This living systematic review and meta-analysis evaluated the agreement between Apple Watch health metrics and criterion measures. Nine databases were searched from inception to 24 September 2025. The review included 82 studies assessing 14 health metrics (430,052 participants), across all Apple Watch models through Series 9 and Ultra 2.

The review is designed as a **living** synthesis. Searches will be updated every 12 months and updates will be shared via the [Open Science Framework](https://osf.io/v5d3k).

## Meta-analysis results

| Metric | Studies (n) | Pooled result |
|---|---|---|
| Heart rate (all conditions) | 22 (1,247) | Mean bias −0.27 bpm (95% CI −0.72 to 0.17); LoA −7.19 to 6.64 bpm |
| Atrial fibrillation detection (ECG app) | 11 (3,144) | Sensitivity 0.79 (95% CI 0.61–0.90); specificity 0.91 (95% CI 0.81–0.96); AUC 0.93 |
| Blood oxygen saturation (SpO₂) | 9 (969) | Mean bias −0.04% (95% CI −0.42 to 0.35); LoA −4.01 to 3.94% |

See the paper for subgroup (rest vs exercise, optical sensor generation, hypoxic range) and sensitivity analyses, and for the narrative synthesis of the remaining metrics.

## Repository contents

| File | Description |
|---|---|
| `hr_tipton&shuster_meta_analysis.r` | Heart rate: pooled mean bias and population limits of agreement using the Tipton & Shuster framework. |
| `spo2_tipton&shuster_meta_analysis.r` | Blood oxygen saturation: as above, with an option to restrict to hypoxic ranges. |
| `afib_meta_analysis.R` | Atrial fibrillation detection: bivariate meta-analysis (Reitsma model, `mada` package) with SROC curve. |
| `tipton_shus_example_MA.R` | Worked example of the Tipton & Shuster limits of agreement meta-analysis. |
| `study_characteristics.xlsx` | Characteristics of included studies. |

## Statistical methods

**Heart rate and SpO₂.** These were meta-analysed with the Tipton & Shuster (2017) framework for Bland–Altman studies. The method uses random-effects models with inverse-variance weighting. Mean bias and the log of the bias-adjusted variance of differences are pooled separately, then combined to give population limits of agreement that account for between-study heterogeneity (τ²). The scripts report 95% CIs for the outer limits using both model-based and robust variance estimation. Only one estimate per study per condition was included, to avoid unit-of-analysis errors.

**Atrial fibrillation detection.** Pooled sensitivity and specificity were estimated using the bivariate random-effects model of Reitsma et al. (2005), implemented in the `mada` package, and a summary ROC curve was produced.

**Sensitivity and subgroup analyses.** Each script contains commented-out blocks that can be enabled to rerun the analysis:
- excluding studies at high risk of bias;
- by optical heart rate sensor generation (first: up to Series 3; second: Series 4–5 and SE; third: Series 6 onwards, including Ultra);
- by population (healthy vs clinical), for heart rate;
- restricted to hypoxic ranges, for SpO₂.

## Requirements

Analyses were run in R 4.5.1.

```r
install.packages(c("readr", "dplyr", "ggplot2", "mada", "meta", "metafor"))
```

## Data availability

The synthesised results data, risk of bias assessments, and study protocol are available on the [Open Science Framework (osf.io/v5d3k)](https://osf.io/v5d3k). The raw datasets are available from the corresponding author on request.

## Registration

The protocol was prospectively registered with PROSPERO: [CRD42023481841](https://www.crd.york.ac.uk/PROSPERO/view/CRD42023481841).

## Citation

```bibtex
@article{lambe2026applewatch,
  title   = {The accuracy of Apple Watch measurements: a living systematic review and meta-analysis},
  author  = {Lambe, Rory and others},
  journal = {npj Digital Medicine},
  volume  = {9},
  pages   = {63},
  year    = {2026},
  doi     = {10.1038/s41746-025-02238-1}
}
```

## Key methodological references

- Tipton E, Shuster J. A framework for the meta-analysis of Bland–Altman studies based on a limits of agreement approach. *Stat Med*. 2017;36:3621–3635.
- Reitsma JB, et al. Bivariate analysis of sensitivity and specificity produces informative summary measures in diagnostic reviews. *J Clin Epidemiol*. 2005;58:982–990.
- Bland JM, Altman DG. Statistical methods for assessing agreement between two methods of clinical measurement. *Lancet*. 1986;327:307–310.

## Contact

Rory Lambe, University College Dublin.
To suggest a study for the next update, please email me.