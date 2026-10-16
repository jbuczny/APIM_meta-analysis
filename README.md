# Meta-analysis of APIM paths across dyadic studies

Pools the paths of an Actor–Partner Interdependence Mediation Model (APIMeM)
across dyadic studies, in two parallel components: an ordinary meta-analysis
of the observed paths, and one in which every path is corrected for the
unreliability of the questionnaires used to measure X, M and Y in both
partners.

The correction is applied to each study's model rather than to the pooled
coefficients. APIM paths are partial regression coefficients, so the
classical disattenuation formula `r / sqrt(rxx * ryy)` does not apply to
them: unreliability in one predictor also biases the coefficients of the
others. Each study is therefore refitted as a single-indicator latent
APIMeM, with the error variance of every scale score fixed at
`(1 − reliability) × variance`. That is equivalent to disattenuating the
6 × 6 correlation matrix and refitting the saturated model, and `lavaan`
returns the correct sampling covariances of the corrected paths.

## What the script does

`apim_meta_analysis.R` runs the whole pipeline:

1. **Per study** — builds scale scores from the item columns and estimates
   each scale's reliability (McDonald's omega total from a one-factor
   congeneric CFA, or Cronbach's alpha).
2. **Per study** — fits the APIMeM in `lavaan` twice, uncorrected and
   corrected, and extracts the 12 standardized paths together with their
   full 12 × 12 sampling covariance matrix (delta method), under normal-theory
   (ML) and robust (MLR sandwich) estimation. An optional dyad bootstrap also
   propagates the sampling error of the reliability coefficients themselves.
3. **Pooling** — `metafor::rma.mv` with the block-diagonal sampling
   covariance matrix, so the dependence among a study's twelve paths is
   carried into the pooling rather than assumed away. The multilevel model
   (partner-specific estimates nested within studies) is primary; a
   common-effect model runs as a sensitivity check. Cluster-robust CR2
   inference with Satterthwaite degrees of freedom comes from `clubSandwich`.
4. **Indirect and total effects** — actor–actor, partner–partner,
   actor–partner and partner–actor products, plus totals, with Monte Carlo
   confidence intervals drawn from the joint sampling distribution of the
   pooled paths (Ledermann, Macho & Kenny, 2011).
5. **Reporting** — APA-formatted tables written to a single `.docx` via
   `flextable` and `officer`.

### Path notation

Partner 1 and Partner 2; `12` means "from Partner 1 to Partner 2".

| | actor | partner |
|---|---|---|
| **a** (X → M) | `a1`, `a2` | `a12`, `a21` |
| **b** (M → Y) | `b1`, `b2` | `b12`, `b21` |
| **c′** (X → Y, direct) | `c1`, `c2` | `c12`, `c21` |

## Running it

```r
source("apim_meta_analysis.R")
```

R ≥ 4.2 with `lavaan`, `metafor`, `clubSandwich`, `MASS`, `flextable` and
`officer`.

Settings live in the `CFG` list at the top of the script: the study files,
the item-naming pattern, omega versus alpha, whether to pool the two
partners' estimates of each effect, which sampling covariance matrix and
pooling model are primary, and whether to run the bootstrap. To change them
without editing the file, define `CFG_OVERRIDE` first:

```r
CFG_OVERRIDE <- list(reliability = "alpha", primary_V = "MLR")
source("apim_meta_analysis.R")
```

### Data format

One file per study, wide, one row per dyad. Items are named
`<construct>_<member>_<item>` — `x_p1_1`, `m_p2_3`, `y_p1_4` and so on. A
study that provides only scale scores names them `<construct>_<member>`
(`x_p1`, `m_p2`) and supplies its reliabilities through
`CFG$reliability_supplied`.

## Files

| File | What it is |
|---|---|
| `apim_meta_analysis.R` | the analysis pipeline |
| `simulate_toy_data.R` | regenerates the toy datasets from known population paths |
| `check_recovery.R` | large-sample check that the corrected model recovers the true paths |
| `toy_study1.csv`, `toy_study2.csv` | toy data: 130 and 174 dyads, 24 Likert items each |
| `toy_population_values.csv` | the population paths behind the toy data |
| `APIM_meta_analysis_tables.docx` | APA tables: descriptives, pooled direct effects, indirect and total effects, heterogeneity, sensitivity |
| `APIM_meta_analysis_results.csv` | long-format results for every component, sampling covariance type and model |
| `APIM_meta_analysis_results.rds` | full results object, including each study's paths and sampling covariances |

## Toy data

The two datasets are simulated, not real dyadic data. They are generated
from a known APIMeM with fixed population paths, so the pipeline can be run
and verified end to end by anyone. Study 2 is given less reliable scales
(omega ≈ .60–.68) than Study 1 (omega ≈ .83–.87), which makes the effect of
the reliability correction visible: the pooled actor *a* path rises from
about .25 to about .34, and the actor–actor indirect effect from about .072
to about .130.

`check_recovery.R` verifies the correction against the truth at N = 20,000:
the corrected model recovers the population paths to within about .02, while
the uncorrected paths stay attenuated.

## A caution about two studies

With only two studies the between-study variance is barely identified, so
the multilevel and common-effect models give almost identical answers, and
the cluster-robust CR2 tests run on roughly one degree of freedom. The
script therefore computes CR2 results but keeps them out of the tables
unless every Satterthwaite *df* reaches 4 (Tipton, 2015). Read the robust
intervals as descriptive until more studies are available.

## Related

The companion repository
[APIM_meta-analysis_sensitivity](https://github.com/jbuczny/APIM_meta-analysis_sensitivity)
asks how large the pooled paths would have to be to reach significance, and
takes `APIM_meta_analysis_results.rds` from this repository as its input.

## References

- Ledermann, T., Macho, S., & Kenny, D. A. (2011). Assessing mediation in
  dyadic data using the actor–partner interdependence model.
  *Structural Equation Modeling, 18*, 595–612.
- Pustejovsky, J. E., & Tipton, E. (2021). Meta-analysis with robust variance
  estimation. *Prevention Science, 23*, 425–438.
- Tipton, E. (2015). Small sample adjustments for robust variance estimation
  with meta-regression. *Psychological Methods, 20*, 375–393.
- Cheung, M. W.-L. (2014). Modeling dependent effect sizes with three-level
  meta-analyses. *Psychological Methods, 19*, 211–229.

## License

Creative Commons Attribution-NonCommercial 4.0 International — see
`LICENSE.md`.
