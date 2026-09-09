<img width="1460" alt="Screenshot 2025-01-13 at 7 20 00 PM" src="https://github.com/user-attachments/assets/cd1cee84-74af-4a35-9f28-6c70a4eedc14" />


`POP-TOOLS` (**PO**st-**P**rediction **TOOLS**) is a Python3-based command line toolkit for conducting valid and powerful machine learning (ML)-assisted genetic association studies. 

The `POP-TOOLS` toolkit can be used to conduct
* **ML-Assisted Genome-Wide Association Studies (GWAS)**.
* **ML-Assisted Rare-Variant Association Studies (RVAS) and gene-based test**.

## Manual

The `POP-TOOLS` and its required modules can be installed via 

```
git clone https://github.com/qlu-lab/POP-TOOLS
cd POP-TOOLS
pip install -r requirements.txt
```

Please see the [TL;DR](https://github.com/qlu-lab/POP-TOOLS/wiki/1.-POP%E2%80%90GWAS#tldr) for **ML-assisted GWAS** using `POP-GWAS`

Please see the [TL;DR](https://github.com/qlu-lab/POP-TOOLS/wiki/2.-POP%E2%80%90GWAS-for-Rare%E2%80%90Variant-Association-Studies#tldr) for **ML-assisted Rare-Variant Association Studies** (single-variant and gene-level burden test) using `POP-GWAS`.

Please see the [wiki](https://github.com/qlu-lab/POP-TOOLS/wiki) for tutorials describing the basic function along with a detailed manual of `POP-TOOLS`. 

Please see the [FAQ](https://github.com/qlu-lab/POP-TOOLS/wiki/FAQ) for frequently asked questions related to `POP-TOOLS`.

## Binary-trait gene-level burden tests

Use `POP-RARE.py --burden` with three tab-separated inputs containing
`GENE`, `MASK`, signed `Z`, and total sample count `N`. This input format
works for quantitative and binary traits; `--bt` is optional when all inputs
already contain `N`. The upstream observed-trait burden test must itself be
appropriate for a binary phenotype. Use the same burden definition and effect
direction across all three inputs, and analyze one mask at a time.

With `--burden --bt`, each input may instead provide `N_case` and `N_control`;
the reader calculates `N = N_case + N_control`. Predicted-phenotype inputs may
still provide `N` alone. If both formats are supplied, their totals must agree.
Counts must be positive finite integers. `N` is the number of analyzed samples,
not the number of carriers or a case/control effective sample size.

```bash
python3 POP-RARE.py \
  --gwas-y-lab observed.Burden.tsv \
  --gwas-yhat-lab predicted_labeled.Burden.tsv \
  --gwas-yhat-unlab predicted_unlabeled.Burden.tsv \
  --burden --bt \
  --out binary_trait
```

`--bt` does not rescale burden test statistics. `Z` and `P` are the association
outputs; their validity depends on calibrated input statistics and the method's
sampling/correlation assumptions. This option does not correct an uncalibrated
upstream test or severe case/control imbalance.

For backward compatibility, burden output retains `BETA` and `SE`, and adds
`EFFECT_SCALE=Z_OVER_SQRT_N`. These are working-scale quantities constructed
from signed Z scores and sample counts, not raw 0/1 phenotype coefficients,
per-allele burden effects, or log odds ratios. `N_eff` is the existing working-scale
precision measure, not an inferred count of cases and controls.

Do **not** apply `_qt_to_bt()` to these burden outputs: case fraction alone
does not recover an OR. For example, multiplying every burden score by a
constant preserves its Z statistic and N but changes its per-unit regression
coefficient. Recovering a meaningful OR requires additional upstream effect
scale/model information and a justified conversion. No burden `OR` or `LOG_OR`
is reported; the helper now rejects gene-level output to prevent this misuse.

Run the regression tests with `python -m pytest test/test_burden.py`
(install `pytest` in addition to `requirements.txt`).

## Power and sample size calculator

We provide a [web interface](https://jmiao24.shinyapps.io/pop-gwas/) for the power and sample size calculation for ML-assisted GWAS.

## Version History
[Version 1.2.0] (Nov 29, 2024): Added POP-GWAS for rare-variant association studies (single-variant and Burden test).

[Version 1.1.0] (May 1, 2024): Added quality control to remove SNPs with duplicate IDs; Added a version of the sample overlap correction; Modified scripts to accommodate the latest version of Polaris.
  
[Version 1.0.0] (Jan 2, 2024): Initial release.

## Reference

If you use `POP-GWAS`, please cite

Miao, J., Wu, Y., Sun, Z. et al. Valid inference for machine learning-assisted genome-wide association studies. Nat Genet, 56(11), 2361-2369 (2024). https://doi.org/10.1038/s41588-024-01934-0

## Contact

For questions and comments, please open a GitHub issue (preferred) or contact Jiacheng Miao at jiacheng.miao@wisc.edu or Qiongshi Lu at qlu@biostat.wisc.edu.
