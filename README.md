<img width="1460" alt="Screenshot 2025-01-13 at 7 20 00 PM" src="https://github.com/user-attachments/assets/cd1cee84-74af-4a35-9f28-6c70a4eedc14" />


`POP-TOOLS` (**PO**st-**P**rediction **TOOLS**) is a Python3-based command line toolkit for conducting valid and powerful machine learning (ML)-assisted genetic association studies. 

The `POP-TOOLS` toolkit can be used to conduct
* **ML-Assisted Genome-Wide Association Studies (GWAS)**.
* **ML-Assisted Rare-Variant Association Studies (RVAS)**.

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

## Power and sample size calculator

We provide a [web interface](https://jmiao24.shinyapps.io/pop-gwas/) for the power and sample size calculation for ML-assisted GWAS.

## Version History
[Version 1.2.0] (Nov 29, 2024): Added POP-GWAS for rare-variant association studies (single-variant and Burden test).

[Version 1.1.0] (May 1, 2024): Added quality control to remove SNPs with duplicate IDs; Added a version of the sample overlap correction; Modified scripts to accommodate the latest version of Polaris.
  
[Version 1.0.0] (Jan 2, 2024): Initial release.

## Reference

If you use `POP-GWAS`, please cite

Miao, J., Wu, Y., Sun, Z. et al. Valid inference for machine learning-assisted genome-wide association studies. Nat Genet (2024). https://doi.org/10.1038/s41588-024-01934-0

## Contact

For questions and comments, please open a GitHub issue (preferred) or contact Jiacheng Miao at jiacheng.miao@wisc.edu or Qiongshi Lu at qlu@biostat.wisc.edu.
