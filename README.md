<h1 align="center">
<p> POP-TOOLS
</h1>

`POP-TOOLS` (**PO**st-**P**rediction **TOOLS**) is a Python3-based command line toolkit for conducting valid and powerful machine learning (ML)-assisted genetic association studies. 

The `POP-TOOLS` toolkit can be used to conduct
* `POP-GWAS` [**PO**st-**P**rediction Genome-wide Association Studies (**GWAS**)]

## Manual

The `POP-TOOLS` and its required modules can be installed via 

```
git clone https://github.com/qlu-lab/POP-TOOLS
cd POP-TOOLS
pip install -r requirements.txt
```

Please see the [TL;DR](https://github.com/qlu-lab/POP-TOOLS/wiki/1.-POP%E2%80%90GWAS#tldr) to conduct `POP-GWAS`

Please see the [wiki](https://github.com/qlu-lab/POP-TOOLS/wiki) for tutorials describing the basic function and along with detailed manual of `POP-TOOLS`. 

Please see the [FAQ](https://github.com/qlu-lab/POP-TOOLS/wiki/FAQ) for frequently asked questions related to `POP-TOOLS`.

## Power and sample size calculator

We provide a [web interface](https://jmiao24.shinyapps.io/pop-gwas/) for the power and sample size calculation for ML-assisted GWAS.

## Version History
[Version 1.1.0] (May 1, 2024): Added quality control to remove SNPs with duplicate IDs; Added a version of the sample overlap correction; Modified scipts to accommodate the latest version of polars.
  
[Version 1.0.0] (Jan 2, 2024): Initial release.

## Reference

[Valid inference for machine learning-assisted genome-wide association studies](https://www.nature.com/articles/s41588-024-01934-0)

[Assumption-Lean and Data-adaptive Post-Prediction Inference](https://arxiv.org/abs/2311.14220)

[Task-agnostic machine-learning-assisted inference](https://arxiv.org/abs/2405.20039)

## Contact

For questions and comments, please open a GitHub issue (preferred) or contact Jiacheng Miao at jiacheng.miao@wisc.edu or Qiongshi Lu at qlu@biostat.wisc.edu.

## "POP" familial links
* [POPInf](https://github.com/qlu-lab/POPinf) (**PO**st-**P**rediction **Inf**erence) is a generic toolkit for conducting valid and powerful post-prediction inference. It is more general (can be applied to a wider range of statistical quantities), but is not optimized for genetic applications.
