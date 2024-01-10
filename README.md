# POP-TOOLS
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

## Power and sample size calculator

We provide a [web interface](https://jmiao24.shinyapps.io/pop-gwas/) for the power and sample size calculator for ML-assisted GWAS.

## Version History
* Jan 2, 2024: Initial release.

## Reference

[Valid inference for machine learning-assisted GWAS](https://www.medrxiv.org/content/10.1101/2024.01.03.24300779v1)

[Assumption-lean and Data-adaptive Post-Prediction Inference](https://arxiv.org/abs/2311.14220)

## Contact

For questions and comments, please open a GitHub issue (preferred) or contact Jiacheng Miao at jiacheng.miao@wisc.edu or Yixuan Wu at wu638@wisc.edu or Qiongshi Lu at qlu@biostat.wisc.edu.

## "POP" familial links
* [POPInf](https://github.com/qlu-lab/POP-TOOLS) (**PO**st-**P**rediction **Inf**erence) is a generic toolkit for conducting valid and powerful post-prediction inference. It is more general (can be applied to a wider range of statistical quantities), but is not optimized for genetic applications.
