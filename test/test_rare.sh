#!/bin/bash
cd ./POP-TOOLS

trait=Head_BMD
# Single-variant test
/s/bin/python3.10 ./POP-RARE.py \
--gwas-yhat-unlab ./test/data/${trait}_yhat_unlab.RVAS.txt.gz \
--gwas-y-lab ./test/data/${trait}_y_lab.RVAS.txt.gz \
--gwas-yhat-lab ./test/data/${trait}_yhat_lab.RVAS.txt.gz \
--r 1 \
--out ./test/results/${trait}_r1

/s/bin/python3.10 ./POP-RARE.py \
--gwas-yhat-unlab ./test/data/${trait}_yhat_unlab.RVAS.txt.gz \
--gwas-y-lab ./test/data/${trait}_y_lab.RVAS.txt.gz \
--gwas-yhat-lab ./test/data/${trait}_yhat_lab.RVAS.txt.gz \
--out ./test/results/${trait}_no_ovp

/s/bin/python3.10 ./POP-RARE.py \
--gwas-yhat-unlab ./test/data/${trait}_yhat_unlab.RVAS.txt.gz \
--gwas-y-lab ./test/data/${trait}_y_lab.RVAS.txt.gz \
--gwas-yhat-lab ./test/data/${trait}_yhat_lab.RVAS.txt.gz \
--sample-overlap \
--out ./test/results/${trait}_ovp

# Burden test
/s/bin/python3.10 ./POP-RARE.py \
--gwas-yhat-unlab ./test/data/${trait}_yhat_unlab.Burden.txt.gz \
--gwas-y-lab ./test/data/${trait}_y_lab.Burden.txt.gz \
--gwas-yhat-lab ./test/data/${trait}_yhat_lab.Burden.txt.gz \
--burden \
--r 1 \
--out ./test/results/${trait}_burden_r1

/s/bin/python3.10 ./POP-RARE.py \
--gwas-yhat-unlab ./test/data/${trait}_yhat_unlab.Burden.txt.gz \
--gwas-y-lab ./test/data/${trait}_y_lab.Burden.txt.gz \
--gwas-yhat-lab ./test/data/${trait}_yhat_lab.Burden.txt.gz \
--burden \
--out ./test/results/${trait}_burden_no_ovp

/s/bin/python3.10 ./POP-RARE.py \
--gwas-yhat-unlab ./test/data/${trait}_yhat_unlab.Burden.txt.gz \
--gwas-y-lab ./test/data/${trait}_y_lab.Burden.txt.gz \
--gwas-yhat-lab ./test/data/${trait}_yhat_lab.Burden.txt.gz \
--burden \
--ovp \
--out ./test/results/${trait}_burden_ovp