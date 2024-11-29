#!/bin/bash
cd ./POP-TOOLS

trait=Head_BMD

/s/bin/python3.10 ./POP-GWAS.py \
--gwas-yhat-unlab ./test/data/${trait}_yhat_unlab.txt.gz \
--gwas-y-lab ./test/data/${trait}_y_lab.txt.gz \
--gwas-yhat-lab ./test/data/${trait}_yhat_lab.txt.gz \
--r 1 \
--out ./test/results/${trait}

/s/bin/python3.10 ./POP-GWAS.py \
--gwas-yhat-unlab ./test/data/${trait}_yhat_unlab.txt.gz \
--gwas-y-lab ./test/data/${trait}_y_lab.txt.gz \
--gwas-yhat-lab ./test/data/${trait}_yhat_lab.txt.gz \
--out ./test/results/${trait}

/s/bin/python3.10 ./POP-GWAS.py \
--gwas-yhat-unlab ./test/data/${trait}_yhat_unlab.txt.gz \
--gwas-y-lab ./test/data/${trait}_y_lab.txt.gz \
--gwas-yhat-lab ./test/data/${trait}_yhat_lab.txt.gz \
--sample-overlap \
--out ./test/results/${trait}

cd ./POP-TOOLS
trait=T2D

/s/bin/python3.10 ./POP-GWAS.py \
--gwas-yhat-unlab ./test/data/${trait}_yhat_unlab.txt.gz \
--gwas-y-lab ./test/data/${trait}_y_lab.txt.gz \
--gwas-yhat-lab ./test/data/${trait}_yhat_lab.txt.gz \
--r 1 \
--bt \
--out ./test/results/${trait}_ovp

/s/bin/python3.10 ./POP-GWAS.py \
--gwas-yhat-unlab ./test/data/${trait}_yhat_unlab.txt.gz \
--gwas-y-lab ./test/data/${trait}_y_lab.txt.gz \
--gwas-yhat-lab ./test/data/${trait}_yhat_lab.txt.gz \
--bt \
--out ./test/results/${trait}_ovp

/s/bin/python3.10 ./POP-GWAS.py \
--gwas-yhat-unlab ./test/data/${trait}_yhat_unlab.txt.gz \
--gwas-y-lab ./test/data/${trait}_y_lab.txt.gz \
--gwas-yhat-lab ./test/data/${trait}_yhat_lab.txt.gz \
--sample-overlap \
--bt \
--out ./test/results/${trait}_ovp
