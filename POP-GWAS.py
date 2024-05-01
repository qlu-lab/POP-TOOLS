#!/usr/bin/env python3
import argparse, time, utils
from compute import estimate_popgwas

TopHEAD = "*********************************************************************\n"
TopHEAD += "* POst-Prediction Genome-Wide Association Studies (POP-GWAS) \n"
TopHEAD += "* Version 1.1.0 \n"
TopHEAD += "* (C) Jiacheng Miao and Yixuan Wu \n"
TopHEAD += "* University of Wisconsin-Madison \n"
TopHEAD += "* https://github.com/qlu-lab/POP-TOOLS \n"
TopHEAD += "* GNU General Public License v3\n"
TopHEAD += "*********************************************************************\n"

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gwas-yhat-unlab", action="store", default=None, type=str, required=True,
        help="Path to GWAS summary statistics file of imputed phenotype in unlabeled data (at least include columns: CHR, BP, SNP, A1, A2, EAF, Z, P and N / [N_case, N_control])")       
    parser.add_argument("--gwas-y-lab", action="store", default=None, type=str, required=True,
        help="Path to GWAS summary statistics file of observed phenotype in labeled data (at least include columns: CHR, BP, SNP, A1, A2, Z, P and N / [N_case, N_control])")
    parser.add_argument("--gwas-yhat-lab", action="store", default=None, type=str, required=True,
        help="Path to GWAS summary statistics file of imputed phenotype in labeled data (at least include columns: CHR, BP, SNP, A1, A2, Z, and P)") 
    parser.add_argument("--out", action="store", default='POP-TOOLS', type=str, required=True,
        help="The prefix of path to output summary statistics")
    parser.add_argument("--bt", "--binary-trait", dest="bt", action="store_true", default=False,
        help="Whether the trait is binary or not")
    parser.add_argument("--ovp", "--sample-overlap", dest="ovp", action="store_true", default=False,
        help="Whether to correct for sample overlap or not")
    parser.add_argument("--r", action="store", default=None, type=float,
        help="The imputation quality (correlation between observed and imputed phenotype after adjusting for GWAS covariates) in labeled dataset")  
    return parser.parse_args()

def main():

    try:
        # parse input argumemts
        args = parse_args()
        if args.out is None:
            raise ValueError('--out is required.')
        
        out_prefix = args.out + "_POP-GWAS"
        log = utils.Logger(out_prefix + '.log')

        header = TopHEAD
        header += 'Call:\n'
        header += './POP-GWAS.py'
        for arg, value in vars(args).items():
            if value:
                if isinstance(value, bool):
                    header += f" \\\n--{arg.replace('_','-')}"
                else:
                    header += f" \\\n--{arg.replace('_','-')} {value}"
        header += '\n'
        log.log(header)
        log.log(f'--- Analysis began at {time.ctime()} ---\n')
        start_time = time.time()

        if not args.ovp:
            r12 = utils.extract_single_r_from_ldsc(args=args, log=log) if args.r is None else args.r
            log.log(f"--- The imputation quality r = {r12}")
        else:
            r12,r13,r23 = utils.extract_multi_r_from_ldsc(args=args, log=log)
            log.log(f"\n--- The correlation between GWAS estimates on y and yhat in labeled data = {r12}")
            log.log(f"--- The correlation between GWAS estimates on y in labeled data and yhat in unlabeled data = {r13}")
            log.log(f"--- The correlation between GWAS estimates on yhat in labeled and unlabeled data = {r23}")
        
        log.log("\n### Reading the GWAS summary statistics ###\n")
        # parse summary statistics -> Z table
        z_df, n_col, N1_col, n_case_col, N2_col, eaf_col = utils.read_z(args=args, log=log)
        
        # calculate
        if not args.ovp:
            log.log("\n### Obatining the POP-GWAS results ###\n")
            df = estimate_popgwas(z_df=z_df, n_col=n_col, N1_col=N1_col, n_case_col=n_case_col, N2_col=N2_col, eaf_col=eaf_col, bt=args.bt, r12=r12, r13 = 0, r23 = 0)
            log.log("\n--- Finish\n")
        else:
            log.log("\n### Obatining the POP-GWAS results with sample overlap correction ###")
            df = estimate_popgwas(z_df=z_df, n_col=n_col, N1_col=N1_col, n_case_col=n_case_col, N2_col=N2_col, eaf_col=eaf_col, bt=args.bt, r12=r12, r13 = r13, r23 = r23)
            log.log("\n--- Finish")


        log.log("\n### Writing the POP-GWAS results ###\n")
        out_fh = utils.save_output(df=df, out_prefix=out_prefix)
        log.log(f"--- POP-GWAS summary statistics has been written to {out_prefix}.txt\n")

        log.log('--- Analysis finished at {T} ---'.format(T=time.ctime()))
        log.log('--- Total time elapsed: {T} ---'.format(T=utils.sec_to_str(time.time()-start_time)))
        
    except Exception as e:
        log.log(f"Error during operation: {e}")
        raise
        
if __name__ == "__main__":
    main()
