#!/usr/bin/env python3
import argparse, time, utils
from compute import estimate_popgwas, estimate_poprare, estimate_popburden

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
        help="Path to rare-variant associataion analysis summary statistics file")       
    parser.add_argument("--gwas-y-lab", action="store", default=None, type=str, required=True,
        help="Path to rare-variant associataion analysis summary statistics")
    parser.add_argument("--gwas-yhat-lab", action="store", default=None, type=str, required=True,
        help="Path to rare-variant associataion analysis file of imputed phenotype in labeled data") 
    parser.add_argument("--burden", dest="burden", action="store_true", default=False,
        help="Whether the input is burden test summary statistics") 
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
        
        if args.burden:
            out_prefix = args.out + "_POP-Burden"
        else:
            out_prefix = args.out + "_POP-Rare-variant"
        log = utils.Logger(out_prefix + '.log')

        header = TopHEAD
        header += 'Call:\n'
        header += './POP-RARE.py'
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
        
        if args.burden:
            log.log("\n### Reading the Burden test summary statistics ###\n")
            # parse summary statistics -> Z table
            z_df, n_col, N1_col, N2_col = utils.read_z_gene(args=args, log=log)
        else:
            log.log("\n### Reading the RVAS summary statistics ###\n")
            # parse summary statistics -> Z table
            z_df, n_col, N1_col, n_case_col, N2_col, eaf_col = utils.read_z(args=args, rare = True, log=log)
        
        # calculate
        if args.burden:
            if args.ovp:
                log.log("\n### Obatining the POP-GWAS results ###\n")
                df = estimate_popburden(z_df=z_df, n_col=n_col, N1_col=N1_col, N2_col=N2_col, ovp = True, log = log)
                log.log("\n--- Finish\n")
            else:
                log.log("\n### Obatining the POP-GWAS results ###\n")
                df = estimate_popburden(z_df=z_df, n_col=n_col, N1_col=N1_col, N2_col=N2_col, ovp = False, log = log, r12=args.r)
                log.log("\n--- Finish\n")
        else:
            if args.ovp:
                log.log("\n### Obatining the POP-RARE results with sample overlap correction ###")
                df = estimate_poprare(z_df=z_df, n_col=n_col, N1_col=N1_col, n_case_col=n_case_col, N2_col=N2_col, eaf_col=eaf_col, bt=args.bt, ovp = True, log = log)
                log.log("\n--- Finish")
            else:
                log.log("\n### Obatining the POP-RARE results ###\n")
                df = estimate_poprare(z_df=z_df, n_col=n_col, N1_col=N1_col, n_case_col=n_case_col, N2_col=N2_col, eaf_col=eaf_col, bt=args.bt, ovp = False, log = log, r12=args.r)
                log.log("\n--- Finish\n")


        log.log("\n### Writing the POP-RARE results ###\n")
        if args.burden:
            out_fh = utils.save_output_gene(df=df, out_prefix=out_prefix)
        else:
            out_fh = utils.save_output(df=df, out_prefix=out_prefix)
        log.log(f"--- POP-RARE summary statistics has been written to {out_prefix}.txt\n")

        log.log('--- Analysis finished at {T} ---'.format(T=time.ctime()))
        log.log('--- Total time elapsed: {T} ---'.format(T=utils.sec_to_str(time.time()-start_time)))
        
    except Exception as e:
        log.log(f"Error during operation: {e}")
        raise
        
if __name__ == "__main__":
    main()
