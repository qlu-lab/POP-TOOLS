import os, sys, subprocess
import polars as pl
from functools import reduce
from argparse import ArgumentTypeError
from typing import List

_self_dir = os.path.dirname(os.path.abspath(__file__))
_self_exec =  sys.executable

const_dict = {
    'CHR': 'CHR', 'BP': 'BP', 'SNP': 'SNP', 'A1': 'A1', 'A2': 'A2', 'Z': 'Z',
    'Z_YHAT_LAB': 'Z_yhat_lab', 'Z_YHAT_UNLAB': 'Z_yhat_unlab', 'Z_Y_LAB': 'Z_y_lab',
    'Z_YHAT_LAB': 'Z_yhat_lab', 'Z_YHAT_UNLAB': 'Z_yhat_unlab', 'Z_Y_LAB': 'Z_y_lab',
    'N_LAB1': 'N_lab1', 'N_LAB2': 'N_lab2','N_LAB_CASE': 'N_lab_case', 'N_LAB_CONTROL': 'N_lab_control',
    'N_UNLAB': 'N_unlab', 'N_UNLAB_CASE': 'N_unlab_case', 'N_UNLAB_CONTROL': 'N_unlab_control',
    'N_EFF': 'N_eff', 'N_EFF_CASE': 'N_eff_case', 'N_EFF_CONTROL': 'N_eff_control',
    'BETA': 'BETA', 'OR': 'OR', 'SE': 'SE', 'P': 'P', 'EAF': 'EAF',
    'LD_score': f'{_self_dir}/ldsc/EUR_1KGphase3/LDscore/LDscore.',
    'LD_weights': f'{_self_dir}/ldsc/EUR_1KGphase3/weights/weights.hm3_noMHC.'
}

def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = reduce(lambda ll, b : divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=round(d))
    if h > 0:
        f += '{H}h:'.format(H=round(h))
    if m > 0:
        f += '{M}m:'.format(M=round(m))

    f += '{S}s'.format(S=round(s))
    return f

class Logger(object):
    '''
    Lightweight logging.
    '''
    def __init__(self, fh):
        self.log_fh = open(fh, 'w')

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.
        '''
        print(msg, file=self.log_fh)
        print(msg)

def extract_single_r_from_ldsc(args, log):
    log.log('### Computing r using LDSC ###\n')

    log.log(f"--- Parsing GWAS summary statistics on y in labeled data: {args.gwas_y_lab}")
    ss_y_lab = _munge_ss_for_ldsc(args.gwas_y_lab, f'{args.out}_y_lab', args.bt)

    log.log(f"--- Parsing GWAS summary statistics on yhat in labeled data: {args.gwas_yhat_lab}")
    ss_yhat_lab = _munge_ss_for_ldsc(args.gwas_yhat_lab, f'{args.out}_yhat_lab', args.bt)

    r = _extract_r_from_ldsc_log(_run_ldsc(y1=ss_y_lab, y2=ss_yhat_lab, out=args.out, binary=args.bt, log=log))
    for f in [f'{args.out}_y_lab.sumstats.gz', f'{args.out}_y_lab.log', f'{args.out}_yhat_lab.sumstats.gz', f'{args.out}_yhat_lab.log', f'{args.out}_ldsc.log']:
        os.remove(f)
    return r

def extract_multi_r_from_ldsc(args, log):
    log.log('### Computing pairwise r using LDSC ###\n')

    log.log(f"--- Parsing GWAS summary statistics on y in labeled data: {args.gwas_y_lab}")
    ss_y_lab = _munge_ss_for_ldsc(args.gwas_y_lab, f'{args.out}_y_lab', args.bt)

    log.log(f"--- Parsing GWAS summary statistics on yhat in labeled data: {args.gwas_yhat_lab}")
    ss_yhat_lab = _munge_ss_for_ldsc(args.gwas_yhat_lab, f'{args.out}_yhat_lab', args.bt)

    log.log(f"--- Parsing GWAS summary statistics on yhat in unlabeled data: {args.gwas_yhat_unlab}")
    ss_yhat_unlab = _munge_ss_for_ldsc(args.gwas_yhat_unlab, f'{args.out}_yhat_unlab', args.bt)

    r12 = _extract_r_from_ldsc_log(_run_ldsc(y1=ss_y_lab, y2=ss_yhat_lab, out=f'{args.out}_12', binary=args.bt, log=log))
    r13 = _extract_r_from_ldsc_log(_run_ldsc(y1=ss_y_lab, y2=ss_yhat_unlab, out=f'{args.out}_13', binary=args.bt, log=log))
    r23 = _extract_r_from_ldsc_log(_run_ldsc(y1=ss_yhat_lab, y2=ss_yhat_unlab, out=f'{args.out}_23', binary=args.bt, log=log))

    for f in [f'{args.out}_y_lab.sumstats.gz', f'{args.out}_y_lab.log', f'{args.out}_yhat_lab.sumstats.gz', f'{args.out}_yhat_lab.log', f'{args.out}_yhat_unlab.sumstats.gz', f'{args.out}_yhat_unlab.log', f'{args.out}_12_ldsc.log', f'{args.out}_13_ldsc.log', f'{args.out}_23_ldsc.log']:
        os.remove(f)
    return r12,r13,r23

def read_z(args, log):
    return _read_z(ss_in_yhat_unlab=args.gwas_yhat_unlab, ss_in_y_lab=args.gwas_y_lab, ss_in_yhat_lab=args.gwas_yhat_lab, binary=args.bt, log=log)

def save_output(df, out_prefix):
    # _format_out(df=df).collect().write_csv(file=out_prefix+'.txt.gz', include_header=True, separator="\t", quote_style="never", null_value='NA', float_precision=5)
    _format_out(df=df).sink_csv(path=out_prefix+'.txt', include_header=True, separator="\t", quote_style="never", null_value='NA', float_precision=5)


"""
Below are helper functions to upper utils
"""
def _munge_ss_for_ldsc(ss_fh, out, binary):
    subprocess.run(
        [
            _self_exec,
            f'{_self_dir}/ldsc/munge_sumstats.py',
            "--sumstats", ss_fh,
            "--out", out,
            "--chunksize", "500000",
            "--signed-sumstats", "Z,0",
            "--merge-alleles", f'{_self_dir}/ldsc/w_hm3.snplist'
        ] + (
            ["--N-cas-col", "N_case", "--N-con-col", "N_control"] if binary else []
        ),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

    return f'{out}.sumstats.gz'
    
def _run_ldsc(y1, y2, out, binary, log):

    subprocess.run(
        [
            _self_exec,
            f'{_self_dir}/ldsc/ldsc.py',
            "--rg", f"{y1},{y2}",
            "--ref-ld-chr", const_dict['LD_score'],
            "--w-ld-chr", const_dict['LD_weights'],
            "--out", f'{out}_ldsc'
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

    return f'{out}_ldsc.log'

def _extract_r_from_ldsc_log(fh):
    with open(fh, 'r') as f:
        ll = f.readlines()
    for i, l in enumerate(ll):
        if "gcov_int" in l:
            header = l.strip().split()
            values = ll[i + 1].strip().split()
            return float(values[header.index("gcov_int")])

def _read_ss(ss_fh, binary, unlab, ylab, log):
    if binary:
        tmp = pl.read_csv(ss_fh, has_header=True, separator="\t", try_parse_dates=False, null_values='NA', n_threads=1, n_rows=1).columns
        if "N_case" not in tmp or "N_control" not in tmp:
            binary = False
    old_cols = ["CHR", "BP", "SNP", "A1", "A2"] + (["EAF"] if ylab else []) + ["Z"] + (["N_case", "N_control"] if binary else ["N"])
    
    z = const_dict['Z_Y_LAB'] if ylab else const_dict['Z']
    new_cols = [const_dict['CHR'], const_dict['BP'], const_dict['SNP'], const_dict['A1'], const_dict['A2']] + ([const_dict['EAF']] if ylab else []) + [z]
    if unlab:
        if binary:
            new_cols.extend([const_dict['N_UNLAB_CASE'], const_dict['N_UNLAB_CONTROL']])
        else:
            new_cols.append(const_dict['N_UNLAB'])
    else:
        if binary:
            new_cols.extend([const_dict['N_LAB_CASE'], const_dict['N_LAB_CONTROL']])
        else:
            if ylab: # If y in labedled data
                new_cols.append(const_dict['N_LAB1'])        
            else:
                new_cols.append(const_dict['N_LAB2'])    
    if isinstance(ss_fh, str) and os.path.exists(ss_fh):
        log.log(f"--- Reading GWAS on {'y' if ylab else 'yhat'} in {'unlabeled' if unlab else 'labeled'} data: {ss_fh}")
        try: 
            ss = pl.read_csv(ss_fh, has_header=True, separator="\t", try_parse_dates=False, null_values='NA') \
                 .select(old_cols) \
                 .rename(dict(zip(old_cols,new_cols)))
        except ValueError as e:
            log.log(f"ValueError error occurred: {e}")
        except pl.exceptions.ComputeError as e:
            # log.log(f"ComputeError error occurred: {e}")
            ss = pl.read_csv(ss_fh, has_header=True, columns=old_cols, separator="\t", try_parse_dates=False, null_values='NA') \
                 .lazy() \
                 .select(old_cols) \
                 .rename(dict(zip(old_cols,new_cols)))
    else:
        raise FileNotFoundError(f"File not found or invalid input: {ss_fh}")
    
    if binary:
        if unlab:
            ss = ss.with_columns((pl.col(const_dict['N_UNLAB_CASE']) + pl.col(const_dict['N_UNLAB_CONTROL'])).alias(const_dict['N_UNLAB'])).drop([const_dict['N_UNLAB_CASE'], const_dict['N_UNLAB_CONTROL']])
        elif ylab:
            ss = ss.with_columns((pl.col(const_dict['N_LAB_CASE']) + pl.col(const_dict['N_LAB_CONTROL'])).alias(const_dict['N_LAB1'])).drop(const_dict['N_LAB_CONTROL'])
        else:
            ss = ss.with_columns((pl.col(const_dict['N_LAB_CASE']) + pl.col(const_dict['N_LAB_CONTROL'])).alias(const_dict['N_LAB2'])).drop(const_dict['N_LAB_CONTROL'], const_dict['N_LAB_CASE'])
    return ss.drop_nulls(subset=[z])

def _merge_match_a1a2(ss1, ss2):
    alleles_list = [const_dict['A1'], const_dict['A2'], const_dict['A1'] + 'x', const_dict['A2'] + 'x']

    ss = ss1.join(ss2, on=const_dict['SNP'], how='inner', suffix='x') \
         .select(~pl.selectors.by_name(const_dict['CHR'] + 'x', const_dict['BP'] + 'x')) \
         .alleles.filter_snps(alleles_list)
         
    if const_dict['Z_YHAT_UNLAB'] not in ss1.columns: # first merge between ss
        ss = ss.alleles.align_alleles_z(alleles_list, const_dict['Z']) \
             .rename({const_dict['Z']: const_dict['Z_YHAT_UNLAB']})
    else:
        ss = ss.alleles.align_alleles_z(alleles_list, const_dict['Z']) \
             .rename({const_dict['Z']: const_dict['Z_YHAT_LAB']})

    return ss.drop([const_dict['A1'] + 'x', const_dict['A2'] + 'x'])

def _read_z(ss_in_yhat_unlab, ss_in_y_lab, ss_in_yhat_lab, binary, log):
    
    @pl.api.register_dataframe_namespace("alleles")
    class AllelesOperations:
        def __init__(self, lf: pl.DataFrame):
            self._lf = lf

        def _invalid_snps(self, alleles_cols: List[str]) -> pl.Expr:
            a1, a2, a1x, a2x = alleles_cols
            con0 = pl.col(a1).is_in(['A', 'T', 'G', 'C']) & pl.col(a1x).is_in(['A', 'T', 'G', 'C']) & pl.col(a2).is_in(['A', 'T', 'G', 'C']) & pl.col(a2x).is_in(['A', 'T', 'G', 'C'])
            con1 = pl.col(a1).is_in(['A', 'T']) & pl.col(a2).is_in(['A', 'T']) & (pl.col(a1x).is_in(['G', 'C']) | pl.col(a2x).is_in(['G', 'C']))
            con2 = pl.col(a1).is_in(['G', 'C']) & pl.col(a2).is_in(['G', 'C']) & (pl.col(a1x).is_in(['A', 'T']) | pl.col(a2x).is_in(['A', 'T']))
            con3 = pl.col(a1x).is_in(['A', 'T']) & pl.col(a2x).is_in(['A', 'T']) & (pl.col(a1).is_in(['G', 'C']) | pl.col(a2).is_in(['G', 'C']))
            con4 = pl.col(a1x).is_in(['G', 'C']) & pl.col(a2x).is_in(['G', 'C']) & (pl.col(a1).is_in(['A', 'T']) | pl.col(a2).is_in(['A', 'T']))
            return (~con0) | (con1) | (con2) | (con3) | (con4)
        
        # Count the number of SNPs with invalid alleles
        def filter_snps(self, alleles_cols: List[str]) -> pl.LazyFrame:
            return self._lf.filter(~self._lf.alleles._invalid_snps(alleles_cols))

        def _match_alleles(self, alleles_cols: List[str]) -> pl.Expr:
            a1, a2, a1x, a2x = alleles_cols
            return (pl.col(a1) == pl.col(a1x)) | (pl.col(a1).is_in(['A', 'T']) & pl.col(a1x).is_in(['A', 'T'])) | (pl.col(a1).is_in(['G', 'C']) & pl.col(a1x).is_in(['G', 'C']))

        def align_alleles_z(self, alleles_cols: List[str], z_col: str) -> pl.LazyFrame:
            return self._lf.with_columns(
                pl.when(self._lf.alleles._match_alleles(alleles_cols)).then(pl.col(z_col)).otherwise(-pl.col(z_col)).alias(z_col)
            )
    
    ss_yhat_unlab = _read_ss(ss_fh=ss_in_yhat_unlab, binary=binary, unlab=True, ylab = False, log=log)
    ss_y_lab = _read_ss(ss_fh=ss_in_y_lab, binary=binary, unlab=False, ylab = True, log=log)
    ss_yhat_lab = _read_ss(ss_fh=ss_in_yhat_lab, binary=binary, unlab=False, ylab = False, log=log)
    
    log.log("\n--- Parsing these three input GWAS")
    # Remove SNPs with any duplicated IDs
    nrow = ss_y_lab.select(pl.col(const_dict['SNP']).count())[0,0]
    unique_filter = (
        ss_y_lab.groupby(pl.col(const_dict['SNP']))
        .agg(pl.count().alias("count"))
        .filter(pl.col("count") == 1)
        .select(const_dict['SNP'])
        )
    ss_y_lab = ss_y_lab.join(unique_filter, on=const_dict['SNP'], how='inner')
    nrow -= ss_y_lab.select(pl.col(const_dict['SNP']).count())[0,0]
    log.log(f"--- Removed {nrow} SNPs with any duplicated IDs")

    # Match A1 and A2; align Z; Remove SNPs with non-biallelic alleles or unmatched alleles across different input GWAS
    nrows = ss_y_lab.shape[0]
    df = reduce(_merge_match_a1a2, [ss_y_lab, ss_yhat_unlab, ss_yhat_lab])
    nrows -= df.shape[0]
    log.log(f"--- Removed {nrows} SNPs with non-biallelic alleles or unmatched alleles across different input GWAS")

    return df, const_dict['N_LAB1'], const_dict['N_LAB2'], const_dict['N_LAB_CASE'] if binary else None, const_dict['N_UNLAB'], const_dict['EAF']

def _format_out(df):
    @pl.api.register_expr_namespace("decimal")
    class FormatDecimalOperations:
        def __init__(self, expr: pl.Expr):
            self._expr = expr

        def to_scientific(self, decimals: int) -> pl.Expr:
            exponent = pl.when(self._expr==0).then(0).otherwise(self._expr.abs().log10().floor().cast(pl.Int32, strict=False))
            mantissa = (self._expr / (10.0 ** exponent)).round(decimals).cast(pl.Utf8).str.ljust(decimals + 2, '0')
            return pl.concat_str([mantissa, pl.lit('e'), pl.when(exponent >= 0).then(pl.lit('+')).otherwise(pl.lit('-')), exponent.abs().cast(pl.Utf8).str.rjust(2, '0')])

        def to_positional(self, decimals: int) -> pl.Expr:
            s = self._expr.round(decimals).cast(pl.Utf8).str.split(by='.')
            return pl.concat_str([s.list.get(0), pl.lit('.'), s.list.get(1).str.ljust(decimals, '0')])
        
        def to_int(self) -> pl.Expr:
            return self._expr.round(0).cast(pl.Int32).cast(pl.Utf8)
        
    df = df.with_columns([pl.col(const_dict['P']).decimal.to_scientific(3).alias(const_dict['P']), pl.col(const_dict['Z']).decimal.to_positional(3).alias(const_dict['Z']), pl.col(const_dict['N_EFF']).decimal.to_int().alias(const_dict['N_EFF'])])
    if const_dict['N_EFF_CASE'] in df.columns:
        df = df.with_columns([pl.col(const_dict['N_EFF_CASE']).decimal.to_int().alias(const_dict['N_EFF_CASE']), pl.col(const_dict['N_EFF_CONTROL']).decimal.to_int().alias(const_dict['N_EFF_CONTROL'])])
    
    # Sorting the DataFrame
    if (const_dict['CHR'] in df.columns and const_dict['BP'] in df.columns):
        return df.sort(by=[const_dict['CHR'], const_dict['BP']])

    return df.sort(by=const_dict['SNP'])
