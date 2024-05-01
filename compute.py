import polars as pl
from scipy.special import chdtrc
from utils import const_dict
import numpy as np

def _z_to_beta(z, n):
    return z / n.sqrt() # return beta

def _beta_se_to_z_p(beta, se):
    try:
        z = beta / se
    except (FloatingPointError, ZeroDivisionError):
        z = float('inf')
    p = chdtrc(1, z ** 2)
    return z, p

def _linear_to_log_or(beta, case_prop, eaf):
    return beta / (case_prop * (1 - case_prop))

def _effect_std_to_allele(effect_std, factor):
    # return effect_allele
    return effect_std / factor

def _case_prop(n_case, n):
    return n_case / n
    

def _cal_qt_wtd(z_df, n, N1, N2, r12, r13, r23, eaf):
    # r:
    ## r12: correlation between y_lab and yhat_lab
    ## r13: correlation between y_lab and yhat_unlab
    ## r23: correlation between yhat_lab and yhat_unlab

    # n:
    ## n: number of samples in y_lab
    ## N1: number of samples in yhat_lab
    ## N2: number of samples in yhat_unlab
    rho12 = r12 * (1/N1).sqrt() * (1/n).sqrt()
    rho13 = r13 * (1/N2).sqrt() * (1/n).sqrt()
    rho23 = r23 * (1/N1).sqrt() * (1/N2).sqrt()

    omega = (rho12 - rho13) / (1/N1 + 1/N2 - 2*rho23)

    # 1 Transform the Z-score into effect size
    beta_yhat_unlab = _z_to_beta(pl.col(const_dict['Z_YHAT_UNLAB']), N2)
    beta_yhat_lab = _z_to_beta(pl.col(const_dict['Z_YHAT_LAB']), N1)
    beta_y_lab = _z_to_beta(pl.col(const_dict['Z_Y_LAB']), n)

    # 2 Calculate the corrected effect size and its correspond standard error
    beta_popgwas = beta_y_lab + omega * (beta_yhat_unlab - beta_yhat_lab)
    se_popgwas = (1/n - (rho12 - rho13) / (1/N1 + 1/N2 - 2 * rho23) * (rho12 - rho13)).sqrt()
    z_popgwas, p_popgwas = _beta_se_to_z_p(beta_popgwas, se_popgwas)
    n_eff = 1 / se_popgwas**2

    # 3 Transform the standardized scale into the allele scale
    factor_eaf = (2*eaf*(1-eaf)).sqrt()
    beta_popgwas = _effect_std_to_allele(beta_popgwas, factor_eaf)
    se_popgwas = _effect_std_to_allele(se_popgwas, factor_eaf)

    return z_df.with_columns([beta_popgwas.alias(const_dict['BETA']), se_popgwas.alias(const_dict['SE']), z_popgwas.alias(const_dict['Z']), p_popgwas.alias(const_dict['P']), n_eff.alias(const_dict['N_EFF'])]).drop([const_dict['Z_YHAT_LAB'], const_dict['Z_YHAT_UNLAB'], const_dict['Z_Y_LAB']])

# d is data calculated by qt
def _qt_to_bt(d, n, n_case, eaf):
    case_prop = _case_prop(n_case, n)
    or_popgwas = _linear_to_log_or(pl.col(const_dict['BETA']), case_prop, eaf)
    se_popgwas = _linear_to_log_or(pl.col(const_dict['SE']), case_prop, eaf)
    n_eff_case = case_prop * pl.col(const_dict['N_EFF'])
    n_eff_control = (1 - case_prop) * pl.col(const_dict['N_EFF'])

    return d.with_columns([or_popgwas.alias(const_dict['BETA']), se_popgwas.alias(const_dict['SE']), n_eff_case.alias(const_dict['N_EFF_CASE']), n_eff_control.alias(const_dict['N_EFF_CONTROL'])]).rename({const_dict['BETA']: const_dict['OR']})


def estimate_popgwas(z_df, n_col, N1_col, n_case_col, N2_col, eaf_col, bt, r12, r13, r23):
    n, N1, n_case, N2, eaf = pl.col(n_col), pl.col(N1_col), pl.col(n_case_col) if bt else None, pl.col(N2_col), pl.col(eaf_col)
    r12 = pl.lit(r12)
    r13 = pl.lit(r13)
    r23 = pl.lit(r23)

    d = _cal_qt_wtd(z_df=z_df.lazy(), n=n, N1=N1, N2=N2, r12=r12, r13 = r13, r23 = r23, eaf=eaf)
    
    return _qt_to_bt(d=d, n=n, n_case=n_case, eaf=eaf).drop([n_col, N1_col, n_case_col, N2_col]) if bt else d.drop([n_col, N1_col, N2_col])