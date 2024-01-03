import polars as pl
from scipy.special import chdtrc
from utils import const_dict

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

def _effect_std_to_allele(effect_std, eaf):
    # return effect_allele
    factor = (2*eaf*(1-eaf)).sqrt()
    return effect_std / factor

def _case_prop(n_case, n):
    return n_case / n
    
def _cal_qt_wtd(z_df, n, N, c, r, eaf):
    # 1 Transform the Z-score into effect size
    beta_yhat_unlab = _z_to_beta(pl.col(const_dict['Z_YHAT_UNLAB']), N)
    beta_yhat_lab = _z_to_beta(pl.col(const_dict['Z_YHAT_LAB']), n)
    beta_y_lab = _z_to_beta(pl.col(const_dict['Z_Y_LAB']), n)
    
    # 2 Calculate the corrected effect size and its correspond standard error
    beta_popgwas = beta_y_lab + c * (beta_yhat_unlab - beta_yhat_lab)
    se_popgwas = ((1 + c**2 - 2*c*r) / n + c**2 / N).sqrt()
    z_popgwas, p_popgwas = _beta_se_to_z_p(beta_popgwas, se_popgwas)
    n_eff = 1 / se_popgwas**2
    
    # 3 Transform the standardized scale into the allele scale
    return z_df.with_columns([_effect_std_to_allele(beta_popgwas, eaf).alias(const_dict['BETA']), _effect_std_to_allele(se_popgwas, eaf).alias(const_dict['SE']), z_popgwas.alias(const_dict['Z']), p_popgwas.alias(const_dict['P']), n_eff.alias(const_dict['N_EFF'])]).drop([const_dict['Z_YHAT_LAB'], const_dict['Z_YHAT_UNLAB'], const_dict['Z_Y_LAB']])
    
# d is data calculated by qt
def _qt_to_bt(d, n, n_case, eaf):
    case_prop = _case_prop(n_case, n)
    or_popgwas = _linear_to_log_or(pl.col(const_dict['BETA']), case_prop, eaf)
    se_popgwas = _linear_to_log_or(pl.col(const_dict['SE']), case_prop, eaf)
    n_eff_case = case_prop * pl.col(const_dict['N_EFF'])
    n_eff_control = (1 - case_prop) * pl.col(const_dict['N_EFF'])

    return d.with_columns([or_popgwas.alias(const_dict['BETA']), se_popgwas.alias(const_dict['SE']), n_eff_case.alias(const_dict['N_EFF_CASE']), n_eff_control.alias(const_dict['N_EFF_CONTROL'])]).rename({const_dict['BETA']: const_dict['OR']})

def estimate_popgwas(z_df, n_col, n_case_col, N_col, eaf_col, bt, r):
    n, n_case, N, eaf = pl.col(n_col), pl.col(n_case_col) if bt else None, pl.col(N_col), pl.col(eaf_col)
    r = pl.lit(r)
    c = r * N / (N + n)
    
    d = _cal_qt_wtd(z_df=z_df, n=n, N=N, c=c, r=r, eaf=eaf)
    
    return _qt_to_bt(d=d, n=n, n_case=n_case, eaf=eaf).drop([n_col, n_case_col, N_col]) if bt else d.drop([n_col, N_col])
