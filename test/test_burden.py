"""Binary burden input, scale, and end-to-end regression tests."""
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import polars as pl
from polars.testing import assert_frame_equal
import pytest
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import compute
import utils


class Log:
    def log(self, message):
        pass


def frame(z=(1., -2., 0.5, 3.)):
    return pl.DataFrame({"GENE": ["A", "B", "C", "D"], "MASK": ["pLoF"] * 4,
                         "Z": list(z), "N": [1000] * 4})


def read(tmp_path, data, binary=True):
    path = tmp_path / "input.tsv"
    data.write_csv(path, separator="\t", null_value="NA")
    return utils._read_ss_gene(str(path), binary, False, True, Log())


def test_total_n_and_case_counts_are_equivalent(tmp_path):
    total = frame()
    counts = total.drop("N").with_columns(pl.lit(100).alias("N_case"), pl.lit(900).alias("N_control"))
    assert_frame_equal(read(tmp_path, total, False), read(tmp_path, total, True))
    assert_frame_equal(read(tmp_path, total), read(tmp_path, counts))
    assert_frame_equal(read(tmp_path, total), read(tmp_path, counts.with_columns(pl.lit(1000).alias("N"))))


@pytest.mark.parametrize("value", [0, -1, 1.5, float("inf"), float("nan"), None])
@pytest.mark.parametrize("column", ["N", "N_case", "N_control"])
def test_invalid_counts_fail(tmp_path, value, column):
    data = frame().with_columns(pl.lit(100).alias("N_case"), pl.lit(900).alias("N_control"))
    data = data.with_columns(pl.lit(value).alias(column))
    with pytest.raises(ValueError, match="positive finite integer"):
        read(tmp_path, data)


def test_inconsistent_or_incomplete_counts_fail(tmp_path):
    with pytest.raises(ValueError, match="both N_case"):
        read(tmp_path, frame().with_columns(pl.lit(100).alias("N_case")))
    with pytest.raises(ValueError, match="N must equal"):
        read(tmp_path, frame().with_columns(pl.lit(100).alias("N_case"), pl.lit(800).alias("N_control")))
    with pytest.raises(ValueError, match="requires N"):
        read(tmp_path, frame().drop("N"))


def test_duplicate_tests_fail_and_na_z_is_dropped(tmp_path):
    with pytest.raises(ValueError, match="duplicate GENE/MASK"):
        read(tmp_path, pl.concat([frame(), frame().head(1)]))
    assert read(tmp_path, frame((None, -2., 0.5, 3.))).height == 3
    with pytest.raises(ValueError, match="Z must be finite"):
        read(tmp_path, frame((float("inf"), -2., 0.5, 3.)))


def inputs(tmp_path, binary=False, mismatched=False):
    paths = []
    for index, z in enumerate([(1., -2., 0.5, 3.), (0.5, -1., 2., -0.3), (1.5, 0.2, -0.4, 2.)]):
        data = frame(z)
        if binary and index == 0:
            data = data.drop("N").with_columns(pl.lit(100).alias("N_case"), pl.lit(900).alias("N_control"))
        if mismatched and index == 1:
            data = data.with_columns(pl.lit("missense").alias("MASK"))
        path = tmp_path / f"{index}.tsv"
        data.write_csv(path, separator="\t")
        paths.append(str(path))
    return paths


def combine(paths, binary=False, r=0.5):
    df, n, n1, n2 = utils._read_z_gene(paths[2], paths[0], paths[1], binary, Log())
    return compute.estimate_popburden(df, n, n1, n2, False, Log(), r12=r).collect()


def test_analytic_equal_sample_size_result(tmp_path):
    result = combine(inputs(tmp_path)).sort("GENE")
    z0 = np.array([1., -2., 0.5, 3.])
    z1 = np.array([0.5, -1., 2., -0.3])
    z2 = np.array([1.5, 0.2, -0.4, 2.])
    # Equal n, no overlap: omega=r/2 and Var(beta)=(1-r^2/2)/n.
    beta = (z0 + 0.25 * (z2 - z1)) / math.sqrt(1000)
    se = math.sqrt(0.875 / 1000)
    np.testing.assert_allclose(result["BETA"], beta)
    np.testing.assert_allclose(result["SE"], se)
    np.testing.assert_allclose(result["Z"], beta / se)
    np.testing.assert_allclose(result["P"], 2 * norm.sf(abs(beta / se)))
    assert result["EFFECT_SCALE"].unique().to_list() == ["Z_OVER_SQRT_N"]
    assert not {"OR", "LOG_OR"}.intersection(result.columns)


def test_zero_correlation_recovers_observed_z(tmp_path):
    result = combine(inputs(tmp_path), r=0).sort("GENE")
    np.testing.assert_allclose(result["Z"], [1., -2., 0.5, 3.])
    np.testing.assert_allclose(result["N_eff"], 1000)


def test_binary_inputs_preserve_results(tmp_path):
    expected = combine(inputs(tmp_path)).sort("GENE")
    actual = combine(inputs(tmp_path, binary=True), binary=True).sort("GENE")
    assert_frame_equal(expected, actual)


def test_mismatched_masks_fail(tmp_path):
    with pytest.raises(ValueError, match="No matching GENE/MASK"):
        combine(inputs(tmp_path, mismatched=True))


def test_merge_matches_gene_and_mask():
    left = pl.DataFrame({"GENE": ["A", "A"], "MASK": ["pLoF", "missense"],
                         "Z_y_lab": [1., 2.], "N_lab1": [1000, 1000]})
    right = pl.DataFrame({"GENE": ["A", "A"], "MASK": ["missense", "pLoF"],
                          "Z": [20., 10.], "N_unlab": [2000, 2000]})
    merged = utils._merge_gene(left, right).sort("MASK")
    assert merged.height == 2
    assert merged["Z_yhat_unlab"].to_list() == [20., 10.]


@pytest.mark.parametrize("overlap,r", [(False, 1.), (False, None), (True, None)])
def test_bundled_burden_data_binary_flag_is_invariant(overlap, r):
    paths = [str(ROOT / "test" / "data" / f"Head_BMD_{part}.Burden.txt.gz")
             for part in ["yhat_unlab", "y_lab", "yhat_lab"]]
    results = []
    for binary in [False, True]:
        df, n, n1, n2 = utils._read_z_gene(*paths, binary, Log())
        result = compute.estimate_popburden(df, n, n1, n2, overlap, Log(), r12=r).collect()
        assert result.height > 0
        assert result["P"].is_finite().all()
        results.append(result.sort("GENE", "MASK"))
    assert_frame_equal(*results)


@pytest.mark.parametrize("lazy", [False, True])
def test_manual_binary_conversion_rejected(tmp_path, lazy):
    result = combine(inputs(tmp_path))
    with pytest.raises(ValueError, match="cannot be converted to OR"):
        compute._qt_to_bt(result.lazy() if lazy else result, pl.lit(1000), pl.lit(100), None)


@pytest.mark.parametrize("binary", [False, True])
@pytest.mark.parametrize("overlap", [False, True])
def test_cli(tmp_path, binary, overlap):
    paths = inputs(tmp_path, binary=binary)
    prefix = tmp_path / "result"
    cmd = [sys.executable, str(ROOT / "POP-RARE.py"), "--burden",
           "--gwas-y-lab", paths[0], "--gwas-yhat-lab", paths[1],
           "--gwas-yhat-unlab", paths[2], "--out", str(prefix)]
    if binary:
        cmd += ["--bt"]
    cmd += ["--ovp"] if overlap else ["--r", "0.5"]
    run = subprocess.run(cmd, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    result = pl.read_csv(str(prefix) + "_POP-Burden.txt", separator="\t")
    assert result.height == 4
    assert result["P"].is_finite().all()
    assert ((result["P"] >= 0) & (result["P"] <= 1)).all()
    assert result["EFFECT_SCALE"].unique().to_list() == ["Z_OVER_SQRT_N"]
    assert "not original-scale effects or log-OR" in run.stdout
    assert "OR" not in result.columns
