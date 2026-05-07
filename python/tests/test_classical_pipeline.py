# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:01:26 2026

@author: Erik
"""
# %% Import
import pandas as pd

from orion.methods.saint.classical_saint import run_classical_pipeline

# %% Test
def test_classical_pipeline_basic_structure():
    # Minimal synthetic wide-format input
    input_data = pd.DataFrame(
        {
            "Protein": ["P1", "P2"],
            "condA_B1_1": [5, 1],
            "condA_B1_2": [3, 0],
            "condA_CTRL_1": [2, 0],
            "condA_CTRL_2": [1, 1],
        }
    )

    results = run_classical_pipeline(
        input_data=input_data,
        max_iter=5,
        seed=1,
        make_plots=False,
    )

    # Top-level structure
    assert set(results.keys()) == {"raw_outputs", "metadata", "results_df"}

    # Raw outputs structure
    raw = results["raw_outputs"]
    assert "em_results" in raw
    assert "tau_info" in raw

    # tau_info must be dict of empty dicts (classical has no tau grid)
    assert isinstance(raw["tau_info"], dict)
    for v in raw["tau_info"].values():
        assert v == {}

    # Metadata structure
    meta = results["metadata"]
    assert hasattr(meta, "inferred_fields")
    assert hasattr(meta.inferred_fields, "baits")
    assert hasattr(meta.inferred_fields, "proteins_by_bait")
    assert hasattr(meta.inferred_fields, "control_baits")
    assert hasattr(meta.inferred_fields, "treatment_baits")

    # Results table
    df = results["results_df"]
    assert isinstance(df, pd.DataFrame)
    assert not df.empty

    # Required columns (2-component classical model)
    expected_cols = {
        "Protein",
        "bait",
        "lambda1",
        "lambda2",
        "pi1",
        "pi2",
        "gamma1",
        "gamma2",
    }
    assert expected_cols.issubset(df.columns)

    # No third-component columns
    forbidden_cols = {"lambda3", "pi3", "gamma3", "tau"}
    assert forbidden_cols.isdisjoint(df.columns)

    # Sorting rule: gamma2 desc, then Protein asc
    df_sorted = df.sort_values(
        by=["gamma2", "Protein"],
        ascending=[False, True],
        ignore_index=True,
    )
    assert df.equals(df_sorted)