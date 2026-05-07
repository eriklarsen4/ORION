# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:03:53 2026

@author: Erik
"""
# %% Import
import numpy as np
import pandas as pd

from orion.methods.saint.classical_saint import run_classical_pipeline

# %% Test
def test_classical_pipeline_integration():
    # Minimal synthetic wide-format input
    input_data = pd.DataFrame(
        {
            "Protein": ["P1", "P2", "P3"],
            "condA_B1_1": [10, 1, 0],
            "condA_B1_2": [5, 0, 0],
            "condA_CTRL_1": [2, 0, 0],
            "condA_CTRL_2": [1, 0, 0],
        }
    )

    results = run_classical_pipeline(
        input_data=input_data,
        max_iter=10,
        seed=1,
        make_plots=False,
    )

    # --- Top-level structure ---
    assert set(results.keys()) == {"raw_outputs", "metadata", "results_df"}

    raw = results["raw_outputs"]
    meta = results["metadata"]
    df = results["results_df"]

    # --- Raw outputs ---
    assert "em_results" in raw
    assert "tau_info" in raw

    # tau_info must be dict of empty dicts
    assert isinstance(raw["tau_info"], dict)
    for v in raw["tau_info"].values():
        assert v == {}

    # EM results must exist for each bait
    assert isinstance(raw["em_results"], dict)
    assert set(raw["em_results"].keys()) == set(meta.inferred_fields.baits)

    # --- Metadata structure ---
    inferred = meta.inferred_fields

    assert hasattr(inferred, "baits")
    assert hasattr(inferred, "proteins_by_bait")
    assert hasattr(inferred, "control_baits")
    assert hasattr(inferred, "treatment_baits")

    assert isinstance(inferred.baits, list)
    assert isinstance(inferred.proteins_by_bait, dict)

    # --- Results table ---
    assert isinstance(df, pd.DataFrame)
    assert not df.empty

    # Required classical columns
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

    # Forbidden third-component columns
    forbidden_cols = {"lambda3", "pi3", "gamma3", "tau"}
    assert forbidden_cols.isdisjoint(df.columns)

    # Sorting rule: gamma2 desc, then Protein asc
    df_sorted = df.sort_values(
        by=["gamma2", "Protein"],
        ascending=[False, True],
        ignore_index=True,
    )
    assert df.equals(df_sorted)

    # --- EM output sanity checks ---
    for bait, em in raw["em_results"].items():
        # Required keys
        assert set(em.keys()) >= {
            "lambda1",
            "lambda2",
            "pi",
            "gamma",
            "loglik_history",
            "convergence_info",
            "iteration_count",
        }

        # Shapes
        n = len(inferred.proteins_by_bait[bait])
        assert em["lambda1"].shape == (n,)
        assert em["lambda2"].shape == (n,)
        assert em["gamma"].shape == (n, 2)

        # pi must be length 2
        assert em["pi"].shape == (2,)

        # gamma rows sum to 1
        assert np.allclose(em["gamma"].sum(axis=1), 1.0, atol=1e-6)