# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 08:30:57 2026

@author: Erik
"""

#%% Import
import pandas as pd

from orion.methods.iris.pipeline import run_iris_pipeline

#%% tests
def test_iris_pipeline_runs_and_returns_expected_structure():
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

    results = run_iris_pipeline(
        input_data,
        hyperparams=None,
        max_iter=5,
        tol_loglik=1e-6,
        tol_params=1e-6,
        seed=1,
        verbose=False,
        make_plots=False,
        show_plots=True,
        save_plots=False,
        plot_dir=None,
        save_results=False,
        results_csv=None,
    )

    # Top-level keys
    assert set(results.keys()) == {"raw_outputs", "metadata", "results_df"}

    # Metadata structure
    meta = results["metadata"]
    assert isinstance(meta, dict)
    for key in ["bait_units", "proteins_by_bait_unit", "control_bait_units", "treatment_bait_units"]:
        assert key in meta

    # Raw outputs structure
    raw = results["raw_outputs"]
    assert isinstance(raw, dict)
    assert "em_results" in raw
    assert "tau_info" in raw

    # EM results must contain all bait-units
    assert set(raw["em_results"].keys()) == set(meta["bait_units"])

    # Results table
    df = results["results_df"]
    assert isinstance(df, pd.DataFrame)
    assert not df.empty

    # Core columns
    for col in ["Protein", "BaitUnit", "gamma3"]:
        assert col in df.columns

    # Replicate columns present
    assert any(col.startswith("rep") for col in df.columns)

    # Control bait-units excluded from results_df
    assert all(bu not in df["BaitUnit"].unique() for bu in meta["control_bait_units"])
