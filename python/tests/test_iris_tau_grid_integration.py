# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 16:39:06 2026

@author: Erik
"""

# %% Import
import pandas as pd

from orion.methods.iris.pipeline import run_iris_pipeline

# %% Tests
def test_tau_grid_integration_structure():
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

    # Run the full IRIS pipeline (tau-grid included)
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

    # Top-level structure
    assert set(results.keys()) == {"raw_outputs", "metadata", "results_df"}

    # Metadata structure
    meta = results["metadata"]
    assert isinstance(meta, dict)
    for key in ["bait_units", "proteins_by_bait_unit", "control_bait_units", "treatment_bait_units"]:
        assert key in meta

    # Tau-grid integration
    raw = results["raw_outputs"]
    assert "tau_info" in raw
    tau_info = raw["tau_info"]
    assert isinstance(tau_info, dict)

    # Each tau_info entry must contain the tau-grid structure
    for bu, info in tau_info.items():
        assert isinstance(info, dict)
        for key in ["best_tau", "best_result", "tau_grid_results", "convergence_info", "iteration_counts", "tau_grid"]:
            assert key in info

    # EM results must exist for each bait-unit
    assert "em_results" in raw
    assert set(raw["em_results"].keys()) == set(meta["bait_units"])

    # Results table
    df = results["results_df"]
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "Protein" in df.columns
    assert "BaitUnit" in df.columns
    assert "gamma3" in df.columns
