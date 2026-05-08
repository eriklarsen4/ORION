# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:03:53 2026

@author: Erik
"""
# %% Import
import numpy as np
import pandas as pd

from orion.methods.saint.pipeline import run_saint_pipeline

# %% Test
def test_saint_pipeline_integration():
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

    results = run_saint_pipeline(
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
    assert hasattr(inferred, "proteins_by_bait_unit")

    # Bait units inferred from input
    assert set(inferred.bait_units) == {"B1", "CTRL"}

    # Proteins by bait_unit
    assert inferred.proteins_by_bait_unit == {
        "B1": ["P1", "P2", "P3"],
        "CTRL": ["P1", "P2", "P3"],
    }

    # SAINT-specific: control/treatment classification
    assert set(inferred.control_bait_units) == {"CTRL"}
    assert set(inferred.treatment_bait_units) == {"B1"}

    # Pipeline-derived fields exist
    assert hasattr(meta, "pipeline_derived_fields")
    assert isinstance(meta.pipeline_derived_fields.hyperparameters, dict)

    # SAINT has no tau grid, convergence, or iteration counts
    assert meta.pipeline_derived_fields.tau_grid == []
    assert meta.pipeline_derived_fields.convergence == {}
    assert meta.pipeline_derived_fields.iteration_counts == {}

    # Results dataframe structure
    assert "Protein" in df.columns
    assert "BaitUnit" in df.columns

    # CTRL is skipped in SAINT scoring
    assert set(df["BaitUnit"].unique()) == {"B1"}