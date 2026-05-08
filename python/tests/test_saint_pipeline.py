# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:01:26 2026

@author: Erik
"""
# %% Import
import pandas as pd

from orion.methods.saint.pipeline import run_saint_pipeline

# %% Test
def test_saint_pipeline_basic_structure():
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

    results = run_saint_pipeline(
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

    # tau_info must be dict of empty dicts (classical SAINT has no tau grid)
    assert isinstance(raw["tau_info"], dict)
    for v in raw["tau_info"].values():
        assert v == {}

    # Metadata structure
    meta = results["metadata"]

    # Unified metadata model: PipelineMetadata
    assert hasattr(meta, "inferred_fields")
    assert hasattr(meta.inferred_fields, "baits")
    assert hasattr(meta.inferred_fields, "proteins_by_bait_unit")

    # Bait units inferred from input
    assert set(meta.inferred_fields.bait_units) == {"B1", "CTRL"}

    # Proteins by bait_unit
    assert meta.inferred_fields.proteins_by_bait_unit == {
        "B1": ["P1", "P2"],
        "CTRL": ["P1", "P2"],
    }

    # SAINT-specific: control/treatment classification
    assert set(meta.inferred_fields.control_bait_units) == {"CTRL"}
    assert set(meta.inferred_fields.treatment_bait_units) == {"B1"}

    # Pipeline-derived fields exist
    assert hasattr(meta, "pipeline_derived_fields")
    assert isinstance(meta.pipeline_derived_fields.hyperparameters, dict)

    # SAINT has no tau grid, convergence, or iteration counts
    assert meta.pipeline_derived_fields.tau_grid == []
    assert meta.pipeline_derived_fields.convergence == {}
    assert meta.pipeline_derived_fields.iteration_counts == {}

    # Results dataframe structure
    df = results["results_df"]
    assert "Protein" in df.columns
    assert "BaitUnit" in df.columns
    assert set(df["BaitUnit"].unique()) == {"B1"}  # CTRL is skipped in SAINT scoring