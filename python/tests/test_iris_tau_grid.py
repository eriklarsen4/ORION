# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 16:37:51 2026

@author: Erik
"""

# %% Imports
import numpy as np
from orion.methods.iris.tau_grid import run_tau_grid


# %% Tests

def test_tau_grid_runs_and_returns_expected_keys():
    # Minimal synthetic data
    X = np.array([[5, 3], [1, 0]])

    # Minimal hyperparameters (only those used by EM)
    hyperparams = {
        "lambda1_init": np.array([1.0, 1.0]),
        "lambda2_init": np.array([2.0, 2.0]),
        "lambda3_init": np.array([4.0, 4.0]),
        "pi_init": np.array([0.5, 0.3, 0.2]),
    }

    # Run tau grid
    results = run_tau_grid(
        X=X,
        hyperparams=hyperparams,
        bait_unit="BAIT",
    )

    expected_keys = {
        "best_tau",
        "best_result",
        "tau_grid_results",
        "convergence_info",
        "iteration_counts",
        "tau_grid",
    }

    assert expected_keys.issubset(results.keys())

    # Structural checks
    assert isinstance(results["best_tau"], float)
    assert isinstance(results["best_result"], dict)
    assert isinstance(results["tau_grid_results"], dict)
    assert set(results["tau_grid_results"].keys()) == set(results["tau_grid"])

    # EM result structure sanity check
    em0 = results["best_result"]
    assert "lambda1" in em0
    assert "lambda2" in em0
    assert "lambda3" in em0
    assert "pi" in em0
    assert "gamma" in em0
    assert "loglik_history" in em0
