# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 08:59:20 2026

@author: Erik
"""

# %%Import
import numpy as np
from orion.methods.saint.classical_em_wrapper import run_em_classical

# %% Test
def test_run_em_classical_basic():
    # Simple synthetic counts
    X = np.array([5, 0, 3, 1], dtype=float)

    hyper = {
        "lambda1_init": np.array([1.0, 1.0, 1.0, 1.0]),
        "lambda2_init": np.array([3.0, 3.0, 3.0, 3.0]),
        "pi_init": np.array([0.6, 0.4]),
    }

    results = run_em_classical(
        X,
        hyper,
        max_iter=10,
        tol_loglik=1e-6,
        tol_params=1e-6,
        seed=1,
        verbose=False,
    )

    # Required keys
    assert set(results.keys()) >= {
        "loglik_history",
        "lambda1_history",
        "lambda2_history",
        "pi_history",
        "gamma_history",
        "lambda1",
        "lambda2",
        "pi",
        "gamma"
    }

    # Shapes
    assert results["lambda1"].shape == X.shape
    assert results["lambda2"].shape == X.shape
    assert results["gamma"].shape == (len(X), 2)

    # pi must be length 2
    assert results["pi"].shape == (2,)

    # gamma rows must sum to 1
    row_sums = results["gamma"].sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-6)

    # loglik must be monotone non-decreasing
    ll = results["loglik_history"]
    assert all(ll[i] <= ll[i + 1] for i in range(len(ll) - 1))