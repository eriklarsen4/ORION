# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 08:35:19 2026

@author: Erik
"""

# %% Import
import numpy as np
import matplotlib
matplotlib.use("Agg")

from orion.methods.iris.diagnostics import make_iris_plots

# %% tests
def test_make_iris_plots_returns_figures():
    # Minimal synthetic hierarchical EM results matching the new architecture
    results_em = {
        "loglik_history": [0.0, 1.0, 2.0],

        "lambda1_history": [
            np.array([1.0, 2.0]),
            np.array([1.5, 2.5])
        ],
        "lambda2_history": [
            np.array([3.0, 4.0]),
            np.array([3.5, 4.5])
        ],
        "lambda3_history": [
            np.array([5.0, 6.0]),
            np.array([5.5, 6.5])
        ],

        "pi_history": [
            np.array([0.6, 0.3, 0.1]),
            np.array([0.5, 0.3, 0.2])
        ],

        "alpha_history": [
            np.array([2.0, 2.0, 2.0]),
            np.array([2.1, 2.0, 1.9])
        ],
        "a_history": [
            np.array([1.0, 1.0, 1.0]),
            np.array([1.1, 1.0, 0.9])
        ],
        "b_history": [
            np.array([0.5, 0.5, 0.5]),
            np.array([0.6, 0.5, 0.4])
        ],

        "gamma_history": [
            np.array([
                [0.8, 0.1, 0.1],
                [0.2, 0.3, 0.5]
            ])
        ],

        # Final gamma matrix for histogram
        "gamma": np.array([
            [0.8, 0.1, 0.1],
            [0.2, 0.3, 0.5]
        ]),
    }

    figs = make_iris_plots(results_em, bait_name="B1")

    # New architecture: only trajectories + gamma
    assert set(figs.keys()) == {"trajectories", "gamma"}

    # Both should be matplotlib Figure objects
    for fig in figs.values():
        # Avoid importing matplotlib in the test; just check duck-typing
        assert hasattr(fig, "savefig")
        assert hasattr(fig, "clf")