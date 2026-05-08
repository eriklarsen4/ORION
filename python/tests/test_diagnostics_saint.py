# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:02:15 2026

@author: Erik
"""
# %% Import
import numpy as np
import matplotlib
matplotlib.use("Agg")

from orion.methods.saint.diagnostics import make_saint_plots

# %% Test
def test_make_saint_plots_returns_figures():
    # Minimal synthetic classical EM results
    results_em = {
        "loglik_history": [0.0, 1.0, 2.0],

        "lambda1_history": [
            np.array([1.0, 2.0]),
            np.array([1.5, 2.5]),
        ],
        "lambda2_history": [
            np.array([3.0, 4.0]),
            np.array([3.5, 4.5]),
        ],

        "pi_history": [
            np.array([0.7, 0.3]),
            np.array([0.6, 0.4]),
        ],

        "gamma_history": [
            np.array([
                [0.8, 0.2],
                [0.3, 0.7],
            ])
        ],

        # Final gamma matrix for histogram
        "gamma": np.array([
            [0.8, 0.2],
            [0.3, 0.7],
        ]),
    }

    figs = make_saint_plots(results_em, bait_name="B1")

    assert set(figs.keys()) == {"loglik", "lambda", "pi", "gamma"}

    # Both should be matplotlib Figure-like objects
    for fig in figs.values():
        assert hasattr(fig, "savefig")
        assert hasattr(fig, "clf")