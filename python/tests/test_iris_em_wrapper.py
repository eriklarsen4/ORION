# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 09:34:40 2026

@author: Erik
"""

# %% Imports
import numpy as np
from orion.methods.iris.em_wrapper import run_em


# %% Tests

def test_iris_runs_and_returns_expected_keys():
    X = np.array([[5, 3], [1, 0]])

    hyperparams = {
        "lambda1_init": np.array([1.0, 1.0]),
        "lambda2_init": np.array([2.0, 2.0]),
        "lambda3_init": np.array([4.0, 4.0]),
        "pi_init": np.array([0.5, 0.3, 0.2]),
    }

    results = run_em(
        X=X,
        hyperparams=hyperparams,
        bait_unit='BAIT',
        max_iter=5,
        tol_loglik=1e-6,
        tol_params=1e-6,
        seed=1,
        verbose=False,
    )

    expected_keys = {
        "loglik_history",
        "lambda1_history",
        "lambda2_history",
        "lambda3_history",
        "pi_history",
        "gamma_history",
        "alpha_history",
        "a_history",
        "b_history",
        "lambda1",
        "lambda2",
        "lambda3",
        "tau",
        "pi",
        "gamma",
    }

    assert expected_keys.issubset(results.keys())

