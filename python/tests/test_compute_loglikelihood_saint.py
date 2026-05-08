# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 13:55:20 2026

@author: Erik
"""

# %% Imports
import numpy as np
from orion.methods.saint.likelihood import compute_loglik


# %% Tests

def test_compute_saint_loglik_runs_and_returns_scalar():
    # X: 2 preys × 2 replicates
    X = np.array([
        [1, 2],
        [0, 1],
    ], dtype=float)

    lambda1 = np.array([1.0, 1.0])
    lambda2 = np.array([2.0, 2.0])
    pi = np.array([0.6, 0.4])

    loglik = compute_loglik(X, lambda1, lambda2, pi)

    # Basic structural checks
    assert isinstance(loglik, float) or np.isscalar(loglik)
    assert np.isfinite(loglik)