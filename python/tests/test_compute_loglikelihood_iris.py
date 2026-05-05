# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 09:14:57 2026

@author: Erik
"""

# %% Imports
import numpy as np
from orion.methods.iris.likelihood import compute_loglik


# %% Tests

def test_compute_iris_loglik_runs_and_returns_scalar():
    X = np.array([
        [1, 2],
        [0, 1],
    ], dtype=float)

    lambda1 = np.array([1.0, 1.0])
    lambda2 = np.array([2.0, 2.0])
    lambda3 = np.array([4.0, 4.0])
    pi = np.array([0.3, 0.3, 0.4])

    loglik = compute_loglik(X, lambda1, lambda2, lambda3, pi)

    assert isinstance(loglik, float) or np.isscalar(loglik)
    assert np.isfinite(loglik)
