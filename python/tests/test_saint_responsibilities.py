# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 09:28:38 2026

@author: Erik
"""

# %% Imports
import numpy as np
from orion.methods.saint.responsibilities import compute_responsibilities


# %% Tests

def test_saint_responsibilities_shape_and_normalization():
    # X: 3 preys × 2 replicates
    X = np.array([
        [5, 3],
        [1, 0],
        [2, 2],
    ], dtype=float)

    lambda1 = np.array([1.0, 1.0, 1.0])
    lambda2 = np.array([4.0, 4.0, 4.0])
    pi = np.array([0.7, 0.3])

    gamma = compute_responsibilities(X, lambda1, lambda2, pi)

    # Correct shape
    assert gamma.shape == (3, 2)

    # Rows must sum to 1
    row_sums = gamma.sum(axis=1)
    assert np.allclose(row_sums, 1.0)


def test_saint_responsibilities_extreme_cases():
    # X: 2 preys × 2 replicates
    X = np.array([
        [10, 10],   # strong signal
        [0, 0],     # strong background
    ], dtype=float)

    lambda1 = np.array([1.0, 1.0])
    lambda2 = np.array([20.0, 20.0])
    pi = np.array([0.5, 0.5])

    gamma = compute_responsibilities(X, lambda1, lambda2, pi)

    # Prey 0 should strongly favor component 2 (signal)
    assert gamma[0, 1] > gamma[0, 0]

    # Prey 1 should strongly favor component 1 (background)
    assert gamma[1, 0] > gamma[1, 1]