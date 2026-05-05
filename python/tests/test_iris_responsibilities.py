# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 09:29:14 2026

@author: Erik
"""

# %% Imports
import numpy as np
from orion.methods.iris.responsibilities import compute_responsibilities


# %% Tests

def test_iris_responsibilities_shape_and_normalization():
    # X: 3 preys × 2 replicates
    X = np.array([
        [5, 3],
        [1, 0],
        [2, 2],
    ], dtype=float)

    lambda1 = np.array([1.0, 1.0, 1.0])
    lambda2 = np.array([3.0, 3.0, 3.0])
    lambda3 = np.array([6.0, 6.0, 6.0])
    pi = np.array([0.4, 0.4, 0.2])

    gamma = compute_responsibilities(X, lambda1, lambda2, lambda3, pi)

    # Correct shape
    assert gamma.shape == (3, 3)

    # Rows must sum to 1
    row_sums = gamma.sum(axis=1)
    assert np.allclose(row_sums, 1.0)


def test_iris_responsibilities_extreme_cases():
    # X: 2 preys × 2 replicates
    X = np.array([
        [10, 10],   # strong signal
        [0, 0],     # strong background
    ], dtype=float)

    lambda1 = np.array([1.0, 1.0])
    lambda2 = np.array([2.0, 2.0])
    lambda3 = np.array([20.0, 20.0])
    pi = np.array([0.3, 0.3, 0.4])

    gamma = compute_responsibilities(X, lambda1, lambda2, lambda3, pi)

    # Prey 0 should strongly favor component 3 (signal)
    assert gamma[0, 2] > gamma[0, 0]
    assert gamma[0, 2] > gamma[0, 1]

    # Prey 1 should strongly favor component 1 (background)
    assert gamma[1, 0] > gamma[1, 1]
    assert gamma[1, 0] > gamma[1, 2]