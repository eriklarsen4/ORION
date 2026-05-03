# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 12:50:02 2026

@author: Erik
"""

# %% Imports

import numpy as np


# %% Validation

def validate_hyperparams(max_iter, tol, lambda_vecs, pi):
    """
    Validate all numeric hyperparameters and parameter shapes.

    Parameters
    max_iter : integer
        Maximum number of EM iterations.

    tol : float
        Convergence threshold.

    lambda_vecs : list of three arrays
        Each array has shape n_preys.

    pi : array with shape three
        Mixture weights.

    Raises
    ValueError if any validation check fails.
    """

    # max_iter must be a positive integer
    if not isinstance(max_iter, int):
        raise ValueError("max_iter must be an integer")
    if max_iter <= 0:
        raise ValueError("max_iter must be positive")

    # tol must be a positive float
    if not isinstance(tol, float):
        raise ValueError("tol must be a float")
    if tol <= 0.0:
        raise ValueError("tol must be positive")

    # lambda_vecs must be a list of three arrays
    if not isinstance(lambda_vecs, list):
        raise ValueError("lambda_vecs must be a list")
    if len(lambda_vecs) != 3:
        raise ValueError("lambda_vecs must contain three arrays")

    n_preys = None
    for lam in lambda_vecs:
        if not isinstance(lam, np.ndarray):
            raise ValueError("each lambda vector must be a numpy array")
        if lam.ndim != 1:
            raise ValueError("each lambda vector must be one dimensional")
        if n_preys is None:
            n_preys = lam.shape[0]
        else:
            if lam.shape[0] != n_preys:
                raise ValueError("all lambda vectors must have the same length")

    # pi must be a length three array that sums to one
    if not isinstance(pi, np.ndarray):
        raise ValueError("pi must be a numpy array")
    if pi.shape != (3,):
        raise ValueError("pi must have shape three")
    if np.any(pi < 0.0):
        raise ValueError("pi must contain non negative values")
    if not np.isclose(pi.sum(), 1.0):
        raise ValueError("pi must sum to one")