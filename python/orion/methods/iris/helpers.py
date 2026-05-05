# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 13:08:00 2026

@author: Erik
"""

# %% Imports

import numpy as np

# %% Hyperparameter builder for IRIS

def build_hyperparams(X):
    """
    Build initial hyperparameters for IRIS.
    
    
    PARAMETERS
    
    X is an array with shape n_preys by n_conditions
    
    
    RETURNS
    
    A dictionary with initial lambda values, pi values, and alpha
    """

    n_preys, n_conditions = X.shape

    row_means = np.mean(X, axis=1)

    lambda1_init = row_means.copy()
    lambda2_init = 0.5 * row_means
    lambda3_init = 0.1 * row_means

    pi_init = np.ones(3, dtype=float) / 3.0

    alpha = np.ones(3, dtype=float)

    hyperparams = {
        "lambda1_init": lambda1_init,
        "lambda2_init": lambda2_init,
        "lambda3_init": lambda3_init,
        "pi_init": pi_init,
        "alpha": alpha
    }

    return hyperparams