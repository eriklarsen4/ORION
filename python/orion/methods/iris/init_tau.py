# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 13:31:37 2026

@author: Erik
"""

# %% Imports

import numpy as np

# %% Tau initialization and update

def init_tau():
    """
    Initialize tau for the hierarchical model.
    
    Parameters
    None
    
    Returns
    A scalar initial value for tau
    """
    return 1.0


def update_tau(lambdas, alpha):
    """
    Update tau using a closed form expression under a Gamma prior.
    
    Parameters
    lambdas is a list of three arrays, one per component
    alpha is an array of length three with Gamma shape parameters
    
    Returns
    A scalar updated value for tau
    """

    lambda1, lambda2, lambda3 = lambdas
    all_lambda = np.concatenate([lambda1, lambda2, lambda3])

    n_preys = lambda1.shape[0]

    num = np.sum(alpha) * float(n_preys)
    denom = np.sum(all_lambda)

    tau = num / np.clip(denom, 1e-12, None)
    return tau