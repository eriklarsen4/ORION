# -*- coding: utf-8 -*-
"""
Created on Sun Feb  8 08:24:09 2026

@author: Erik
"""

# %% Import

import numpy as np


# %% Pi updates

def update_pi(gamma):
    """
    Update mixture weights using responsibilities.

    pi is the vector of mixture proportions for the three components ("background", "noise", "signal").
    Each pi value is the average posterior probability of its component
    across all preys.

    Parameters
    gamma : array with shape (n_preys, 3)
        Responsibilities for each component.

    Returns
    pi : array with shape (3,)
        Updated mixture weights.
    """

    n_preys = gamma.shape[0]

    pi = gamma.sum(axis=0) / n_preys

    return pi



