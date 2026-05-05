# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 21:06:40 2026

@author: Erik
"""

# %% Imports

import numpy as np


# %% Responsibilities for three component Poisson mixture

def compute_responsibilities(
    X,
    lambda1,
    lambda2,
    lambda3,
    pi
):
    """
    Compute responsibilities for a three-component Poisson mixture ("background", "noise", "signal").


    PARAMETERS
    
    X: prey-level count matrix by bait_unit (treatment/condition/tag + replicate)
    lambda1: prey-level Poisson rate for background component
    lambda2: prey-level Poisson rate for noise/ambiguous component
    lambda3: prey-level Poisson rate for signal/true interactor component
    pi: vector of mixture proportions for the three components.
    
    
    RETURNS
    
    gamma: array represents the posterior probability that each prey belongs
    to each component (gamma1: pp for background component
                       gamma2: pp for noise/ambiguous component
                       gamma3: pp for signal/true interactor component).
    """

    rate1 = lambda1[:, None]
    rate2 = lambda2[:, None]
    rate3 = lambda3[:, None]

    loglik1 = X * np.log(rate1) - rate1
    loglik2 = X * np.log(rate2) - rate2
    loglik3 = X * np.log(rate3) - rate3

    loglik1 = loglik1.sum(axis=1) + np.log(pi[0])
    loglik2 = loglik2.sum(axis=1) + np.log(pi[1])
    loglik3 = loglik3.sum(axis=1) + np.log(pi[2])

    maxlog = np.maximum(np.maximum(loglik1, loglik2), loglik3)

    w1 = np.exp(loglik1 - maxlog)
    w2 = np.exp(loglik2 - maxlog)
    w3 = np.exp(loglik3 - maxlog)

    total = w1 + w2 + w3

    gamma = np.column_stack([w1 / total, w2 / total, w3 / total])

    return gamma


