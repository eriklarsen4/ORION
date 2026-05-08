# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 09:08:01 2026

@author: Erik
"""
# %% Import

import numpy as np

# %% Classical SAINT compute responsibilities

def compute_responsibilities(X, lambda1, lambda2, pi):
    """
    Classical responsibilities for a two-component Poisson mixture.

    lambda1 is the background rate.
    lambda2 is the signal rate.
    pi is the vector of mixture proportions.
    gamma is the posterior probability for each component.
    """

    rate1 = lambda1[:, None]
    rate2 = lambda2[:, None]

    loglik1 = X * np.log(rate1) - rate1
    loglik2 = X * np.log(rate2) - rate2

    loglik1 = loglik1.sum(axis=1) + np.log(pi[0])
    loglik2 = loglik2.sum(axis=1) + np.log(pi[1])

    maxlog = np.maximum(loglik1, loglik2)

    w1 = np.exp(loglik1 - maxlog)
    w2 = np.exp(loglik2 - maxlog)

    total = w1 + w2

    gamma = np.column_stack([w1 / total, w2 / total])

    return gamma