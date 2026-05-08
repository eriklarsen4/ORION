# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 10:25:15 2026

@author: Erik
"""

# %% Import

import numpy as np

# %% Likelihood functions

def poisson_logpmf(x, lam):
    """
    Compute the log of the Poisson probability mass function.

    x is the observed count for each prey (protein).
    lam is the Poisson rate parameter for each prey (protein).
    """

    eps = 1e-12

    out = x * np.log(lam + eps) - lam

    return out


def compute_loglik(
    X,
    lambda1,
    lambda2,
    pi
):
    """
    Compute the observed-data log likelihood for the two-component mixture model.

    Each lambda is a component-specific Poisson rate for each prey (protein).
    pi is the vector of mixture proportions for the two components
    ("background", "signal").
    The log likelihood is the sum over preys (proteins) of the log of the
    mixture-weighted Poisson likelihood.

    X has shape n_preys by n_conditions.
    lambda1 and lambda2 are prey-specific rate vectors.
    pi is a vector (length of two) of mixture weights.
    """

    x_sum = X.sum(axis=1)

    ll1 = poisson_logpmf(x_sum, lambda1) + np.log(pi[0])
    ll2 = poisson_logpmf(x_sum, lambda2) + np.log(pi[1])

    ll = np.vstack([ll1, ll2]).T

    ll_max = np.max(ll, axis=1, keepdims=True)
    ll_stable = ll - ll_max

    exp_ll = np.exp(ll_stable)
    row_sum = exp_ll.sum(axis=1)

    loglik = np.log(row_sum).sum() + ll_max.sum()

    return loglik

