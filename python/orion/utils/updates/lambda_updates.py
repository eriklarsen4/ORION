# -*- coding: utf-8 -*-
"""
Created on Sun Feb  8 08:16:25 2026

@author: Erik
"""

# %% Import

import numpy as np


# %% Classical lambda updates

def update_lambda1(X, gamma):
    """
    Classical SAINT lambda1 update.

    lambda1 is the background Poisson rate per prey.
    No priors are used.
    alpha and tau do not apply here.
    """
    x_sum = np.sum(X, axis=1)
    gamma1 = gamma[:, 0]
    eps = 1e-12
    lambda1 = x_sum / np.clip(gamma1, eps, None)
    return np.clip(lambda1, 1e-6, 1e6)


def update_lambda2(X, gamma):
    """
    Classical SAINT lambda2 update.

    lambda2 is the signal Poisson rate per prey.
    No priors are used.
    alpha and tau do not apply here.
    """
    x_sum = np.sum(X, axis=1)
    gamma2 = gamma[:, 1]
    eps = 1e-12
    lambda2 = x_sum / np.clip(gamma2, eps, None)
    return np.clip(lambda2, 1e-6, 1e6)


# %% Hierarchical lambda updates

def update_lambda_hierarchical(X, gamma, alpha, tau):
    """
    Hierarchical lambda updates for all three components.

    Each lambda is a component-specific (background, noise, and signal) Poisson rate per prey.
    alpha is the shape parameter of the Gamma prior on each lambda.
    tau is the shared rate parameter that controls global shrinkage of lambda.
    """

    n_preys, n_conditions = X.shape
    x_sum = np.sum(X, axis=1)

    gamma1 = gamma[:, 0]
    gamma2 = gamma[:, 1]
    gamma3 = gamma[:, 2]

    eps = 1e-12

    denom1 = n_conditions * gamma1 + tau
    denom2 = n_conditions * gamma2 + tau
    denom3 = n_conditions * gamma3 + tau

    lambda1 = (x_sum + (alpha[0] - 1.0)) / np.clip(denom1, eps, None)
    lambda2 = (x_sum + (alpha[1] - 1.0)) / np.clip(denom2, eps, None)
    lambda3 = (x_sum + (alpha[2] - 1.0)) / np.clip(denom3, eps, None)

    lambda1 = np.clip(lambda1, 1e-6, 1e6)
    lambda2 = np.clip(lambda2, 1e-6, 1e6)
    lambda3 = np.clip(lambda3, 1e-6, 1e6)

    return [lambda1, lambda2, lambda3]


