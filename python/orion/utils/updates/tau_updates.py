# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 12:47:28 2026

@author: Erik
"""

# %% Imports

# import numpy as np


# %% Tau update for hierarchical SAINT

def update_tau(
    X,
    gamma,
    lambda1,
    lambda2,
    lambda3
):
    """
    Update tau for hierarchical SAINT.

    tau is the shared rate parameter of the Gamma prior on each lambda (Poisson rate parameter for background, noise, and signal).
    Larger tau increases shrinkage of all lambda values toward zero.

    X has shape n_preys by n_conditions.
    gamma has shape n_preys by three (three latent variables).
    lambda1, lambda2, lambda3, are prey-specific Poisson rate vectors (again, 1, 2, and 3 represent latent "background", "noise", and "signal" components of the mixture).

    Returns a prey-specific tau vector.
    """

    n_preys, n_conditions = X.shape

    gamma1 = gamma[:, 0]
    gamma2 = gamma[:, 1]
    gamma3 = gamma[:, 2]

    # Expected counts under each component
    exp1 = n_conditions * lambda1
    exp2 = n_conditions * lambda2
    exp3 = n_conditions * lambda3

    # Weighted expected counts
    weighted = gamma1 * exp1 + gamma2 * exp2 + gamma3 * exp3

    # Tau update
    eps = 1e-12
    tau_new = weighted / (gamma1 + gamma2 + gamma3 + eps)

    return tau_new



