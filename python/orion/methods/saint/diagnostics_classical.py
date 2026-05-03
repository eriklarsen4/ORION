# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 18:52:18 2026

@author: Erik
"""

# %% Imports

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# %% Classical diagnostics

def make_classical_plots(results_em, bait_name):
    """
    Create diagnostic plots for the classical SAINT EM results.
    This function returns a dictionary of matplotlib Figure objects.
    Figures will automatically appear in IDEs such as Spyder or VS Code.

    The diagnostics include:
        - Log-likelihood trajectory
        - Lambda1 and lambda2 trajectories
        - Mixture weight (pi1, pi2) trajectories
        - Gamma1 and gamma2 distributions

    All plot titles explicitly include the bait name.
    """

    sns.set_style("whitegrid")

    figs = {}

    # %% Log likelihood trajectory
    fig_loglik = plt.figure()
    loglik_history = results_em["loglik_history"]
    plt.plot(loglik_history, marker="o")
    plt.xlabel("Iteration")
    plt.ylabel("Log likelihood")
    plt.title(f"Classical SAINT log likelihood trajectory for {bait_name}")
    figs["loglik"] = fig_loglik

    # %% Lambda trajectories
    fig_lambda = plt.figure()
    lambda1_history = np.array(results_em["lambda1_history"])
    lambda2_history = np.array(results_em["lambda2_history"])

    mean_lambda1_per_iteration = lambda1_history.mean(axis=1)
    mean_lambda2_per_iteration = lambda2_history.mean(axis=1)

    plt.plot(mean_lambda1_per_iteration, label="lambda1 mean")
    plt.plot(mean_lambda2_per_iteration, label="lambda2 mean")
    plt.xlabel("Iteration")
    plt.ylabel("Mean lambda value")
    plt.title(f"Classical SAINT lambda trajectories for {bait_name}")
    plt.legend()
    figs["lambda"] = fig_lambda

    # %% Pi trajectory
    fig_pi = plt.figure()
    pi_history = np.array(results_em["pi_history"])

    mixture_weight_background = pi_history[:, 0]
    mixture_weight_signal = pi_history[:, 1]

    plt.plot(mixture_weight_background, label="pi1 (background)")
    plt.plot(mixture_weight_signal, label="pi2 (signal)")
    plt.xlabel("Iteration")
    plt.ylabel("Mixture weight")
    plt.title(f"Classical SAINT pi trajectories for {bait_name}")
    plt.legend()
    figs["pi"] = fig_pi

    # %% Gamma distributions
    fig_gamma = plt.figure()
    posterior_membership_final = results_em["gamma"]

    plt.hist(posterior_membership_final[:, 0], bins=30, alpha=0.5, label="gamma1 (background)")
    plt.hist(posterior_membership_final[:, 1], bins=30, alpha=0.5, label="gamma2 (signal)")
    plt.xlabel("Posterior membership probability")
    plt.ylabel("Count")
    plt.title(f"Classical SAINT gamma distributions for {bait_name}")
    plt.legend()
    figs["gamma"] = fig_gamma

    return figs

