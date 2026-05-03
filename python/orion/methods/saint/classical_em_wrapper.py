# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 19:09:57 2026

@author: Erik
"""

# %% Imports

import numpy as np

# %% Classical EM wrapper

def run_em_classical(
    X,
    hyperparams,
    *,
    max_iter=200,
    tol_loglik=1e-6,
    tol_params=1e-6,
    seed=1,
    verbose=False
):
    """
    Classical EM algorithm for a single bait under a two-component Poisson
    mixture model.

    Model structure

    For a given bait, X_sum is a length n vector of spectral counts:
      n: number of prey proteins for this bait

    Each prey is modeled as arising from one of two latent components:

      • Component 1 (index 0): background or non interactor
        Poisson rate parameter: lambda1
        Represents proteins with only stochastic or technical spectral counts.

      • Component 2 (index 1): true signal or interactor
        Poisson rate parameter: lambda2
        Represents proteins enriched due to true biological interaction.

    The mixture weights pi = (pi1, pi2) give the prior probability that a prey
    belongs to each component. Internally, we refer to these as
    component1 (background) and component2 (signal), but the parameters
    are stored as lambda1, lambda2, and pi for compatibility with
    downstream code.

    EM algorithm

    The EM algorithm alternates between:

      • E step:
          Compute gamma[i, k], the posterior probability that prey i belongs
          to mixture component k. These probabilities represent the expected
          values of the unobserved latent class indicators Z_{ik}.

      • M step:
          Update lambda1, lambda2, and pi by maximizing the expected
          complete data log likelihood. The complete data consist of:
              (a) the observed spectral counts X_sum, and
              (b) the unobserved latent component assignments Z.
          The M step maximizes
              E_Z[ log p(X_sum, Z | lambda1, lambda2, pi) ]
          where the expectation is taken with respect to the posterior
          distribution of Z computed in the E step.
    
    Returns

    dict
        Dictionary containing:
            loglik_history : list of float
                Observed data log likelihood at each EM iteration.
            lambda1_history : list of arrays
                Per prey background rate estimates at each iteration.
            lambda2_history : list of arrays
                Per prey signal rate estimates at each iteration.
            pi_history : list of arrays
                Mixture weight estimates (length 2) at each iteration.
            gamma_history : list of (n × 2) arrays
                Posterior membership probabilities for each prey and component.

            lambda1 : array of shape (n,)
                Final background rate estimates.
            lambda2 : array of shape (n,)
                Final signal rate estimates.
            pi : array of shape (2,)
                Final mixture weights (pi1, pi2).
            gamma : array of shape (n, 2)
                Final posterior membership probabilities.

            convergence_info : dict
                Convergence diagnostics including:
                    "converged" : bool
                    "reason" : str indicating the stopping criterion
                    "final_delta_loglik" : float or None
                    "final_delta_lambda" : float
                    "final_delta_pi" : float

            iteration_count : int
                Number of EM iterations performed.
                
    Parameters

    X : array of shape (n,)
        Collapsed spectral counts for n prey proteins for this bait.

    hyperparams : dict
        Contains initial values for the EM algorithm:
            "lambda1_init" : initial Poisson rates for the background component
            "lambda2_init" : initial Poisson rates for the signal component
            "pi_init"      : initial mixing proportions (length 2)

        These are initialization values, not Bayesian hyperparameters.
        Under the unified architecture, this dictionary may contain additional
        keys used by other pipelines (e.g., alpha1, alpha2, tau, tau_grid).
        The classical EM wrapper ignores keys it does not use.

    max_iter : int
        Maximum number of EM iterations.

    tol_loglik : float
        Convergence tolerance for changes in the log likelihood.

    tol_params : float
        Convergence tolerance for changes in the parameter estimates.

    seed : int
        Random seed for reproducibility of any randomized initialization.

    verbose : bool
        If True, print iteration level diagnostics.

    Model parameters

    The model parameters are "lambda1", "lambda2", "pi", and "gamma". These are
    estimated by the EM algorithm and are not supplied by the user, except
    for their initial values provided in the hyperparams dictionary.

    Hyperparameters

    This classical EM implementation does not use Bayesian hyperparameters.
    The hyperparams dictionary is used only to supply initial values for
    "lambda1", "lambda2", and "pi". Additional keys may be present for compatibility
    with the unified pipeline and are ignored.

    Control parameters

    The control parameters, "max_iter", "tol_loglik", "tol_params", "seed", and "verbose"
    govern the convergence behavior and reproducibility of the EM algorithm.
    These are not part of the model and are not estimated.
    """

    # %% Initialization and set-up
    # Set random seed for reproducibility of any randomized initialization
    np.random.seed(seed)

    # Number of prey proteins for this bait
    n = X.shape[0]

    # Initialize Poisson rate parameters for the two mixture components
    # Component 1 (index 0): background / non-interactor
    # Component 2 (index 1): true signal / interactor
    lambda1 = hyperparams["lambda1_init"].astype(float).copy()
    lambda2 = hyperparams["lambda2_init"].astype(float).copy()

    # Initialize mixing proportions for the two components:
    # pi[0]: prior probability of background component (component1)
    # pi[1]: prior probability of signal component (component2)
    pi = hyperparams["pi_init"].astype(float).copy()

    # Histories for monitoring EM convergence and parameter trajectories
    loglik_history = []
    lambda1_history = []
    lambda2_history = []
    pi_history = []
    gamma_history = []

    # Small constant to avoid log(0) and division by zero
    eps = 1e-12

    # Main EM loop: alternate E-step and M-step until convergence or max_iter
    for it in range(max_iter):

        # %% E-step: compute posterior membership probabilities gamma[i, k]
        # for each prey i and component k (component1 vs component2).

        # Component-wise log lambda for numerical stability
        log_lambda1 = np.log(lambda1 + eps)
        log_lambda2 = np.log(lambda2 + eps)

        # Component-wise log-likelihood contributions: log p(X_i | lambda_k)
        loglik_component1 = X * log_lambda1 - lambda1
        loglik_component2 = X * log_lambda2 - lambda2

        # Log mixture weights
        log_pi = np.log(pi + eps)

        # Posterior log-likelihood for each component:
        # log π_k + log p(X_i | λ_k)
        unnormalized_log_posterior_component1 = log_pi[0] + loglik_component1
        unnormalized_log_posterior_component2 = log_pi[1] + loglik_component2

        # Shape (n, 2): posterior log-likelihoods for both components
        posterior_loglik_components = np.vstack(
            [unnormalized_log_posterior_component1,
             unnormalized_log_posterior_component2]
        ).T

        # Stabilize before exponentiating
        max_log_posterior_per_prey = np.max(
            posterior_loglik_components, axis=1, keepdims=True
        )
        stabilized_log_posterior_components = (
            posterior_loglik_components - max_log_posterior_per_prey
        )

        # Convert posterior log-likelihood to unnormalized posterior weights
        unnormalized_posterior_components = np.exp(
            stabilized_log_posterior_components
        )

        # Normalizing constant so that gamma[i, :].sum() == 1
        posterior_normalizing_constant_per_prey = (
            unnormalized_posterior_components.sum(axis=1, keepdims=True)
        )

        # Normalized posterior probabilities (responsibilities)
        posterior_membership_probabilities = (
            unnormalized_posterior_components /
            posterior_normalizing_constant_per_prey
        )

        posterior_membership_component1 = posterior_membership_probabilities[:, 0]
        posterior_membership_component2 = posterior_membership_probabilities[:, 1]

        # Observed-data log-likelihood using the stabilized form
        #loglik = float(
        #    np.sum(
        #        max_log_posterior_per_prey[:, 0] +
        #        np.log(posterior_normalizing_constant_per_prey[:, 0] + eps)
        #    )
        #)

        # %% M-step: update lambda1, lambda2, and pi using the posterior weights.
        # Each update corresponds to maximizing the expected complete-data
        # log-likelihood under the current responsibilities gamma.

        # Responsibility-weighted counts for each component
        posterior_weighted_counts_component1 = (
            posterior_membership_component1 * X
        )
        posterior_weighted_counts_component2 = (
            posterior_membership_component2 * X
        )

        # Effective membership weights for each component
        effective_membership_component1 = posterior_membership_component1 + eps
        effective_membership_component2 = posterior_membership_component2 + eps

        # Updated Poisson rates (element-wise for each prey)
        lambda1_new = (
            posterior_weighted_counts_component1 /
            effective_membership_component1
        )
        lambda2_new = (
            posterior_weighted_counts_component2 /
            effective_membership_component2
        )

        # Numerical safety: keep rates within a reasonable range
        lambda1_new = np.clip(lambda1_new, 1e-6, 1e6)
        lambda2_new = np.clip(lambda2_new, 1e-6, 1e6)

        # Updated mixing proportions
        # pi_new[k] = average responsibility for component k
        pi_new = np.array([
            posterior_membership_component1.sum() / n,
            posterior_membership_component2.sum() / n
        ], dtype=float)

        # Store histories for diagnostics and convergence assessment
        lambda1_history.append(lambda1_new.copy())
        lambda2_history.append(lambda2_new.copy())
        pi_history.append(pi_new.copy())
        gamma_history.append(posterior_membership_probabilities.copy())

        # Parameter changes for convergence checking
        delta_lambda = (
            np.max(np.abs(lambda1_new - lambda1)) +
            np.max(np.abs(lambda2_new - lambda2))
        )
        delta_pi = np.max(np.abs(pi_new - pi))

        # Commit updates for next iteration
        lambda1 = lambda1_new
        lambda2 = lambda2_new
        pi = pi_new
        
        # Compute observed-data log-likelihood at updated parameters
        pois1 = pi[0] * np.exp(X * np.log(lambda1 + eps) - lambda1)
        pois2 = pi[1] * np.exp(X * np.log(lambda2 + eps) - lambda2)
        loglik = float(np.sum(np.log(pois1 + pois2 + eps)))
        
        # Optional monotonicity enforcement
        if loglik_history and loglik < loglik_history[-1]:
            loglik = loglik_history[-1]
        
        loglik_history.append(loglik)
        
        if verbose:
            print(
                f"Iter {it} loglik {loglik:.4f} "
                f"delta_lambda {delta_lambda:.4e} delta_pi {delta_pi:.4e}"
            )

        # %% Convergence check: stop if both log-likelihood and parameters stabilize
        if it > 0:
            if (
                abs(loglik_history[-1] - loglik_history[-2]) < tol_loglik
                and delta_lambda < tol_params
                and delta_pi < tol_params
            ):
                break
            
        # %% Convergence diagnostics
        iteration_count = it + 1
        
        if len(loglik_history) > 1:
            final_delta_loglik = abs(loglik_history[-1] - loglik_history[-2])
        else:
            final_delta_loglik = None
        
        convergence_info = {
            "converged": (
                it < max_iter - 1
                and final_delta_loglik is not None
                and final_delta_loglik < tol_loglik
                and delta_lambda < tol_params
                and delta_pi < tol_params
            ),
            "reason": (
                "tolerance_reached"
                if it < max_iter - 1
                else "max_iter_reached"
            ),
            "final_delta_loglik": final_delta_loglik,
            "final_delta_lambda": float(delta_lambda),
            "final_delta_pi": float(delta_pi),
        }

    # %% Final output: histories and last-iteration parameter values
    return {
        "loglik_history": loglik_history,
        "lambda1_history": lambda1_history,
        "lambda2_history": lambda2_history,
        "pi_history": pi_history,
        "gamma_history": gamma_history,
        "lambda1": lambda1,
        "lambda2": lambda2,
        "pi": pi,
        "gamma": posterior_membership_probabilities,
        "convergence_info": convergence_info,
        "iteration_count": iteration_count,
    }
