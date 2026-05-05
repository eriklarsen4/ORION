# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 12:48:36 2026

@author: Erik
"""

# %% Imports
import numpy as np
from orion.methods.iris.em_wrapper import run_em

# %% Tau grid parameters
COARSE_MIN = -1.0
COARSE_MAX = 1.0
COARSE_POINTS = 7

REFINE_HALF_RANGE = 0.5
REFINE_POINTS = 7

# %% Tau grid search
def run_tau_grid(X, hyperparams, bait_unit):
    """
    Perform a two-stage grid search over the global shrinkage hyperparameter, tau,
    for a single bait_unit. The function evaluates a coarse log-spaced grid followed
    by a refinement grid centered on the best coarse-grid tau. For each tau value,
    the IRIS algorithm is run to convergence and all model outputs and
    diagnostics are recorded.
    
    Tau selection is performed only after both grid stages have been evaluated.
    Instead of choosing the tau with the single highest final log-likelihood
    (argmax), the function identifies all tau values whose final log-likelihood is
    within a small tolerance of the maximum, and selects the smallest tau among
    these near-optimal candidates. This prevents runaway selection of extreme 
    tau values that can arise from the interaction between hierarchical 
    shrinkage and the collapsed Poisson likelihood.

    
    Parameters
    
    X : array-like, shape (n_proteins, n_replicates)
        The data matrix for the current bait_unit. Each row corresponds to a protein
        and each column corresponds to a replicate measurement.
    
    hyperparams : dict
        Hyperparameters for the IRIS model. Any missing entries are
        filled internally by the EM wrapper. Only the "tau" entry is modified
        during the grid search; all other hyperparameters remain unchanged.
    
    bait_unit : str
        Unique identifier for the current modeling unit (condition + bait_name). Used only for labeling and diagnostics.
    
    Returns
    
    dict
        A dictionary containing the following entries:
    
        best_tau : float
            The tau value selected by the near-optimal rule: the smallest tau
            whose final log-likelihood lies within a fixed tolerance of the
            maximum log-likelihood across all evaluated tau values.
    
        best_result : dict
            The full EM output corresponding to best_tau. This includes:
                - lambda1, lambda2, lambda3 : component means
                - pi : mixture proportions
                - gamma : posterior component probabilities
                - loglik_history : list of log-likelihood values per EM iteration
                - convergence_info : dict describing EM convergence diagnostics
                - iteration_count : number of EM iterations performed
    
        tau_grid_results : dict
            A mapping from each tau value (float) to the full EM result for that
            tau. Each entry contains the same fields as best_result.
    
        convergence_info : dict
            A mapping from each tau value to the EM convergence diagnostics
            returned by the EM wrapper. This allows inspection of EM behavior
            across the entire grid.
    
        iteration_counts : dict
            A mapping from each tau value to the number of EM iterations
            performed for that tau.
    
        tau_grid : list of float
            The complete set of tau values evaluated across both the coarse and
            refinement grids, sorted in ascending order.
    """
    
    # %% Stage 1: coarse grid
    coarse_grid = np.logspace(COARSE_MIN, COARSE_MAX, num=COARSE_POINTS)
    
    tau_grid_results = {}
    convergence_info = {}
    iteration_counts = {}
    
    # --- OLD incremental argmax selection (commented out) ---
    # best_tau = None
    # best_loglik = -np.inf
    # best_result = None
    
    for tau in coarse_grid:
        hyper = dict(hyperparams)
        hyper["tau"] = float(tau)
    
        result = run_em(
            X=X,
            hyperparams=hyper,
            bait_unit=bait_unit,
            max_iter=hyper.get("max_iter", 100),
            tol_loglik=hyper.get("tol_loglik", 1e-6),
            tol_params=hyper.get("tol_params", 1e-6),
            seed=hyper.get("seed", None),
        )
        
        #final_loglik = result["loglik_history"][-1]
        #mean_gamma3  = result["gamma"][:, 2].mean()
        #print(f"tau={tau:.4g}  loglik={final_loglik:.3f}  mean_gamma3={mean_gamma3:.3f}")
    
        tau_grid_results[tau] = result
    
        convergence_info[tau] = result.get("convergence_info", {})
        iteration_counts[tau] = result.get("iteration_count", None)
    
        # --- OLD incremental argmax logic (commented out) ---
        # final_loglik = result["loglik_history"][-1]
        # if final_loglik > best_loglik:
        #     best_loglik = final_loglik
        #     best_tau = tau
        #     best_result = result
        
    # %% Stage 2: refinement grid
    # Compute coarse-grid best tau safely
    coarse_ll = {}
    for tau, res in tau_grid_results.items():
        llhist = res.get("loglik_history", [])
        if llhist:
            coarse_ll[tau] = llhist[-1]
    
    # If no valid log-likelihoods exist, abort cleanly
    if not coarse_ll:
        return {
            "best_tau": None,
            "best_result": None,
            "tau_grid_results": tau_grid_results,
            "convergence_info": convergence_info,
            "iteration_counts": iteration_counts,
            "tau_grid": sorted(tau_grid_results.keys()),
        }
    
    coarse_best_tau = max(coarse_ll, key=coarse_ll.get)
    center = np.log10(coarse_best_tau)
    
    refine_grid = np.logspace(
        center - REFINE_HALF_RANGE,
        center + REFINE_HALF_RANGE,
        num=REFINE_POINTS,
    )
    
    for tau in refine_grid:
        hyper = dict(hyperparams)
        hyper["tau"] = float(tau)
    
        result = run_em(
            X=X,
            hyperparams=hyper,
            bait_unit=bait_unit,
            max_iter=hyper.get("max_iter", 100),
            tol_loglik=hyper.get("tol_loglik", 1e-6),
            tol_params=hyper.get("tol_params", 1e-6),
            seed=hyper.get("seed", None),
        )
    
        tau_grid_results[tau] = result
        convergence_info[tau] = result.get("convergence_info", {})
        iteration_counts[tau] = result.get("iteration_count", None)
    
        # --- OLD incremental argmax logic (commented out) ---
        # final_loglik = result["loglik_history"][-1]
        # if final_loglik > best_loglik:
        #     best_loglik = final_loglik
        #     best_tau = tau
        #     best_result = result
        
    # %% NEW best-tau selection rule (local, no imports)
    taus = []
    logliks = []
    for tau, res in tau_grid_results.items():
        llhist = res.get("loglik_history", [])
        if llhist:
            taus.append(tau)
            logliks.append(llhist[-1])
    
    taus = np.array(taus)
    logliks = np.array(logliks)
    
    # If no valid log-likelihoods exist, abort cleanly
    if len(taus) == 0:
        return {
            "best_tau": None,
            "best_result": None,
            "tau_grid_results": tau_grid_results,
            "convergence_info": convergence_info,
            "iteration_counts": iteration_counts,
            "tau_grid": sorted(tau_grid_results.keys()),
        }
    
    max_ll = logliks.max()
    eps = 1.0  # tolerance for near-optimal taus
    
    # mask of taus within eps of max loglik
    mask = logliks >= max_ll - eps
    
    # choose smallest tau among near-optimal candidates
    best_tau = float(taus[mask].min())
    best_result = tau_grid_results[best_tau]
    
    tau_grid = sorted(tau_grid_results.keys())
    
    # %% Return results
    return {
        "best_tau": best_tau,
        "best_result": best_result,
        "tau_grid_results": tau_grid_results,
        "convergence_info": convergence_info,
        "iteration_counts": iteration_counts,
        "tau_grid": tau_grid,
    }