# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 18:20:10 2026

@author: Erik
"""

# %% Imports
import os
import pandas as pd

from orion.utils.io.data_input import load_bait_data
from orion.methods.iris.hierarchical_tau_grid import run_tau_grid
from orion.methods.iris.hierarchical_em_wrapper import run_em_hierarchical

from orion.methods.iris.diagnostics_hierarchical import (
    make_hierarchical_plots,
    plot_gamma3_density,
)
from orion.methods.iris.diagnostics_tau_grid import (
    diagnostics_tau_grid,
)


def run_hierarchical_pipeline(
    input_data,
    hyperparams=None,
    max_iter=200,
    tol_loglik=1e-6,
    tol_params=1e-6,
    seed=1,
    verbose=False,
    make_plots=False,
    show_plots=True,
    save_plots=False,
    plot_dir=None,
    save_results=False,
    results_csv=None,
):
    """
    Run the IRIS pipeline for all bait_units using a per-bait_unit tau grid
    search. This function accepts raw input_data from the user, converts it into
    per-bait_unit replicate matrices using the data loader, applies the hierarchical
    EM model with a two-stage tau grid for each bait_unit, selects the best-fitting
    tau, and assembles a unified prey-level results table.
    
    Hyperparameters may be omitted entirely. If the user supplies no dictionary,
    or supplies a dictionary missing some entries, all missing hyperparameters
    are filled with internal defaults before any access occurs. Only the `tau`
    entry is modified during the grid search; all other hyperparameters remain
    unchanged unless explicitly provided by the user. EM control parameters
    (max_iter, tol_loglik, tol_params, seed) are pipeline-level arguments and
    are not part of the hyperparameter dictionary.
    
    Parameters
    
    input_data : dict or DataFrame
        Raw IPMS data provided by the user. The pipeline calls the data loader
        to convert this into per-bait_unit replicate matrices (X_by_bait_unit) and to
        construct metadata including the per-bait_unit protein ordering.
    
    hyperparams : dict, optional
        Optional hyperparameters for the hierarchical EM model. Users may supply
        overrides, but any missing entries are filled with internal defaults
        before the model is run. The `tau` entry is overridden during the grid
        search; all other entries remain unchanged unless explicitly provided.
    
    max_iter : int
        Maximum number of EM iterations for each tau evaluation.
    
    tol_loglik : float
        Convergence tolerance for changes in the log-likelihood.
    
    tol_params : float
        Convergence tolerance for changes in the parameter estimates.
    
    seed : int
        Random seed used for reproducibility.
    
    verbose : bool
        If True, print progress information during the pipeline.
    
    make_plots : bool
        If True, generate diagnostic plots for each bait (EM diagnostics,
        tau-grid diagnostics, and gamma3 density).
    
    show_plots : bool
        If True and make_plots is True, display plots (e.g., in an IDE plot
        pane). Has no effect if make_plots is False.
    
    save_plots : bool
        If True and make_plots is True, save plots as .png files to plot_dir.
        Has no effect if make_plots is False.
    
    plot_dir : str or Path, optional
        Directory in which to save .png files when save_plots is True. If
        save_plots is True and plot_dir is None, a ValueError is raised.
    
    Returns
    
    dict
        Unified results object containing:
    
            - raw_outputs :
                em_results : dict mapping each bait to the best-fitting EM result
                tau_info   : dict mapping each bait to tau-grid diagnostics
                             (full tau grid, refinement grid, log-likelihoods,
                              convergence info, iteration counts)
    
            - metadata :
                Experiment metadata produced by the data loader, including
                protein ordering, replicate structure, and pipeline-level
                hyperparameters and control settings.
    
            - results_df :
                Combined prey-level results table with one row per prey–bait_unit
                pair, containing:
                    Protein, BaitUnit,
                    lambda1, lambda2, lambda3,
                    tau,
                    pi1, pi2, pi3,
                    gamma1, gamma2, gamma3
    
    Model parameters
    
    The model parameters lambda1, lambda2, lambda3, pi, and gamma are estimated
    by the EM algorithm and are not supplied by the user.
    
    Hyperparameters
    
    The hyperparameters alpha, a_k, b_k, and tau may be supplied in the
    hyperparameter dictionary but are not required. Any missing entries are
    filled with internal defaults before any access occurs. The `tau` entry is
    overridden by the grid search. All other hyperparameters remain unchanged
    unless explicitly provided.
    
    Tau grid tuning parameters
    
    The "coarse" and "refinement" grid ranges and resolutions define the search
    over tau. These are tuning parameters for the grid search procedure and are
    not model hyperparameters. The final evaluated tau grid is the union of the
    coarse and refinement grids, and the best tau is selected as the smallest
    tau whose final log-likelihood lies within a fixed tolerance of the maximum.
    """
    
    # %% Load data and metadata
    bait_unit_list, X_by_bait_unit, metadata = load_bait_data(input_data)
    
    # Fill missing hyperparameters with defaults
    if hyperparams is None:
        hyperparams = {}
    #hyperparams = metadata.fill_missing_hyperparams(hyperparams)
    
    em_results = {}
    tau_info = {}
    rows = []
    
    # %% Per-bait tau grid + EM refinement
    for bait_unit, X in X_by_bait_unit.items():
        
        # Skip negative-control baits in the results table
        if bait_unit in metadata["control_bait_units"]:
            continue

    
        # 1. Run tau grid
        tau_grid_result = run_tau_grid(
            X=X,
            hyperparams=hyperparams,
            bait_unit=bait_unit,
        )
        tau_info[bait_unit] = tau_grid_result
    
        # 2. Select best tau
        #taus = tau_grid_result["taus"]
        #logliks = tau_grid_result["logliks"]
        #best_tau = float(taus[int(np.argmax(logliks))])
        #best_tau = metadata.select_best_tau(taus, logliks)
        best_tau = tau_grid_result["best_tau"]
        best_em  = tau_grid_result["best_result"]
    
        # 3. Run EM wrapper at best tau
        hyperparams_best = hyperparams.copy()
        hyperparams_best["tau"] = best_tau
    
        best_em = run_em_hierarchical(
            X=X,
            hyperparams=hyperparams_best,
            bait_unit=bait_unit,
            max_iter=max_iter,
            tol_loglik=tol_loglik,
            tol_params=tol_params,
            seed=seed,
            verbose=verbose,
        )
        em_results[bait_unit] = best_em
    
        # 4. Optional per-bait diagnostics
        if make_plots:
            if save_plots and plot_dir is None:
                raise ValueError("plot_dir must be provided when save_plots=True.")
    
            # Hierarchical diagnostics (returns dict of figures)
            figs_hier = make_hierarchical_plots(best_em, bait_unit)
    
            # Tau-grid diagnostics (returns dict with figure)
            tau_results = tau_grid_result["tau_grid_results"]
            
            taus = []
            logliks = []
            em_results_tau = {}

            for tau, res in tau_results.items():
                llhist = res.get("loglik_history", [])
                if llhist:
                    taus.append(tau)
                    logliks.append(llhist[-1])
                    em_results_tau[tau] = res
            
            diag_input = {
                "taus": taus,
                "logliks": logliks,
                "em_results": em_results_tau,
            }
            
            diag_tau = diagnostics_tau_grid(diag_input, bait_unit, best_tau)
            fig_tau = diag_tau["figure"]

            # Show/save hierarchical diagnostics
            for name, fig in figs_hier.items():
                if show_plots:
                    fig.show()
                if save_plots:
                    os.makedirs(plot_dir, exist_ok=True)
                    fig.savefig(
                        os.path.join(plot_dir, f"{bait_unit}_hier_{name}.png"),
                        bbox_inches="tight",
                    )
    
            # Show/save tau-grid diagnostics
            if show_plots:
                fig_tau.show()
            if save_plots:
                fig_tau.savefig(
                    os.path.join(plot_dir, f"{bait_unit}_tau_grid.png"),
                    bbox_inches="tight",
                )
    
        # 5. Assemble rows for results_df
        #proteins = metadata.inferred_fields.proteins_by_bait[bait]
        proteins = metadata["proteins_by_bait_unit"][bait_unit]
        for i, protein in enumerate(proteins):
            rows.append(
                {
                    "Protein": protein,
                    "BaitUnit": bait_unit,
                    **{f"rep{j+1}": float(X[i, j]) for j in range(X.shape[1])},
                    "lambda1": best_em["lambda1"][i],
                    "lambda2": best_em["lambda2"][i],
                    "lambda3": best_em["lambda3"][i],
                    "tau": best_em["tau"],
                    "pi1": best_em["pi"][0],
                    "pi2": best_em["pi"][1],
                    "pi3": best_em["pi"][2],
                    "gamma1": best_em["gamma"][i, 0],
                    "gamma2": best_em["gamma"][i, 1],
                    "gamma3": best_em["gamma"][i, 2],
                }
            )
            
    # %% Assemble results_df
    results_df = pd.DataFrame(rows)
    
    # Enforce deterministic column order
    results_df = results_df[sorted(results_df.columns)]
    
    results_df = results_df.sort_values(
                                        ["Protein", "gamma3"],
                                        ascending=[True, False],
                                        ignore_index=True
                                    )
    
    # %% Global gamma3 density plot (only once, after results_df)
    if make_plots:
        ax_gamma3 = plot_gamma3_density(results_df)
    
        if show_plots:
            ax_gamma3.figure.show()
    
        if save_plots:
            os.makedirs(plot_dir, exist_ok=True)
            ax_gamma3.figure.savefig(
                os.path.join(plot_dir, "global_gamma3_density.png"),
                bbox_inches="tight",
            )
    
    # %% Unified return object
    output = {
        "raw_outputs": {
            "em_results": em_results,
            "tau_info": tau_info,
        },
        "metadata": metadata,
        "results_df": results_df,
    }
    
    if save_results and results_csv is not None:
        results_df.to_csv(results_csv, index=False)
    
    return output