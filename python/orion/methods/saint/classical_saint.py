# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 08:51:47 2026

@author: Erik
"""

# %% Imports

import numpy as np
import pandas as pd

from orion.methods.saint.classical_em_wrapper import run_em_classical
from orion.utils.io.data_input import load_bait_data
from orion.methods.saint.diagnostics_classical import make_classical_plots

# %% Classical SAINT pipeline

def run_classical_pipeline(
    input_data,
    hyperparams=None,
    make_plots=False,
    plot_dir=None,
    max_iter=200,
    tol_loglik=1e-6,
    tol_params=1e-6,
    seed=1,
    verbose=False,
):
    """
    Run the classical SAINT pipeline. This function collapses replicate-level counts
    to a single total count per prey, runs classical EM for each bait_unit, and returns
    a unified results object containing raw EM outputs, typed metadata, and a
    prey-level results dataframe aligned with the hierarchical pipeline.

    Parameters

    input_data : dict or DataFrame
        Raw APMS data in the format expected by the data loader. The loader
        converts the wide-format APMS table into per-bait_unit replicate matrices
        and constructs metadata including the per-bait protein ordering.

    hyperparams : dict, optional
        Optional hyperparameters and initialization values for the classical EM
        model. The classical EM wrapper uses:
            "lambda1_init" : initial Poisson rates for the background component
            "lambda2_init" : initial Poisson rates for the signal component
            "pi_init"      : initial mixing proportions (length 2)
        Any additional entries are ignored by the classical EM wrapper but are
        preserved in the metadata for architectural symmetry.

    make_plots : bool
        If True, generate diagnostic plots for each bait.

    plot_dir : str or None
        Directory in which to save diagnostic plots. If None and make_plots is
        True, a default directory may be used by the plotting function.

    max_iter : int
        Maximum number of EM iterations.

    tol_loglik : float
        Convergence tolerance for changes in the log likelihood.

    tol_params : float
        Convergence tolerance for changes in the parameter estimates.

    seed : int
        Random seed for reproducibility of any randomized initialization.

    verbose : bool
        If True, print iteration-level diagnostics.

    Returns

    dict
        Unified results object containing:
            - raw_outputs:
                em_results : dict mapping each bait to the classical EM result
                tau_info   : dict mapping each bait to an empty tau info dict
            - metadata : typed metadata object mirroring the hierarchical pipeline
            - results_df : combined prey-level results table with one row per
              prey–bait_unit pair, containing:
                  Protein, bait,
                  lambda1, lambda2
                  tau,
                  pi1, pi2
                  gamma1, gamma2
    """

    # %% Load data via unified loader
    bait_unit_list, X_by_bait_unit, loader_metadata = load_bait_data(input_data)

    # %% Storage
    all_results = {}
    all_tau_info = {}
    all_convergence_info = {}
    all_iteration_counts = {}

    # Classical pipeline has no tau grid; keep empty for symmetry
    tau_grid = []

    # Initialize hyperparams container
    if hyperparams is None:
        hyperparams = {}

    # %% Per-bait classical EM
    rows = []

    for bait_unit in bait_unit_list:
        X = X_by_bait_unit[bait_unit]
        X_sum = X.sum(axis=1).astype(float)

        # Classical pipeline does not use biological_bait_names
        bait_unit = bait_unit

        # Mean level for initialization
        mean_level = max(X_sum.mean(), 1.0)

        # Per-bait_unit hyperparameters (copy to avoid cross-bait_unit mutation)
        hyperparams_bait_unit = dict(hyperparams)

        # Initialization values (only set if not already provided)
        hyperparams_bait_unit.setdefault(
            "lambda1_init",
            np.full(X_sum.shape[0], 0.5 * mean_level),
        )
        hyperparams_bait_unit.setdefault(
            "lambda2_init",
            np.full(X_sum.shape[0], 1.5 * mean_level),
        )
        hyperparams_bait_unit.setdefault(
            "pi_init",
            np.array([0.7, 0.3], dtype=float),
        )

        # Run classical EM
        results_em = run_em_classical(
            X=X_sum,
            hyperparams=hyperparams_bait_unit,
            max_iter=max_iter,
            tol_loglik=tol_loglik,
            tol_params=tol_params,
            seed=seed,
            verbose=verbose,
        )
        
        # Extract convergence diagnostics
        convergence_info = results_em.get("convergence_info", {})
        iteration_count = results_em.get("iteration_count", None)
        
        all_results[bait_unit] = results_em
        all_tau_info[bait_unit] = {}  # no tau grid in classical pipeline
        all_convergence_info[bait_unit] = convergence_info
        all_iteration_counts[bait_unit] = iteration_count
        
        # Optional plotting
        if make_plots:
            figs = make_classical_plots(results_em, bait_unit, plot_dir=plot_dir)
            all_results[bait_unit]["figures"] = figs
        
        # Build prey-level rows for this bait
        proteins = loader_metadata["proteins_by_bait_unit"][bait_unit]
        
        lambda1 = results_em["lambda1"]
        lambda2 = results_em["lambda2"]
        pi = results_em["pi"]
        gamma = results_em["gamma"]
        
        # Classical model has only two components
        pi1, pi2 = pi[0], pi[1]
        
        gamma1 = gamma[:, 0]
        gamma2 = np.full_like(gamma1, np.nan, dtype=float)
        
        df_bait_unit = pd.DataFrame({
            "Protein": proteins,
            "BaitUnit": bait_unit,
            "lambda1": lambda1,
            "lambda2": lambda2,
            "pi1": pi1,
            "pi2": pi2,
            "gamma1": gamma1,
            "gamma2": gamma2,
        })
        
        rows.append(df_bait_unit)
    
    # %% Build results_df
    results_df = pd.concat(rows, ignore_index=True)
    
    # Sort by Protein then gamma2 descending (analogous to hierarchical convention)
    results_df = results_df.sort_values(
        by=["Protein", "gamma2"],
        ascending=[True, False],
        ignore_index=True,
    )
    
    # %% Build the typed metadata
    
    from orion.utils.metadata_types import (
        HierarchicalMetadata,
        UserProvidedFields,
        InferredFields,
        PipelineDerivedFields,
    )
    
    # Construct typed metadata
    metadata_obj = HierarchicalMetadata(
        user_provided_fields=UserProvidedFields(
            biological_bait_names={},
            AN={},
            MW={},
            extra_fields={},
        ),
        inferred_fields=InferredFields(
            bait_units=loader_metadata.get("bait_units", []),
            proteins_by_bait_unit=loader_metadata.get("proteins_by_bait_unit", {}),
            control_bait_units=loader_metadata.get("control_bait_units", []),
            extra_fields={
                k: v
                for k, v in loader_metadata.items()
                if k not in {
                    "bait_units",
                    "proteins_by_bait_unit",
                    "control_bait_units",
                }
            },
        ),
        pipeline_derived_fields=PipelineDerivedFields(
            bait_units=loader_metadata.get("bait_units", []),
            hyperparameters=hyperparams,
            tau_grid=tau_grid,
            convergence=all_convergence_info,
            iteration_counts=all_iteration_counts,
        ),
    )
    
    # %% Final unified output
    return {
        "raw_outputs": {
            "em_results": all_results,
            "tau_info": all_tau_info,
        },
        "metadata": metadata_obj,
        "results_df": results_df.sort_values(
            by=["gamma2", "Protein"],
            ascending=[False, True],
            ignore_index=True,
        ),
    }