# SAINT (Python)

Statistical Analysis of Interactomes with Explicit, Reproducible EM Models.  
A modular, research‑grade Python implementation of classical and hierarchical EM models for AP‑MS data, with tau‑grid tuning, diagnostics, and a fully explicit architecture.

---

## Package Structure

<details>
<summary><strong>Click to expand</strong></summary>

The Python package lives under ---python/saint/--- and is organized into explicit, modular components:

```text
saint/
    __init__.py

    diagnostics/
        diagnostics.py
        diagnostics_classical.py
        diagnostics_hierarchical.py
        diagnostics_pipeline.py
        diagnostics_tau_grid.py
        parameter_trajectories.py
        __init__.py

    io/
        data_input.py
        __init__.py

    model/
        classical_em_wrapper.py
        classical_likelihood.py
        classical_responsibilities.py
        hierarchical_em_wrapper.py
        hierarchical_likelihood.py
        hierarchical_responsibilities.py
        hierarchical_tau_grid.py
        init_tau.py
        __init__.py

        updates/
            lambda_updates.py
            pi_updates.py
            tau_updates.py
            __init__.py

    pipeline/
        classical_saint.py
        helpers_hierarchical.py
        hierarchical_saint.py
        __init__.py

    tests/
        test_classical_em_wrapper.py
        test_classical_pipeline.py
        test_classical_responsibilities.py
        test_compute_loglikelihood_classical.py
        test_compute_loglikelihood_hierarchical.py
        test_data_input.py
        test_diagnostics_classical.py
        test_diagnostics_hierarchical.py
        test_extract_bait_matrix.py
        test_hierarchical_em_wrapper.py
        test_hierarchical_pipeline.py
        test_hierarchical_responsibilities.py
        test_package_integrity.py
        test_validate_hyperparams.py
        _roadmap_future_tests.md
        __init__.py

    validation/
        validate_hyperparams.py
        __init__.py
```

</details>

## Classical Model

<details>
<summary><strong>Click to expand</strong></summary>

### Two-component Poisson mixture

The classical model assumes each prey belongs to one of two latent components:

- Component 1 (background) — Poisson rate lambda1  
- Component 2 (signal) — Poisson rate lambda2  

Mixing proportions:

- pi1, pi2

Posterior membership probabilities:

- gamma1, gamma2 per prey

The classical EM wrapper estimates:

- lambda1, lambda2  
- pi1, pi2  
- gamma1, gamma2  

All intermediate variables use explicit, component-indexed names for clarity and symmetry with the hierarchical model.

### Classical diagnostics include:

- log-likelihood trajectory  
- lambda1 and lambda2 trajectories  
- pi1 and pi2 trajectories  
- gamma1 and gamma2 distributions  
- parameter trajectory plots via parameter_trajectories.py  

</details>

## Hierarchical Model

<details>
<summary><strong>Click to expand</strong></summary>

### Three-component Poisson–Gamma mixture with empirical Bayes updates

The hierarchical model extends the classical model by introducing a third latent component and Gamma priors on the Poisson rate parameters. Each prey belongs to one of three components:

- Component 1 (background) — Poisson rate lambda1  
- Component 2 (contaminant) — Poisson rate lambda2  
- Component 3 (signal) — Poisson rate lambda3  

Each component has a Gamma(a_k, b_k) prior on lambda_k, with hyperparameters updated via empirical Bayes.

Mixing proportions:

- pi1, pi2, pi3

Posterior membership probabilities:

- gamma1, gamma2, gamma3 per prey

The hierarchical EM wrapper estimates:

- lambda1, lambda2, lambda3  
- pi1, pi2, pi3  
- gamma1, gamma2, gamma3  
- alpha1, alpha2, alpha3  
- a1, a2, a3 and b1, b2, b3  

The global shrinkage strength is controlled by:

- tau (a scalar hyperparameter)

### Hierarchical diagnostics include:

- log-likelihood trajectory  
- lambda1–3 trajectories  
- pi1–3 trajectories  
- alpha1–3, a1–3, b1–3 trajectories  
- gamma1–3 distributions  
- shrinkage behavior under different tau values  
- parameter trajectory visualizations  

These outputs help assess whether the hierarchical model is stable, well-calibrated, and biologically interpretable.

</details>

## Tau-grid Tuning

<details>
<summary><strong>Click to expand</strong></summary>

### Selecting the global shrinkage strength

The hierarchical model includes a global hyperparameter:

- tau (controls the strength of Gamma shrinkage on lambda1–3)

Choosing tau affects:

- the degree of shrinkage applied to lambda estimates  
- separation between background, contaminant, and signal components  
- stability of gamma3 (signal occupancy)  
- overall robustness under sparse or noisy data  

### The tau-grid tuner performs:

1. Running the hierarchical EM model across a user-defined grid of tau values.  
2. Recording final log-likelihood, lambda1–3, pi1–3, and mean gamma3 for each tau.  
3. Identifying the tau that optimizes model fit and stability.  
4. Producing a four-panel diagnostic figure summarizing the grid search.

### Diagnostics include:

- tau vs final log-likelihood  
- tau vs lambda1–3  
- tau vs pi1–3  
- tau vs mean gamma3  

These diagnostics help determine whether the model is under-shrunk, over-shrunk, or well-calibrated.

</details>

## Diagnostics Suite

<details>
<summary><strong>Click to expand</strong></summary>

### Three layers of model evaluation and visualization

The diagnostics subpackage provides a comprehensive set of tools for evaluating model behavior, convergence, and biological interpretability. Diagnostics are organized into three conceptual layers:

---

### Classical diagnostics

These diagnostics evaluate the two-component Poisson mixture model and include:

- log-likelihood trajectory  
- lambda1 and lambda2 trajectories  
- pi1 and pi2 trajectories  
- gamma1 and gamma2 distributions  
- parameter trajectory visualizations via parameter_trajectories.py  

These outputs help verify that the classical model converges cleanly and separates background from signal.

---

### Hierarchical diagnostics

These diagnostics evaluate the three-component Poisson–Gamma mixture model and include:

- log-likelihood trajectory  
- lambda1–3 trajectories  
- pi1–3 trajectories  
- alpha1–3, a1–3, b1–3 hyperparameter trajectories  
- gamma1–3 distributions  
- shrinkage behavior under different tau values  
- parameter trajectory visualizations  

These outputs help assess whether the hierarchical model is stable, well-calibrated, and biologically interpretable.

---

### Pipeline-level diagnostics

These diagnostics summarize results across all baits and include:

- density of mean gamma3 across baits  
- optional highlight of a positive control bait  
- summary tables of mean gamma3 per bait  
- cross-bait comparisons for signal occupancy  

These outputs provide a global view of model performance across the entire interactome.

---

### Diagnostic philosophy

All diagnostics follow the same architectural principles:

- explicit, component-indexed naming  
- bait-explicit titles and labels  
- no hidden defaults  
- reproducible, transparent visualizations  
- separation of classical, hierarchical, tau-grid, and pipeline layers  

This ensures that every diagnostic figure or table is interpretable by collaborators and future maintainers.

</details>

## Naming Conventions

<details>
<summary><strong>Click to expand</strong></summary>

### Explicit, symmetric, algebraic naming

All parameters and hyperparameters use explicit, component-indexed names to ensure clarity, symmetry, and interpretability:

- lambda1, lambda2, lambda3  
- pi1, pi2, pi3  
- gamma1, gamma2, gamma3  
- alpha1, alpha2, alpha3  
- a1, a2, a3  
- b1, b2, b3  

Intermediate variables also follow explicit naming:

- posterior_loglik_components  
- stabilized_log_posterior_components  
- posterior_membership_probabilities  
- posterior_weighted_counts_componentk  
- effective_membership_componentk  

No shorthand is used anywhere in the package. Every variable name reflects the underlying algebraic quantity.

</details>

---

## Reproducibility and Architecture

<details>
<summary><strong>Click to expand</strong></summary>

### Transparent, modular, and future-proof

The SAINT Python package enforces a strict architectural philosophy:

- explicit separation of model parameters, hyperparameters, and control variables  
- explicit component-wise naming across all layers  
- bait-explicit labels and titles in all diagnostics  
- no hidden defaults or implicit behavior  
- reproducible initialization and deterministic EM updates  
- clear segmentation using ---# %%--- markers for IDE navigation  
- modular separation of classical, hierarchical, tau-grid, diagnostics, and pipeline layers  
- no cross-layer leakage of logic or state  

This architecture ensures that:

- collaborators can understand and extend the code without inference  
- future maintainers can trace every computation step  
- scientific results are reproducible and transparent  
- the package remains stable as new models or diagnostics are added  

</details>

## Workflow and Usage

<details>
<summary><strong>Click to expand</strong></summary>

### End-to-end workflow

A standard analysis using the SAINT Python package proceeds through the following stages:

1. Prepare input data using the io layer (data_input.py).  
2. Run the classical or hierarchical EM model for each bait using the pipeline layer.  
3. If using the hierarchical model, perform tau-grid tuning to select the optimal shrinkage strength.  
4. Generate diagnostics at the classical, hierarchical, tau-grid, and pipeline levels.  
5. Aggregate results into a final scoring table or ranked interactome.

Each stage is modular and can be run independently or integrated into a full pipeline.

---

### Minimal usage example

The following pseudocode illustrates a typical usage pattern.

```python
from saint.io.data_input import load_data
from saint.pipeline.classical_saint import run_classical_saint
from saint.pipeline.hierarchical_saint import run_hierarchical_saint
from saint.model.hierarchical_tau_grid import run_tau_grid
from saint.diagnostics.diagnostics_pipeline import plot_gamma3_density

# Load and prepare data
bait_matrix = load_data("path/to/input.csv")

# Classical model
classical_results = run_classical_saint(bait_matrix)

# Hierarchical model
hierarchical_results = run_hierarchical_saint(bait_matrix)

# Tau-grid tuning
tau_grid_results = run_tau_grid(input_data=bait_matrix, tau_grid=[0.1, 0.5, 1.0, 2.0])

# Pipeline-level diagnostics
plot_gamma3_density(hierarchical_results)
```

This example demonstrates the modularity of the package: each layer (io, model, pipeline, diagnostics) is cleanly separated and can be used independently.

---

### Notes on data formatting

The io layer expects:

- a bait-by-prey count matrix  
- integer-valued spectral counts  
- consistent bait naming  
- no missing values  

The data_input.py module includes helper functions for:

- validating input structure  
- extracting bait-specific matrices  
- mapping identifiers  
- preparing data for classical or hierarchical models  

</details>

## Development Status

<details>
<summary><strong>Click to expand</strong></summary>

### Current state of the Python package

The SAINT Python package is under active development and includes:

- complete classical and hierarchical EM wrappers  
- full likelihood, responsibility, and update modules  
- tau-grid tuning functionality  
- a comprehensive diagnostics suite  
- a modular pipeline layer for end-to-end analysis  
- a robust test suite covering all major components  
- explicit naming conventions and reproducible architecture  

The package is suitable for research-grade analysis and is designed for long-term maintainability.

</details>

---

## Future Improvements

<details>
<summary><strong>Click to expand</strong></summary>

Planned enhancements include:

- renaming validate_hyperparams.py to a more explicit and descriptive name  
- expanding the diagnostics suite with additional cross-bait summaries  
- adding more integration tests for pipeline-level workflows  
- refining the tau-grid tuner for adaptive grid selection  
- improving documentation and adding example notebooks  
- expanding the io layer for additional data formats  

These improvements will further strengthen reproducibility, clarity, and usability.

</details>

---

## Contributing

<details>
<summary><strong>Click to expand</strong></summary>

Contributions are welcome. The repository includes:

- a full test suite under ---saint/tests/---  
- modular code organized by function and model layer  
- explicit naming conventions and architectural guidelines  

Contributors should follow the existing structure and naming discipline to maintain consistency across the package.

</details>

---

## License

<details>
<summary><strong>Click to expand</strong></summary>

This project is distributed under an open-source license. See the top-level repository for details.

</details>
