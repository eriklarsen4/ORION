# ORION (Observed Relationship Inference of Organized Networks) Python package

The ORION Python package implements two related probabilistic models for analyzing protein-protein interaction signal in bait‑specific experiments:

- SAINT — Significance Analysis of INTeractome  
- IRIS — Inference of Regularized Interaction Signals  

SAINT provides a classical two‑component Poisson mixture baseline.  
IRIS is a hierarchical, regularized three‑component Poisson mixture model that produces stabilized, biologically meaningful per‑prey posterior probabilities.

This README gives both a modeler‑level mathematical description and a biologist/user‑level interpretation.

---

## Package layout
```text
    python/
    └── orion/
        ├── methods/
        │   ├── saint/
        │   └── iris/
        └── utils/
```
You typically interact with:

- `orion.methods.saint` — SAINT baseline model  
- `orion.methods.iris` — IRIS regularized three‑component model  
- `orion.utils` — shared utilities (I/O, summaries, etc.)

---

## Installation

From the repo root:

```python
    pip install -e python
```

Then in Python:

```python
    from orion.methods.saint import saint
    from orion.methods.iris import iris
```

---

## SAINT: Significance Analysis of INTeractome

SAINT models observed spectral counts for a given bait‑specific experiment as arising from a two‑component Poisson mixture:


SAINT models observed spectral counts for a given bait‑specific experiment as arising from a two‑component Poisson mixture:

$$
p(x_i \mid \Theta) = \pi_1 \,\text{Pois}(x_i \mid \lambda_1)
                   + \pi_2 \,\text{Pois}(x_i \mid \lambda_2)
$$




where:

- $$\(\lambda_1\)$$ is the background mean  
- $$\(\lambda_2\)$$ is the enriched mean  
- $$\(\pi_1, \pi_2\)$$ are mixture weights with $$\(\pi_1 + \pi_2 = 1\)$$

SAINT produces continuous posterior probabilities for each prey belonging to background or enriched.  
However, the **two‑component structure, itself,** introduces several limitations that affect interpretability and stability:

- **Ambiguous signal distorts both components**  
  Even though SAINT outputs continuous probabilities, it must explain all observations using only two latent ("background"/non-interaction & "enriched"/signal/true-interaction) states.  
  Borderline or inconsistent counts therefore pull the background mean upward and the enriched mean downward, reducing the separation between them.  
  This makes the posterior probability of belonging to the signal component less stable and less interpretable.

- **No hierarchical shrinkage**  
  SAINT estimates $$\(\lambda_1\)$$ and $$\(\lambda_2\)$$ independently for each experiment.  
  Experiments with sparse counts, replicate imbalance, or outliers can produce extreme or biologically implausible component means.

- **Unregularized mixture weights**  
  The mixture weights $$\(\pi_1, \pi_2\)$$ are updated purely from the data.  
  In low‑information settings, one component can collapse to nearly zero weight, producing degenerate fits and unstable posteriors.

- **Sensitivity to outliers and replicate imbalance**  
  A single high count can inflate $$\(\lambda_2\)$$, distorting responsibilities for all other preys.  
  Without shrinkage or an intermediate component, SAINT cannot buffer against these distortions.

IRIS is designed specifically to address these structural limitations.

---

## IRIS: Inference of Regularized Interaction Signals

IRIS is a three‑component hierarchical Poisson mixture model fit independently to each bait‑unit*.
Its primary purpose is to produce **stabilized, biologically meaningful per‑prey posterior probabilities** that are less distorted by the behavior of other preys, serving to more clearly prioritize follow-up candidate interactors.  
To achieve this, IRIS uses hierarchical priors (Gamma and Dirichlet) as regularizers within a maximum-likelihood EM-style alternating procedure to produce point estimates for parameters rather than full Bayesian posterior distributions.
IRIS achieves this with:

- fixed global Gamma shrinkage on component means  
- Dirichlet‑stabilized mixture weights  
- three components: background, noise/ambiguous, signal  

IRIS generalizes SAINT by:

- providing stabilized, biologically meaningful per‑prey posterior probabilities  
- adding an intermediate component to prevent ambiguous proteins from distorting background and signal  
- introducing hierarchical regularization via a global Gamma prior to stabilize component means  
- stabilizing mixture weights with a Dirichlet prior to prevent component collapse  
- improving robustness to sparse counts, outliers, and replicate imbalance  


*NOTE: A *bait‑unit* is the specific antibody-protein immunoprecipitation setup
whose replicates share the same underlying interaction distribution.

In other words, a bait-unit corresponds to a particular antibody pulling down a 
particular tagged or endogenous protein under a specific construct/condition, 
with all its technical replicates grouped together.

e.g. 

     TP53-myc + all replicates = a bait-unit;
     
     TP53 (endogenous) + all replicates = another bait-unit;
     
     Control-protein-myc + all replicates = another bait-unit;
     
     TP53-KO (endogenous) + all replicates = another bait-unit
     
IRIS fits one independent mixture model per bait-unit
---

## IRIS model architecture (modeler view)

For a given bait‑specific experiment, IRIS models each collapsed prey count $$\(x_i = \sum_j X_{ij}\)$$ as:


$$
p(x_i \mid \Theta) = \sum_{k=1}^{3} \pi_k \, \text{Pois}(x_i \mid \lambda_k)
$$


where:

- $$\(k = 1,2,3\)$$ correspond to background, noise/ambiguous, signal  
- $$\(\lambda_k\)$$ are component-specific Poisson means for collapsed counts
- $$\(\pi_k\)$$ are mixture weights with $$\(\sum_k \pi_k = 1\)$$

### E‑step: responsibilities


$$
\gamma_{ik} =
\frac{\pi_k \, \text{Pois}(x_i \mid \lambda_k)}
{\sum_{j=1}^{3} \pi_j \, \text{Pois}(x_i \mid \lambda_j)}
$$


### M‑step: component means with global Gamma shrinkage

IRIS applies a $$Gamma(\(\alpha, \tau\))$$ prior shared across experiments
(These priors regularize the collapsed Poisson means, not replicate-specific rates):


$$
\lambda_k^{\text{new}} =
\frac{\alpha + \sum_i \gamma_{ik} x_i}
{\tau + \sum_i \gamma_{ik}}
$$


Biologically, this reflects the assumption that background, ambiguous, and enriched signal levels should fall within broadly similar ranges across experiments-- or at least that, with more replicates, the estimates approach true mean signal levels.  
Shrinkage prevents sparse or noisy experiments from producing extreme or biologically implausible component means.

### M‑step: mixture weights with Dirichlet stabilization

Mixture weights represent the proportion of preys that are background, ambiguous, or enriched in an experiment.  
IRIS updates them using a $$Dirichlet(\(\beta_1, \beta_2, \beta_3\))$$ prior:


$$
\pi_k^{\text{new}} =
\frac{\beta_k + \sum_i \gamma_{ik}}
{\sum_j \left(\beta_j + \sum_i \gamma_{ij}\right)}
$$


Stabilization prevents any component from collapsing to zero prevalence, ensuring that IRIS maintains meaningful component structure and can express uncertainty even in low‑information settings.

### Per‑bait-unit independence

**IRIS fits this model independently for each bait‑unit**, allowing:

- bait-unit‑specific signal distributions  
- bait-unit‑specific mixture weights  
- shared global hyperparameters $\alpha, \tau, \beta$ across experiments

---

## Why IRIS improves upon SAINT

IRIS introduces several structural enhancements that directly address SAINT’s limitations:

- **Stabilized, biologically meaningful per‑prey posterior probabilities**  
  By generating per-prey estimates, IRIS prevents ambiguous or noisy proteins that, in a non- per-prey paradigm, would disproportionately influence the component parameters used to evaluate all other preys.
  A disproprtionate regression across preys also shifts component mixtures together. Thus, per-prey independence yields clearer separation between components and more reliable signal probabilities.

- **Intermediate component for ambiguous signal**  
  The intermediate component absorbs borderline counts so that they do not distort background or signal counts.

- **Global hierarchical shrinkage**  
  Regressing all prey estimates toward a global mean stabilizes extreme estimates.
  This also implies the biological assumption that background, ambiguous, and enriched signal levels should fall within broadly similar ranges across experiments (sampling variability should stabilize as replication number increases).

- **Dirichlet‑stabilized mixture weights**  
  Ensures that all components remain viable and prevents collapse in low‑information settings.

---

## Per‑prey independence and why it matters

To reiterate: *the most important improvement of IRIS over SAINT is*: **IRIS provides per‑prey, per‑bait-unit posterior probabilities** that are more stable and interpretable because they are less distorted by the behavior of other preys.

In SAINT, all preys jointly determine only two global components (background and signal).  
This means that background, noisy, and signal proteins all influence the same two component means and mixture weights.  

As a result, ambiguous or inconsistent proteins can distort the estimates used for every other prey.
In other words, in a non- per-prey paradigm, prey proteins (potentially) completely unrelated to another prey directly effect that prey protein's estimates, masking its empirical evidence and fundamentally mis-representing the underlying biology.

IRIS' independent, per-prey estimates **prevent this undue influence by noisy, ambiguous, and/or unrelated proteins on component parameters** that determine the posterior probabilities.
However, IRIS' prey proteins are intentionally linked by:

- bait-unit-wide hierarchical shrinkage-- regressing per-prey estimates toward bait-unit global mean
- mixture-weight stabilization

These enforce consistency and acknowledge the underlying biological context of a given sample/replicate.
This shrinkage intentionally has more influence on ambiguous or high-variability prey, allowing more consistent or high-confidence prey estimates to retain empirical signal.

The drives a clearer separation between components with more stable signal probabilities due to the shrinkage, and therefore, results in more biologically meaningful per‑prey estimates.

---

## IRIS interpretation (user & biologist view)

IRIS separates observed spectral counts into three probabilistic signal categories:

- background — noise, contaminants, non‑specific binders  
- noise/ambiguous — borderline, inconsistent, or low‑confidence signal  
- signal — strong, reproducible interaction signal  

The algorithm iteratively:

1. Estimates how likely each prey belongs to each signal category  
2. Updates the definitions of those categories based on the data  
3. Repeats until prey assignments to each component and with their given parameters stabilize  

### Per‑prey (per‑protein) estimates

For each prey $i$, IRIS provides:

- a posterior probability vector  

$$
\gamma_i = \left(
  \gamma_{i,\text{background}},
  \gamma_{i,\text{ambiguous}},
  \gamma_{i,\text{signal}}
\right)
$$


- an optional classification by taking the most probable component  
- a stabilized, biologically interpretable view of its interaction behavior that accounts for bait-unit-wide structure

Downstream analyses focus on the posterior probability of belonging to the **signal** component, using the background and ambiguous components as context for uncertainty and noise:

- posterior probability belonging to the signal component for all target bait-units (i.e. not negative control) are averaged to provide the first, primary metric for candidate scoring.

- in other words, **for target protein X, all $\gamma_3$ values are averaged together** 

- (e.g., for the TP53 example above, all non‑negative‑control TP53 condition $\gamma_{3,\text{signal}}$ values are averaged together).

---

## Usage examples

### 1. Fitting IRIS to a single bait-unit (minimal example)

```python
    import numpy as np
    from orion.methods.iris.em_wrapper import run_em
    
    # Spectral counts for a single bait-specific experiment
    # X should be an (n_prey × r) matrix; if your data is 1D, reshape or load accordingly
    X = np.load("counts.npy")  # shape: (n_prey, r)
    
    # Run the IRIS EM algorithm directly
    fit = run_em(
        X=X,
        max_iter=200,
        tol_loglik=1e-6,
        tol_params=1e-6,
        seed=1,
        verbose=True,
    )
    
    # Posterior probabilities for each component (background, ambiguous, signal)
    gamma = fit["responsibilities"]  # shape: (n_prey, 3)
    
    # Component means and mixture weights
    lambdas = fit["lambdas"]  # length 3
    pis = fit["pis"]          # length 3
```

### 2. Fitting IRIS to all bait-units

```python
    from orion.methods.iris.pipeline import run_iris_pipeline
    from orion.utils.io.data_input import load_bait_data
    
    # Load bait-specific spectral count data
    # Returns a dictionary: bait → 1D array of prey counts
    input_data = load_bait_data("path/to/input_directory")
    
    # Run the full IRIS pipeline
    results = run_iris_pipeline(
        input_data=input_data,
        max_iter=200,
        tol_loglik=1e-6,
        tol_params=1e-6,
        seed=1,
        verbose=True,
        make_plots=False,
        save_results=True,
        results_csv="iris_results.csv",
    )
    
    # 'results' is a dictionary keyed by bait, containing:
    # - posterior probabilities (n_prey × 3)
    # - component means
    # - mixture weights
    # - tau-grid diagnostics (if enabled)

```

For comparison, a SAINT‑style two‑component fit might look like:

```python
    from orion.methods.saint.pipeline import run_classical_pipeline
    from orion.utils.io.data_input import load_bait_data
    
    # Load bait-specific spectral count data
    # Returns a dictionary: bait → 1D array of prey counts
    input_data = load_bait_data("path/to/input_directory")
    
    # Run the classical SAINT pipeline
    results = run_classical_pipeline(
        input_data=input_data,
        hyperparams=None,      # optional initialization values
        make_plots=False,
        plot_dir=None,
        max_iter=200,
        tol_loglik=1e-6,
        tol_params=1e-6,
        seed=1,
        verbose=True,
    )
    
    # 'results' is a dictionary keyed by bait, containing:
    # - posterior probabilities (n_prey × 2)
    # - component means (background, signal)
    # - mixture weights
    # - optional diagnostic plots (if enabled)
```

---

## Audience guide

- If you care about implementation and math → focus on the model architecture and EM update equations.  
- If you care about biological interpretation → focus on the interpretation sections and the three signal categories.  
- If you are comparing SAINT vs IRIS → note that IRIS:
  - provides stabilized per‑prey posterior probabilities,  
  - adds an intermediate component,  
  - uses global shrinkage and Dirichlet stabilization,  
  - and is designed as a theoretical and practical improvement over the SAINT two‑component mixture.

---

## License

MIT License (see repository root).


