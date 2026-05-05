# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 13:08:34 2026

@author: Erik
"""

# %% Imports
import numpy as np

# %% IRIS EM wrapper
def run_em(
    X,
    hyperparams,
    bait_unit,
    max_iter=200,
    tol_loglik=1e-6,
    tol_params=1e-6,
    seed=1,
    verbose=False,
):
    """
    Fit a three-component Poisson mixture model for a single bait_unit (across 
    baits) using an EM-style alternating procedure with fixed Gamma shrinkage
    and Dirichlet-based stabilization.
    
    
    BAIT UNIT
    
    A bait_unit is defined by a specific immunoprecipitated protein (the bait),
    either endogenous or tagged, under a specific experimental condition. All
    replicate IP samples for that bait under that condition belong to the same
    bait_unit. For example, with target protein, TP53, tagged with myc (TP53-myc),
    the bait_unit is TP53-myc, which includes all replicates; however, an 
    endogenous target, TP53, in the same study, would be considered a different 
    bait_unit.
    
    Each bait_unit is modeled independently, with its own λ, π, γ, α,
    and (a_k, b_k). Global hyperparameters such as τ and the model structure are
    shared across bait_units. Parameters are defined below.
    
    
    DATA STRUCTURE
    
    For a given bait_unit, X is an n × r matrix of non-negative spectral counts:
        n : number of prey proteins
        r : number of technical or biological replicates (samples)
    Each row X_i contains the replicate counts for prey i. All replicates in X
    belong to the same bait_unit. Thus the entire data frame for a study contains
    all X's for all bait_units, joined together.
    
    
    MIXTURE MODEL
    
    As defined in the DATA STRUCTURE section, r is the number of replicate 
    samples for a given bait_unit. Each replicate sample contains one abundance 
    measurement for each prey protein. Therefore, for prey i, X_i is a vector 
    of r measurements, one measurement per replicate sample of a given 
    bait_unit.

    A mixture implies different amounts of multiple entities. Here, the mixture
    is the different proportions of 3 components. Therefore, mixture 
    proportions are:
    
        π = (π_1, π_2, π_3)
    
    with π_k ≥ 0 and π_1 + π_2 + π_3 = 1. These proportions specify how 
    frequently each component occurs across prey proteins.
    
    
    LATENT COMPONENTS
    
    Each component k has its own prey-specific Poisson rate vector:

        λ_k = (λ_{k1}, …, λ_{kn})

    where λ_{ki} is the Poisson rate associated with prey i under component k.
    The three components differ only in their rate vectors and represent:

        • Component 1 (index 0): background or non-interactor
        • Component 2 (index 1): contaminant or high background
        • Component 3 (index 2): true signal or interactor
    
    
    GENERATIVE MODEL

    This section describes the probabilistic model IRIS fits. It defines the
    likelihood used in the EM-like algorithm. The process does not draw samples
    from this model; it evaluates the likelihood implied by it.

    IRIS uses a collapsed Poisson formulation. For each prey i:

        1. compute the collapsed count X_i• = ∑_j X_{ij} across replicates
        2. draw a component label z_i ∈ {1,2,3} with probability π_k
        3. generate X_i• from Poisson(r * λ_{k,i})

    This collapsed formulation is equivalent to assuming replicate-specific
    Poisson rates that are equal within a component and integrating out the
    replicate dimension. The mixture proportions enter the E-step through
    log(π_k).
    
    
    RESPONSIBILITIES
    
    For each prey i and component k, the algorithm maintains γ_{ik}, the 
    conditional probability that prey i belongs to component k given the 
    current λ_k and π_k.
    
    These γ_{ik} are posterior probabilities over latent labels only.
    
    SHRINKAGE TERMS
    
    Each λ_k has associated constants a_k and b_k defining a Gamma(a_k, b_k)
    density. These act as fixed shrinkage strengths in the "M-step" updates, 
    serving as regularization of iterative estimates. The hyperparameters 
    (a_k, b_k) are fixed, estimated initially from the initial λ_k values using
    empirical mean/variance, and remain fixed throughout the EM iterations.
    
    A scalar τ is a global shrinkage strength chosen outside the EM loop (for
    example, by grid search) and applied uniformly across bait_units.
    
    Mixture proportions π are stabilized using parameters α derived from the
    empirical mean and variance of γ. These α values are updated by a variance-based
    heuristic for numerical stability.
    
    EM-like ALGORITHM
    
    E-step:
        For each prey i and component k, compute:
            log(π_k) + log Poisson(X_i | λ_k)
        Apply row-wise log-sum-exp stabilization, exponentiate, and normalize to
        obtain γ_{ik}.
    
    M-step:
        Update λ_k using γ-weighted Poisson means with fixed Gamma(a_k, b_k)
        shrinkage controlled by τ. Update π using normalized component
        responsibilities and the stabilized α parameters. These updates maximize the
        standard EM surrogate objective, not the true observed-data likelihood.
    
    Convergence:
        The algorithm stops when changes in the log-likelihood and changes in λ_k
        and π fall below specified tolerances.
    
    PARAMETERS
    
    X : array
        An (n × r) matrix of non-negative counts for a single bait_unit.
    
    hyperparams : dict
        Contains:
            - alpha1, alpha2, alpha3 : initial Dirichlet-like stabilization parameters
            - pi_init                : initial mixture proportions
            - tau                    : global shrinkage strength
        Missing keys are assigned defaults.
    
    bait_unit : str
        Identifier for the current modeling unit (condition + bait_name).
    
    max_iter : int
        Maximum number of EM iterations.
    
    tol_loglik : float
        Convergence tolerance for changes in log-likelihood.
    
    tol_params : float
        Convergence tolerance for changes in λ and π.
    
    seed : int
        Random seed for initialization.
    
    verbose : bool
        If True, prints iteration-level diagnostics.
    
    RETURNS
    
    A dictionary containing:
    
        loglik_history :
            List of log-likelihood values, one per EM iteration.
    
        lambda1_history, lambda2_history, lambda3_history :
            Lists of component-specific rate vectors at each iteration.
    
        tau_history :
            List containing the τ value used at each iteration.
    
        pi_history :
            List of mixture proportion vectors at each iteration.
    
        gamma_history :
            List of responsibility matrices γ_{ik} at each iteration.
    
        alpha_history :
            List of α vectors used at each iteration.
    
        a_history, b_history :
            Lists containing the fixed Gamma(a_k, b_k) hyperparameters.
    
        lambda1, lambda2, lambda3 :
            Final component-specific rate vectors.
    
        tau :
            The τ value used for shrinkage.
    
        pi :
            Final mixture proportion vector.
    
        gamma :
            Final responsibility matrix of shape (n, 3).
    """

    # %% Initialization and set-up
    if hyperparams is None:
        hyperparams = {}

    # Dirichlet prior defaults
    hyperparams.setdefault("alpha1", 2.0)
    hyperparams.setdefault("alpha2", 2.0)
    hyperparams.setdefault("alpha3", 2.0)

    # Mixing proportions initialization (uniform)
    hyperparams.setdefault("pi_init", np.ones(3, dtype=float) / 3.0)

    # Global shrinkage parameter
    hyperparams.setdefault("tau", 1.0)

    # Tau grid (not used here but kept for API stability)
    hyperparams.setdefault("tau_grid", [0.1, 0.5, 1.0, 2.0])

    np.random.seed(seed)

    # Dimensions
    n, r = X.shape

    # Prey-dimensional lambda initializations
    global_mean = X.mean() / r

    lambda1 = np.full(n, 0.8 * global_mean)
    lambda2 = np.full(n, 1.0 * global_mean)
    lambda3 = np.full(n, 1.2 * global_mean)


    # Global shrinkage strength
    tau = float(hyperparams.get("tau", 1.0))

    # Mixing proportions
    pi = hyperparams["pi_init"].astype(float).copy()

    # Histories
    loglik_history = []
    lambda1_history = []
    lambda2_history = []
    lambda3_history = []
    tau_history = []
    pi_history = []
    gamma_history = []

    alpha_history = []
    a_history = []
    b_history = []

    # Precompute row sums
    X_sum = X.sum(axis=1)

    eps = 1e-12

    # Dirichlet prior for pi
    alpha = np.array(
        [hyperparams["alpha1"], hyperparams["alpha2"], hyperparams["alpha3"]],
        dtype=float,
    )

    # Initialize Gamma(a_k, b_k) shrinkage hyperparams from empirical lambda_k statistics
    def init_gamma_hyperparams(lam):
        m = float(np.mean(lam))
        v = float(np.var(lam))
        v = max(v, eps)
        a = max(m * m / v, 1.0)
        b = max(m / v, eps)
        return a, b

    # Initial Gamma hyperparameters
    a1, b1 = init_gamma_hyperparams(lambda1)
    a2, b2 = init_gamma_hyperparams(lambda2)
    a3, b3 = init_gamma_hyperparams(lambda3)
    
    # %% EM iterations
    for it in range(max_iter):
    
        # %% E-step
    
        # Component-wise log lambda
        log_lambda1 = np.log(lambda1 + eps)
        log_lambda2 = np.log(lambda2 + eps)
        log_lambda3 = np.log(lambda3 + eps)
    
        # Component-wise log-likelihood contributions: log p(X_i | lambda_k)
        loglik1 = X_sum * log_lambda1 - r * lambda1
        loglik2 = X_sum * log_lambda2 - r * lambda2
        loglik3 = X_sum * log_lambda3 - r * lambda3
    
        # Log mixture weights
        log_pi = np.log(pi + eps)
    
        # Posterior log-likelihood for each component
        log_post1 = log_pi[0] + loglik1
        log_post2 = log_pi[1] + loglik2
        log_post3 = log_pi[2] + loglik3
    
        # Shape (n, 3)
        posterior_loglik = np.vstack([log_post1, log_post2, log_post3]).T
    
        # Stabilize before exponentiating (prevent Inf's upstream of exponentiating)
        max_log = np.max(posterior_loglik, axis=1, keepdims=True)
        stabilized_log_posterior = posterior_loglik - max_log
    
        # Convert posterior log-likelihood to posterior weights
        posterior_weights = np.exp(stabilized_log_posterior)
    
        # Normalizing constant
        normalizing_constant = posterior_weights.sum(axis=1, keepdims=True)
    
        # Posterior responsibilities
        gamma = posterior_weights / (normalizing_constant + eps)
    
        gamma1 = gamma[:, 0]
        gamma2 = gamma[:, 1]
        gamma3 = gamma[:, 2]
    
        # Log-likelihood for this iteration
        loglik = float(
            np.sum(max_log[:, 0] + np.log(normalizing_constant[:, 0] + eps))
        )
    
        # %% Hyperparameter stabilization (Dirichlet alpha update; Gamma shrinkage fixed)
    
        # Stabilize Dirichlet alpha using empirical mean/var of gamma
        m = gamma.mean(axis=0)
        v = gamma.var(axis=0) + eps
        alpha0 = (m * (1.0 - m) / v - 1.0).clip(min=1.0)
        alpha = np.maximum(m * alpha0, 1.0)
    
        # Gamma(a_k, b_k) shrinkage hyperparameters
        # (initialized once from empirical lambda_k stats; not updated in EM)
        # NOTE: These are intentionally *not* updated each iteration
        # unless explicitly desired. The original code commented them out.
        # We preserve that behavior for stability.
        # a1, b1 = init_gamma_hyperparams(lambda1)
        # a2, b2 = init_gamma_hyperparams(lambda2)
        # a3, b3 = init_gamma_hyperparams(lambda3)
    
        # %% M-step (regularized update of lambda_k)
    
        # Component 1
        # gamma-weighted counts
        data_contrib_1 = gamma1 * X_sum
        # pull lambda1 toward (a1 - 1) / b1
        prior_contrib_1 = tau * (a1 - 1.0)
        # gamma-weighted replicate count
        data_weight_1 = gamma1 * r
        # shrinkage weight: tau * b1
        prior_weight_1 = tau * b1
    
        lambda1_new = np.clip(
            (data_contrib_1 + prior_contrib_1) /
            (data_weight_1 + prior_weight_1 + eps),
            1e-6, 1e6
        )
    
        # Component 2
        data_contrib_2 = gamma2 * X_sum
        prior_contrib_2 = tau * (a2 - 1.0)
        data_weight_2 = gamma2 * r
        prior_weight_2 = tau * b2
    
        lambda2_new = np.clip(
            (data_contrib_2 + prior_contrib_2) /
            (data_weight_2 + prior_weight_2 + eps),
            1e-6, 1e6
        )
    
        # Component 3
        data_contrib_3 = gamma3 * X_sum
        prior_contrib_3 = tau * (a3 - 1.0)
        data_weight_3 = gamma3 * r
        prior_weight_3 = tau * b3
    
        lambda3_new = np.clip(
            (data_contrib_3 + prior_contrib_3) /
            (data_weight_3 + prior_weight_3 + eps),
            1e-6, 1e6
        )
    
        # Updated mixing proportions
        pi_new = gamma.sum(axis=0) + alpha - 1.0
        pi_new = pi_new / np.sum(pi_new)
    
        # Store histories
        loglik_history.append(loglik)
        lambda1_history.append(lambda1_new.copy())
        lambda2_history.append(lambda2_new.copy())
        lambda3_history.append(lambda3_new.copy())
        tau_history.append(tau)
        pi_history.append(pi_new.copy())
        gamma_history.append(gamma.copy())
    
        alpha_history.append(alpha.copy())
        a_history.append(np.array([a1, a2, a3], dtype=float))
        b_history.append(np.array([b1, b2, b3], dtype=float))
    
        # Convergence diagnostics
        delta_lambda = (
            np.max(np.abs(lambda1_new - lambda1)) +
            np.max(np.abs(lambda2_new - lambda2)) +
            np.max(np.abs(lambda3_new - lambda3))
        )
        delta_pi = np.max(np.abs(pi_new - pi))
    
        if it > 0:
            if (
                abs(loglik_history[-1] - loglik_history[-2]) < tol_loglik
                and delta_lambda < tol_params
                and delta_pi < tol_params
            ):
                break
    
        # Commit updates
        lambda1 = lambda1_new
        lambda2 = lambda2_new
        lambda3 = lambda3_new
        pi = pi_new
        
        # %% Final output
    return {
        "loglik_history": loglik_history,
        "lambda1_history": lambda1_history,
        "lambda2_history": lambda2_history,
        "lambda3_history": lambda3_history,
        "tau_history": tau_history,
        "pi_history": pi_history,
        "gamma_history": gamma_history,
        "alpha_history": alpha_history,
        "a_history": a_history,
        "b_history": b_history,
        "lambda1": lambda1.copy(),
        "lambda2": lambda2.copy(),
        "lambda3": lambda3.copy(),
        "tau": float(tau),
        "pi": pi.copy(),
        "gamma": gamma.copy(),
    }