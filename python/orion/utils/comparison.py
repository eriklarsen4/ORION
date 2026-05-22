# -*- coding: utf-8 -*-
"""
Created on Tue May 19 09:54:33 2026

@author: Erik
"""
# %% Import
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %% Benchmark function
def compare_posteriors(saint_df, iris_df):
    """
    Compare SAINT and IRIS posteriors for treatment bait units only.
    
    SAINT:
        - One row per Protein × BaitUnit
        - gamma1 = background posterior
        - signal posterior = 1 - gamma1

    IRIS:
        - One row per Protein × BaitUnit
        - BaitUnit includes prefixes: Treat_*
        - gamma3 = signal posterior

    Returns:
        DataFrame with per-protein averaged posteriors:
            Protein
            SAINT_gamma
            IRIS_gamma
    """

    # --- IRIS: keep only treatment bait units ---
    iris_treat = iris_df[iris_df["BaitUnit"].str.startswith("Treat_")].copy()

    # strip prefix so it matches SAINT's bait names
    iris_treat["BaseBait"] = iris_treat["BaitUnit"].str.replace("^Treat_", "", regex=True)

    # keep only what we need
    iris_treat = iris_treat[["Protein", "BaseBait", "gamma3"]]


    # --- SAINT: compute signal posterior ---
    saint_df = saint_df.copy()
    saint_df["BaseBait"] = saint_df["BaitUnit"]          # already no prefix
    saint_df["SAINT_gamma"] = 1 - saint_df["gamma1"]      # gamma1 = background posterior

    saint_df = saint_df[["Protein", "BaseBait", "SAINT_gamma"]]


    # --- MERGE ON Protein + BaseBait ---
    merged = saint_df.merge(
        iris_treat,
        on=["Protein", "BaseBait"],
        how="inner"
    )


    # --- AVERAGE PER PROTEIN ---
    posterior_avg = (
        merged.groupby("Protein")[["SAINT_gamma", "gamma3"]]
              .mean()
              .reset_index()
              .rename(columns={"gamma3": "IRIS_gamma"})
    )

    return posterior_avg

# %% Visualizing Model Outputs
def diagnostic_panel(df, saint_col="SAINT_gamma", iris_col="IRIS_gamma",
                                protein_col="Protein", controls=None, figsize=(12, 10)):
    """
    Diagnostic panel for comparing SAINT and IRIS posteriors.

    df          : DataFrame with SAINT and IRIS columns.
    saint_col   : column name for SAINT posterior.
    iris_col    : column name for IRIS posterior.
    protein_col : column name for protein IDs (only used if controls is not None).
    controls    : optional iterable of protein IDs to highlight (no hard-coding).
    """

    # compute delta = IRIS - SAINT
    delta = df[iris_col] - df[saint_col]

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    ax_scatter, ax_delta, ax_saint, ax_iris = axes.ravel()

    # ---------- Panel A: SAINT vs IRIS scatter ----------
    sns.scatterplot(
        data=df,
        x=saint_col,
        y=iris_col,
        s=12,
        alpha=0.4,
        color="gray",
        edgecolor=None,
        ax=ax_scatter
    )

    # optional control highlighting (user-specified, not hard-coded)
    if controls is not None and protein_col in df.columns:
        ctrl_df = df[df[protein_col].isin(controls)]
        for _, row in ctrl_df.iterrows():
            ax_scatter.scatter(
                row[saint_col],
                row[iris_col],
                s=120,
                edgecolor="black",
                facecolor="none",
                linewidth=1.2,
                zorder=10
            )
            ax_scatter.text(
                row[saint_col] + 0.01,
                row[iris_col] + 0.01,
                str(row[protein_col]),
                fontsize=9,
                weight="bold",
                zorder=11
            )

    ax_scatter.plot([0, 1], [0, 1], "k--", linewidth=1)
    ax_scatter.set_xlabel("SAINT signal posterior")
    ax_scatter.set_ylabel("IRIS gamma3 posterior")
    ax_scatter.set_title("A. SAINT vs IRIS posteriors")

    # ---------- Panel B: Δ = IRIS − SAINT ----------
    sns.histplot(delta, bins=80, color="steelblue", alpha=0.7, ax=ax_delta)
    ax_delta.axvline(0, color="k", linestyle="--", linewidth=1)
    ax_delta.set_xlabel("delta = IRIS − SAINT")
    ax_delta.set_ylabel("Count")
    ax_delta.set_title("B. delta distribution (model disagreement)")

    # ---------- Panel C: SAINT marginal ----------
    sns.histplot(df[saint_col], bins=80, color="darkgray", alpha=0.8, ax=ax_saint)
    ax_saint.set_xlabel("SAINT signal posterior")
    ax_saint.set_ylabel("Count")
    ax_saint.set_title("C. SAINT posterior marginal")

    # ---------- Panel D: IRIS marginal ----------
    sns.histplot(df[iris_col], bins=80, color="darkslateblue", alpha=0.8, ax=ax_iris)
    ax_iris.set_xlabel("IRIS gamma3 posterior")
    ax_iris.set_ylabel("Count")
    ax_iris.set_title("D. IRIS posterior marginal")

    plt.tight_layout()
    plt.show()

    # numeric summary
    corr = df[saint_col].corr(df[iris_col])
    print(f"Correlation(SAINT, IRIS) = {corr:.3f}")
    print(f"delta > 0 (IRIS > SAINT): {(delta > 0).mean():.3f}")
    print(f"delta < 0 (IRIS < SAINT): {(delta < 0).mean():.3f}")
