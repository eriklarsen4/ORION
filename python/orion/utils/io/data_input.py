# -*- coding: utf-8 -*-
"""
Created on Sun Feb  8 23:36:08 2026

@author: Erik
"""

# %% Imports

import pandas as pd
from pathlib import Path

# %%
def load_bait_data(user_provided_df):
    """
    Load IPMS data and construct per-bait replicate matrices for the hierarchical
    SAINT pipeline. This function accepts either a wide-format DataFrame or a
    path to a CSV file, converts the data into the long format required by the
    model, and assembles the per-bait numeric matrices used by the hierarchical
    EM algorithm and tau grid search.

    The loader preserves all replicates, performs no collapsing, and enforces a
    deterministic protein ordering within each bait. The resulting metadata
    includes the bait list and the per-bait protein lists, which define the
    canonical row ordering for all downstream model outputs.

    Parameters
    
    user_provided_df : DataFrame or str
        Raw IPMS data in wide format (Protein name column plus replicate columns
        of the form <condition>_<bait>_<rep>), or a path to a .CSV file
        containing such a table.

    Returns
    
    bait_unit_list : list of str
        Sorted list of bait names present in the dataset.

    X_by_bait_unit : dict
        Dictionary mapping each bait_unit to its numeric replicate matrix
        (rows = proteins in deterministic order, columns = replicates).

    metadata : dict
        Dictionary containing:
            - "bait_units": the sorted list of bait_units
            - "proteins_by_bait_unit": dict mapping each bait_unit to its ordered protein list
            - "control_bait_units": bait_units with Condition == 'Control'
            - "treatment_bait_units": bait_units with Condition != 'Control'

    Notes
    
    This function performs only data loading and restructuring. It does not
    apply the EM model, tau grid search, or any statistical transformations.
    All model parameters (lambda1, lambda2, lambda3, pi, gamma) are estimated
    downstream by the hierarchical EM algorithm and are not supplied here.

    The replicate structure is preserved exactly as provided in the input.
    No averaging, collapsing, or normalization is performed at this stage.
    """

    # Case 1: user passed a DataFrame
    if isinstance(user_provided_df, pd.DataFrame):
        user_provided_df = user_provided_df.copy()
    else:
        path = Path(user_provided_df)
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {user_provided_df}")
        user_provided_df = pd.read_csv(path)

    # Convert wide → long
    long_df = extract_bait_matrix(user_provided_df)
    
    # Identify control vs treatment baits using exact match on Condition
    control_mask = long_df["Condition"] == "Control"
    control_baits = sorted(long_df.loc[control_mask, "BaitUnit"].unique())
    treatment_baits = sorted(long_df.loc[~control_mask, "BaitUnit"].unique())
    
    # Identify all baits
    bait_list = sorted(long_df["BaitUnit"].unique())

    X_by_bait = {}
    proteins_by_bait = {}

    for bait in bait_list:
        df_b = long_df[long_df["BaitUnit"] == bait].copy()

        # Deterministic row order
        df_b = df_b.sort_values("Protein").reset_index(drop=True)

        # Store per-bait protein list
        proteins_by_bait[bait] = df_b["Protein"].tolist()

        # Extract numeric replicate columns
        numeric_cols = [c for c in df_b.columns if c not in ("Protein", "Bait", "Condition")]
        X_by_bait[bait] = df_b[numeric_cols].astype(float).to_numpy()

    metadata = {
        "bait_units": bait_unit_list,
        "proteins_by_bait_unit": proteins_by_bait,
        "control_bait_units": control_bait_units,
        "treatment_bait_units": treatment_bait_units,
    }

    return bait_list, X_by_bait, metadata


# %% Extract bait matrix (wide → long with replicate preservation)

def extract_bait_matrix(user_provided_df):
    """
    Convert a wide-format IPMS table into the long-format structure required by
    the hierarchical SAINT pipeline. The input table must contain a Protein
    column and replicate measurement columns of the form:

        <condition>_<bait>_<rep>

    where <bait> identifies the bait name and <rep> is an integer replicate
    index. Replicates are preserved exactly as provided; no collapsing or
    averaging is performed.

    This function parses the wide-format replicate columns, groups them by bait,
    and reconstructs a long-format table in which each row corresponds to a
    (Protein, Bait) pair and each replicate is placed in its own column
    (rep1, rep2, rep3, ...). Missing replicates for a bait are filled with zeros
    to maintain a consistent replicate structure.

    Parameters
    
    user_provided_df : DataFrame
        Wide-format IPMS data containing a Protein column and replicate columns
        named according to the <condition>_<bait>_<rep> convention.

    Returns
    
    long_df : DataFrame
        Long-format table with columns:
        Protein   : str
        BaitUnit  : str   (unique modeling unit = condition × bait_name)
        BaitName  : str   (biological bait name)
        Condition : str
        rep1, rep2, rep3, ... : float

    Each BaitUnit appears once per protein, and replicate columns are ordered
    by increasing replicate index.

    Notes
    
    This function performs only structural reshaping. It does not normalize,
    filter, or transform the data. The resulting long-format table is consumed
    by `load_bait_data`, which enforces deterministic protein ordering and
    constructs the per-bait numeric matrices used by the hierarchical EM model
    and tau grid search.

    """

    protein = user_provided_df["Protein"].tolist()
    value_cols = [c for c in user_provided_df.columns if c != "Protein"]
    
    parsed = {}
    
    for col in value_cols:
        parts = col.split("_")
        if len(parts) < 3:
            continue
    
        # NAMING SCHEME:
        # <condition>_<bait>_<rep>
        condition = parts[0]
        bait_name = parts[1]
        rep_label = parts[2]
    
        try:
            rep_index = int(rep_label)
        except ValueError:
            continue
    
        counts = pd.to_numeric(user_provided_df[col], errors="coerce").fillna(0).to_numpy()
    
        if bait_name not in parsed:
            parsed[bait_name] = {
                "condition": condition,
                "reps": {}
            }
    
        parsed[bait_name]["reps"][rep_index] = counts
    
    # Build long-format dataframe
    records = []
    
    for bait_name, info in parsed.items():
        rep_dict = info["reps"]
        condition = info["condition"]
        max_rep = max(rep_dict.keys())
    
        for i, p in enumerate(protein):
            #row = {
            #    "Protein": p,
            #    "Bait": bait_name,
            #    "Condition": condition
            #}
            
            bait_unit = f"{condition}_{bait_name}"
            
            row = {"Protein": p,
                   "BaitUnit": bait_unit,
                   "BaitName": bait_name,
                   "Condition": condition}
    
            for r in range(1, max_rep + 1):
                row[f"rep{r}"] = float(rep_dict.get(r, pd.Series([0]*len(protein)))[i])
    
            records.append(row)
    
    long_df = pd.DataFrame(records)
    return long_df


