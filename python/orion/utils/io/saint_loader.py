# -*- coding: utf-8 -*-
"""
Created on Thu May  7 15:32:00 2026

@author: Erik
"""

# %% Import
from pathlib import Path
from typing import Dict, List, Tuple, Union

import pandas as pd
import numpy as np

# %% SAINT data import function

def load_bait_data_saint(
    user_provided_df: Union[pd.DataFrame, str]
) -> Tuple[List[str], Dict[str, np.ndarray], Dict]:
    """
    Load IPMS data and construct per-bait replicate matrices for the
    classical SAINT pipeline.

    This loader accepts either a wide-format DataFrame or a path to a CSV file,
    converts the data into a long format, and assembles per-bait numeric
    matrices. Unlike the IRIS pipeline, the classical SAINT model treats the
    bait name itself as the modeling unit; condition is not part of the
    bait_unit identity.

    Replicates are preserved (no collapsing), and a deterministic protein
    ordering is enforced within each bait. The resulting metadata includes the
    bait list and per-bait protein lists, which define the canonical row
    ordering for downstream classical SAINT outputs.


    PARAMETERS
    
    user_provided_df : DataFrame or str
        Raw IPMS data in wide format (Protein column plus replicate columns
        of the form <condition>_<bait>_<rep>), or a path to a .CSV file
        containing such a table.


    RETURNS
    
    bait_unit_list : list of str
        Sorted list of bait names present in the dataset (bait-only units).

    X_by_bait_unit : dict
        Dictionary mapping each bait (bait name) to its numeric replicate matrix
        (rows = proteins in deterministic order, columns = replicates).

    metadata : dict
        Dictionary containing:
            - "bait_units": sorted list of bait names
            - "proteins_by_bait_unit": dict mapping bait → ordered protein list
            - "control_bait_units": bait names with Condition == 'Control'
            - "treatment_bait_units": bait names with Condition != 'Control'
    """

    # Case 1: user passed a DataFrame
    if isinstance(user_provided_df, pd.DataFrame):
        user_provided_df = user_provided_df.copy()
    else:
        path = Path(user_provided_df)
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {user_provided_df}")
        user_provided_df = pd.read_csv(path)

    # Convert wide → long with bait = bait_name (no condition in bait_unit)
    long_df = _extract_bait_matrix_classical(user_provided_df)

    # Identify control vs treatment baits using Condition, but group by bait name
    control_mask = long_df["Condition"] == "Control"
    control_bait_units = sorted(long_df.loc[control_mask, "BaitUnit"].unique())
    treatment_bait_units = sorted(long_df.loc[~control_mask, "BaitUnit"].unique())

    # Identify all baits (bait names only)
    bait_unit_list = sorted(long_df["BaitUnit"].unique())
    
    # Classical SAINT: control bait is named "CTRL"
    control_bait_units = sorted([b for b in bait_unit_list if b.upper() == "CTRL"])
    treatment_bait_units = sorted([b for b in bait_unit_list if b.upper() != "CTRL"])

    X_by_bait_unit: Dict[str, np.ndarray] = {}
    proteins_by_bait_unit: Dict[str, List[str]] = {}

    for bait in bait_unit_list:
        df_b = long_df[long_df["BaitUnit"] == bait].copy()

        # Deterministic row order
        df_b = df_b.sort_values("Protein").reset_index(drop=True)

        # Store per-bait protein list
        proteins_by_bait_unit[bait] = df_b["Protein"].tolist()

        # Extract numeric replicate columns
        numeric_cols = [c for c in df_b.columns if c.startswith("rep")]
        X_by_bait_unit[bait] = df_b[numeric_cols].astype(float).to_numpy()

    metadata = {
        "bait_units": bait_unit_list,
        "proteins_by_bait_unit": proteins_by_bait_unit,
        "control_bait_units": control_bait_units,
        "treatment_bait_units": treatment_bait_units,
    }

    return bait_unit_list, X_by_bait_unit, metadata


def _extract_bait_matrix_classical(user_provided_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert a wide-format IPMS table into the long-format structure required by
    the classical SAINT pipeline.

    The input table must contain a Protein column and replicate measurement
    columns of the form:

        <condition>_<bait>_<rep>

    where <bait> identifies the bait name and <rep> is an integer replicate
    index. Replicates are preserved exactly as provided; no collapsing or
    averaging is performed.

    For the classical pipeline, the modeling unit is the bait name alone.
    Condition is tracked but is not part of the bait_unit identity.


    RETURNS
    
    long_df : DataFrame
        Long-format table with columns:
            Protein   : str
            BaitUnit  : str   (bait name only)
            BaitName  : str   (same as BaitUnit)
            Condition : str
            rep1, rep2, rep3, ... : float
    """

    protein = user_provided_df["Protein"].tolist()
    value_cols = [c for c in user_provided_df.columns if c != "Protein"]

    parsed: Dict[str, Dict] = {}

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
                "reps": {},
            }

        parsed[bait_name]["reps"][rep_index] = counts

    # Build long-format dataframe
    records = []

    for bait_name, info in parsed.items():
        rep_dict = info["reps"]
        condition = info["condition"]
        max_rep = max(rep_dict.keys())

        for i, p in enumerate(protein):
            bait_unit = bait_name  # classical: bait_unit = bait only

            row = {
                "Protein": p,
                "BaitUnit": bait_unit,
                "BaitName": bait_name,
                "Condition": condition,
            }

            for r in range(1, max_rep + 1):
                row[f"rep{r}"] = float(
                    rep_dict.get(r, pd.Series([0] * len(protein)))[i]
                )

            records.append(row)

    long_df = pd.DataFrame(records)
    return long_df
