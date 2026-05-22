# -*- coding: utf-8 -*-
"""
Created on Thu May  7 15:36:39 2026

@author: Erik
"""
# %% Import
from pathlib import Path
from typing import Dict, List, Tuple, Union

import pandas as pd
import numpy as np

# %% IRIS function: load the data
def load_bait_data_iris(
    user_provided_df: Union[pd.DataFrame, str],
    metadata: Dict
) -> Tuple[List[str], Dict[str, np.ndarray], Dict]:
    """
    Load IPMS data and construct per-bait_unit replicate matrices for the IRIS
    pipeline.

    This loader accepts either a wide-format DataFrame or a path to a CSV file,
    converts the data into the long format required by IRIS, and assembles
    per-bait_unit numeric matrices. For IRIS, the modeling unit is the
    condition × bait combination, so bait_unit is defined as:

        bait_unit = f"{condition}_{bait_name}"

    Replicates are preserved (no collapsing), and a deterministic protein
    ordering is enforced within each bait_unit. The resulting metadata includes
    the bait_unit list, per-bait_unit protein lists, and now the expression
    system associated with each bait_unit.

    PARAMETERS

    user_provided_df : DataFrame or str
        Raw IPMS data in wide format (Protein column plus replicate columns
        of the form <condition>_<bait>_<rep>), or a path to a .CSV file
        containing such a table.

    metadata : dict
        User-supplied metadata including:
            - "expression_system": dict mapping bait_name → endogenous/non_endogenous

    RETURNS

    bait_unit_list : list of str
        Sorted list of bait_units (condition × bait) present in the dataset.

    X_by_bait_unit : dict
        Dictionary mapping each bait_unit to its numeric replicate matrix
        (rows = proteins in deterministic order, columns = replicates).

    metadata_out : dict
        Dictionary containing:
            - "bait_units": sorted list of bait_units
            - "proteins_by_bait_unit": dict mapping bait_unit → ordered protein list
            - "control_bait_units": bait_units with Condition == 'Control'
            - "treatment_bait_units": bait_units with Condition != 'Control'
            - "expression_system_by_bait_unit": dict mapping bait_unit → expression system
    """

    # Case 1: user passed a DataFrame
    if isinstance(user_provided_df, pd.DataFrame):
        user_provided_df = user_provided_df.copy()
    else:
        path = Path(user_provided_df)
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {user_provided_df}")
        user_provided_df = pd.read_csv(path)

    # Convert wide → long with bait_unit = condition_bait
    long_df = _extract_bait_matrix_iris(
        user_provided_df,
        metadata["expression_system"]
    )

    # Identify control vs treatment bait_units using Condition
    control_mask = long_df["Condition"] == "Control"
    control_bait_units = sorted(long_df.loc[control_mask, "BaitUnit"].unique())
    treatment_bait_units = sorted(long_df.loc[~control_mask, "BaitUnit"].unique())

    # Identify all bait_units
    bait_unit_list = sorted(long_df["BaitUnit"].unique())

    X_by_bait_unit: Dict[str, np.ndarray] = {}
    proteins_by_bait_unit: Dict[str, List[str]] = {}

    for bait_unit in bait_unit_list:
        df_b = long_df[long_df["BaitUnit"] == bait_unit].copy()

        # Deterministic row order
        df_b = df_b.sort_values("Protein").reset_index(drop=True)

        # Store per-bait_unit protein list
        proteins_by_bait_unit[bait_unit] = df_b["Protein"].tolist()

        # Extract numeric replicate columns (drop all-NaN reps for this bait_unit)
        numeric_cols = sorted(
            [
                c
                for c in df_b.columns
                if c.startswith("rep") and df_b[c].notna().any()
            ]
        )
        X_by_bait_unit[bait_unit] = df_b[numeric_cols].astype(float).to_numpy()

    # expression system per bait_unit
    expression_system_by_bait_unit = {
        bu: long_df.loc[long_df["BaitUnit"] == bu, "ExpressionSystem"].iloc[0]
        for bu in bait_unit_list
    }

    metadata_out = {
            "bait_units": bait_unit_list,
            "proteins_by_bait_unit": proteins_by_bait_unit,
            "control_bait_units": control_bait_units,
            "treatment_bait_units": treatment_bait_units,
            "expression_system_by_bait_unit": expression_system_by_bait_unit,
        }
    
    # IRIS-specific metadata for PipelineMetadata.extra_fields
    bait_name_by_bait_unit = {
        bu: long_df.loc[long_df["BaitUnit"] == bu, "BaitName"].iloc[0]
        for bu in bait_unit_list
    }

    condition_by_bait_unit = {
        bu: long_df.loc[long_df["BaitUnit"] == bu, "Condition"].iloc[0]
        for bu in bait_unit_list
    }

    metadata_out["extra_fields"] = {
        "bait_name_by_bait_unit": bait_name_by_bait_unit,
        "condition_by_bait_unit": condition_by_bait_unit,
        "expression_system_by_bait_unit": expression_system_by_bait_unit,
    }

    return bait_unit_list, X_by_bait_unit, metadata_out


# %% Extract the bait matrix, having parsed the conditions, samples, and replicates
def _extract_bait_matrix_iris(
    user_provided_df: pd.DataFrame,
    expression_system: Dict[str, str]
) -> pd.DataFrame:

    protein = user_provided_df["Protein"].tolist()
    value_cols = [c for c in user_provided_df.columns if c != "Protein"]

    parsed = {}

    # Parse columns into bait_unit → rep_index → vector
    for col in value_cols:
        parts = col.split("_")
        if len(parts) < 3:
            continue

        condition = parts[0]
        bait_name = "_".join(parts[1:-1])
        rep_label = parts[-1]

        try:
            rep_index = int(rep_label)
        except ValueError:
            continue

        bait_unit = f"{condition}_{bait_name}"

        if bait_unit not in parsed:
            parsed[bait_unit] = {
                "condition": condition,
                "bait_name": bait_name,
                "reps": {}
            }

        parsed[bait_unit]["reps"][rep_index] = (
            pd.to_numeric(user_provided_df[col], errors="coerce")
            .fillna(0)
            .to_numpy()
        )

    # Build a separate DataFrame per bait_unit
    dfs = []

    for bait_unit, info in parsed.items():
        condition = info["condition"]
        bait_name = info["bait_name"]
        rep_dict = info["reps"]

        # lookup expression system
        expr = expression_system.get(bait_name, None)

        rows = []
        for i, p in enumerate(protein):
            row = {
                "Protein": p,
                "BaitUnit": bait_unit,
                "BaitName": bait_name,
                "Condition": condition,
                "ExpressionSystem": expr,
            }
            for r in sorted(rep_dict.keys()):
                row[f"rep{r}"] = float(rep_dict[r][i])
            rows.append(row)

        dfs.append(pd.DataFrame(rows))

    return pd.concat(dfs, ignore_index=True, sort=False)
