# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 12:06:37 2026

@author: Erik
"""

# %% Imports

import numpy as np
import pandas as pd
from orion.utils.io.data_input import load_bait_data


# %% Tests

def test_load_bait_data_basic_structure():
    # Minimal synthetic wide-format input
    df = pd.DataFrame({
        "Protein": ["P1", "P2"],
        "condX_BaitA_1": [5, 1],
        "condX_BaitA_2": [7, 0],
        "condX_BaitB_1": [3, 0],
        "condX_BaitB_2": [4, 1],
    })

    bait_unit_list, X_by_bait_unit, metadata = load_bait_data(df)

    # Bait list must contain both baits
    assert set(bait_unit_list) == {"BaitA", "BaitB"}

    # Replicate matrices must exist for each bait
    assert "BaitA" in X_by_bait_unit
    assert "BaitB" in X_by_bait_unit

    XA = X_by_bait_unit["BaitA"]
    XB = X_by_bait_unit["BaitB"]

    # Each matrix must be (n_proteins × n_replicates)
    assert XA.shape == (2, 2)
    assert XB.shape == (2, 2)

    # Values must match input
    assert np.allclose(XA[:, 0], [5, 1])
    assert np.allclose(XA[:, 1], [7, 0])

    assert np.allclose(XB[:, 0], [3, 0])
    assert np.allclose(XB[:, 1], [4, 1])

    # Metadata invariants
    assert "proteins_by_bait_unit" in metadata
    assert metadata["proteins_by_bait_unit"]["BaitA"] == ["P1", "P2"]
    assert metadata["proteins_by_bait_unit"]["BaitB"] == ["P1", "P2"]

    assert "bait_units" in metadata
    assert "proteins_by_bait_unit" in metadata
    assert "control_bait_units" in metadata
    assert "treatment_bait_units" in metadata

