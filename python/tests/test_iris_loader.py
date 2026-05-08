# -*- coding: utf-8 -*-
"""
Created on Thu May  7 16:26:36 2026

@author: Erik
"""
# %% Import
import numpy as np
import pandas as pd
from orion.utils.io.iris_loader import load_bait_data_iris

# %% test function for IRIS loader
def test_load_bait_data_iris_basic_structure():
    df = pd.DataFrame({
        "Protein": ["P1", "P2"],
        "condX_BaitA_1": [5, 1],
        "condX_BaitA_2": [7, 0],
        "condY_BaitA_1": [9, 2],
        "condX_BaitB_1": [3, 0],
    })

    bait_units, X_by_bait_unit, metadata = load_bait_data_iris(df)

    # IRIS bait_units are condition × bait
    assert set(bait_units) == {"condX_BaitA", "condY_BaitA", "condX_BaitB"}

    # Replicate matrices exist for each bait_unit
    assert "condX_BaitA" in X_by_bait_unit
    assert "condY_BaitA" in X_by_bait_unit
    assert "condX_BaitB" in X_by_bait_unit

    XA = X_by_bait_unit["condX_BaitA"]
    YA = X_by_bait_unit["condY_BaitA"]
    XB = X_by_bait_unit["condX_BaitB"]

    # Shapes (IRIS uses per-bait_unit replicate structure)
    assert XA.shape == (2, 2)   # two replicates for condX_BaitA
    assert YA.shape == (2, 1)   # only one replicate for condY_BaitA
    assert XB.shape == (2, 1)   # only one replicate for condX_BaitB

    # Values for condX_BaitA
    assert np.allclose(XA[:, 0], [5, 1])
    assert np.allclose(XA[:, 1], [7, 0])

    # Values for condY_BaitA
    assert np.allclose(YA[:, 0], [9, 2])

    # Values for condX_BaitB
    assert np.allclose(XB[:, 0], [3, 0])

    # Metadata invariants
    assert metadata["proteins_by_bait_unit"]["condX_BaitA"] == ["P1", "P2"]
    assert metadata["proteins_by_bait_unit"]["condY_BaitA"] == ["P1", "P2"]
    assert metadata["proteins_by_bait_unit"]["condX_BaitB"] == ["P1", "P2"]

    assert "bait_units" in metadata
    assert "proteins_by_bait_unit" in metadata
    assert "control_bait_units" in metadata
    assert "treatment_bait_units" in metadata
