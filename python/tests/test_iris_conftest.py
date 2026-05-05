# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 08:34:21 2026

@author: Erik
"""
# %% Import
import pytest
import pandas as pd

from orion.methods.iris.pipeline import run_iris_pipeline

# %% fixtures
@pytest.fixture
def minimal_input_data():
    return pd.DataFrame(
        {
            "Protein": ["P1", "P2"],
            "condA_B1_1": [5, 1],
            "condA_B1_2": [3, 0],
            "condA_CTRL_1": [2, 0],
            "condA_CTRL_2": [1, 1],
        }
    )


@pytest.fixture
def iris_results(minimal_input_data):
    return run_iris_pipeline(
        input_data=minimal_input_data,
        max_iter=5,
        seed=1,
        make_plots=False,
    )


@pytest.fixture
def iris_metadata(iris_results):
    return iris_results["metadata"]


@pytest.fixture
def iris_raw(iris_results):
    return iris_results["raw_outputs"]


@pytest.fixture
def iris_df(iris_results):
    return iris_results["results_df"]