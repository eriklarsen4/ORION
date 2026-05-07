# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:20:41 2026

@author: Erik
"""

# %% Imports
import numpy as np
import pytest

from orion.utils.validation.validate_hyperparams import validate_hyperparams


# %% Tests

def test_validate_hyperparams_accepts_valid_input():
    max_iter = 10
    tol = 1e-6

    lambda_vecs = [
        np.array([1.0, 2.0]),
        np.array([1.5, 2.5]),
        np.array([2.0, 3.0]),
    ]

    pi = np.array([0.3, 0.4, 0.3])

    # Should not raise
    validate_hyperparams(max_iter, tol, lambda_vecs, pi)


def test_validate_hyperparams_rejects_non_integer_max_iter():
    with pytest.raises(ValueError):
        validate_hyperparams(
            max_iter=5.5,
            tol=1e-6,
            lambda_vecs=[np.array([1]), np.array([1]), np.array([1])],
            pi=np.array([1/3, 1/3, 1/3]),
        )


def test_validate_hyperparams_rejects_non_positive_max_iter():
    with pytest.raises(ValueError):
        validate_hyperparams(
            max_iter=0,
            tol=1e-6,
            lambda_vecs=[np.array([1]), np.array([1]), np.array([1])],
            pi=np.array([1/3, 1/3, 1/3]),
        )


def test_validate_hyperparams_rejects_non_float_tol():
    with pytest.raises(ValueError):
        validate_hyperparams(
            max_iter=10,
            tol="not-a-float",
            lambda_vecs=[np.array([1]), np.array([1]), np.array([1])],
            pi=np.array([1/3, 1/3, 1/3]),
        )


def test_validate_hyperparams_rejects_non_positive_tol():
    with pytest.raises(ValueError):
        validate_hyperparams(
            max_iter=10,
            tol=0.0,
            lambda_vecs=[np.array([1]), np.array([1]), np.array([1])],
            pi=np.array([1/3, 1/3, 1/3]),
        )


def test_validate_hyperparams_rejects_wrong_number_of_lambda_vectors():
    with pytest.raises(ValueError):
        validate_hyperparams(
            max_iter=10,
            tol=1e-6,
            lambda_vecs=[np.array([1]), np.array([1])],  # only 2
            pi=np.array([1/3, 1/3, 1/3]),
        )


def test_validate_hyperparams_rejects_mismatched_lambda_lengths():
    with pytest.raises(ValueError):
        validate_hyperparams(
            max_iter=10,
            tol=1e-6,
            lambda_vecs=[
                np.array([1.0, 2.0]),
                np.array([1.0, 2.0]),
                np.array([1.0]),  # wrong length
            ],
            pi=np.array([1/3, 1/3, 1/3]),
        )


def test_validate_hyperparams_rejects_non_numpy_lambda():
    with pytest.raises(ValueError):
        validate_hyperparams(
            max_iter=10,
            tol=1e-6,
            lambda_vecs=[1, 2, 3],  # not arrays
            pi=np.array([1/3, 1/3, 1/3]),
        )


def test_validate_hyperparams_rejects_wrong_pi_shape():
    with pytest.raises(ValueError):
        validate_hyperparams(
            max_iter=10,
            tol=1e-6,
            lambda_vecs=[np.array([1]), np.array([1]), np.array([1])],
            pi=np.array([0.5, 0.5]),  # wrong length
        )


def test_validate_hyperparams_rejects_negative_pi():
    with pytest.raises(ValueError):
        validate_hyperparams(
            max_iter=10,
            tol=1e-6,
            lambda_vecs=[np.array([1]), np.array([1]), np.array([1])],
            pi=np.array([0.5, -0.2, 0.7]),
        )


def test_validate_hyperparams_rejects_pi_not_summing_to_one():
    with pytest.raises(ValueError):
        validate_hyperparams(
            max_iter=10,
            tol=1e-6,
            lambda_vecs=[np.array([1]), np.array([1]), np.array([1])],
            pi=np.array([0.5, 0.5, 0.2]),  # sums to 1.2
        )