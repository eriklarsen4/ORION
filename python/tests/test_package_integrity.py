# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 09:37:10 2026

@author: Erik
"""

# %% Imports

import importlib
import os
import pkgutil
import orion


# %% Helpers

def _module_exists(module_path):
    try:
        importlib.import_module(module_path)
        return True
    except Exception:
        return False


# %% Tests

def test_all_expected_modules_importable():
    """
    Ensure that all core modules in the ORION package import cleanly.
    This catches missing files, missing __init__.py files, and broken imports.
    """

    expected_modules = [
        "orion.methods.iris.em_wrapper",
        "orion.methods.iris.responsibilities",
        "orion.methods.iris.likelihood",
        "orion.methods.iris.helpers",
        "orion.methods.iris.init_tau",
        "orion.methods.iris.tau_grid",
        "orion.methods.iris.pipeline",
        "orion.methods.iris.diagnostics",
        "orion.methods.iris.diagnostics_tau_grid",
        "orion.methods.saint.em_wrapper",
        "orion.methods.saint.responsibilities",
        "orion.methods.saint.likelihood",
        "orion.methods.saint.pipeline",
        "orion.methods.saint.diagnostics",
        "orion.utils.metadata_types",
        "orion.utils.io.iris_loader",
        "orion.utils.io.saint_loader",
        "orion.utils.updates.lambda_updates",
        "orion.utils.updates.pi_updates",
        "orion.utils.updates.tau_updates",
        "orion.utils.metadata_types",
        "orion.utils.validation.validate_hyperparams"
    ]

    for module in expected_modules:
        assert _module_exists(module), f"Module failed to import: {module}"


def test_all_directories_have_init_files():
    """
    Ensure every directory in the orion package contains an __init__.py file.
    This enforces stable imports.
    """

    orion_path = os.path.dirname(orion.__file__)

    for root, dirs, files in os.walk(orion_path):
        # Skip hidden directories like __pycache__
        dirs[:] = [d for d in dirs if not d.startswith("__")]

        if root.endswith("orion"):
            continue

        assert "__init__.py" in files, f"Missing __init__.py in {root}"


def test_no_broken_submodules():
    """
    Ensure that pkgutil can walk the entire orion package without errors.
    This catches namespace issues and missing modules.
    """

    for module_info in pkgutil.walk_packages(orion.__path__, prefix="orion."):
        module_name = module_info.name
        assert _module_exists(module_name), f"Broken submodule: {module_name}"
