# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 10:59:07 2026

@author: Erik
"""
# %% Import
from dataclasses import dataclass, field
from typing import Dict, List, Any

# %% Dataclasses
@dataclass
class UserProvidedFields:
    biological_bait_names: Dict[str, str] = field(default_factory=dict)
    AN: Dict[str, str] = field(default_factory=dict)   # Accession Number
    MW: Dict[str, float] = field(default_factory=dict)
    extra_fields: Dict[str, Any] = field(default_factory=dict)


@dataclass
class InferredFields:
    bait_units: List[str]
    proteins_by_bait_unit: Dict[str, List[str]]
    replicate_map: Dict[str, Any] = field(default_factory=dict)
    conditions: Dict[str, Any] = field(default_factory=dict)
    negative_controls_inferred: List[str] = field(default_factory=list)
    control_bait_units: List[str] = field(default_factory=list)
    treatment_bait_units: List[str] = field(default_factory=list)
    extra_fields: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PipelineDerivedFields:
    bait_units: List[str] = field(default_factory=list)
    hyperparameters: Dict[str, Any] = field(default_factory=dict)
    tau_grid: List[float] = field(default_factory=list)
    convergence: Dict[str, Any] = field(default_factory=dict)
    iteration_counts: Dict[str, int] = field(default_factory=dict)
    extra_fields: Dict[str, Any] = field(default_factory=dict)


@dataclass
class HierarchicalMetadata:
    user_provided_fields: UserProvidedFields
    inferred_fields: InferredFields
    pipeline_derived_fields: PipelineDerivedFields