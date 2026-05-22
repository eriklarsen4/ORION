# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 10:59:07 2026

@author: Erik
"""
# %% Import
from dataclasses import dataclass, field
from typing import Dict, List, Any

# %% Define data classes
@dataclass
class UserProvidedFields:
    biological_bait_names: Dict[str, str] = field(default_factory=dict)
    AN: Dict[str, str] = field(default_factory=dict)
    MW: Dict[str, float] = field(default_factory=dict)
    extra_fields: Dict[str, Any] = field(default_factory=dict)


@dataclass
class InferredFields:
    bait_units: list[str]
    proteins_by_bait_unit: dict[str, list[str]]
    control_bait_units: list[str]
    treatment_bait_units: list[str]

    # IRIS‑specific metadata here; SAINT ignores it
    extra_fields: dict = field(default_factory=dict)

    @property
    def baits(self):
        return self.bait_units


@dataclass
class PipelineDerivedFields:
    bait_units: List[str] = field(default_factory=list)
    hyperparameters: Dict[str, Any] = field(default_factory=dict)
    tau_grid: List[float] = field(default_factory=list)
    convergence: Dict[str, Any] = field(default_factory=dict)
    iteration_counts: Dict[str, int] = field(default_factory=dict)
    extra_fields: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PipelineMetadata:
    user_provided_fields: UserProvidedFields
    inferred_fields: InferredFields
    pipeline_derived_fields: PipelineDerivedFields
    
    @property
    def baits(self) -> List[str]:
        return self.inferred_fields.bait_units
