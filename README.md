# ORION

A unified repository for the ORION project, containing:

- The ORION R package (package name: orion), stored under the Rpkg/ directory
- The ORION Python package (python/orion), which contains the methods/ submodule with SAINT and IRIS
- Shared datasets, documentation, and analysis workflows

This repo provides a coherent structure for scientific modeling, data management, and reproducible analysis across R and Python.

---

## Repository structure

ORION/
├── Rpkg/                       # R package source (package name: orion)
│   ├── R/
│   ├── man/
│   ├── data/
│   ├── DESCRIPTION
│   ├── NAMESPACE
│   └── LICENSE
│
├── python/
│   └── orion/                  # Python ORION package
│       ├── methods/            # Methods subpackage
│       │   ├── iris/           # IRIS method implementation
│       │   └── saint/          # SAINT method implementation
│       ├── tests/
│       ├── pyproject.toml
│       └── README.md
│
├── data/                       # Raw or shared datasets (not package-internal)
├── analysis/                   # Rmds, notebooks, exploratory work
└── README.md                   # This file

---

## ORION R package (R)

The R side of ORION (package name: orion) is a lightweight container for:

- Experiment metadata
- Processed datasets
- Documentation used in downstream analyses

It intentionally contains no modeling code.

Install locally from the repo root with:

```
devtools::install("Rpkg")
```

---

## ORION Python package (methods: SAINT and IRIS)

The Python side lives under python/orion and currently provides the methods:

- SAINT (python/orion/methods/saint)
- IRIS (python/orion/methods/iris)

Install locally in editable mode from the repo root with:

```
pip install -e python
```

---

## Development workflow

- R package development happens inside Rpkg/ (package: orion)
- Python method development happens inside python/orion/methods/ (subpackages: saint, iris)
- Shared data lives in data/
- Analysis notebooks live in analysis/

Both sides share the same project identity: ORION.

---

## License

The project is licensed under the terms specified in the R and Python components.
