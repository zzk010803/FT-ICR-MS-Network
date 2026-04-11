# Molecular Transformation Network Analysis for FT-ICR-MS Data

This repository contains a custom Python script for constructing and analyzing molecular transformation networks from FT-ICR-MS datasets. The workflow integrates molecular formula information with a PageRank-based network algorithm to identify key molecular species and possible transformation pathways in dissolved organic matter (DOM).

---

## Overview

This script takes FT-ICR-MS molecular formula data as input and performs the following steps:

1. Reads molecular formula and peak intensity data from a CSV file
2. Calculates monoisotopic masses for all assigned formulas
3. Defines a set of possible reaction transformations based on exact mass differences
4. Constructs a molecular transformation network
5. Applies the PageRank algorithm to rank formulas by their network connectivity
6. Computes correlations between PageRank scores, peak intensity, and molecular mass
7. Exports the network in `.gexf` format for visualization in **Gephi**

---

## Main Features

- Reaction-rule-based molecular transformation network construction
- Exact-mass matching between molecular formulas
- PageRank analysis for identifying influential molecular formulas
- Exportable network files for Gephi visualization
- Annotation of nodes with PageRank and peak intensity values
- Simple and reproducible workflow for DOM transformation analysis

---

## Dependencies

Please install the following Python packages before running the script:

- `pandas`
- `numpy`
- `scipy`
- `matplotlib`
- `pykrev`

Example:

```bash
pip install pandas numpy scipy matplotlib

