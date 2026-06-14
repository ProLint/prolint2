**ProLint** is a Python library for analyzing biomolecular interactions from molecular dynamics simulations. Built on [MDAnalysis](https://www.mdanalysis.org/), it offers great efficiency in contact calculations and trajectory manipulation, proposing a simple four-step workflow: load your simulation, compute contacts between user-defined atom groups, analyze the results, and generate publication-ready plots. In addition, it provides a React-based web dashboard for interactive exploration of interaction results.

#### Basic Example

```python
from prolint import Universe

# Load simulation
universe = Universe("topology.gro", "trajectory.xtc")

# Compute contacts (default: protein vs non-protein)
contacts = universe.compute_contacts(cutoff=7.0)

# Run analysis
result = contacts.analyze("timeseries", database_type="CHOL")

# Visualize
from prolint.plotting import plot
fig, ax = plot("heatmap", result, colorscheme="viridis")
fig.savefig('heatmap.png', dpi=150, bbox_inches="tight")
```

:::{tip}
By default, ProLint analyzes interactions between protein and non-protein atoms. Customize this by setting `universe.query` and `universe.database` to any valid MDAnalysis selection.
:::

#### Available Analysis Types

```{list-table}
:header-rows: 1
:align: left

* - Analysis
  - Use when
  - Plotters
* - `timeseries`
  - Track contact counts over trajectory time
  - `heatmap`, `timeseries`
* - `metrics`
  - Compare residues by occupancy, mean, max, or sum
  - `residue_metrics`, `logo_grid`
* - `density_map`
  - Visualize spatial distribution around query
  - `density_map`
* - `radial_density`
  - Measure density as function of distance from query
  - `radial_density`
* - `shared_contacts`
  - Find residues contacting the same molecules
  - `heatmap`, `network`
* - `database_contacts`
  - Track which molecules contact a specific residue
  - `database_contacts_heatmap`
* - `distances`
  - Track distance between residue pairs over time
  - `distance_timeseries`
* - `atom_distances`
  - Get atom-atom distance matrix at a frame
  - `distance_heatmap`
* - `kinetics`
  - Measure binding dynamics and residence times
  - `survival_curve`, `residence_distribution`, `contact_events`
```

See the [Analysis](autoapi/prolint/analysis/index) and [Plotting](autoapi/prolint/plotting/index) APIs for details.

#### Web Dashboard

For users who prefer a graphical interface, ProLint includes an interactive web dashboard. It allows you to upload simulation files, configure analysis parameters, and explore results, all without writing code. See the [Web Dashboard Tutorial](tutorials/web-dashboard) for an introduction to its usage.

#### Citation

If ProLint contributes to your research, please cite:

```bibtex
@article{prolint_2024,
  title={ProLint v. 2: An optimized tool for the analysis and visualization of lipid-protein interactions},
  author={Ramirez-Echemendia, Daniel P.; Sejdiu, Besian I.; Tieleman, D. Peter},
  journal={Biophysical Journal},
  year={2024}
}
```

#### License

ProLint is released under the [MIT License](https://opensource.org/license/MIT).
