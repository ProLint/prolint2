## Installation

ProLint can be easily installed via `pip`. To ensure a smooth and isolated installation process, we recommend creating a dedicated conda environment. Please follow the steps below:

1. Create a new conda environment:
```bash
conda create -n prolint python=3.10
```

2. Activate the newly created environment:
```bash
conda activate prolint
```

2. Install `prolint` via `pip`:
```bash
python -m pip install prolint
```

That’s it! You have successfully installed `prolint` within a controlled and reproducible environment for your biomolecular interaction analyses.

#### Troubleshooting

:::{warning}
ProLint requires MDAnalysis >= 2.4.0. Older versions may have compatibility issues with certain trajectory formats.
:::

<!-- #### Web Dashboard

The web dashboard requires additional dependencies:

```bash
# Backend
pip install fastapi uvicorn python-multipart websockets

# Frontend (requires Node.js)
cd prolint/web/frontend
npm install -->
<!-- ``` -->

