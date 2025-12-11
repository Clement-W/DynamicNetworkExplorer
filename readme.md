# Dynamic Network Explorer (Dash + Cytoscape)

Interactive network explorer for gene–gene interactions across time points.  

## Features

- Visualize time-varying gene–gene networks
- Search and highlight a gene with its directed connections
- Sign-aware edge coloring (positive/negative)
- TF highlighting (pre-defined mouse/human lists)

## 0) Get the code

Clone the repository:

**SSH**

```bash
git clone git@github.com:Clement-W/DynamicNetworkExplorer.git
cd DynamicNetworkExplorer
```

**HTTPS**

```bash
git clone https://github.com/Clement-W/DynamicNetworkExplorer.git
cd DynamicNetworkExplorer
```

> Alternatively, download the repo as a ZIP, extract it, then open a terminal in the extracted folder.

## 1) Prerequisites

- **Python 3.11**
- One of:
  - **Conda/Miniconda**, or
  - **venv** + **pip**

## 2) Installation (choose Conda **or** venv)

### A) Conda (recommended)

```bash
# from the project root
conda env create -f environment.yml          
conda activate dynamic-network-explorer
python app.py -h
```

### B) venv + pip

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -U pip wheel
pip install -r requirements.txt
python app.py -h
```

## 3) Data layout & required input

Place CSVs under `data/` (or set `--path`). Filenames matched by `--pattern` become time-point labels.

Each CSV must contain:

- `gene1` (str)
- `gene2` (str)
- `strength` (float; absolute value used for ranking)
- `sign` (int; **1** = positive, **-1** = negative)

The pre-defined lists of TFs are here for highlighting purpose in the visualization. TFs for human and mouse are already included in `assets/`:

- Mus_musculus_TF downloaded July 9th 2025 (https://guolab.wchscu.cn/AnimalTFDB4/#/)
- human_TF_names_v_1.01.txt downloaded July 9th 2025 (https://humantfs.ccbr.utoronto.ca/download.php)

## 4) Running the app

Available parameters:

```bash
python app.py -h
```

```
options:
  -h, --help              Show this help message and exit
  --topn TOPN             Number of top interactions to display per time-point (default: 5000)
  --path PATH             Path to the directory containing the CSV files (default: data)
  --pattern PATTERN       File pattern to match CSV files (default: ranked_interactions_*.csv)
  --scale SCALE           Scale factor for node layout (default: 1200). Increase/Decrease if you increase/decrease topn.
  --organism ORGANISM     Organism for TF highlighting (default: mouse). Options: mouse, human, none
  --subnetwork SUBNETWORK Comma-separated list of source genes (gene1) to restrict the network to up to the topn interactions for each source (e.g. --subnetwork Tcf21,Zeb2. If empty, the full network is used)
```

Examples:

```bash
# Default data folder; 6000 edges per time point
python app.py --topn 6000
```

The server starts at `http://127.0.0.1:8050/`.

## 5) Controls

- **Slider:** switch time points
- **Dropdown:** search & highlight a gene
- **Click node:** highlight neighborhood for current time point
- **Click background:** clear highlight
- **Zoom/Pan:** mouse/trackpad

## 6) Additional information

To explore all the links without restricting to top N links, run `explore_csv.ipynb`.
