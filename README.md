# CMML_mini2_final

## Project Structure

```
CMML_mini2_final/
├── bench/: run under metric.yml
│   │── Horizontal_calc.ipynb: benchmark file for horizontal integration
│   └── Vertical_calc.ipynb: benchmark file for vertical integration
├── bench_res/: store bench result
├── data/: store data in 10x Genomics format, barcodes.csv, matrix.mtx, features.tsv
├── envs/
│   │── metric.yml: yml for benchmark
│   └── mira.yml: yml for MIRA
├── run/: run integration
│   └── horizontal
│       │── RunMIRA.ipynb: run MIRA in horizontal, jupyter notebook
│       │── RunMIRA.py: run MIRA in horizontal, terminal
│       │── RunPCA.R: run PCA in horizontal
│       └── RunSeurat.R: run Seurat in horizontal
│   └── Vertical
│       │── MIRA_run.ipynb: run MIRA in vertical, jupyter notebook
│       │── MIRA_run.py: run MIRA in vertical, terminal
│       │── PCA_run.R: run PCA in vertical
│       └── Seurat_run.R: run Seurat in vertical
├── run_res/ : store result
└── README.md
```

## Getting Started

1. Clone the repository.
2. Install environment and register
     ```bash
     conda env create -f envs/mira.yml
     conda activate mira-env
     python -m ipykernel install --user --name=mira-env --display-name="mira-env"
     ```
3. Run the main script (guidance in code)

