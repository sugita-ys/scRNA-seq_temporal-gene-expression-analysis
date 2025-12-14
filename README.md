# scRNA-seq_temporal-gene-expression-analysis

# DevVIS Gene Viewer (Streamlit)

A tiny Streamlit app to visualize gene expression across developmental stages (SynchronizedAge) comparing **L6a vs L6b** neurons.

## Data
This app uses the open processed h5ad dataset:
Gao, Y., van Velthoven, C.T.J., Lee, C. et al. Continuous cell-type diversification in mouse visual cortex development. Nature 647, 127–142 (2025). 
https://doi.org/10.1038/s41586-025-09644-1

https://allen-developmental-mouse-atlas.s3.amazonaws.com/scRNA/DevVIS_scRNA_processed.h5ad

## Setup
### Create environment (optional)
conda create -n devvis-gene-viewer python=3.12 -y
conda activate devvis-gene-viewer
### Install deps
pip install -r requirements.txt

## Run
### Download within the app (recommended)
python -m streamlit run app.py
