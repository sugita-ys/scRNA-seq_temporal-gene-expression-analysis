import os
from pathlib import Path
import time

import streamlit as st
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import issparse
import requests


# =========================
# Config
# =========================
st.set_page_config(page_title="DevVIS Gene Viewer (L6a vs L6b)", layout="centered")

DATA_URL = "https://allen-developmental-mouse-atlas.s3.amazonaws.com/scRNA/DevVIS_scRNA_processed.h5ad"
DEFAULT_DATA_DIR = Path("data")
DEFAULT_H5_NAME = "DevVIS_scRNA_processed.h5ad"
DEFAULT_H5_PATH = DEFAULT_DATA_DIR / DEFAULT_H5_NAME


# =========================
# Helpers: download
# =========================
def _download_file(url: str, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)

    # Stream download with progress
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        chunk_size = 1024 * 1024  # 1MB

        progress = st.progress(0)
        status = st.empty()

        downloaded = 0
        t0 = time.time()

        with open(dst, "wb") as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total > 0:
                        progress.progress(min(downloaded / total, 1.0))
                    elapsed = max(time.time() - t0, 1e-6)
                    mb = downloaded / (1024 * 1024)
                    speed = mb / elapsed
                    if total > 0:
                        status.info(f"Downloading... {mb:.1f} MB / {total/(1024*1024):.1f} MB  ({speed:.1f} MB/s)")
                    else:
                        status.info(f"Downloading... {mb:.1f} MB  ({speed:.1f} MB/s)")

        progress.progress(1.0)
        status.success(f"Downloaded to: {dst}")


# =========================
# Data load (cached)
# =========================
@st.cache_resource(show_spinner=True)
def load_and_split(h5_path: str):
    adata_l6_all = sc.read_h5ad(h5_path)

    # L6b 
    l6b_labels = ["L6b/CT ENT Glut", "L6b CTX Glut"]
    adata_l6b = adata_l6_all[adata_l6_all.obs["subclass_label"].isin(l6b_labels)].copy()

    # L6（CT/IT）
    l6_labels = ["L6 CT CTX Glut", "L6 IT CTX Glut"]
    adata_l6 = adata_l6_all[adata_l6_all.obs["subclass_label"].isin(l6_labels)].copy()

    return adata_l6_all, adata_l6, adata_l6b


# =========================
# Your analysis functions
# =========================
def _gene_age_stats(adata, gene, birth_e=19.5):
    adata_plot = adata[~adata.obs["SynchronizedAge"].isna()].copy()
    if gene not in adata_plot.var_names:
        raise ValueError(f"{gene} isn't found in var_names")

    age_str = adata_plot.obs["SynchronizedAge"].astype(str)

    def age_to_num(s: str, birth_e=birth_e):
        if s.startswith("E"):
            val = float(s[1:])
            return val - birth_e      # E19.5 = 0
        elif s.startswith("P"):
            val = float(s[1:])
            return val
        else:
            return np.nan

    age_num = age_str.map(age_to_num)
    adata_plot.obs["age_label"] = age_str
    adata_plot.obs["age_num"] = age_num

    if issparse(adata_plot.X):
        expr_gene = adata_plot[:, gene].X.toarray().flatten()
    else:
        expr_gene = np.asarray(adata_plot[:, gene].X).flatten()

    df = pd.DataFrame({
        "age_label": adata_plot.obs["age_label"].values,
        "age_num": adata_plot.obs["age_num"].values,
        "expr": expr_gene,
    })

    stats = (
        df.groupby("age_label")
          .agg(mean=("expr", "mean"),
               std=("expr", "std"),
               count=("expr", "count"),
               age_num=("age_num", "first"))
          .sort_values("age_num")
    )
    stats["sem"] = stats["std"] / np.sqrt(stats["count"])
    return stats


def _make_plot(stats_l6: pd.DataFrame, stats_l6b: pd.DataFrame, gene: str):
    age_union = (
        pd.concat([
            stats_l6[["age_num"]].reset_index(),
            stats_l6b[["age_num"]].reset_index()
        ])
        .drop_duplicates(subset="age_label")
        .sort_values("age_num")
    )

    x = age_union["age_num"].values
    x_labels = age_union["age_label"].tolist()

    mean_l6  = [stats_l6.loc[l, "mean"]  if l in stats_l6.index  else np.nan for l in x_labels]
    sem_l6   = [stats_l6.loc[l, "sem"]   if l in stats_l6.index  else np.nan for l in x_labels]
    mean_l6b = [stats_l6b.loc[l, "mean"] if l in stats_l6b.index else np.nan for l in x_labels]
    sem_l6b  = [stats_l6b.loc[l, "sem"]  if l in stats_l6b.index else np.nan for l in x_labels]

    fig = plt.figure(figsize=(6, 4))
    plt.fill_between(x, np.array(mean_l6b) - np.array(sem_l6b), np.array(mean_l6b) + np.array(sem_l6b), alpha=0.2)
    plt.plot(x, mean_l6b, "-o", label="L6b")

    plt.fill_between(x, np.array(mean_l6) - np.array(sem_l6), np.array(mean_l6) + np.array(sem_l6), alpha=0.2)
    plt.plot(x, mean_l6, "-o", label="L6a")

    plt.xticks(x, x_labels, rotation=45, fontsize=7)
    plt.xlabel("Age")
    plt.ylabel(f"Mean expression of {gene}")
    plt.title(f"{gene} expression across developmental stages (L6a vs L6b)")
    plt.legend()
    plt.tight_layout()
    return fig


# =========================
# UI
# =========================
st.title("DevVIS Gene Viewer (L6a vs L6b)")
st.caption("Search Gene, you can see the average expression (mean±SEM) difference between L6a and L6b by SynchronizedAge.")

with st.expander("Data source (open data)"):
    st.code(DATA_URL, language="text")
    st.write("Paper (Nature):")
    st.code("https://doi.org/10.1038/s41586-025-09644-1", language="text")

# Determine data path
env_path = os.environ.get("DEVVIS_H5AD", "").strip()
candidate_paths = []
if env_path:
    candidate_paths.append(Path(env_path))
candidate_paths.append(DEFAULT_H5_PATH)

existing = next((p for p in candidate_paths if p.exists()), None)

colA, colB = st.columns([2, 1])
with colA:
    h5_path = st.text_input(
        "h5ad path (DEVVIS_H5AD env var supported)",
        value=str(existing) if existing else str(DEFAULT_H5_PATH),
    )
with colB:
    st.write("")
    st.write("")
    if st.button("Download h5ad"):
        dst = Path(h5_path)
        if dst.is_dir():
            dst = dst / DEFAULT_H5_NAME
        st.info("Starting download (this file is large)...")
        try:
            _download_file(DATA_URL, dst)
        except Exception as e:
            st.error(f"Download failed: {e}")

# Load
if not Path(h5_path).exists():
    st.warning("h5ad file not found yet. If you haven't downloaded it, click **Download h5ad**.")
    st.stop()

with st.spinner("Loading h5ad (this can take a while on first run)..."):
    adata_all, adata_l6, adata_l6b = load_and_split(h5_path)

st.success(f"Loaded: all={adata_all.n_obs:,} cells | L6={adata_l6.n_obs:,} | L6b={adata_l6b.n_obs:,}")

# Gene input
gene = st.selectbox("Gene (from var_names)", options=list(adata_all.var_names), index=0)
manual = st.checkbox("Write gene name manually", value=False)
if manual:
    gene = st.text_input("Gene symbol", value=gene).strip()

birth_e = st.number_input("birth_e (E?? = 0 base: prefer E19.5)", value=19.5, step=0.5)
run = st.button("Plot")

if run:
    try:
        stats_l6  = _gene_age_stats(adata_l6, gene, birth_e=birth_e)
        stats_l6b = _gene_age_stats(adata_l6b, gene, birth_e=birth_e)

        fig = _make_plot(stats_l6, stats_l6b, gene)
        st.pyplot(fig)

        with st.expander("stats table (L6)"):
            st.dataframe(stats_l6)

        with st.expander("stats table (L6b)"):
            st.dataframe(stats_l6b)

    except Exception as e:
        st.error(str(e))
