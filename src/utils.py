# src/utils.py

"""Utility functions: data loading, metrics, and plotting helpers."""

import os
import numpy as np
import pandas as pd


def load_pkpd_data(path="data/sample_pkpd.csv"):
    """Load PK/PD IC50 dataset."""
    df = pd.read_csv(path)
    required = {"cell_line", "drug", "IC50"}
    if not required.issubset(df.columns):
        raise ValueError(f"CSV must contain columns {required}, got {set(df.columns)}")
    return df


def load_spatial_data(path="data/sample_spatial.csv"):
    """Load 3-D spatial concentration dataset."""
    df = pd.read_csv(path)
    required = {"x", "y", "z", "concentration"}
    if not required.issubset(df.columns):
        raise ValueError(f"CSV must contain columns {required}, got {set(df.columns)}")
    return df


def rmse(y_true, y_pred):
    """Root Mean Squared Error."""
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    return np.sqrt(np.mean((y_true - y_pred) ** 2))


def r_squared(y_true, y_pred):
    """Coefficient of determination (R²)."""
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    if ss_tot == 0:
        return 1.0
    return 1 - ss_res / ss_tot


def plot_concentration_profile(x, c, title="Drug Concentration Profile",
                                xlabel="Distance (mm)", ylabel="Concentration (µM)",
                                save_path=None):
    """Plot a 1-D concentration profile."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, c, "b-", linewidth=2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fig
