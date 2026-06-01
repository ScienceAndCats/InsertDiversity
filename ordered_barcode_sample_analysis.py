#!/usr/bin/env python3
"""
ordered_barcode_sample_analysis.py

Analyze barcode count CSVs across an ordered set of samples/treatments/timepoints.

Input CSVs should contain at least:
    barcode,count

Behavior:
    - Files are analyzed in alphabetical order unless explicit_files is provided.
    - A barcode is retained if it has >= min_count_to_call_barcode_present reads in
      at least min_samples_present_to_retain_barcode samples.
    - If a retained barcode is absent from another sample, it is kept as 0 reads.
    - Outputs ranked tables for largest increase, largest decrease, highest before,
      and highest after.
    - Outputs presence/absence differences.
    - Outputs all-barcode trajectories and a gain/loss-ratio heatmap.

Run:
    python ordered_barcode_sample_analysis.py pipeline_config.json
"""

from __future__ import annotations

import argparse
import concurrent.futures
import math
import re
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config_helpers import available_threads, configured_threads, load_script_config


def load_config(config_path: Path) -> dict:
    cfg = load_script_config(config_path, "ordered_barcode_sample_analysis")

    defaults = {
        "input_dir": ".",
        "output_dir": "ordered_barcode_analysis_outputs",
        "csv_pattern": "*.csv",
        "explicit_files": [],
        "sample_name_remove_regex": "_barcode_counts\\.csv$",
        "barcode_column": "barcode",
        "count_column": "count",
        "existing_real_column": "is_real_barcode",
        "len_minus_expected_column": "len_minus_expected",
        "expected_len_column": "expected_barcode_len",
        "barcode_len_column": "barcode_len",
        "use_existing_real_column": False,
        "require_expected_length": True,
        "expected_barcode_length": 25,
        "min_count_to_call_barcode_present": 3,
        "min_samples_present_to_retain_barcode": 1,
        "min_total_count_across_all_samples": 0,
        "frequency_denominator": "all_barcodes",
        "before_sample_indices": [],
        "after_sample_indices": [],
        "before_n_samples": 1,
        "after_n_samples": 1,
        "pseudocount_frequency": 1e-9,
        "top_n_for_each_ranked_output": 500,
        "top_n_to_plot": 50,
        "max_barcodes_in_heatmap": 200,
        "make_plots": True,
        "make_heatmaps": True,
        "make_first_vs_last_scatter": True,
        "make_top_trajectory_plots": True,
        "make_all_barcode_trajectory_plot": True,
        "all_trajectory_min_count_in_any_sample": 3,
        "all_trajectory_max_barcodes_to_plot": 5000,
        "make_gain_loss_ratio_heatmap": True,
        "gain_loss_heatmap_min_count_in_any_sample": 3,
        "gain_loss_heatmap_max_barcodes": 200,
        "gain_loss_heatmap_sort": "last_vs_first",
        "gain_loss_heatmap_show_numbers": True,
        "gain_loss_heatmap_number_decimals": 1,
        "log_scale_y_axis": True,
        "pseudocount_for_log_plots": 1e-9,
        "dpi": 300,
    }

    for k, v in defaults.items():
        cfg.setdefault(k, v)

    if cfg["frequency_denominator"] not in {"all_barcodes", "retained_barcodes"}:
        raise ValueError("frequency_denominator must be 'all_barcodes' or 'retained_barcodes'")

    if cfg["gain_loss_heatmap_sort"] not in {"last_vs_first", "largest_absolute_change", "max_count"}:
        raise ValueError(
            "gain_loss_heatmap_sort must be last_vs_first, largest_absolute_change, or max_count"
        )

    return cfg


def parse_bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)

    return s.astype(str).str.strip().str.lower().isin(
        ["true", "t", "1", "yes", "y"]
    )


def sample_name_from_file(path: Path, remove_regex: str) -> str:
    name = path.name
    if remove_regex:
        name = re.sub(remove_regex, "", name)
    return re.sub(r"\.csv$", "", name)


def discover_files(cfg: dict) -> List[Path]:
    if cfg.get("explicit_files"):
        files = [Path(p).expanduser().resolve() for p in cfg["explicit_files"]]
    else:
        input_dir = Path(cfg["input_dir"]).expanduser().resolve()
        files = sorted(input_dir.glob(cfg["csv_pattern"]))

    if len(files) < 2:
        raise SystemExit(f"Found {len(files)} file(s). Need at least 2 CSV files.")

    missing = [str(p) for p in files if not p.exists()]
    if missing:
        raise FileNotFoundError("These input files do not exist:\n" + "\n".join(missing))

    return files


def read_one_csv(path: Path, cfg: dict) -> Tuple[str, pd.DataFrame, dict]:
    barcode_col = cfg["barcode_column"]
    count_col = cfg["count_column"]

    df = pd.read_csv(path)

    if barcode_col not in df.columns:
        raise ValueError(f"{path.name} does not contain barcode column '{barcode_col}'")
    if count_col not in df.columns:
        raise ValueError(f"{path.name} does not contain count column '{count_col}'")

    df = df.copy()
    df[barcode_col] = df[barcode_col].astype(str)
    df[count_col] = pd.to_numeric(df[count_col], errors="coerce").fillna(0).astype(int)

    grouped = df.groupby(barcode_col, as_index=False)[count_col].sum()

    meta_cols = [c for c in df.columns if c != count_col]
    meta = df[meta_cols].drop_duplicates(subset=[barcode_col], keep="first")
    df = grouped.merge(meta, on=barcode_col, how="left")

    sample_name = sample_name_from_file(path, cfg["sample_name_remove_regex"])

    count_present = df[count_col] >= int(cfg["min_count_to_call_barcode_present"])
    usable = count_present.copy()

    if cfg["require_expected_length"]:
        len_minus_col = cfg["len_minus_expected_column"]
        barcode_len_col = cfg["barcode_len_column"]
        expected_len_col = cfg["expected_len_column"]
        expected_len_config = cfg.get("expected_barcode_length")

        if len_minus_col in df.columns:
            length_ok = pd.to_numeric(
                df[len_minus_col], errors="coerce"
            ).fillna(999999).eq(0)

        elif barcode_len_col in df.columns and expected_len_col in df.columns:
            length_ok = pd.to_numeric(df[barcode_len_col], errors="coerce").eq(
                pd.to_numeric(df[expected_len_col], errors="coerce")
            )

        elif expected_len_config is not None:
            length_ok = df[barcode_col].str.len().eq(int(expected_len_config))

        else:
            raise ValueError(
                f"{path.name}: no usable length column and require_expected_length=True"
            )

        usable &= length_ok

    if cfg["use_existing_real_column"]:
        real_col = cfg["existing_real_column"]
        if real_col not in df.columns:
            raise ValueError(f"{path.name}: missing real barcode column '{real_col}'")
        usable &= parse_bool_series(df[real_col])

    df["passes_per_sample_filter"] = usable

    summary = {
        "sample": sample_name,
        "file": str(path),
        "total_raw_reads": int(df[count_col].sum()),
        "unique_barcodes_raw_count_gt_0": int((df[count_col] > 0).sum()),
        "unique_barcodes_count_ge_present_cutoff": int(count_present.sum()),
        "unique_barcodes_passing_per_sample_filter": int(usable.sum()),
        "reads_in_barcodes_passing_per_sample_filter": int(
            df.loc[usable, count_col].sum()
        ),
    }

    return sample_name, df[[barcode_col, count_col, "passes_per_sample_filter"]], summary


def build_matrices(files: List[Path], cfg: dict):
    barcode_col = cfg["barcode_column"]
    count_col = cfg["count_column"]

    sample_order: List[str] = []
    per_sample: Dict[str, pd.DataFrame] = {}
    summaries: List[dict] = []

    max_workers = configured_threads(cfg, workload_size=len(files))
    print(f"Using {max_workers} thread(s) for {len(files)} ordered-analysis CSV file(s) (available: {available_threads()}).")

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_index = {
            executor.submit(read_one_csv, path, cfg): index
            for index, path in enumerate(files)
        }
        ordered_results = [None] * len(files)
        for future in concurrent.futures.as_completed(future_to_index):
            ordered_results[future_to_index[future]] = future.result()

    for sample, df, summary in ordered_results:
        if sample in per_sample:
            raise ValueError(f"Duplicate sample name after regex cleanup: {sample}")

        sample_order.append(sample)
        per_sample[sample] = df
        summaries.append(summary)

    all_barcodes = sorted(set().union(*[set(df[barcode_col]) for df in per_sample.values()]))

    raw_counts = pd.DataFrame(index=all_barcodes)
    present = pd.DataFrame(index=all_barcodes)

    for sample in sample_order:
        df = per_sample[sample]

        s_counts = df.set_index(barcode_col)[count_col]
        s_present = df.set_index(barcode_col)["passes_per_sample_filter"]

        raw_counts[sample] = s_counts.reindex(all_barcodes).fillna(0).astype(int)
        present[sample] = s_present.reindex(all_barcodes).fillna(False).astype(bool)

    n_present = present.sum(axis=1)
    total_count = raw_counts.sum(axis=1)

    keep = (
        (n_present >= int(cfg["min_samples_present_to_retain_barcode"]))
        & (total_count >= int(cfg["min_total_count_across_all_samples"]))
    )

    barcode_filter_summary = pd.DataFrame(
        {
            "barcode": raw_counts.index,
            "total_count_across_all_samples": total_count.values,
            "n_samples_present_after_filter": n_present.values,
            "kept_for_analysis": keep.values,
        }
    ).sort_values(
        [
            "kept_for_analysis",
            "n_samples_present_after_filter",
            "total_count_across_all_samples",
        ],
        ascending=[False, False, False],
    )

    counts = raw_counts.loc[keep].copy()
    present = present.loc[keep].copy()

    if cfg["frequency_denominator"] == "all_barcodes":
        denominators = raw_counts.sum(axis=0).replace(0, np.nan)
    else:
        denominators = counts.sum(axis=0).replace(0, np.nan)

    freqs = counts.div(denominators, axis=1).fillna(0)

    counts.index.name = "barcode"
    freqs.index.name = "barcode"
    present.index.name = "barcode"

    return sample_order, counts, freqs, present, barcode_filter_summary, pd.DataFrame(summaries)


def get_group_samples(sample_order: List[str], cfg: dict, group: str) -> List[str]:
    if group == "before":
        indices = cfg.get("before_sample_indices", [])
        n = int(cfg["before_n_samples"])
        default = list(range(n))

    elif group == "after":
        indices = cfg.get("after_sample_indices", [])
        n = int(cfg["after_n_samples"])
        default = list(range(len(sample_order) - n, len(sample_order)))

    else:
        raise ValueError(group)

    if not indices:
        indices = default

    bad = [i for i in indices if i < 0 or i >= len(sample_order)]
    if bad:
        raise ValueError(
            f"Bad {group}_sample_indices: {bad}. There are {len(sample_order)} samples."
        )

    return [sample_order[i] for i in indices]


def compute_stats(
    sample_order: List[str],
    counts: pd.DataFrame,
    freqs: pd.DataFrame,
    present: pd.DataFrame,
    cfg: dict,
) -> pd.DataFrame:
    before_samples = get_group_samples(sample_order, cfg, "before")
    after_samples = get_group_samples(sample_order, cfg, "after")

    pseudo = float(cfg["pseudocount_frequency"])

    first_sample = sample_order[0]
    last_sample = sample_order[-1]

    rows = []

    for barcode in freqs.index:
        f = freqs.loc[barcode]
        c = counts.loc[barcode]
        p = present.loc[barcode]

        first_freq = float(f[first_sample])
        last_freq = float(f[last_sample])

        first_count = int(c[first_sample])
        last_count = int(c[last_sample])

        diff = last_freq - first_freq
        log2fc = math.log2((last_freq + pseudo) / (first_freq + pseudo))

        max_freq_any = float(f.max())
        max_count_any = int(c.max())

        n_present = int(p.sum())
        n_observed = int((c > 0).sum())

        rows.append(
            {
                "barcode": barcode,
                "first_sample": first_sample,
                "last_sample": last_sample,
                "first_count": first_count,
                "last_count": last_count,
                "first_frequency": first_freq,
                "last_frequency": last_freq,
                "frequency_difference_first_to_last": diff,
                "log2_fold_change_first_to_last": log2fc,
                "abs_log2_fold_change_first_to_last": abs(log2fc),
                "max_frequency_any_sample": max_freq_any,
                "sample_at_max_frequency": str(f.idxmax()),
                "max_count_any_sample": max_count_any,
                "sample_at_max_count": str(c.idxmax()),
                "max_frequency_before": float(f[before_samples].max()),
                "sample_at_max_frequency_before": str(f[before_samples].idxmax()),
                "max_frequency_after": float(f[after_samples].max()),
                "sample_at_max_frequency_after": str(f[after_samples].idxmax()),
                "mean_frequency": float(f.mean()),
                "median_frequency": float(f.median()),
                "total_count_across_samples": int(c.sum()),
                "n_samples_present": n_present,
                "n_samples_observed_count_gt_0": n_observed,
                "present_in_all_samples": bool(n_present == len(sample_order)),
                "absent_or_below_cutoff_in_at_least_one_sample": bool(
                    n_present < len(sample_order)
                ),
                "variable_presence_absence": bool(
                    n_present > 0 and n_present < len(sample_order)
                ),
            }
        )

    return pd.DataFrame(rows)


def save_matrix_subset(
    outdir: Path,
    barcodes: List[str],
    counts: pd.DataFrame,
    freqs: pd.DataFrame,
    prefix: str,
):
    if barcodes:
        counts.loc[barcodes].to_csv(outdir / f"{prefix}_count_matrix.csv")
        freqs.loc[barcodes].to_csv(outdir / f"{prefix}_frequency_matrix.csv")
    else:
        pd.DataFrame().to_csv(outdir / f"{prefix}_count_matrix.csv")
        pd.DataFrame().to_csv(outdir / f"{prefix}_frequency_matrix.csv")


def save_trajectory_plot(
    barcodes: List[str],
    freqs: pd.DataFrame,
    outpath: Path,
    title: str,
    cfg: dict,
):
    if not barcodes:
        return

    samples = list(freqs.columns)
    x = np.arange(len(samples))
    pseudo = float(cfg["pseudocount_for_log_plots"])

    fig, ax = plt.subplots(figsize=(max(8, 0.45 * len(samples) + 4), 6))

    for bc in barcodes:
        ax.plot(
            x,
            freqs.loc[bc].values.astype(float) + pseudo,
            linewidth=1,
            alpha=0.7,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=90)
    ax.set_ylabel("Barcode frequency")

    if cfg.get("log_scale_y_axis", True):
        ax.set_yscale("log")
        ax.set_ylabel("Barcode frequency, log scale")

    ax.set_title(title)

    fig.tight_layout()
    fig.savefig(outpath, dpi=int(cfg["dpi"]))
    plt.close(fig)


def save_frequency_heatmap(
    barcodes: List[str],
    freqs: pd.DataFrame,
    outpath: Path,
    title: str,
    cfg: dict,
):
    if not barcodes:
        return

    barcodes = barcodes[: int(cfg["max_barcodes_in_heatmap"])]

    mat = freqs.loc[barcodes]
    pseudo = float(cfg["pseudocount_for_log_plots"])
    values = np.log10(mat.values.astype(float) + pseudo)

    fig_h = max(5, min(24, 0.18 * len(barcodes) + 3))
    fig_w = max(8, min(24, 0.45 * len(mat.columns) + 5))

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(values, aspect="auto")

    ax.set_title(title + "\nvalues = log10(frequency + pseudocount)")
    ax.set_xticks(np.arange(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=90)
    ax.set_yticks(np.arange(len(barcodes)))
    ax.set_yticklabels(barcodes, fontsize=5 if len(barcodes) > 50 else 7)

    fig.colorbar(im, ax=ax, shrink=0.85, label="log10 frequency")

    fig.tight_layout()
    fig.savefig(outpath, dpi=int(cfg["dpi"]))
    plt.close(fig)


def save_ranked_output(
    output_dir: Path,
    folder_name: str,
    stats: pd.DataFrame,
    sort_col: str,
    ascending: bool,
    counts: pd.DataFrame,
    freqs: pd.DataFrame,
    cfg: dict,
):
    outdir = output_dir / folder_name
    plotdir = outdir / "plots"

    outdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)

    ranked = stats.sort_values(sort_col, ascending=ascending).copy()

    top_n = int(cfg["top_n_for_each_ranked_output"])
    ranked_top = ranked.head(top_n)

    ranked.to_csv(outdir / f"{folder_name}_all_ranked_barcodes.csv", index=False)
    ranked_top.to_csv(outdir / f"{folder_name}_top_{top_n}.csv", index=False)

    barcodes = list(ranked_top["barcode"])

    save_matrix_subset(outdir, barcodes, counts, freqs, folder_name)

    if cfg.get("make_plots", True) and cfg.get("make_top_trajectory_plots", True):
        plot_barcodes = barcodes[: int(cfg["top_n_to_plot"])]
        save_trajectory_plot(
            plot_barcodes,
            freqs,
            plotdir / f"{folder_name}_top_trajectories.png",
            f"Top trajectories: {folder_name}",
            cfg,
        )

    if cfg.get("make_plots", True) and cfg.get("make_heatmaps", True):
        save_frequency_heatmap(
            barcodes,
            freqs,
            plotdir / f"{folder_name}_frequency_heatmap.png",
            f"Frequency heatmap: {folder_name}",
            cfg,
        )


def save_all_barcode_trajectory_plot(
    stats: pd.DataFrame,
    freqs: pd.DataFrame,
    counts: pd.DataFrame,
    output_dir: Path,
    cfg: dict,
):
    plotdir = output_dir / "plots"
    plotdir.mkdir(parents=True, exist_ok=True)

    min_count = int(cfg["all_trajectory_min_count_in_any_sample"])
    max_barcodes = int(cfg["all_trajectory_max_barcodes_to_plot"])

    selected = stats[stats["max_count_any_sample"] >= min_count].copy()

    selected = selected.sort_values(
        ["max_count_any_sample", "abs_log2_fold_change_first_to_last"],
        ascending=[False, False],
    )

    selected = selected.head(max_barcodes)

    barcodes = list(selected["barcode"])

    selected.to_csv(
        plotdir / "all_barcode_trajectories_plotted_barcodes.csv",
        index=False,
    )

    save_matrix_subset(plotdir, barcodes, counts, freqs, "all_barcode_trajectories")

    save_trajectory_plot(
        barcodes,
        freqs,
        plotdir / "all_barcode_trajectories.png",
        f"All barcode trajectories with max count >= {min_count}",
        cfg,
    )


def save_gain_loss_ratio_heatmap(
    stats: pd.DataFrame,
    freqs: pd.DataFrame,
    counts: pd.DataFrame,
    output_dir: Path,
    cfg: dict,
):
    plotdir = output_dir / "plots"
    plotdir.mkdir(parents=True, exist_ok=True)

    min_count = int(cfg["gain_loss_heatmap_min_count_in_any_sample"])
    max_barcodes = int(cfg["gain_loss_heatmap_max_barcodes"])
    pseudo = float(cfg["pseudocount_frequency"])

    first_sample = freqs.columns[0]

    selected = stats[stats["max_count_any_sample"] >= min_count].copy()

    if cfg["gain_loss_heatmap_sort"] == "last_vs_first":
        selected = selected.sort_values("log2_fold_change_first_to_last", ascending=False)

    elif cfg["gain_loss_heatmap_sort"] == "largest_absolute_change":
        selected = selected.sort_values(
            "abs_log2_fold_change_first_to_last",
            ascending=False,
        )

    elif cfg["gain_loss_heatmap_sort"] == "max_count":
        selected = selected.sort_values("max_count_any_sample", ascending=False)

    selected = selected.head(max_barcodes)

    barcodes = list(selected["barcode"])

    selected.to_csv(plotdir / "gain_loss_ratio_heatmap_barcodes.csv", index=False)

    if not barcodes:
        return

    sub_freq = freqs.loc[barcodes]
    baseline = sub_freq[first_sample].replace(0, 0.0)

    raw_ratio = sub_freq.add(pseudo).div(baseline + pseudo, axis=0)
    log2_ratio = np.log2(raw_ratio)

    raw_ratio.to_csv(plotdir / "gain_loss_raw_ratio_matrix.csv")
    log2_ratio.to_csv(plotdir / "gain_loss_log2_ratio_matrix.csv")
    counts.loc[barcodes].to_csv(plotdir / "gain_loss_heatmap_count_matrix.csv")

    n_rows, n_cols = log2_ratio.shape

    fig_w = max(8, min(28, 0.75 * n_cols + 5))
    fig_h = max(5, min(30, 0.24 * n_rows + 3))

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(log2_ratio.values.astype(float), aspect="auto")

    ax.set_title(
        f"Gain/loss ratio heatmap\n"
        f"log2((sample freq + pseudo) / ({first_sample} freq + pseudo))"
    )

    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(log2_ratio.columns, rotation=90)

    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(
        log2_ratio.index,
        fontsize=5 if n_rows > 50 else 7,
    )

    if cfg.get("gain_loss_heatmap_show_numbers", True) and n_rows * n_cols <= 2500:
        decimals = int(cfg["gain_loss_heatmap_number_decimals"])

        for i in range(n_rows):
            for j in range(n_cols):
                ax.text(
                    j,
                    i,
                    f"{log2_ratio.iat[i, j]:.{decimals}f}",
                    ha="center",
                    va="center",
                    fontsize=5,
                )

    fig.colorbar(im, ax=ax, shrink=0.85, label="log2 gain/loss ratio vs first sample")

    fig.tight_layout()
    fig.savefig(plotdir / "gain_loss_ratio_heatmap.png", dpi=int(cfg["dpi"]))
    plt.close(fig)


def save_first_vs_last_scatter(stats: pd.DataFrame, output_dir: Path, cfg: dict):
    if stats.empty:
        return

    plotdir = output_dir / "plots"
    plotdir.mkdir(parents=True, exist_ok=True)

    pseudo = float(cfg["pseudocount_for_log_plots"])

    fig, ax = plt.subplots(figsize=(7, 7))

    ax.scatter(
        stats["first_frequency"] + pseudo,
        stats["last_frequency"] + pseudo,
        s=10,
        alpha=0.5,
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("First sample frequency")
    ax.set_ylabel("Last sample frequency")
    ax.set_title("First vs last barcode frequency")

    fig.tight_layout()
    fig.savefig(plotdir / "first_vs_last_frequency_scatter.png", dpi=int(cfg["dpi"]))
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(
        stats["log2_fold_change_first_to_last"],
        stats["max_frequency_any_sample"] + pseudo,
        s=10,
        alpha=0.5,
    )

    ax.set_yscale("log")
    ax.set_xlabel("log2 fold-change, last vs first")
    ax.set_ylabel("Max frequency in any sample")
    ax.set_title("Barcode change vs maximum abundance")

    fig.tight_layout()
    fig.savefig(
        plotdir / "log2fc_first_to_last_vs_max_frequency.png",
        dpi=int(cfg["dpi"]),
    )
    plt.close(fig)


def write_summary(
    outpath: Path,
    cfg: dict,
    files: List[Path],
    sample_order: List[str],
    before_samples: List[str],
    after_samples: List[str],
    sample_summary: pd.DataFrame,
    stats: pd.DataFrame,
):
    n_variable = int(stats["variable_presence_absence"].sum()) if not stats.empty else 0

    with open(outpath, "w") as f:
        f.write("Ordered barcode sample analysis\n")
        f.write("===============================\n\n")

        f.write("Sample order used:\n")
        for i, (sample, path) in enumerate(zip(sample_order, files), start=0):
            f.write(f"  {i}: {sample}  ({path.name})\n")
        f.write("\n")

        f.write("Before/after groups used for max-frequency summaries:\n")
        f.write(f"  Before samples: {', '.join(before_samples)}\n")
        f.write(f"  After samples: {', '.join(after_samples)}\n\n")

        f.write("Barcode retention/filtering:\n")
        f.write(
            f"  - Barcode is called present in a sample if count >= "
            f"{cfg['min_count_to_call_barcode_present']}\n"
        )
        f.write(
            f"  - Barcode is retained if present in >= "
            f"{cfg['min_samples_present_to_retain_barcode']} sample(s)\n"
        )
        f.write("  - Missing barcode rows are treated as 0 reads in that sample.\n")
        f.write(f"  - Frequency denominator: {cfg['frequency_denominator']}\n\n")

        f.write("Main outputs:\n")
        f.write("  - all_barcode_ordered_sample_stats.csv\n")
        f.write("  - sample_barcode_counts_summary.csv\n")
        f.write("  - barcodes_with_presence_absence_differences.csv\n")
        f.write("  - largest_increase_first_to_last/\n")
        f.write("  - largest_decrease_first_to_last/\n")
        f.write("  - highest_frequency_before/\n")
        f.write("  - highest_frequency_after/\n")
        f.write("  - plots/all_barcode_trajectories.png\n")
        f.write("  - plots/gain_loss_ratio_heatmap.png\n\n")

        f.write(f"Retained barcodes analyzed: {len(stats):,}\n")
        f.write(f"Barcodes with presence/absence differences: {n_variable:,}\n")
        f.write(
            f"All-trajectory min count threshold: "
            f"{cfg['all_trajectory_min_count_in_any_sample']}\n"
        )
        f.write(
            f"Gain/loss heatmap min count threshold: "
            f"{cfg['gain_loss_heatmap_min_count_in_any_sample']}\n"
        )
        f.write(f"Gain/loss heatmap sort mode: {cfg['gain_loss_heatmap_sort']}\n\n")

        f.write("Per-sample barcode counts:\n")
        f.write(sample_summary.to_string(index=False))
        f.write("\n\n")

        if not stats.empty:
            cols = [
                "barcode",
                "first_frequency",
                "last_frequency",
                "frequency_difference_first_to_last",
                "log2_fold_change_first_to_last",
            ]

            f.write("Top 10 largest increases, first to last:\n")
            f.write(
                stats.sort_values(
                    "frequency_difference_first_to_last",
                    ascending=False,
                )
                .head(10)[cols]
                .to_string(index=False)
            )

            f.write("\n\nTop 10 largest decreases, first to last:\n")
            f.write(
                stats.sort_values(
                    "frequency_difference_first_to_last",
                    ascending=True,
                )
                .head(10)[cols]
                .to_string(index=False)
            )

            f.write("\n")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze barcode counts across ordered samples or treatments."
    )
    parser.add_argument("config", help="Path to config JSON file")
    args = parser.parse_args()

    cfg = load_config(Path(args.config))

    output_dir = Path(cfg["output_dir"]).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    files = discover_files(cfg)

    (
        sample_order,
        counts,
        freqs,
        present,
        barcode_filter_summary,
        sample_summary,
    ) = build_matrices(files, cfg)

    before_samples = get_group_samples(sample_order, cfg, "before")
    after_samples = get_group_samples(sample_order, cfg, "after")

    stats = compute_stats(sample_order, counts, freqs, present, cfg)

    counts.to_csv(output_dir / "retained_counts_matrix.csv")
    freqs.to_csv(output_dir / "retained_frequency_matrix.csv")
    present.to_csv(output_dir / "retained_presence_matrix.csv")
    barcode_filter_summary.to_csv(output_dir / "barcode_filter_summary.csv", index=False)
    sample_summary.to_csv(output_dir / "sample_barcode_counts_summary.csv", index=False)
    stats.to_csv(output_dir / "all_barcode_ordered_sample_stats.csv", index=False)

    variable = stats[stats["variable_presence_absence"]].copy()
    variable = variable.sort_values(
        ["n_samples_present", "max_count_any_sample"],
        ascending=[True, False],
    )
    variable.to_csv(
        output_dir / "barcodes_with_presence_absence_differences.csv",
        index=False,
    )

    save_ranked_output(
        output_dir,
        "largest_increase_first_to_last",
        stats,
        "frequency_difference_first_to_last",
        False,
        counts,
        freqs,
        cfg,
    )

    save_ranked_output(
        output_dir,
        "largest_decrease_first_to_last",
        stats,
        "frequency_difference_first_to_last",
        True,
        counts,
        freqs,
        cfg,
    )

    save_ranked_output(
        output_dir,
        "highest_frequency_before",
        stats,
        "max_frequency_before",
        False,
        counts,
        freqs,
        cfg,
    )

    save_ranked_output(
        output_dir,
        "highest_frequency_after",
        stats,
        "max_frequency_after",
        False,
        counts,
        freqs,
        cfg,
    )

    if cfg.get("make_plots", True) and cfg.get("make_first_vs_last_scatter", True):
        save_first_vs_last_scatter(stats, output_dir, cfg)

    if cfg.get("make_plots", True) and cfg.get("make_all_barcode_trajectory_plot", True):
        save_all_barcode_trajectory_plot(stats, freqs, counts, output_dir, cfg)

    if cfg.get("make_plots", True) and cfg.get("make_gain_loss_ratio_heatmap", True):
        save_gain_loss_ratio_heatmap(stats, freqs, counts, output_dir, cfg)

    write_summary(
        output_dir / "analysis_summary.txt",
        cfg,
        files,
        sample_order,
        before_samples,
        after_samples,
        sample_summary,
        stats,
    )

    print(f"Done. Results written to: {output_dir}")
    print(f"Samples analyzed in order: {len(sample_order)}")
    print(f"Retained barcodes: {len(freqs):,}")
    print(f"Barcodes with presence/absence differences: {len(variable):,}")
    print("Key outputs:")
    print(f"  {output_dir / 'all_barcode_ordered_sample_stats.csv'}")
    print(f"  {output_dir / 'sample_barcode_counts_summary.csv'}")
    print(f"  {output_dir / 'barcodes_with_presence_absence_differences.csv'}")
    print(f"  {output_dir / 'analysis_summary.txt'}")


if __name__ == "__main__":
    main()