#!/usr/bin/env python3
"""Generate QC summary plots from the qc.db SQLite database."""
import argparse
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db", required=True)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def fetch_row(db_path, sample_id):
    con = sqlite3.connect(db_path)
    con.row_factory = sqlite3.Row
    cur = con.cursor()
    cur.execute("SELECT * FROM qc_metrics WHERE sample_id = ?", (sample_id,))
    row = cur.fetchone()
    con.close()
    if row is None:
        raise SystemExit(f"No qc_metrics row found for sample_id={sample_id}")
    return dict(row)


def plot_variant_composition(row, outdir):
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.bar(["SNPs", "Indels"], [row["n_snps"], row["n_indels"]], color=["#4C72B0", "#DD8452"])
    ax.set_ylabel("Variant count")
    ax.set_title(f"{row['sample_id']}: variant composition")
    fig.tight_layout()
    fig.savefig(Path(outdir) / f"{row['sample_id']}.variant_composition.png", dpi=150)
    plt.close(fig)


def plot_qc_summary(row, outdir):
    metrics = ["mean_depth", "pct_bases_q30", "pct_variants_pass", "pct_duplicates", "titv_ratio"]
    labels = ["Mean depth (x)", "% bases Q30", "% variants PASS", "% duplicates", "Ti/Tv ratio"]
    values = [row[m] for m in metrics]

    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(labels, values, color="#55A868")
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{v:.1f}",
                ha="center", va="bottom", fontsize=9)
    ax.set_title(f"{row['sample_id']}: QC metrics summary ({row['overall_determination']})")
    plt.xticks(rotation=20, ha="right")
    fig.tight_layout()
    fig.savefig(Path(outdir) / f"{row['sample_id']}.qc_summary.png", dpi=150)
    plt.close(fig)


def main():
    args = parse_args()
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    row = fetch_row(args.db, args.sample_id)
    plot_variant_composition(row, args.outdir)
    plot_qc_summary(row, args.outdir)
    print(f"Wrote plots for {args.sample_id} to {args.outdir}")


if __name__ == "__main__":
    main()
