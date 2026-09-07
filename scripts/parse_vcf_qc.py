#!/usr/bin/env python3
"""Compute QC metrics from an aligned BAM + called VCF and write them to SQLite + JSON."""
import argparse
import json
import re
import sqlite3
from pathlib import Path

import pysam


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--bam", required=True)
    p.add_argument("--vcf", required=True)
    p.add_argument("--dup-metrics", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--db", required=True)
    p.add_argument("--json", required=True)
    p.add_argument("--min-mean-depth", type=float, required=True)
    p.add_argument("--min-q30-pct", type=float, required=True)
    p.add_argument("--titv-min", type=float, required=True)
    p.add_argument("--titv-max", type=float, required=True)
    p.add_argument("--min-variant-pass-pct", type=float, required=True)
    p.add_argument("--max-dup-rate-pct", type=float, required=True)
    return p.parse_args()


TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def compute_depth_and_q30(bam_path):
    """Per-base depth (via pileup) and fraction of bases with Q>=30, restricted to
    the region(s) covered by the BAM (which is already sliced to our target region)."""
    bam = pysam.AlignmentFile(bam_path, "rb")

    depths = []
    for pileup_col in bam.pileup(stepper="samtools", min_base_quality=0):
        depths.append(pileup_col.n)

    total_bases = 0
    q30_bases = 0
    for read in bam.fetch():
        if read.is_unmapped or read.query_qualities is None:
            continue
        quals = read.query_qualities
        total_bases += len(quals)
        q30_bases += sum(1 for q in quals if q >= 30)

    bam.close()

    mean_depth = sum(depths) / len(depths) if depths else 0.0
    min_depth = min(depths) if depths else 0
    max_depth = max(depths) if depths else 0
    q30_pct = (100.0 * q30_bases / total_bases) if total_bases else 0.0

    return {
        "mean_depth": mean_depth,
        "min_depth": min_depth,
        "max_depth": max_depth,
        "bases_covered": len(depths),
        "pct_bases_q30": q30_pct,
    }


def compute_variant_metrics(vcf_path):
    vcf = pysam.VariantFile(vcf_path)

    n_total = 0
    n_pass = 0
    n_snps = 0
    n_indels = 0
    n_ts = 0
    n_tv = 0

    for rec in vcf:
        n_total += 1
        is_pass = (rec.filter.keys() == ["PASS"]) or (len(rec.filter.keys()) == 0)
        if is_pass:
            n_pass += 1

        ref = rec.ref
        for alt in rec.alts or []:
            if len(ref) == 1 and len(alt) == 1:
                n_snps += 1
                if is_pass:
                    if (ref, alt) in TRANSITIONS:
                        n_ts += 1
                    else:
                        n_tv += 1
            else:
                n_indels += 1

    vcf.close()

    titv = (n_ts / n_tv) if n_tv else 0.0
    pct_pass = (100.0 * n_pass / n_total) if n_total else 0.0

    return {
        "n_variants_total": n_total,
        "n_variants_pass": n_pass,
        "pct_variants_pass": pct_pass,
        "n_snps": n_snps,
        "n_indels": n_indels,
        "titv_ratio": titv,
    }


def parse_dup_metrics(path):
    """samtools markdup -s writes a small text summary to stderr-captured file, e.g.:
    READ COUNT: 244604, WRITTEN: 244604, EXCLUDED: 0, EXAMINED: 244604,
    PAIRED: 244000, SINGLE: 604, DUPLICATE PAIR: 12000, DUPLICATE SINGLE: 10,
    DUPLICATE TOTAL: 12010
    """
    text = Path(path).read_text()
    examined = _extract_int(text, r"EXAMINED:\s*(\d+)")
    dup_total = _extract_int(text, r"DUPLICATE\s+TOTAL:\s*(\d+)")
    dup_rate_pct = (100.0 * dup_total / examined) if examined else 0.0
    return {
        "reads_examined": examined,
        "duplicate_reads": dup_total,
        "pct_duplicates": dup_rate_pct,
    }


def _extract_int(text, pattern):
    m = re.search(pattern, text)
    return int(m.group(1)) if m else 0


def determine_result(metrics, thresholds):
    flags = []
    if metrics["mean_depth"] < thresholds["min_mean_depth"]:
        flags.append(
            f"Mean depth {metrics['mean_depth']:.1f}x is below the "
            f"{thresholds['min_mean_depth']}x threshold."
        )
    if metrics["pct_bases_q30"] < thresholds["min_q30_pct"]:
        flags.append(
            f"% bases Q>=30 ({metrics['pct_bases_q30']:.1f}%) is below "
            f"{thresholds['min_q30_pct']}% threshold."
        )
    if not (thresholds["titv_min"] <= metrics["titv_ratio"] <= thresholds["titv_max"]):
        flags.append(
            f"Ti/Tv ratio {metrics['titv_ratio']:.2f} is outside the expected "
            f"{thresholds['titv_min']}-{thresholds['titv_max']} range for this assay type, "
            "which can indicate elevated false-positive SNP calls."
        )
    if metrics["pct_variants_pass"] < thresholds["min_variant_pass_pct"]:
        flags.append(
            f"Only {metrics['pct_variants_pass']:.1f}% of called variants pass filters "
            f"(threshold {thresholds['min_variant_pass_pct']}%)."
        )
    if metrics["pct_duplicates"] > thresholds["max_dup_rate_pct"]:
        flags.append(
            f"Duplicate read rate {metrics['pct_duplicates']:.1f}% exceeds "
            f"{thresholds['max_dup_rate_pct']}% threshold."
        )

    if not flags:
        overall = "PASS"
    elif len(flags) <= 1 and metrics["mean_depth"] >= thresholds["min_mean_depth"] * 0.5:
        overall = "PASS WITH FLAGS"
    else:
        overall = "FAIL"

    return overall, flags


def main():
    args = parse_args()
    thresholds = {
        "min_mean_depth": args.min_mean_depth,
        "min_q30_pct": args.min_q30_pct,
        "titv_min": args.titv_min,
        "titv_max": args.titv_max,
        "min_variant_pass_pct": args.min_variant_pass_pct,
        "max_dup_rate_pct": args.max_dup_rate_pct,
    }

    depth_metrics = compute_depth_and_q30(args.bam)
    variant_metrics = compute_variant_metrics(args.vcf)
    dup_metrics = parse_dup_metrics(args.dup_metrics)

    metrics = {**depth_metrics, **variant_metrics, **dup_metrics}
    overall, flags = determine_result(metrics, thresholds)

    result = {
        "sample_id": args.sample_id,
        **metrics,
        "overall_determination": overall,
        "flags": flags,
        "thresholds": thresholds,
    }

    with open(args.json, "w") as f:
        json.dump(result, f, indent=2)

    con = sqlite3.connect(args.db)
    cur = con.cursor()
    cur.execute("""
        CREATE TABLE IF NOT EXISTS qc_metrics (
            sample_id TEXT PRIMARY KEY,
            mean_depth REAL,
            min_depth INTEGER,
            max_depth INTEGER,
            bases_covered INTEGER,
            pct_bases_q30 REAL,
            n_variants_total INTEGER,
            n_variants_pass INTEGER,
            pct_variants_pass REAL,
            n_snps INTEGER,
            n_indels INTEGER,
            titv_ratio REAL,
            reads_examined INTEGER,
            duplicate_reads INTEGER,
            pct_duplicates REAL,
            overall_determination TEXT,
            flags TEXT
        )
    """)
    cur.execute("""
        INSERT OR REPLACE INTO qc_metrics VALUES
        (:sample_id, :mean_depth, :min_depth, :max_depth, :bases_covered, :pct_bases_q30,
         :n_variants_total, :n_variants_pass, :pct_variants_pass, :n_snps, :n_indels,
         :titv_ratio, :reads_examined, :duplicate_reads, :pct_duplicates,
         :overall_determination, :flags)
    """, {
        **metrics,
        "sample_id": args.sample_id,
        "overall_determination": overall,
        "flags": json.dumps(flags),
    })
    con.commit()
    con.close()

    print(f"[{args.sample_id}] overall determination: {overall}")
    for flag in flags:
        print(f"  FLAG: {flag}")


if __name__ == "__main__":
    main()
