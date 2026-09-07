# Validation Report: NGS Variant-Calling QC Pipeline

**Run date:** 2026-09-07
**Pipeline version:** initial run (Nextflow 26.04.6, no tagged release yet)
**Analyst:** Vedika Harnathka

## 1. Purpose

This report documents the quality control checks performed on a variant-calling pipeline run and states whether the output meets the pre-defined QC thresholds for downstream use. It follows the structure of a lightweight production validation report: what was run, what was measured, what the pass/fail criteria were, and what the outcome was.

## 2. Input

| Field | Value |
|---|---|
| Sample ID | NA12878_chr7 |
| Data source | Genome in a Bottle (GIAB) reference material HG001 / NA12878; 121,325 real Illumina HiSeq 2500 read pairs extracted region-restricted from NIST's public `RMNISTHS_30xdownsample.bam` |
| Reference genome | GRCh37/hg19, chr7:117,000,000–118,000,000 |
| Region of interest | CFTR locus (chr7:117,120,016–117,308,718 in hg19 falls within this 1 Mb window) |

## 3. Pipeline steps executed

| Step | Tool | Purpose |
|---|---|---|
| Raw read QC | FastQC 0.12.1 | Flag adapter contamination, low-quality bases |
| Alignment | BWA-MEM 0.7.19 | Align reads to reference |
| Post-processing | Samtools 1.24 | Sort, mark duplicates, index |
| Variant calling | bcftools 1.24 (mpileup / call) | Call SNPs/indels |
| QC metrics | Python (pysam) | Compute depth, Q30%, Ti/Tv, duplicate rate, filter pass rate |

## 4. QC metrics and thresholds

| Metric | Threshold | Observed | Pass/Fail |
|---|---|---|---|
| Mean read depth | ≥ 20x | 35.5x | [x] Pass |
| % bases with Q ≥ 30 | ≥ 90% | 50.2% | [ ] **Fail** |
| Ti/Tv ratio | 1.9–2.3 | 1.62 | [ ] **Fail** |
| % variants passing filters | ≥ 90% | 98.6% (756/767) | [x] Pass |
| Duplicate read rate | ≤ 20% | 0.13% | [x] Pass |

> **Note on thresholds:** Mean depth ≥20x is a common lower bound for confident germline SNP/indel calling. Q30% ≥90% is a standard modern-Illumina expectation for base-calling accuracy — this run's HiSeq 2500 (2013–2014, v1 rapid-mode chemistry) data falls well short of that, which reflects the sequencer generation rather than a processing error (see §6). Ti/Tv 1.9–2.3 is the expected genome-wide range for human variant calls; a lower observed ratio typically indicates an elevated false-positive SNP rate. Duplicate rate ≤20% is a standard PCR/optical duplicate ceiling for WGS libraries.

## 5. Result

**Overall determination:** FAIL

Depth, variant-filter pass rate, and duplicate rate are all comfortably within spec. However, both the % bases ≥Q30 and the Ti/Tv ratio fall outside their thresholds, and together they point to the same root cause rather than two independent problems: the input reads are older-chemistry HiSeq 2500 data with a lower baseline quality-score ceiling (most bases sit in the Q28–33 range rather than well above Q30), and this pipeline calls variants with `bcftools mpileup/call` directly off raw alignments — no base quality score recalibration (BQSR) or indel realignment. With more failure-prone bases (Q rather than Q30) with an uncorrected quality distribution, false-positive SNP calls (which skew Ti/Tv toward 1) are expected. This run should **not** be used as-is for a downstream diagnostic claim; it should be re-run with a modern-chemistry dataset and/or a variant caller step that includes quality recalibration before being treated as representative.

## 6. Flagged issues (if any)

- **% bases ≥Q30 (50.2%) below 90% threshold.** Root cause: this HiSeq 2500 rapid-mode-v1 dataset (2013–2014 GIAB NA12878 sequencing) predates the quality-score chemistry improvements of later Illumina platforms; this is a property of the input data, not a pipeline defect. Recommend re-running with newer chemistry data (e.g. GIAB's NovaSeq NA12878 datasets) before drawing conclusions from Q30% on its own.
- **Ti/Tv ratio (1.62) below the 1.9–2.3 expected range.** Likely compounded by (a) the Q30 issue above feeding lower-confidence bases into calling, and (b) the absence of BQSR/indel realignment before `bcftools call`, and (c) the small sample size (767 total variants in a 1 Mb window), which makes the ratio more sensitive to a handful of false-positive calls than a whole-genome Ti/Tv estimate would be. Recommend adding a BQSR step and comparing Ti/Tv on a larger region before treating this as a hard fail on the calling logic itself.

## 7. Reproducibility

```bash
nextflow run main.nf -c nextflow.config \
  --reads 'data/*_R{1,2}.fastq.gz' \
  --ref reference/chr7_slice.fa \
  --outdir results/
```

Pipeline code and QC scripts: see repository README. All QC thresholds are defined in `nextflow.config`; tool versions are recorded in this report (§3) and in the BAM/VCF header `@PG` lines under `results/`.

## 8. Limitations of this analysis

- Single sample, restricted to a 1 Mb region — not representative of full-exome or full-panel QC behavior, and the small variant count (n=767) makes ratio-based metrics like Ti/Tv noisier than they would be genome-wide.
- No orthogonal confirmation (e.g. GIAB high-confidence truth-set comparison, or Sanger validation) of called variants — the flagged Ti/Tv issue is diagnosed from read/quality properties, not confirmed against a ground truth.
- Variant calling used `bcftools mpileup/call` with a basic QUAL/DP filter only — no BQSR, no indel realignment, no joint genotyping. A production pipeline (e.g. GATK HaplotypeCaller with BQSR) would very plausibly produce a materially better Ti/Tv ratio on the same input reads.
