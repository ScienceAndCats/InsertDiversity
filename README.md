## InsertDiversity

**A lightweight FASTQ analysis toolkit for quantifying library diversity between known flanking sequences.**

InsertDiversity scans Illumina FASTQ files to identify reads containing user-defined upstream and downstream flanking sequences (“buns”), extracts the sequence between them (the library region of interest or barcode), summarizes barcode diversity, and can compare barcode abundance across an ordered set of samples.

---

## How It Works

### Step 1 — Flank Detection

`bun_extract.py`:

* Searches reads in forward and reverse-complement orientations.
* Finds exact partial matches at least `min_match` bases long for the configured upstream and downstream buns.
* Determines the best read orientation.
* Writes one `{name}_bun_matches.csv` file per primary FASTA/FASTQ input file.
* Skips paired-end R2 files as primary inputs, but uses an R2 partner when present to verify the barcode extracted from R1.
* Writes a cross-file barcode membership CSV.

### Step 2 — Insert / Barcode Extraction

`barcode_from_bun_csv.py`:

* Reads one bun-match CSV or a folder of `*_bun_matches.csv` files.
* Extracts the sequence between upstream and downstream matches.
* Adds `barcode` and `barcode_len` columns to augmented CSV outputs.
* Produces per-file barcode count reports with:

  * Unique barcode count across all extracted barcodes.
  * Real unique barcode count using a configurable minimum count threshold (default 3).
  * Frequency per barcode.
  * Percent abundance.
  * A per-barcode flag showing whether it meets the real-barcode threshold.
  * Length deviation from the expected barcode length.

### Step 3 — Ordered Barcode Sample Analysis

`ordered_barcode_sample_analysis.py`:

* Reads barcode count CSVs across samples, treatments, or timepoints in either alphabetical order or an explicitly configured order.
* Calls barcode presence per sample using a minimum count threshold.
* Retains barcodes that pass the configured sample/count filters, fills missing sample/barcode combinations with zero reads, and builds count, frequency, and presence matrices.
* Reports first-to-last changes, largest increases/decreases, highest before/after abundance, and presence/absence differences.
* Optionally writes trajectory plots, frequency heatmaps, first-vs-last scatter plots, and a gain/loss ratio heatmap.

---

## Dependencies

Core FASTA/FASTQ parsing uses Biopython. Ordered sample analysis also uses pandas, NumPy, and matplotlib.

```bash
conda install -c bioconda -c conda-forge biopython pandas numpy matplotlib
```

---

## Configuration Files

The repository keeps individual config files for running each script directly:

* `config.json` configures `bun_extract.py`.
* `barcode_config.json` configures `barcode_from_bun_csv.py`.
* `config_ordered_barcode_analysis.json` configures `ordered_barcode_sample_analysis.py`.

It also includes a consolidated pipeline config:

* `pipeline_config.json` contains the config blocks for all scripts in one file.
* Each entry under `scripts` has an `enabled` flag. Set that flag to `true` to run the step or `false` to skip it when using `run_pipeline.py`.
* `run_order` controls the order used by `run_pipeline.py`.

Example toggle:

```json
{
  "scripts": {
    "bun_extract": { "enabled": true },
    "barcode_from_bun_csv": { "enabled": true },
    "ordered_barcode_sample_analysis": { "enabled": false }
  }
}
```

The example above is abbreviated; keep the full `config` objects in `pipeline_config.json`.

---

## Running the Full or Partial Pipeline

Edit `pipeline_config.json`, then run:

```bash
python run_pipeline.py --config pipeline_config.json
```

To check which enabled commands would run without executing them:

```bash
python run_pipeline.py --config pipeline_config.json --dry-run
```

`run_pipeline.py` writes temporary per-step JSON files from the consolidated config and passes them to each script. This preserves compatibility with the original scripts while letting you manage all settings and run toggles in one place.

---

## Running Scripts Individually

### 1. Detect buns/flanking matches

Configure `config.json`, then run:

```bash
python bun_extract.py --config config.json
```

In the default `config.json`, `bun_extract.py` scans `data/s17-bc` and writes one `{name}_bun_matches.csv` per primary input file into `outputs/bun_matches/`, plus a cross-file barcode membership CSV.

### 2. Extract and summarize barcodes

Configure `barcode_config.json`, then run:

```bash
python barcode_from_bun_csv.py --config barcode_config.json
```

In the default `barcode_config.json`, this processes every `*_bun_matches.csv` in `outputs/bun_matches/` and writes per-file augmented CSVs, per-file barcode count reports, and a multi-file barcode membership CSV under `outputs/barcodes/`.

It also writes a per-file summary CSV (`outputs/barcodes/barcode_file_summary.csv`) containing file-level metrics such as:

* total reads in each input CSV
* reads with extracted barcodes
* total barcode observations
* unique and real-unique barcode counts
* counts of barcodes unique to that file vs shared across files
* min/max barcode abundance and mean count per unique barcode
* percent of reads with barcodes and percent unique/shared barcodes

### 3. Compare barcodes across ordered samples

Configure `config_ordered_barcode_analysis.json`, then run:

```bash
python ordered_barcode_sample_analysis.py config_ordered_barcode_analysis.json
```

The default ordered-analysis config reads barcode count CSVs from `outputs/barcodes/per_file_reports/`. If you need a specific sample/timepoint order, set `explicit_files` to the exact file list in the desired order. Otherwise files are analyzed alphabetically using `csv_pattern`.

Important ordered-analysis settings include:

* `min_count_to_call_barcode_present`: minimum count needed to call a barcode present in a sample.
* `min_samples_present_to_retain_barcode`: number of samples in which a barcode must be present to keep it for analysis.
* `frequency_denominator`: choose `all_barcodes` or `retained_barcodes` for frequency normalization.
* `before_sample_indices` / `after_sample_indices`: optional zero-based sample groups used for before/after summaries.
* `make_plots` and plot-specific flags: toggle trajectory plots, heatmaps, and scatter plots.

Primary outputs are written under `ordered_barcode_analysis_outputs/` and include retained matrices, all-barcode stats, presence/absence differences, ranked barcode tables, plots, and `analysis_summary.txt`.
