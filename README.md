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
* Uses file-level threading when multiple primary input files are available.

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
* Uses file-level threading when multiple bun-match CSVs are available.

### Step 3 — Ordered Barcode Sample Analysis

`ordered_barcode_sample_analysis.py`:

* Reads barcode count CSVs across samples, treatments, or timepoints in either alphabetical order or an explicitly configured order.
* Calls barcode presence per sample using a minimum count threshold.
* Retains barcodes that pass the configured sample/count filters, fills missing sample/barcode combinations with zero reads, and builds count, frequency, and presence matrices.
* Reports first-to-last changes, largest increases/decreases, highest before/after abundance, and presence/absence differences.
* Optionally writes trajectory plots, frequency heatmaps, first-vs-last scatter plots, and a gain/loss ratio heatmap.
* Uses threading to read and preprocess per-sample barcode count CSVs.

---

## Dependencies

Core FASTA/FASTQ parsing uses Biopython. Ordered sample analysis also uses pandas, NumPy, and matplotlib.

```bash
conda install -c bioconda -c conda-forge biopython pandas numpy matplotlib
```

---

## Configuration

InsertDiversity now uses **one consolidated config file**:

* `pipeline_config.json` controls every script.
* `run_order` controls the order used by `run_pipeline.py`.
* Each entry under `scripts` has an `enabled` flag. Set that flag to `true` to run the step or `false` to skip it.
* Each entry under `scripts` has a `config` object containing that step's settings.
* Top-level `threads` controls threading for all steps. The default value is `"auto"`, which uses every hardware thread reported by the system. Set it to a positive integer to limit thread usage globally.
* If needed, a specific step can override the global thread count by adding `"threads": <number>` inside that step's `config` object.

Example toggle and thread limit:

```json
{
  "threads": 8,
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

`run_pipeline.py` passes the same consolidated config file to each enabled script. Each script reads its own block from `pipeline_config.json`, so no temporary or per-script config files are required.

---

## Running Scripts Individually

You can still run any script directly, but use the same consolidated config file:

### 1. Detect buns/flanking matches

```bash
python bun_extract.py --config pipeline_config.json
```

### 2. Extract and summarize barcodes

```bash
python barcode_from_bun_csv.py --config pipeline_config.json
```

This writes per-file augmented CSVs, per-file barcode count reports, and a multi-file barcode membership CSV under the locations configured in `pipeline_config.json`.

It also writes a per-file summary CSV containing file-level metrics such as:

* total reads in each input CSV
* reads with extracted barcodes
* total barcode observations
* unique and real-unique barcode counts
* counts of barcodes unique to that file vs shared across files
* min/max barcode abundance and mean count per unique barcode
* percent of reads with barcodes and percent unique/shared barcodes

### 3. Compare barcodes across ordered samples

```bash
python ordered_barcode_sample_analysis.py pipeline_config.json
```

Important ordered-analysis settings include:

* `min_count_to_call_barcode_present`: minimum count needed to call a barcode present in a sample.
* `min_samples_present_to_retain_barcode`: number of samples in which a barcode must pass the per-sample filter to be retained.
* `require_expected_length`: whether to require barcode length to match `expected_barcode_length`.
* `frequency_denominator`: whether frequencies use all barcodes or only retained barcodes.
* `before_sample_indices` / `after_sample_indices`: optional explicit sample groups for before-vs-after comparisons.

---

## Threading Notes

The scripts parallelize across independent input files, which is the safest concurrency boundary for this workflow:

* `bun_extract.py` processes primary read files concurrently.
* `barcode_from_bun_csv.py` processes bun-match CSV files concurrently.
* `ordered_barcode_sample_analysis.py` reads and preprocesses barcode count CSV files concurrently before building ordered matrices.

By default, `threads` is `"auto"`, so the scripts use the maximum hardware threads available to the process. To limit CPU usage, set top-level `threads` in `pipeline_config.json` to a positive integer.
