## InsertDiversity

**A lightweight FASTQ analysis toolkit for quantifying library diversity between known flanking sequences.**

InsertDiversity scans Illumina FASTQ files to identify reads containing user-defined upstream and downstream flanking sequences (“buns”), extracts the sequence between them (the library region of interest or barcode), and reports the diversity of the insert.

 
---

## How It Works

### Step 1 — Flank Detection

`bun_extract.py`:

* Searches reads (forward and reverse-complement)
* Finds exact partial matches ≥ configurable length
* Determines orientation
* Outputs a CSV with flanking match positions and trimmed regions

### Step 2 — Insert / Barcode Extraction

`barcode_from_bun_csv.py`:

* Extracts the sequence between upstream and downstream matches
* Adds barcode + barcode length to CSV
* Produces a summary report with:

  * Unique barcode count across all extracted barcodes
  * Real unique barcode count using a configurable minimum count threshold (default 3)
  * Frequency per barcode
  * Percent abundance
  * A per-barcode flag showing whether it meets the real-barcode threshold
  * Length deviation from expected

---

## Dependencies

```bash
conda install -c bioconda -c conda-forge biopython
```

---

## How to run it:

Just configure the config files, then run the below in the terminal (only tested on Linux):

```bash
python bun_extract.py --config config.json
```

In the default `config.json`, `bun_extract.py` scans a folder and writes one `{name}_bun_matches.csv` per primary input file into `outputs/bun_matches/`, plus a cross-file barcode membership CSV.

Then run:

```bash
python barcode_from_bun_csv.py --config barcode_config.json
```

In the default `barcode_config.json`, this processes every `*_bun_matches.csv` in that folder and writes per-file augmented CSVs, per-file barcode count reports, and a multi-file barcode membership CSV under `outputs/barcodes/`.
It also writes a per-file summary CSV (`outputs/barcodes/barcode_file_summary.csv`) containing file-level metrics such as:

* total reads in each input CSV
* reads with extracted barcodes
* total barcode observations
* unique and real-unique barcode counts
* counts of barcodes unique to that file vs shared across files
* min/max barcode abundance and mean count per unique barcode
* percent of reads with barcodes and percent unique/shared barcodes
