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

This outputs bun_matches.csv. 

Then run:

```bash
python barcode_from_bun_csv.py --config barcode_config.json
```

This outputs bun_matches_with_barcode.csv and barcode_counts.csv.
