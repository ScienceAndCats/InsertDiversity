#!/usr/bin/env python3
"""
barcode_from_bun_csv.py

Input: CSV produced by bun_partial_extract_to_csv.py (the script that outputs oriented_read_seq,
       upstream_start/upstream_end, downstream_start/downstream_end, match_type, etc.)

Outputs:
  1) An "augmented" CSV with two new columns:
        - barcode     (sequence between upstream_end and downstream_start in oriented_read_seq)
        - barcode_len (length of barcode, blank if not available)

  2) A barcode report CSV with:
        - barcode
        - barcode_len
        - count
        - percent   (percent of reads that had a valid barcode)

Barcode definition:
  If BOTH upstream and downstream hits exist AND are ordered (upstream_end < downstream_start),
  then:
      barcode = oriented_read_seq[upstream_end + 1 : downstream_start]
  (upstream_end is inclusive; downstream_start is the first base of the downstream match)

Config (JSON) required fields:
  - in_csv
  - out_csv
  - out_barcode_report

Optional config fields:
  - expected_barcode_len (default 25)  # used only for an extra QC column in report
  - min_barcode_len (default 1)        # only barcodes >= this length are counted in report
  - include_empty_barcodes (default false)  # if true, count empty barcode strings too (not recommended)

Dependencies:
  - standard library only (csv, json)
"""

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Optional


def safe_int(x: str) -> Optional[int]:
    if x is None:
        return None
    s = str(x).strip()
    if s == "":
        return None
    try:
        return int(s)
    except ValueError:
        return None


def main():
    ap = argparse.ArgumentParser(description="Extract barcode (between buns) from bun-match CSV and summarize counts.")
    ap.add_argument("--config", required=True, help="Path to JSON config file")
    args = ap.parse_args()

    cfg_path = Path(args.config)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config not found: {cfg_path}")

    cfg = json.loads(cfg_path.read_text())

    in_csv = cfg.get("in_csv")
    out_csv = cfg.get("out_csv")
    out_barcode_report = cfg.get("out_barcode_report")

    expected_len = int(cfg.get("expected_barcode_len", 25))
    min_barcode_len = int(cfg.get("min_barcode_len", 1))
    include_empty = bool(cfg.get("include_empty_barcodes", False))

    if not in_csv or not out_csv or not out_barcode_report:
        raise ValueError("Config must include: in_csv, out_csv, out_barcode_report")

    in_csv = Path(in_csv)
    out_csv = Path(out_csv)
    out_barcode_report = Path(out_barcode_report)

    if not in_csv.exists():
        raise FileNotFoundError(f"Input CSV not found: {in_csv}")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_barcode_report.parent.mkdir(parents=True, exist_ok=True)

    barcode_counts = Counter()
    total_rows = 0
    rows_with_valid_barcode = 0

    with open(in_csv, "r", newline="") as f_in:
        reader = csv.DictReader(f_in)
        if reader.fieldnames is None:
            raise ValueError("Input CSV appears to have no header row.")

        # Add new columns (avoid duplicates if re-running)
        fieldnames = list(reader.fieldnames)
        for new_col in ["barcode", "barcode_len"]:
            if new_col not in fieldnames:
                fieldnames.append(new_col)

        with open(out_csv, "w", newline="") as f_out:
            writer = csv.DictWriter(f_out, fieldnames=fieldnames)
            writer.writeheader()

            for row in reader:
                total_rows += 1

                oriented = (row.get("oriented_read_seq") or "").strip().upper()
                up_end = safe_int(row.get("upstream_end", ""))
                dn_start = safe_int(row.get("downstream_start", ""))

                barcode = ""
                barcode_len = ""

                # Only definable if both coordinates exist and are ordered
                if oriented and up_end is not None and dn_start is not None and up_end < dn_start:
                    # upstream_end is inclusive; downstream_start is start of downstream match
                    barcode = oriented[up_end + 1 : dn_start]
                    barcode_len = str(len(barcode))

                    # Decide whether to count it in the summary
                    if include_empty:
                        ok_to_count = (len(barcode) >= 0)
                    else:
                        ok_to_count = (len(barcode) >= min_barcode_len)

                    if ok_to_count:
                        barcode_counts[barcode] += 1
                        rows_with_valid_barcode += 1

                row["barcode"] = barcode
                row["barcode_len"] = barcode_len

                writer.writerow(row)

    # Write barcode report
    # Percent is relative to rows_with_valid_barcode (i.e., reads where we successfully extracted a barcode)
    with open(out_barcode_report, "w", newline="") as f_rep:
        rep_fields = ["barcode", "barcode_len", "count", "percent", "expected_barcode_len", "len_minus_expected"]
        rep = csv.DictWriter(f_rep, fieldnames=rep_fields)
        rep.writeheader()

        denom = rows_with_valid_barcode if rows_with_valid_barcode > 0 else 1

        for barcode, count in barcode_counts.most_common():
            blen = len(barcode)
            rep.writerow({
                "barcode": barcode,
                "barcode_len": blen,
                "count": count,
                "percent": round(100.0 * count / denom, 6),
                "expected_barcode_len": expected_len,
                "len_minus_expected": blen - expected_len,
            })

    unique_barcodes = len(barcode_counts)

    print("Done.")
    print(f"Input rows (reads in CSV):                {total_rows:,}")
    print(f"Rows with extracted barcodes (counted):   {rows_with_valid_barcode:,}")
    print(f"Unique barcodes (in counted set):         {unique_barcodes:,}")
    print(f"Augmented CSV written to:                 {out_csv}")
    print(f"Barcode report written to:                {out_barcode_report}")


if __name__ == "__main__":
    main()
