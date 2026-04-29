#!/usr/bin/env python3
"""
barcode_from_bun_csv.py

Extract barcode sequences from bun-match CSV(s) and produce reports.

Supports:
  - single input CSV (legacy mode)
  - folder mode: all *_bun_matches.csv files in a directory

In folder mode, also reports:
  - barcodes found in multiple files
  - barcodes unique to only one file
  - per-file unique barcode counts
"""

import argparse
import concurrent.futures
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Optional


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


def strip_suffix(name: str, suffix: str) -> str:
    if name.endswith(suffix):
        return name[: -len(suffix)]
    return Path(name).stem


def discover_input_csvs(cfg: dict) -> Iterable[Path]:
    in_csv = cfg.get("in_csv")
    in_csv_dir = cfg.get("in_csv_dir")
    in_csv_glob = cfg.get("in_csv_glob", "*_bun_matches.csv")

    if in_csv:
        p = Path(in_csv)
        if not p.exists():
            raise FileNotFoundError(f"Input CSV not found: {p}")
        return [p]

    if in_csv_dir:
        d = Path(in_csv_dir)
        if not d.is_dir():
            raise FileNotFoundError(f"in_csv_dir not found or not a directory: {d}")
        files = sorted(p for p in d.glob(in_csv_glob) if p.is_file())
        if not files:
            raise ValueError(f"No CSV files matched {in_csv_glob} in {d}")
        return files

    raise ValueError("Config must include one of: in_csv or in_csv_dir")


def process_one_csv(
    in_csv: Path,
    out_csv: Path,
    out_barcode_report: Path,
    expected_len: int,
    min_barcode_len: int,
    include_empty: bool,
    min_real_barcode_count: int,
):
    barcode_counts = Counter()
    total_rows = 0
    rows_with_extracted_barcode = 0

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_barcode_report.parent.mkdir(parents=True, exist_ok=True)

    with open(in_csv, "r", newline="") as f_in:
        reader = csv.DictReader(f_in)
        if reader.fieldnames is None:
            raise ValueError(f"Input CSV appears to have no header row: {in_csv}")

        fieldnames = list(reader.fieldnames)
        for new_col in ["barcode", "barcode_len"]:
            if new_col not in fieldnames:
                fieldnames.append(new_col)

        with open(out_csv, "w", newline="") as f_out:
            writer = csv.DictWriter(f_out, fieldnames=fieldnames)
            writer.writeheader()

            for row in reader:
                total_rows += 1

                barcode = (row.get("barcode") or "").strip().upper()
                barcode_len = row.get("barcode_len", "")

                # Backward compatibility for older bun CSVs that do not include barcode columns.
                if not barcode:
                    oriented = (row.get("oriented_read_seq") or "").strip().upper()
                    up_end = safe_int(row.get("upstream_end", ""))
                    dn_start = safe_int(row.get("downstream_start", ""))

                    if oriented and up_end is not None and dn_start is not None and up_end < dn_start:
                        barcode = oriented[up_end + 1 : dn_start]
                        barcode_len = str(len(barcode))

                if include_empty:
                    ok_to_count = barcode_len != ""
                else:
                    ok_to_count = len(barcode) >= min_barcode_len

                if ok_to_count and barcode_len != "":
                    barcode_counts[barcode] += 1
                    rows_with_extracted_barcode += 1

                row["barcode"] = barcode
                row["barcode_len"] = barcode_len
                writer.writerow(row)

    with open(out_barcode_report, "w", newline="") as f_rep:
        rep_fields = [
            "barcode",
            "barcode_len",
            "count",
            "percent",
            "is_real_barcode",
            "min_real_barcode_count",
            "expected_barcode_len",
            "len_minus_expected",
        ]
        rep = csv.DictWriter(f_rep, fieldnames=rep_fields)
        rep.writeheader()

        denom = rows_with_extracted_barcode if rows_with_extracted_barcode > 0 else 1

        for barcode, count in barcode_counts.most_common():
            blen = len(barcode)
            rep.writerow(
                {
                    "barcode": barcode,
                    "barcode_len": blen,
                    "count": count,
                    "percent": round(100.0 * count / denom, 6),
                    "is_real_barcode": "true" if count >= min_real_barcode_count else "false",
                    "min_real_barcode_count": min_real_barcode_count,
                    "expected_barcode_len": expected_len,
                    "len_minus_expected": blen - expected_len,
                }
            )

    return {
        "total_rows": total_rows,
        "rows_with_extracted_barcode": rows_with_extracted_barcode,
        "barcode_counts": barcode_counts,
        "unique_barcodes": len(barcode_counts),
        "real_unique_barcodes": sum(1 for c in barcode_counts.values() if c >= min_real_barcode_count),
    }


def main():
    ap = argparse.ArgumentParser(description="Extract barcode(s) from bun-match CSV(s) and summarize counts.")
    ap.add_argument("--config", required=True, help="Path to JSON config file")
    args = ap.parse_args()

    cfg_path = Path(args.config)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config not found: {cfg_path}")
    cfg = json.loads(cfg_path.read_text())

    expected_len = int(cfg.get("expected_barcode_len", 25))
    min_barcode_len = int(cfg.get("min_barcode_len", 1))
    include_empty = bool(cfg.get("include_empty_barcodes", False))
    min_real_barcode_count = int(cfg.get("min_real_barcode_count", 3))

    input_csvs = [Path(p) for p in discover_input_csvs(cfg)]

    out_csv = cfg.get("out_csv")
    out_barcode_report = cfg.get("out_barcode_report")
    out_csv_dir = Path(cfg.get("out_csv_dir", "."))
    out_report_dir = Path(cfg.get("out_barcode_report_dir", "."))
    out_multi = Path(cfg.get("out_multi_file_barcode_report", "barcode_multi_file_membership.csv"))

    if len(input_csvs) == 1:
        if not out_csv or not out_barcode_report:
            raise ValueError("Single-file mode requires: out_csv and out_barcode_report")
        output_plan = [
            {
                "in": input_csvs[0],
                "out_csv": Path(out_csv),
                "out_report": Path(out_barcode_report),
                "label": strip_suffix(input_csvs[0].name, "_bun_matches.csv"),
            }
        ]
    else:
        output_plan = []
        for in_file in input_csvs:
            label = strip_suffix(in_file.name, "_bun_matches.csv")
            output_plan.append(
                {
                    "in": in_file,
                    "out_csv": out_csv_dir / f"{label}_bun_matches_with_barcode.csv",
                    "out_report": out_report_dir / f"{label}_barcode_counts.csv",
                    "label": label,
                }
            )

    barcode_presence = defaultdict(set)
    per_file_unique_counts: Dict[str, int] = {}

    max_workers = int(cfg.get("max_workers", len(output_plan)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_item = {
            executor.submit(
                process_one_csv,
                in_csv=item["in"],
                out_csv=item["out_csv"],
                out_barcode_report=item["out_report"],
                expected_len=expected_len,
                min_barcode_len=min_barcode_len,
                include_empty=include_empty,
                min_real_barcode_count=min_real_barcode_count,
            ): item
            for item in output_plan
        }

        for future in concurrent.futures.as_completed(future_to_item):
            item = future_to_item[future]
            stats = future.result()
            print(f"Done: {item['in']}")
            print(f"  Input rows (reads in CSV):                {stats['total_rows']:,}")
            print(f"  Rows with extracted barcodes:             {stats['rows_with_extracted_barcode']:,}")
            print(f"  Unique barcodes (all extracted):          {stats['unique_barcodes']:,}")
            print(f"  Real unique barcodes (count>={min_real_barcode_count}): {stats['real_unique_barcodes']:,}")
            print(f"  Augmented CSV written to:                 {item['out_csv']}")
            print(f"  Barcode report written to:                {item['out_report']}")

            for bc in stats["barcode_counts"].keys():
                barcode_presence[bc].add(item["label"])

    if len(output_plan) > 1:
        out_multi.parent.mkdir(parents=True, exist_ok=True)
        with open(out_multi, "w", newline="") as f_out:
            writer = csv.DictWriter(
                f_out,
                fieldnames=["barcode", "num_files", "files", "classification"],
            )
            writer.writeheader()
            for barcode, files in sorted(barcode_presence.items(), key=lambda kv: (-len(kv[1]), kv[0])):
                file_list = sorted(files)
                writer.writerow(
                    {
                        "barcode": barcode,
                        "num_files": len(file_list),
                        "files": ";".join(file_list),
                        "classification": "multiple_files" if len(file_list) > 1 else f"unique_to_{file_list[0]}",
                    }
                )

        shared = sum(1 for files in barcode_presence.values() if len(files) > 1)
        unique_total = sum(1 for files in barcode_presence.values() if len(files) == 1)

        for label in [item["label"] for item in output_plan]:
            per_file_unique_counts[label] = sum(1 for files in barcode_presence.values() if files == {label})

        print(f"Multi-file barcode membership report:       {out_multi}")
        print(f"Barcodes found in multiple files:           {shared:,}")
        print(f"Barcodes unique to single files:            {unique_total:,}")
        for label in sorted(per_file_unique_counts):
            print(f"  Unique to {label}: {per_file_unique_counts[label]:,}")


if __name__ == "__main__":
    main()
