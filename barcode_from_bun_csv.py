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
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

from config_helpers import available_threads, configured_threads, load_script_config, run_cli


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


def parse_bool(value, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def discover_input_csvs(cfg: dict) -> Iterable[Path]:
    in_csv = cfg.get("in_csv")
    in_csv_dir = cfg.get("in_csv_dir")
    in_csv_glob = cfg.get("in_csv_glob", "*_bun_matches.csv")

    if in_csv:
        p = Path(in_csv)
        if not p.is_file():
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


def validate_starcode_config(max_distance: int, cluster_ratio: float, threads: int) -> None:
    if max_distance < 0:
        raise ValueError("collapse_max_distance/starcode_distance must be >= 0")
    if cluster_ratio < 1:
        raise ValueError("starcode_cluster_ratio must be >= 1")
    if threads < 1:
        raise ValueError("starcode_threads must be >= 1")


def resolve_starcode_path(configured_path: str) -> str:
    path = Path(configured_path).expanduser()
    if path.parent != Path(".") or path.is_absolute():
        if path.is_file():
            return str(path)
        raise FileNotFoundError(f"starcode executable not found: {path}")

    resolved = shutil.which(configured_path)
    if resolved:
        return resolved
    raise FileNotFoundError(
        f"starcode executable not found on PATH: {configured_path}. "
        "Install gui11aume/starcode or set starcode_path."
    )


def starcode_algorithm_args(algorithm: str) -> list:
    normalized = algorithm.strip().lower().replace("-", "_")
    if normalized in {"message_passing", "mp", "default"}:
        return []
    if normalized in {"spheres", "sphere", "s"}:
        return ["--spheres"]
    if normalized in {"connected_components", "connected_comp", "connected", "cc", "c"}:
        return ["--connected-comp"]
    raise ValueError(
        "starcode_clustering_algorithm must be one of: message_passing, spheres, connected_components"
    )


def is_starcode_compatible_barcode(barcode: str) -> bool:
    return bool(barcode) and set(barcode) <= {"A", "C", "G", "T"}


def run_starcode_collapse(
    barcode_counts: Counter,
    starcode_path: str,
    max_distance: int,
    cluster_ratio: float,
    algorithm: str,
    threads: int,
) -> Tuple[Counter, Dict[str, dict]]:
    """Collapse barcode counts with Starcode and return canonical cluster counts."""
    validate_starcode_config(max_distance, cluster_ratio, threads)

    compatible_counts = Counter(
        {barcode: count for barcode, count in barcode_counts.items() if is_starcode_compatible_barcode(barcode)}
    )
    unclustered_counts = Counter(
        {barcode: count for barcode, count in barcode_counts.items() if not is_starcode_compatible_barcode(barcode)}
    )

    cluster_info: Dict[str, dict] = {}

    if compatible_counts:
        with tempfile.TemporaryDirectory(prefix="insertdiversity_starcode_") as tmpdir:
            input_path = Path(tmpdir) / "barcodes.tsv"
            output_path = Path(tmpdir) / "starcode.tsv"

            with open(input_path, "w", newline="") as handle:
                for barcode, count in compatible_counts.most_common():
                    handle.write(f"{barcode}\t{count}\n")

            command = [
                starcode_path,
                "-i",
                str(input_path),
                "-o",
                str(output_path),
                "-d",
                str(max_distance),
                "-r",
                str(cluster_ratio),
                "--print-clusters",
                "-t",
                str(threads),
                "-q",
                *starcode_algorithm_args(algorithm),
            ]

            try:
                subprocess.run(command, check=True, capture_output=True, text=True)
            except FileNotFoundError as exc:
                raise FileNotFoundError(f"Could not start starcode executable: {starcode_path}") from exc
            except subprocess.CalledProcessError as exc:
                detail = (exc.stderr or exc.stdout or "").strip()
                message = f"starcode failed with exit code {exc.returncode}: {' '.join(command)}"
                if detail:
                    message += f"\n{detail}"
                raise RuntimeError(message) from exc

            with open(output_path, "r", newline="") as handle:
                for line_number, line in enumerate(handle, start=1):
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 2:
                        raise ValueError(f"Unexpected starcode output at line {line_number}: {line.rstrip()}")

                    canonical = fields[0].strip().upper()
                    if canonical in cluster_info:
                        raise ValueError(f"Duplicate starcode canonical barcode at line {line_number}: {canonical}")
                    try:
                        collapsed_count = int(fields[1])
                    except ValueError as exc:
                        raise ValueError(
                            f"Unexpected starcode cluster size at line {line_number}: {fields[1]}"
                        ) from exc

                    members = [m.strip().upper() for m in fields[2].split(",") if m.strip()] if len(fields) > 2 else [canonical]
                    if canonical not in members:
                        members.insert(0, canonical)
                    unknown_members = sorted(set(members) - set(compatible_counts))
                    if unknown_members:
                        raise ValueError(
                            f"starcode output line {line_number} included barcode(s) not present in input: "
                            + ", ".join(unknown_members[:10])
                        )

                    member_exact_counts = {member: compatible_counts[member] for member in members}
                    expected_cluster_size = sum(member_exact_counts.values())
                    if collapsed_count != expected_cluster_size:
                        raise ValueError(
                            f"Unexpected starcode cluster size at line {line_number}: "
                            f"{collapsed_count} != summed member counts {expected_cluster_size}"
                        )
                    cluster_info[canonical] = {
                        "exact_count": compatible_counts[canonical],
                        "collapsed_count": collapsed_count,
                        "members": members,
                        "member_exact_counts": member_exact_counts,
                        "max_distance": "",
                        "is_real": True,
                        "not_clustered_reason": "",
                    }

        clustered_members = {
            member
            for cluster in cluster_info.values()
            for member in cluster["members"]
        }
        missing_members = sorted(set(compatible_counts) - clustered_members)
        if missing_members:
            raise RuntimeError(
                "starcode output did not include all input barcodes. Missing examples: "
                + ", ".join(missing_members[:10])
            )

    for barcode, count in unclustered_counts.items():
        reason = "empty_barcode" if not barcode else "non_acgt_barcode"
        cluster_info[barcode] = {
            "exact_count": count,
            "collapsed_count": count,
            "members": [barcode],
            "member_exact_counts": {barcode: count},
            "max_distance": "",
            "is_real": False,
            "not_clustered_reason": reason,
        }

    collapsed_counts = Counter(
        {barcode: cluster["collapsed_count"] for barcode, cluster in cluster_info.items()}
    )
    return collapsed_counts, cluster_info


def process_one_csv(
    in_csv: Path,
    out_csv: Path,
    out_barcode_report: Path,
    expected_len: int,
    min_barcode_len: int,
    max_barcode_len: int,
    include_empty: bool,
    min_real_barcode_count: int,
    collapse_barcodes: bool,
    collapse_max_distance: int,
    starcode_path: str,
    starcode_cluster_ratio: float,
    starcode_clustering_algorithm: str,
    starcode_threads: int,
):
    barcode_counts = Counter()
    total_rows = 0
    rows_with_extracted_barcode = 0
    rows_discarded_long_barcode = 0

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

                parsed_barcode_len = safe_int(barcode_len)
                if ok_to_count and parsed_barcode_len is not None and parsed_barcode_len > max_barcode_len:
                    ok_to_count = False
                    rows_discarded_long_barcode += 1

                if ok_to_count and barcode_len != "":
                    barcode_counts[barcode] += 1
                    rows_with_extracted_barcode += 1

                row["barcode"] = barcode
                row["barcode_len"] = barcode_len
                writer.writerow(row)

    if collapse_barcodes:
        report_counts, cluster_info = run_starcode_collapse(
            barcode_counts=barcode_counts,
            starcode_path=starcode_path,
            max_distance=collapse_max_distance,
            cluster_ratio=starcode_cluster_ratio,
            algorithm=starcode_clustering_algorithm,
            threads=starcode_threads,
        )
    else:
        report_counts = Counter(barcode_counts)
        cluster_info = {
            barcode: {
                "exact_count": count,
                "collapsed_count": count,
                "members": [barcode],
                "member_exact_counts": {barcode: count},
                "max_distance": 0,
                "is_real": count >= min_real_barcode_count,
                "not_clustered_reason": "",
            }
            for barcode, count in barcode_counts.items()
        }

    with open(out_barcode_report, "w", newline="") as f_rep:
        rep_fields = [
            "barcode",
            "barcode_len",
            "count",
            "exact_count",
            "percent",
            "is_real_barcode",
            "real_barcode_method",
            "min_real_barcode_count",
            "expected_barcode_len",
            "len_minus_expected",
            "collapse_barcodes",
            "collapse_max_distance",
            "starcode_path",
            "starcode_cluster_ratio",
            "starcode_clustering_algorithm",
            "starcode_threads",
            "num_collapsed_variants",
            "max_collapse_distance",
            "not_clustered_reason",
            "collapsed_from_barcodes",
            "collapsed_from_exact_counts",
        ]
        rep = csv.DictWriter(f_rep, fieldnames=rep_fields)
        rep.writeheader()

        denom = rows_with_extracted_barcode if rows_with_extracted_barcode > 0 else 1

        for barcode, count in report_counts.most_common():
            cluster = cluster_info[barcode]
            members = sorted(
                cluster["members"],
                key=lambda member: (-cluster["member_exact_counts"][member], member),
            )
            blen = len(barcode)
            rep.writerow(
                {
                    "barcode": barcode,
                    "barcode_len": blen,
                    "count": count,
                    "exact_count": cluster["exact_count"],
                    "percent": round(100.0 * count / denom, 6),
                    "is_real_barcode": "true" if cluster["is_real"] else "false",
                    "real_barcode_method": "starcode_canonical_cluster" if collapse_barcodes else "min_real_barcode_count",
                    "min_real_barcode_count": min_real_barcode_count,
                    "expected_barcode_len": expected_len,
                    "len_minus_expected": blen - expected_len,
                    "collapse_barcodes": "true" if collapse_barcodes else "false",
                    "collapse_max_distance": collapse_max_distance,
                    "starcode_path": starcode_path if collapse_barcodes else "",
                    "starcode_cluster_ratio": starcode_cluster_ratio if collapse_barcodes else "",
                    "starcode_clustering_algorithm": starcode_clustering_algorithm if collapse_barcodes else "",
                    "starcode_threads": starcode_threads if collapse_barcodes else "",
                    "num_collapsed_variants": len(members) - 1,
                    "max_collapse_distance": cluster["max_distance"],
                    "not_clustered_reason": cluster["not_clustered_reason"],
                    "collapsed_from_barcodes": ";".join(members[1:]),
                    "collapsed_from_exact_counts": ";".join(
                        str(cluster["member_exact_counts"][member]) for member in members[1:]
                    ),
                }
            )

    return {
        "total_rows": total_rows,
        "rows_with_extracted_barcode": rows_with_extracted_barcode,
        "rows_discarded_long_barcode": rows_discarded_long_barcode,
        "barcode_counts": report_counts,
        "exact_barcode_counts": barcode_counts,
        "unique_barcodes": len(barcode_counts),
        "collapsed_unique_barcodes": len(report_counts),
        "real_unique_barcodes": sum(1 for cluster in cluster_info.values() if cluster["is_real"]),
        "total_barcode_observations": sum(barcode_counts.values()),
        "total_collapsed_barcode_observations": sum(report_counts.values()),
        "max_barcode_count": max(barcode_counts.values()) if barcode_counts else 0,
        "min_barcode_count": min(barcode_counts.values()) if barcode_counts else 0,
        "max_collapsed_barcode_count": max(report_counts.values()) if report_counts else 0,
        "min_collapsed_barcode_count": min(report_counts.values()) if report_counts else 0,
        "collapsed_variants": len(barcode_counts) - len(report_counts),
        "starcode_unclustered_barcodes": sum(1 for cluster in cluster_info.values() if cluster["not_clustered_reason"]),
    }


def main():
    ap = argparse.ArgumentParser(description="Extract barcode(s) from bun-match CSV(s) and summarize counts.")
    ap.add_argument("--config", required=True, help="Path to JSON config file")
    args = ap.parse_args()

    cfg_path = Path(args.config)
    if not cfg_path.is_file():
        raise FileNotFoundError(f"Config not found: {cfg_path}")
    cfg = load_script_config(cfg_path, "barcode_from_bun_csv")

    expected_len = int(cfg.get("expected_barcode_len", 25))
    min_barcode_len = int(cfg.get("min_barcode_len", 1))
    max_barcode_len = int(cfg.get("max_barcode_len", 100))
    if max_barcode_len < min_barcode_len:
        raise ValueError("max_barcode_len must be >= min_barcode_len")
    include_empty = parse_bool(cfg.get("include_empty_barcodes", False))
    min_real_barcode_count = int(cfg.get("min_real_barcode_count", 3))
    collapse_barcodes = parse_bool(cfg.get("collapse_barcodes", False))
    collapse_max_distance = int(cfg.get("starcode_distance", cfg.get("collapse_max_distance", 1)))
    starcode_cluster_ratio = float(
        cfg.get("starcode_cluster_ratio", cfg.get("collapse_min_abundance_ratio", 5.0))
    )
    starcode_clustering_algorithm = str(cfg.get("starcode_clustering_algorithm", "message_passing"))
    starcode_threads = int(cfg.get("starcode_threads", 1))
    starcode_path = str(cfg.get("starcode_path", "starcode"))
    if collapse_barcodes:
        validate_starcode_config(collapse_max_distance, starcode_cluster_ratio, starcode_threads)
        starcode_path = resolve_starcode_path(starcode_path)

    input_csvs = [Path(p) for p in discover_input_csvs(cfg)]

    out_csv = cfg.get("out_csv")
    out_barcode_report = cfg.get("out_barcode_report")
    out_csv_dir = Path(cfg.get("out_csv_dir", "."))
    out_report_dir = Path(cfg.get("out_barcode_report_dir", "."))
    out_multi = Path(cfg.get("out_multi_file_barcode_report", "barcode_multi_file_membership.csv"))
    out_file_summary = Path(cfg.get("out_file_summary_report", "barcode_file_summary.csv"))

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
    per_file_stats: Dict[str, dict] = {}

    max_workers = configured_threads(cfg, workload_size=len(output_plan))
    print(f"Using {max_workers} thread(s) for {len(output_plan)} CSV file(s) (available: {available_threads()}).")
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_item = {
            executor.submit(
                process_one_csv,
                in_csv=item["in"],
                out_csv=item["out_csv"],
                out_barcode_report=item["out_report"],
                expected_len=expected_len,
                min_barcode_len=min_barcode_len,
                max_barcode_len=max_barcode_len,
                include_empty=include_empty,
                min_real_barcode_count=min_real_barcode_count,
                collapse_barcodes=collapse_barcodes,
                collapse_max_distance=collapse_max_distance,
                starcode_path=starcode_path,
                starcode_cluster_ratio=starcode_cluster_ratio,
                starcode_clustering_algorithm=starcode_clustering_algorithm,
                starcode_threads=starcode_threads,
            ): item
            for item in output_plan
        }

        for future in concurrent.futures.as_completed(future_to_item):
            item = future_to_item[future]
            stats = future.result()
            print(f"Done: {item['in']}")
            print(f"  Input rows (reads in CSV):                {stats['total_rows']:,}")
            print(f"  Rows with extracted barcodes:             {stats['rows_with_extracted_barcode']:,}")
            print(f"  Rows discarded over {max_barcode_len} bp:              {stats['rows_discarded_long_barcode']:,}")
            print(f"  Unique barcodes (all extracted):          {stats['unique_barcodes']:,}")
            print(f"  Unique barcodes after collapsing:         {stats['collapsed_unique_barcodes']:,}")
            print(f"  Real unique barcodes:                     {stats['real_unique_barcodes']:,}")
            if collapse_barcodes and stats["starcode_unclustered_barcodes"]:
                print(f"  Barcodes not sent to starcode:            {stats['starcode_unclustered_barcodes']:,}")
            print(f"  Augmented CSV written to:                 {item['out_csv']}")
            print(f"  Barcode report written to:                {item['out_report']}")

            for bc in stats["barcode_counts"].keys():
                barcode_presence[bc].add(item["label"])
            per_file_stats[item["label"]] = stats

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

    out_file_summary.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file_summary, "w", newline="") as f_summary:
        fields = [
            "file_label",
            "input_csv",
            "total_reads_in_csv",
            "reads_with_barcodes",
            "reads_discarded_long_barcode",
            "max_barcode_len",
            "barcode_observations_total",
            "unique_barcodes",
            "collapsed_unique_barcodes",
            "collapsed_variants",
            "real_unique_barcodes",
            "barcodes_unique_to_this_file",
            "barcodes_shared_with_other_files",
            "max_barcode_count",
            "min_barcode_count",
            "max_collapsed_barcode_count",
            "min_collapsed_barcode_count",
            "mean_count_per_unique_barcode",
            "mean_count_per_collapsed_unique_barcode",
            "percent_reads_with_barcodes",
            "percent_barcodes_unique_to_file",
            "percent_barcodes_shared",
            "min_real_barcode_count",
            "collapse_barcodes",
            "collapse_max_distance",
            "starcode_path",
            "starcode_cluster_ratio",
            "starcode_clustering_algorithm",
            "starcode_threads",
            "starcode_unclustered_barcodes",
        ]
        writer = csv.DictWriter(f_summary, fieldnames=fields)
        writer.writeheader()
        for item in sorted(output_plan, key=lambda x: x["label"]):
            label = item["label"]
            stats = per_file_stats[label]
            bc_counts = stats["barcode_counts"]
            unique_in_file = sum(1 for bc in bc_counts if len(barcode_presence[bc]) == 1)
            shared_in_file = sum(1 for bc in bc_counts if len(barcode_presence[bc]) > 1)
            per_file_unique_counts[label] = unique_in_file
            total_unique = stats["unique_barcodes"] or 1
            collapsed_total_unique = stats["collapsed_unique_barcodes"] or 1
            writer.writerow(
                {
                    "file_label": label,
                    "input_csv": item["in"],
                    "total_reads_in_csv": stats["total_rows"],
                    "reads_with_barcodes": stats["rows_with_extracted_barcode"],
                    "reads_discarded_long_barcode": stats["rows_discarded_long_barcode"],
                    "max_barcode_len": max_barcode_len,
                    "barcode_observations_total": stats["total_barcode_observations"],
                    "unique_barcodes": stats["unique_barcodes"],
                    "collapsed_unique_barcodes": stats["collapsed_unique_barcodes"],
                    "collapsed_variants": stats["collapsed_variants"],
                    "real_unique_barcodes": stats["real_unique_barcodes"],
                    "barcodes_unique_to_this_file": unique_in_file,
                    "barcodes_shared_with_other_files": shared_in_file,
                    "max_barcode_count": stats["max_barcode_count"],
                    "min_barcode_count": stats["min_barcode_count"],
                    "max_collapsed_barcode_count": stats["max_collapsed_barcode_count"],
                    "min_collapsed_barcode_count": stats["min_collapsed_barcode_count"],
                    "mean_count_per_unique_barcode": round(stats["total_barcode_observations"] / total_unique, 6),
                    "mean_count_per_collapsed_unique_barcode": round(
                        stats["total_collapsed_barcode_observations"] / collapsed_total_unique, 6
                    ),
                    "percent_reads_with_barcodes": round(100.0 * stats["rows_with_extracted_barcode"] / (stats["total_rows"] or 1), 6),
                    "percent_barcodes_unique_to_file": round(100.0 * unique_in_file / collapsed_total_unique, 6),
                    "percent_barcodes_shared": round(100.0 * shared_in_file / collapsed_total_unique, 6),
                    "min_real_barcode_count": min_real_barcode_count,
                    "collapse_barcodes": "true" if collapse_barcodes else "false",
                    "collapse_max_distance": collapse_max_distance,
                    "starcode_path": starcode_path if collapse_barcodes else "",
                    "starcode_cluster_ratio": starcode_cluster_ratio if collapse_barcodes else "",
                    "starcode_clustering_algorithm": starcode_clustering_algorithm if collapse_barcodes else "",
                    "starcode_threads": starcode_threads if collapse_barcodes else "",
                    "starcode_unclustered_barcodes": stats["starcode_unclustered_barcodes"],
                }
            )
    print(f"Per-file summary report:                    {out_file_summary}")


if __name__ == "__main__":
    run_cli(main)
