#!/usr/bin/env python3
"""
bun_extract.py

Scan one FASTA/FASTQ file OR a folder of FASTA/FASTQ(.gz) files for upstream/downstream bun matches.

Features:
  - Finds partial exact matches (>= min_match) against upstream/downstream buns.
  - Searches both forward and reverse-complement orientations.
  - Writes one CSV per input read file: {name}_bun_matches.csv.
  - For paired-end files, R2 files are skipped as primary inputs. If an R1/R2 pair exists,
    R2 is used only to verify the extracted barcode sequence from R1 reads.
  - Writes an additional cross-file barcode membership CSV with barcodes found in multiple
    files and barcodes unique to each file.
"""

import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

from Bio import SeqIO


READ_EXTENSIONS = (".fasta", ".fa", ".fna", ".fastq", ".fq")


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def open_text_maybe_gzip(path: Path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def detect_seqio_format(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    if any(name.endswith(ext) for ext in (".fasta", ".fa", ".fna")):
        return "fasta"
    if any(name.endswith(ext) for ext in (".fastq", ".fq")):
        return "fastq"
    raise ValueError(f"Unsupported file format for {path}. Expected fasta/fa/fna/fastq/fq (+ optional .gz).")


def normalize_seq(s: str) -> str:
    return "".join(str(s).split()).upper()


def strip_seq_extensions(name: str) -> str:
    n = name
    if n.lower().endswith(".gz"):
        n = n[:-3]
    for ext in READ_EXTENSIONS:
        if n.lower().endswith(ext):
            return n[: -len(ext)]
    return n


def canonical_read_id(read_id: str) -> str:
    rid = read_id.split()[0]
    rid = rid.replace("/1", "").replace("/2", "")
    return rid


def swap_read_marker_for_r2(stem: str) -> str:
    patterns = [
        (r"(_R)1(?=_|$)", r"\g<1>2"),
        (r"(-R)1(?=-|_|$)", r"\g<1>2"),
        (r"(\.R)1(?=\.|_|$)", r"\g<1>2"),
        (r"(_read)1(?=_|$)", r"\g<1>2"),
    ]
    for pattern, repl in patterns:
        swapped = re.sub(pattern, repl, stem, flags=re.IGNORECASE)
        if swapped != stem:
            return swapped
    return stem


def find_r2_partner(r1_path: Path, all_candidates: Dict[str, Path]) -> Optional[Path]:
    stem = strip_seq_extensions(r1_path.name)
    swapped = swap_read_marker_for_r2(stem)
    if swapped == stem:
        return None

    possible = []
    for ext in READ_EXTENSIONS:
        possible.append(swapped + ext)
        possible.append(swapped + ext + ".gz")

    for p in possible:
        if p in all_candidates:
            return all_candidates[p]
    return None


def is_r2_file(path: Path) -> bool:
    stem = strip_seq_extensions(path.name)
    return re.search(r"(?:^|[_\-.])R2(?:$|[_\-.])", stem, flags=re.IGNORECASE) is not None


def collect_input_files(cfg: dict) -> Iterable[Path]:
    if cfg.get("input_dir"):
        input_dir = Path(cfg["input_dir"])
        if not input_dir.is_dir():
            raise FileNotFoundError(f"input_dir not found or not a directory: {input_dir}")

        files = []
        for p in sorted(input_dir.iterdir()):
            if not p.is_file():
                continue
            lower = p.name.lower()
            if lower.endswith(".gz"):
                lower = lower[:-3]
            if any(lower.endswith(ext) for ext in READ_EXTENSIONS):
                files.append(p)
        if not files:
            raise ValueError(f"No supported FASTA/FASTQ files found in {input_dir}")
        return files

    single = cfg.get("fastq") or cfg.get("input_file")
    if not single:
        raise ValueError("Config must include one of: fastq, input_file, input_dir")
    single_path = Path(single)
    if not single_path.exists():
        raise FileNotFoundError(f"Input file not found: {single_path}")
    return [single_path]


def find_best_partial_hit(read_seq: str, bun_seq: str, min_match: int) -> Optional[Tuple[int, int, str]]:
    bun_seq = bun_seq.upper()
    read_seq = read_seq.upper()

    n = len(bun_seq)
    if min_match <= 0:
        raise ValueError("min_match must be >= 1")
    if n < min_match:
        return None

    for length in range(n, min_match - 1, -1):
        for i in range(0, n - length + 1):
            sub = bun_seq[i : i + length]
            pos = read_seq.find(sub)
            if pos != -1:
                return pos, pos + length - 1, sub

    return None


def choose_best_orientation(raw_seq: str, upstream: str, downstream: str, min_match: int):
    fwd = raw_seq
    rc = revcomp(raw_seq)

    def score(oriented: str):
        up_hit = find_best_partial_hit(oriented, upstream, min_match)
        dn_hit = find_best_partial_hit(oriented, downstream, min_match)
        hits = int(up_hit is not None) + int(dn_hit is not None)
        total_len = 0
        if up_hit:
            total_len += up_hit[1] - up_hit[0] + 1
        if dn_hit:
            total_len += dn_hit[1] - dn_hit[0] + 1
        return hits, total_len, up_hit, dn_hit

    f_hits, f_len, f_up, f_dn = score(fwd)
    r_hits, r_len, r_up, r_dn = score(rc)

    if (r_hits > f_hits) or (r_hits == f_hits and r_len > f_len):
        return "revcomp", rc, r_up, r_dn
    return "forward", fwd, f_up, f_dn


def load_r2_lookup(r2_path: Path) -> Dict[str, str]:
    lookup: Dict[str, str] = {}
    fmt = detect_seqio_format(r2_path)
    with open_text_maybe_gzip(r2_path) as r2_handle:
        for rec in SeqIO.parse(r2_handle, fmt):
            lookup[canonical_read_id(rec.id)] = normalize_seq(str(rec.seq))
    return lookup


def verify_barcode_with_r2(barcode: str, r2_seq: str) -> Tuple[bool, str]:
    if not barcode or not r2_seq:
        return False, ""
    if barcode in r2_seq:
        return True, barcode
    rc_barcode = revcomp(barcode)
    if rc_barcode in r2_seq:
        return True, rc_barcode
    return False, ""


def main():
    ap = argparse.ArgumentParser(description="Extract flanking sequence relative to bun matches across one file or a folder.")
    ap.add_argument("--config", required=True, help="Path to JSON config file")
    args = ap.parse_args()

    cfg_path = Path(args.config)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config not found: {cfg_path}")
    cfg = json.loads(cfg_path.read_text())

    upstream = normalize_seq(cfg.get("upstream", ""))
    downstream = normalize_seq(cfg.get("downstream", ""))
    min_match = int(cfg.get("min_match", 10))
    out_dir = Path(cfg.get("out_dir", "."))
    cross_file_csv = Path(cfg.get("cross_file_barcode_csv", out_dir / "barcode_file_membership.csv"))
    max_reads = cfg.get("max_reads", None)
    progress_every = int(cfg.get("progress_every", 200000))

    if not upstream:
        raise ValueError("Config must include: upstream")
    if not downstream:
        raise ValueError("Config must include: downstream")
    if min_match < 1:
        raise ValueError("min_match must be >= 1")
    if min_match > len(upstream) and min_match > len(downstream):
        raise ValueError(
            f"min_match ({min_match}) is longer than BOTH buns (upstream={len(upstream)}, downstream={len(downstream)})."
        )

    input_files = [Path(p) for p in collect_input_files(cfg)]
    out_dir.mkdir(parents=True, exist_ok=True)
    cross_file_csv.parent.mkdir(parents=True, exist_ok=True)

    file_map = {p.name: p for p in input_files}
    primary_files = []
    pair_map: Dict[Path, Optional[Path]] = {}
    for f in input_files:
        if is_r2_file(f):
            continue
        r2 = find_r2_partner(f, file_map)
        primary_files.append(f)
        pair_map[f] = r2

    if not primary_files:
        raise ValueError("No primary (non-R2) input files found to process.")

    barcode_to_files = defaultdict(set)

    for in_path in primary_files:
        file_label = strip_seq_extensions(in_path.name)
        out_path = out_dir / f"{file_label}_bun_matches.csv"
        fmt = detect_seqio_format(in_path)
        r2_path = pair_map.get(in_path)
        r2_lookup = load_r2_lookup(r2_path) if r2_path else {}

        total = 0
        matched = 0

        with open_text_maybe_gzip(in_path) as in_f, open(out_path, "w", newline="") as out_f:
            reader = SeqIO.parse(in_f, fmt)
            w = csv.DictWriter(
                out_f,
                fieldnames=[
                    "source_file",
                    "read_id",
                    "orientation",
                    "match_type",
                    "raw_read_seq",
                    "oriented_read_seq",
                    "upstream_start",
                    "upstream_end",
                    "upstream_match_len",
                    "upstream_matched_seq",
                    "downstream_start",
                    "downstream_end",
                    "downstream_match_len",
                    "downstream_matched_seq",
                    "after_upstream",
                    "before_downstream",
                    "barcode",
                    "barcode_len",
                    "r2_barcode_verified",
                    "r2_verified_barcode_seq",
                    "r2_partner_file",
                ],
            )
            w.writeheader()

            for rec in reader:
                total += 1
                if max_reads is not None and total > int(max_reads):
                    break

                raw_seq = normalize_seq(str(rec.seq))
                if not raw_seq:
                    continue

                orientation, oriented, up_hit, dn_hit = choose_best_orientation(raw_seq, upstream, downstream, min_match)
                has_up = up_hit is not None
                has_dn = dn_hit is not None

                if not (has_up or has_dn):
                    continue

                matched += 1

                if has_up:
                    up_start, up_end, up_sub = up_hit
                    up_len = up_end - up_start + 1
                    after_up = oriented[up_end + 1 :] if up_end + 1 <= len(oriented) else ""
                else:
                    up_start = up_end = up_len = ""
                    up_sub = ""
                    after_up = ""

                if has_dn:
                    dn_start, dn_end, dn_sub = dn_hit
                    dn_len = dn_end - dn_start + 1
                    before_dn = oriented[:dn_start] if dn_start >= 0 else ""
                else:
                    dn_start = dn_end = dn_len = ""
                    dn_sub = ""
                    before_dn = ""

                if has_up and has_dn:
                    match_type = "both"
                elif has_up:
                    match_type = "upstream"
                else:
                    match_type = "downstream"

                barcode = ""
                barcode_len = ""
                r2_verified = ""
                r2_verified_seq = ""

                if has_up and has_dn and up_end < dn_start:
                    barcode = oriented[up_end + 1 : dn_start]
                    barcode_len = len(barcode)

                    if r2_lookup:
                        r2_seq = r2_lookup.get(canonical_read_id(rec.id), "")
                        ok, verified_seq = verify_barcode_with_r2(barcode, r2_seq)
                        r2_verified = "true" if ok else "false"
                        r2_verified_seq = verified_seq
                        if ok:
                            barcode_to_files[barcode].add(file_label)
                    else:
                        barcode_to_files[barcode].add(file_label)

                w.writerow(
                    {
                        "source_file": file_label,
                        "read_id": rec.id,
                        "orientation": orientation,
                        "match_type": match_type,
                        "raw_read_seq": raw_seq,
                        "oriented_read_seq": oriented,
                        "upstream_start": up_start,
                        "upstream_end": up_end,
                        "upstream_match_len": up_len,
                        "upstream_matched_seq": up_sub,
                        "downstream_start": dn_start,
                        "downstream_end": dn_end,
                        "downstream_match_len": dn_len,
                        "downstream_matched_seq": dn_sub,
                        "after_upstream": after_up,
                        "before_downstream": before_dn,
                        "barcode": barcode,
                        "barcode_len": barcode_len,
                        "r2_barcode_verified": r2_verified,
                        "r2_verified_barcode_seq": r2_verified_seq,
                        "r2_partner_file": r2_path.name if r2_path else "",
                    }
                )

                if progress_every and total % progress_every == 0:
                    print(f"[{file_label}] processed={total:,} matched={matched:,}", flush=True)

        paired_note = f" (R2 verification: {r2_path.name})" if r2_path else ""
        print(f"Done {in_path.name}{paired_note}")
        print(f"  Reads processed: {total:,}")
        print(f"  Reads matched:   {matched:,}")
        print(f"  Output CSV:      {out_path}")

    with open(cross_file_csv, "w", newline="") as out_cross:
        w = csv.DictWriter(
            out_cross,
            fieldnames=["barcode", "num_files", "files", "classification"],
        )
        w.writeheader()
        for barcode, files in sorted(barcode_to_files.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            file_list = sorted(files)
            w.writerow(
                {
                    "barcode": barcode,
                    "num_files": len(file_list),
                    "files": ";".join(file_list),
                    "classification": "multiple_files" if len(file_list) > 1 else f"unique_to_{file_list[0]}",
                }
            )

    print(f"Cross-file barcode membership CSV: {cross_file_csv}")


if __name__ == "__main__":
    main()
