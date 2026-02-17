#!/usr/bin/env python3
"""
bun_extract.py

For each read in a FASTQ/FASTQ.GZ:
  - Look for partial EXACT matches (>= min_match bp) to upstream and/or downstream bun
  - Search is done on both the read and its reverse-complement
  - Choose an "oriented" version of the read (forward or revcomp) such that matches
    correspond to the buns as provided in config (5'->3' plasmid orientation).

For matched reads, write a CSV including:
  - read_id
  - raw_read_seq
  - oriented_read_seq
  - match_type: upstream / downstream / both
  - upstream match details (start/end/len/matched_seq) if present
  - downstream match details (start/end/len/matched_seq) if present
  - after_upstream: everything after the upstream match (3' side of match in oriented read)
  - before_downstream: everything before the downstream match (5' side of match in oriented read)

Config:
  JSON file (see example below).

Dependencies:
  conda install -c bioconda -c conda-forge biopython

Example config.json:
{
  "fastq": "reads.fastq.gz",
  "upstream": "ACGTACGTACGTACGTACGTACGTA",
  "downstream": "TTGGAATTCCGGAACTTCCGGAATT",
  "min_match": 10,
  "out_csv": "bun_matches.csv",
  "max_reads": null,
  "progress_every": 200000
}
"""

import argparse
import csv
import gzip
import json
from pathlib import Path
from typing import Optional, Tuple

from Bio import SeqIO


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def open_text_maybe_gzip(path: str):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def normalize_seq(s: str) -> str:
    # Remove whitespace/newlines and uppercase
    return "".join(str(s).split()).upper()


def find_best_partial_hit(read_seq: str, bun_seq: str, min_match: int) -> Optional[Tuple[int, int, str]]:
    """
    Find the best (longest) EXACT partial match of bun_seq within read_seq.
    Returns (start, end_inclusive, matched_subseq) or None if no match >= min_match.

    Strategy (buns are short, <= ~25bp typical):
      - Try substring lengths from len(bun) down to min_match
      - For each length, try all bun substrings of that length
      - Return the FIRST occurrence in the read for the first length that matches
        (i.e., longest match wins; ties broken by earliest found in bun scanning).
    """
    bun_seq = bun_seq.upper()
    read_seq = read_seq.upper()

    n = len(bun_seq)
    if min_match <= 0:
        raise ValueError("min_match must be >= 1")
    if n < min_match:
        return None

    for L in range(n, min_match - 1, -1):
        for i in range(0, n - L + 1):
            sub = bun_seq[i : i + L]
            pos = read_seq.find(sub)
            if pos != -1:
                start = pos
                end = pos + L - 1  # inclusive
                return start, end, sub

    return None


def choose_best_orientation(
    raw_seq: str,
    upstream: str,
    downstream: str,
    min_match: int,
) -> Tuple[str, str, Optional[Tuple[int, int, str]], Optional[Tuple[int, int, str]]]:
    """
    Evaluate forward and reverse-complement orientations.
    For each orientation, detect best partial hits to upstream/downstream.
    Choose orientation by:
      1) more hits (both > one > none)
      2) larger total matched length
      3) prefer forward if tie

    Returns:
      (orientation, oriented_seq, up_hit, dn_hit)
    where each hit is (start, end_inclusive, matched_subseq) or None.
    """
    fwd = raw_seq
    rc = revcomp(raw_seq)

    def score(oriented: str):
        up_hit = find_best_partial_hit(oriented, upstream, min_match)
        dn_hit = find_best_partial_hit(oriented, downstream, min_match)

        hits = int(up_hit is not None) + int(dn_hit is not None)
        total_len = 0
        if up_hit:
            total_len += (up_hit[1] - up_hit[0] + 1)
        if dn_hit:
            total_len += (dn_hit[1] - dn_hit[0] + 1)

        return hits, total_len, up_hit, dn_hit

    f_hits, f_len, f_up, f_dn = score(fwd)
    r_hits, r_len, r_up, r_dn = score(rc)

    if (r_hits > f_hits) or (r_hits == f_hits and r_len > f_len):
        return "revcomp", rc, r_up, r_dn
    else:
        return "forward", fwd, f_up, f_dn


def main():
    ap = argparse.ArgumentParser(description="Extract flanking sequence relative to partial upstream/downstream bun matches.")
    ap.add_argument("--config", required=True, help="Path to JSON config file")
    args = ap.parse_args()

    cfg_path = Path(args.config)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config not found: {cfg_path}")

    cfg = json.loads(cfg_path.read_text())

    fastq = cfg.get("fastq")
    upstream = normalize_seq(cfg.get("upstream", ""))
    downstream = normalize_seq(cfg.get("downstream", ""))
    min_match = int(cfg.get("min_match", 10))
    out_csv = cfg.get("out_csv", "bun_matches.csv")
    max_reads = cfg.get("max_reads", None)
    progress_every = int(cfg.get("progress_every", 200000))

    if not fastq:
        raise ValueError("Config must include: fastq")
    if not upstream:
        raise ValueError("Config must include: upstream (sequence string)")
    if not downstream:
        raise ValueError("Config must include: downstream (sequence string)")
    if min_match < 1:
        raise ValueError("min_match must be >= 1")

    # optional sanity: warn if min_match exceeds bun lengths
    if min_match > len(upstream) and min_match > len(downstream):
        raise ValueError(
            f"min_match ({min_match}) is longer than BOTH buns "
            f"(upstream={len(upstream)}, downstream={len(downstream)})."
        )

    out_path = Path(out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    matched = 0

    with open_text_maybe_gzip(fastq) as in_f, open(out_path, "w", newline="") as out_f:
        reader = SeqIO.parse(in_f, "fastq")
        w = csv.DictWriter(
            out_f,
            fieldnames=[
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

            orientation, oriented, up_hit, dn_hit = choose_best_orientation(
                raw_seq, upstream, downstream, min_match
            )

            has_up = up_hit is not None
            has_dn = dn_hit is not None

            if not (has_up or has_dn):
                continue  # throw out

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

            w.writerow(
                {
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
                }
            )

            if progress_every and total % progress_every == 0:
                print(f"[progress] processed={total:,} matched={matched:,}", flush=True)

    print("Done.")
    print(f"Reads processed: {total:,}")
    print(f"Reads matched:   {matched:,}")
    print(f"Output CSV:      {out_path}")


if __name__ == "__main__":
    main()
