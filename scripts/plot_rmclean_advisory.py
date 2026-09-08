#!/usr/bin/env python3
"""Plot rmclean_cubes' T15 Pass-0 build-time advisory (per-block
predicted RMSF table-GENERATION time), optionally overlaid against a
run's own measured per-block TOTAL wall time (table-gen + CLEAN).

The predicted curve comes from <outfile>.advisory.csv, written by
report_build_time_advisory (src/rmclean_cubes.f90, T15) right after
Pass 0 completes -- one row per block: hits/admitted/evicted/declined
counts from the full Belady simulation, the measured ms/channel rate,
and the predicted table-generation minutes. See
docs/dev/RMCLEAN_INTEGRATION_PLAN.md T15 for the full design.

The actual curve (optional) comes from a run's own debug log
('thread_timing stage=clean event=start|done ... block=N' lines,
log_level=debug) -- the SAME per-block wall time reported to the user
as "Block N of 32 processed."; this script computes it directly from
the start/done timestamps rather than trusting that summary line's own
formatting.

--csv is optional if --actual-log is given: the same run's log already
contains the line report_build_time_advisory prints right after
writing the CSV ('Wrote <path>.advisory.csv (...)'), so it can be
found there automatically -- this only works for a log produced by the
CURRENT binary; a log from before this CSV-writing feature existed has
no such line, and --csv must be given explicitly.

Example:
  scripts/plot_rmclean_advisory.py \\
      --csv out_cleaned.advisory.csv \\
      --out scratch/rmclean_advisory.png
  scripts/plot_rmclean_advisory.py \\
      --actual-log rmclean.run.log \\
      --out scratch/rmclean_advisory_vs_actual.png
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

THREAD_TIMING_RE = re.compile(
    r"^(?P<ts>\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}).*"
    r"stage=clean\s+event=(?P<event>start|done)\s+.*block=(?P<block>\d+)"
)

WROTE_CSV_RE = re.compile(r"^Wrote (?P<path>.+\.advisory\.csv) \(")


def find_advisory_csv_in_log(log_path: Path) -> Optional[Path]:
    """Find the CSV path from this SAME run's own 'Wrote ...' line."""
    with log_path.open() as f:
        for line in f:
            m = WROTE_CSV_RE.match(line)
            if not m:
                continue
            candidate = Path(m.group("path"))
            if candidate.is_file():
                return candidate
            # outfile= in the run's own cfg may have been relative to
            # wherever rmclean_cubes was launched from, not this
            # script's own cwd -- try relative to the log file itself,
            # a common layout (both under the same provenance dir).
            candidate = log_path.parent / candidate.name
            if candidate.is_file():
                return candidate
            return None
    return None


def read_advisory_csv(path: Path) -> Tuple[List[int], List[float]]:
    blocks: List[int] = []
    predicted_min: List[float] = []
    with path.open() as f:
        for row in csv.DictReader(f):
            blocks.append(int(row["block"]))
            predicted_min.append(float(row["predicted_build_min"]))
    return blocks, predicted_min


def read_actual_block_times(path: Path) -> Tuple[List[int], List[float]]:
    starts: Dict[int, datetime] = {}
    ends: Dict[int, datetime] = {}
    with path.open() as f:
        for line in f:
            m = THREAD_TIMING_RE.match(line)
            if not m:
                continue
            ts = datetime.strptime(m.group("ts"), "%Y-%m-%dT%H:%M:%S.%f")
            blk = int(m.group("block"))
            if m.group("event") == "start":
                starts[blk] = min(starts.get(blk, ts), ts)
            else:
                ends[blk] = max(ends.get(blk, ts), ts)
    blocks = sorted(b for b in starts if b in ends)
    minutes = [(ends[b] - starts[b]).total_seconds() / 60.0 for b in blocks]
    return blocks, minutes


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot rmclean_cubes' T15 build-time advisory "
        "(predicted table-generation time per block), optionally "
        "overlaid against a run's own measured total block time."
    )
    parser.add_argument(
        "--csv",
        default=None,
        help="<outfile>.advisory.csv from a rmclean_cubes run. Optional if "
        "--actual-log is given AND that log's own run already printed "
        "'Wrote ...advisory.csv' (i.e. it wasn't cache_eviction_policy=hitcount, "
        "which skips Pass 0 entirely) -- auto-detected from the log in that case.",
    )
    parser.add_argument("--out", required=True, help="Output PNG path")
    parser.add_argument(
        "--actual-log",
        default=None,
        help="Optional: a rmclean_cubes run log (log_level=debug) to overlay "
        "measured per-block TOTAL wall time against. Also used to auto-detect "
        "--csv when --csv is omitted.",
    )
    args = parser.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if args.csv:
        csv_path = Path(args.csv)
    else:
        if not args.actual_log:
            print(
                "ERROR: --csv is required unless --actual-log is given "
                "(the CSV path is then auto-detected from that log).",
                file=sys.stderr,
            )
            return 1
        log_path = Path(args.actual_log)
        if not log_path.is_file():
            print(f"ERROR: no such file: {log_path}", file=sys.stderr)
            return 1
        found = find_advisory_csv_in_log(log_path)
        if not found:
            print(
                f"ERROR: no 'Wrote ...advisory.csv' line found in {log_path} "
                "(older run, predates this feature, or cache_eviction_policy="
                "hitcount) -- pass --csv explicitly.",
                file=sys.stderr,
            )
            return 1
        csv_path = found
        print(f"auto-detected --csv {csv_path} from {log_path}")

    if not csv_path.is_file():
        print(f"ERROR: no such file: {csv_path}", file=sys.stderr)
        return 1
    pred_blocks, pred_min = read_advisory_csv(csv_path)
    if not pred_blocks:
        print(f"ERROR: {csv_path} has no rows", file=sys.stderr)
        return 1
    n_blocks = max(pred_blocks)

    actual_blocks: List[int] = []
    actual_min: List[float] = []
    if args.actual_log:
        log_path = Path(args.actual_log)
        if not log_path.is_file():
            print(f"ERROR: no such file: {log_path}", file=sys.stderr)
            return 1
        actual_blocks, actual_min = read_actual_block_times(log_path)
        if not actual_blocks:
            print(
                f"WARNING: {log_path} has no 'stage=clean' events -- "
                "was it run with log_level=debug? Plotting predicted only.",
                file=sys.stderr,
            )

    fig, ax = plt.subplots(figsize=(10, 5.5))
    if actual_blocks:
        ax.plot(
            actual_blocks,
            actual_min,
            marker="o",
            color="#1f77b4",
            label="Actual TOTAL block wall time (table-gen + CLEAN, measured)",
        )
    ax.plot(
        pred_blocks,
        pred_min,
        marker="s",
        color="#d62728",
        linestyle="--",
        label="Advisory-predicted table-GENERATION time only (Belady simulation)",
    )
    ax.set_xlim(0.5, n_blocks + 0.5)
    ax.set_xlabel(f"Block number (of {n_blocks})")
    ax.set_ylabel("Time (minutes)")
    title = "rmclean_cubes T15 advisory: predicted table-generation time per block"
    if actual_blocks:
        title += "\nvs actual total block wall time"
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.out, dpi=130)
    print(f"saved {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
