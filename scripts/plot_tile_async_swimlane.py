#!/usr/bin/env python3
"""Generate a clean 3-lane swim-lane plot from RM-synthesis debug logs.

Style:
- Lanes: GPU, CPU, I/O
- Colors: GPU compute, CPU prep, CPU scatter, I/O read/write
- X-axis: seconds since selected run start

Example:
  /home/wasim/venv/rmtool/bin/python scripts/plot_tile_async_swimlane.py \
      --log MY_CASA_RMSYNTH_FULLIM_TEST.run.log \
      --out scratch/tile_async_swimlane_clean_latest.png \
      --run latest
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple


LINE_RE = re.compile(
    r"^(?P<ts>\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}).*"
    r"\[(?P<cat>[^\]]+)\]\s+\[tid=\d+\]\s+(?P<msg>.*)$"
)

ASYNC_RE = re.compile(
    r"async\s+(?P<kind>enqueue|start|done)\s+"
    r"(?P<stage>compute|scatter)\s+slot=(?P<slot>\d+)\s+sub=(?P<sub>\d+)"
)

SUBBLOCK_RE = re.compile(
    r"(?P<label>prep|send|start|done)\s+(?P<sub>\d+)/(?P<tot>\d+)\s+y:\[(?P<yb>\d+),(?P<ye>\d+)\]"
)

STAGE_BOUNDS_RE = re.compile(r"\b(?P<ev>start|done)\b.*x:\[")


@dataclass(frozen=True)
class Event:
    ts: datetime
    category: str
    message: str
    label: Optional[str] = None
    sub: Optional[int] = None
    kind: Optional[str] = None
    stage: Optional[str] = None
    slot: Optional[int] = None


@dataclass(frozen=True)
class Interval:
    lane: str
    kind: str
    label: str
    start: datetime
    end: datetime
    sub: Optional[int] = None
    slot: Optional[int] = None


@dataclass(frozen=True)
class RunWindow:
    start: datetime
    end: Optional[datetime]


@dataclass(frozen=True)
class PhaseRow:
    sub: int
    h2d_proxy_s: float
    kernel_s: float
    d2h_proxy_s: float


def parse_events(log_path: Path) -> Tuple[List[Event], List[datetime]]:
    events: List[Event] = []
    run_starts: List[datetime] = []

    with log_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            m = LINE_RE.search(line)
            if not m:
                continue

            ts = datetime.strptime(m.group("ts"), "%Y-%m-%dT%H:%M:%S.%f")
            cat = m.group("cat")
            msg = m.group("msg")

            if cat == "startup" and "run started" in msg:
                run_starts.append(ts)

            if cat == "tile_async":
                ma = ASYNC_RE.search(msg)
                if ma:
                    events.append(
                        Event(
                            ts=ts,
                            category=cat,
                            message=msg,
                            kind=ma.group("kind"),
                            stage=ma.group("stage"),
                            slot=int(ma.group("slot")),
                            sub=int(ma.group("sub")),
                        )
                    )
                continue

            if cat in {"tile_prep", "tile_compute", "tile_scatter"}:
                ms = SUBBLOCK_RE.search(msg)
                if ms:
                    events.append(
                        Event(
                            ts=ts,
                            category=cat,
                            message=msg,
                            label=ms.group("label"),
                            sub=int(ms.group("sub")),
                        )
                    )
                continue

            if cat in {"tile_read", "tile_write"}:
                mb = STAGE_BOUNDS_RE.search(msg)
                if mb:
                    events.append(
                        Event(
                            ts=ts,
                            category=cat,
                            message=msg,
                            label=mb.group("ev"),
                        )
                    )

    events.sort(key=lambda e: e.ts)
    run_starts.sort()
    return events, run_starts


def split_runs(
    events: List[Event], run_starts: List[datetime], gap_seconds: float
) -> List[List[Event]]:
    if not events:
        return []

    if run_starts:
        windows: List[RunWindow] = []
        for i, st in enumerate(run_starts):
            en = run_starts[i + 1] if i + 1 < len(run_starts) else None
            windows.append(RunWindow(st, en))

        by_window: List[List[Event]] = []
        for w in windows:
            run_events = [
                ev for ev in events if ev.ts >= w.start and (w.end is None or ev.ts < w.end)
            ]
            if run_events:
                by_window.append(run_events)
        if by_window:
            return by_window

    runs: List[List[Event]] = [[events[0]]]
    for ev in events[1:]:
        dt = (ev.ts - runs[-1][-1].ts).total_seconds()
        if dt > gap_seconds:
            runs.append([ev])
        else:
            runs[-1].append(ev)
    return runs


def select_run(runs: List[List[Event]], run_selector: str) -> List[Event]:
    if not runs:
        return []
    if run_selector == "latest":
        return runs[-1]
    if run_selector == "first":
        return runs[0]

    idx = int(run_selector)
    if idx < 0:
        idx = len(runs) + idx
    if idx < 0 or idx >= len(runs):
        raise IndexError(f"Run index {run_selector} out of range (0..{len(runs)-1})")
    return runs[idx]


def _pair_intervals(starts: List[datetime], ends: List[datetime]) -> List[Tuple[datetime, datetime]]:
    pairs: List[Tuple[datetime, datetime]] = []
    n = min(len(starts), len(ends))
    for i in range(n):
        if ends[i] > starts[i]:
            pairs.append((starts[i], ends[i]))
    return pairs


def build_intervals(events: List[Event]) -> List[Interval]:
    if not events:
        return []

    prep_open_ts: Dict[int, Optional[datetime]] = {}
    scatter_start_ts: Dict[int, List[datetime]] = {}
    scatter_done_ts: Dict[int, List[datetime]] = {}
    gpu_compute_start_ts: Dict[Tuple[int, int], List[datetime]] = {}
    gpu_compute_done_ts: Dict[Tuple[int, int], List[datetime]] = {}

    io_read_start: List[datetime] = []
    io_read_done: List[datetime] = []
    io_write_start: List[datetime] = []
    io_write_done: List[datetime] = []

    intervals: List[Interval] = []

    for ev in events:
        if ev.category == "tile_prep" and ev.label == "prep" and ev.sub is not None:
            if prep_open_ts.get(ev.sub) is None:
                prep_open_ts[ev.sub] = ev.ts
        elif ev.category == "tile_compute" and ev.label == "send" and ev.sub is not None:
            prep_start = prep_open_ts.get(ev.sub)
            if prep_start is not None and ev.ts > prep_start:
                intervals.append(
                    Interval(
                        "CPU", "cpu_prep", f"prep {ev.sub}", prep_start, ev.ts, sub=ev.sub
                    )
                )
            prep_open_ts[ev.sub] = None
        elif ev.category == "tile_scatter" and ev.label == "start" and ev.sub is not None:
            scatter_start_ts.setdefault(ev.sub, []).append(ev.ts)
        elif ev.category == "tile_scatter" and ev.label == "done" and ev.sub is not None:
            scatter_done_ts.setdefault(ev.sub, []).append(ev.ts)
        elif (
            ev.category == "tile_async"
            and ev.stage == "compute"
            and ev.sub is not None
            and ev.slot is not None
            and ev.kind is not None
        ):
            key = (ev.sub, ev.slot)
            if ev.kind == "start":
                gpu_compute_start_ts.setdefault(key, []).append(ev.ts)
            elif ev.kind == "done":
                gpu_compute_done_ts.setdefault(key, []).append(ev.ts)
        elif ev.category == "tile_read" and ev.label in {"start", "done"}:
            if ev.label == "start":
                io_read_start.append(ev.ts)
            else:
                io_read_done.append(ev.ts)
        elif ev.category == "tile_write" and ev.label in {"start", "done"}:
            if ev.label == "start":
                io_write_start.append(ev.ts)
            else:
                io_write_done.append(ev.ts)

    for s, e in _pair_intervals(io_read_start, io_read_done):
        intervals.append(Interval("I/O", "io_read", "read", s, e, sub=None))
    for s, e in _pair_intervals(io_write_start, io_write_done):
        intervals.append(Interval("I/O", "io_write", "write", s, e, sub=None))

    for sub in sorted(scatter_start_ts.keys()):
        starts = sorted(scatter_start_ts.get(sub, []))
        ends = sorted(scatter_done_ts.get(sub, []))
        n = min(len(starts), len(ends))
        for i in range(n):
            if ends[i] > starts[i]:
                intervals.append(
                    Interval(
                        "CPU",
                        "cpu_scatter",
                        f"scatter {sub}",
                        starts[i],
                        ends[i],
                        sub=sub,
                    )
                )

    for key in sorted(gpu_compute_start_ts.keys()):
        sub, slot = key
        starts = sorted(gpu_compute_start_ts.get(key, []))
        ends = sorted(gpu_compute_done_ts.get(key, []))
        n = min(len(starts), len(ends))
        for i in range(n):
            if ends[i] > starts[i]:
                intervals.append(
                    Interval(
                        "GPU",
                        f"gpu_compute_slot{slot}",
                        f"compute {sub} (slot {slot})",
                        starts[i],
                        ends[i],
                        sub=sub,
                        slot=slot,
                    )
                )

    intervals.sort(key=lambda x: (x.start, x.lane, x.kind))
    return intervals


def build_phase_rows(events: List[Event]) -> List[PhaseRow]:
    send_ts: Dict[int, datetime] = {}
    compute_start_ts: Dict[int, datetime] = {}
    compute_done_ts: Dict[int, datetime] = {}
    scatter_start_ts: Dict[int, datetime] = {}
    scatter_done_ts: Dict[int, datetime] = {}

    for ev in events:
        if ev.category == "tile_compute" and ev.label == "send" and ev.sub is not None:
            send_ts[ev.sub] = ev.ts
        elif (
            ev.category == "tile_async"
            and ev.stage == "compute"
            and ev.sub is not None
            and ev.kind is not None
        ):
            if ev.kind == "start":
                compute_start_ts[ev.sub] = ev.ts
            elif ev.kind == "done":
                compute_done_ts[ev.sub] = ev.ts
        elif (
            ev.category == "tile_async"
            and ev.stage == "scatter"
            and ev.sub is not None
            and ev.kind is not None
        ):
            if ev.kind == "start":
                scatter_start_ts[ev.sub] = ev.ts
            elif ev.kind == "done":
                scatter_done_ts[ev.sub] = ev.ts

    rows: List[PhaseRow] = []
    all_subs = sorted(
        set(send_ts.keys())
        | set(compute_start_ts.keys())
        | set(compute_done_ts.keys())
        | set(scatter_start_ts.keys())
        | set(scatter_done_ts.keys())
    )

    for sub in all_subs:
        h2d_proxy = 0.0
        kernel = 0.0
        d2h_proxy = 0.0

        if sub in send_ts and sub in compute_start_ts:
            h2d_proxy = max(0.0, (compute_start_ts[sub] - send_ts[sub]).total_seconds())
        if sub in compute_start_ts and sub in compute_done_ts:
            kernel = max(0.0, (compute_done_ts[sub] - compute_start_ts[sub]).total_seconds())
        if sub in scatter_start_ts and sub in scatter_done_ts:
            d2h_proxy = max(0.0, (scatter_done_ts[sub] - scatter_start_ts[sub]).total_seconds())

        rows.append(
            PhaseRow(
                sub=sub,
                h2d_proxy_s=h2d_proxy,
                kernel_s=kernel,
                d2h_proxy_s=d2h_proxy,
            )
        )

    return rows


def _self_overlap_seconds(intervals: List[Tuple[float, float]]) -> float:
    points: List[Tuple[float, int]] = []
    for s, e in intervals:
        points.append((s, 1))
        points.append((e, -1))
    points.sort(key=lambda x: (x[0], -x[1]))

    total = 0.0
    active = 0
    prev_t: Optional[float] = None
    for t, delta in points:
        if prev_t is not None and active >= 2 and t > prev_t:
            total += t - prev_t
        active += delta
        prev_t = t
    return total


def _cross_overlap_seconds(
    a_intervals: List[Tuple[float, float]], b_intervals: List[Tuple[float, float]]
) -> float:
    points: List[Tuple[float, str, int]] = []
    for s, e in a_intervals:
        points.append((s, "a", 1))
        points.append((e, "a", -1))
    for s, e in b_intervals:
        points.append((s, "b", 1))
        points.append((e, "b", -1))
    points.sort(key=lambda x: (x[0], -x[2]))

    total = 0.0
    active_a = 0
    active_b = 0
    prev_t: Optional[float] = None
    for t, tag, delta in points:
        if prev_t is not None and active_a > 0 and active_b > 0 and t > prev_t:
            total += t - prev_t
        if tag == "a":
            active_a += delta
        else:
            active_b += delta
        prev_t = t
    return total


def plot_clean_swimlane(
    intervals: List[Interval],
    out_path: Path,
    title: str,
    time_axis: str,
    right_info_lines: Optional[List[str]] = None,
) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Liberation Sans", "DejaVu Sans", "Noto Sans"],
            "text.antialiased": True,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
        }
    )

    if not intervals:
        raise ValueError("No intervals available to plot")

    lane_order = ["GPU", "CPU", "I/O"]
    y_map = {lane: idx for idx, lane in enumerate(lane_order)}

    colors = {
        "gpu_compute_slot1": "#4a3aa7",
        "gpu_compute_slot2": "#7d6fd1",
        "cpu_prep_odd": "#eda100",
        "cpu_prep_even": "#f4c65a",
        "cpu_scatter": "#e34948",
        "io_read": "#2a78d6",
        "io_write": "#1f4f99",
    }

    t0 = min(iv.start for iv in intervals)
    t1 = max(iv.end for iv in intervals)
    wall = max(0.001, (t1 - t0).total_seconds())

    fig, ax = plt.subplots(figsize=(14, 3.8))

    for iv in intervals:
        if time_axis == "absolute":
            left = mdates.date2num(iv.start)
            width = max(0.02, (iv.end - iv.start).total_seconds()) / 86400.0
            width_s = (iv.end - iv.start).total_seconds()
        else:
            left = (iv.start - t0).total_seconds()
            width = max(0.02, (iv.end - iv.start).total_seconds())
            width_s = width

        # Keep 3 primary lanes but add slight vertical offsets so overlapping
        # sub-block intervals remain visible and individually labelable.
        y_base = y_map[iv.lane]
        y_off = 0.0
        if iv.sub is not None and iv.lane == "CPU":
            parity_center = -0.14 if (iv.sub % 2 == 1) else 0.14
            if iv.kind == "cpu_prep":
                y_off = parity_center - 0.08
            elif iv.kind == "cpu_scatter":
                y_off = parity_center + 0.08
            else:
                y_off = parity_center
        elif iv.sub is not None and iv.lane == "GPU":
            y_off = -0.14 if (iv.sub % 2 == 1) else 0.14
        y = y_base + y_off

        bar_color = colors.get(iv.kind, "#777777")
        bar_hatch = None
        if iv.kind == "cpu_prep" and iv.sub is not None:
            bar_color = colors["cpu_prep_odd"] if (iv.sub % 2 == 1) else colors["cpu_prep_even"]
            if iv.sub % 2 == 0:
                bar_hatch = "//"
        if iv.kind == "gpu_compute_slot2":
            bar_hatch = "//"

        ax.barh(
            y=y,
            left=left,
            width=width,
            height=0.34,
            color=bar_color,
            edgecolor="black",
            linewidth=0.5,
            alpha=0.95,
            hatch=bar_hatch,
        )

    ax.set_yticks([y_map[k] for k in lane_order])
    ax.set_yticklabels(lane_order)
    # Add extra bottom margin so GPU lane does not crowd the x-axis labels.
    ax.set_ylim(-0.8, len(lane_order) - 0.2)
    if time_axis == "absolute":
        ax.set_xlabel("absolute time")
        ax.set_xlim(mdates.date2num(t0), mdates.date2num(t1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    else:
        ax.set_xlabel("seconds since band start")
        ax.set_xlim(0.0, wall)
    ax.set_ylabel("lane")
    ax.set_title(title)
    ax.grid(axis="x", linestyle="--", alpha=0.35)

    legend_items = [
        ("I/O read", "#2a78d6"),
        ("I/O write", "#1f4f99"),
        ("CPU prep (odd blocks)", "#eda100"),
        ("CPU prep (even blocks)", "#f4c65a"),
        ("CPU scatter", "#e34948"),
        ("GPU compute slot 1", "#4a3aa7"),
        ("GPU compute slot 2", "#7d6fd1"),
    ]
    handles = []
    for name, color in legend_items:
        hatch = "//" if ("slot 2" in name or "even blocks" in name) else None
        handles.append(plt.Rectangle((0, 0), 1, 1, color=color, ec="black", lw=0.5, hatch=hatch))
    labels = [n for n, _ in legend_items]
    ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(1.02, 0.98),
        ncol=1,
        fontsize=9,
        borderaxespad=0.0,
        framealpha=0.95,
    )

    if right_info_lines:
        info_text = "\n".join(right_info_lines)
        ax.text(
            1.02,
            0.42,
            info_text,
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            clip_on=False,
            bbox={
                "facecolor": "white",
                "edgecolor": "#666666",
                "alpha": 0.9,
                "boxstyle": "round,pad=0.35",
            },
        )

    # Keep room on the right for the out-of-axes legend.
    fig.tight_layout(rect=(0.0, 0.0, 0.82, 1.0))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=170)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Clean 3-lane timing swim-lane")
    parser.add_argument("--log", required=True, help="Path to run log file")
    parser.add_argument("--out", required=True, help="Output PNG path")
    parser.add_argument(
        "--run",
        default="latest",
        help="Run selector: latest, first, or numeric index (supports negative)",
    )
    parser.add_argument(
        "--gap-sec",
        type=float,
        default=300.0,
        help="Gap threshold in seconds when startup markers are unavailable",
    )
    parser.add_argument(
        "--time-axis",
        choices=["absolute", "relative"],
        default="absolute",
        help="X-axis mode: absolute wall-clock or relative seconds since run start",
    )
    args = parser.parse_args()

    log_path = Path(args.log)
    out_path = Path(args.out)

    events, run_starts = parse_events(log_path)
    if not events:
        raise SystemExit(f"No supported timing events found in {log_path}")

    runs = split_runs(events, run_starts, gap_seconds=args.gap_sec)
    run_events = select_run(runs, args.run)
    intervals = build_intervals(run_events)
    phase_rows = build_phase_rows(run_events)
    if not intervals:
        raise SystemExit("No intervals produced for selected run")

    t_start = min(iv.start for iv in intervals)
    t_end = max(iv.end for iv in intervals)
    wall = (t_end - t_start).total_seconds()

    gpu_intervals = [
        ((iv.start - t_start).total_seconds(), (iv.end - t_start).total_seconds())
        for iv in intervals
        if iv.kind.startswith("gpu_compute_slot")
    ]
    cpu_intervals = [
        ((iv.start - t_start).total_seconds(), (iv.end - t_start).total_seconds())
        for iv in intervals
        if iv.kind in {"cpu_prep", "cpu_scatter"}
    ]

    gpu_gpu_overlap = _self_overlap_seconds(gpu_intervals) if gpu_intervals else 0.0
    cpu_gpu_overlap = (
        _cross_overlap_seconds(cpu_intervals, gpu_intervals)
        if cpu_intervals and gpu_intervals
        else 0.0
    )

    title = "Process timeline"
    gpu_slot1_count = sum(1 for iv in intervals if iv.kind == "gpu_compute_slot1")
    gpu_slot2_count = sum(1 for iv in intervals if iv.kind == "gpu_compute_slot2")
    cpu_prep_count = sum(1 for iv in intervals if iv.kind == "cpu_prep")
    cpu_scatter_count = sum(1 for iv in intervals if iv.kind == "cpu_scatter")
    io_read_count = sum(1 for iv in intervals if iv.kind == "io_read")
    io_write_count = sum(1 for iv in intervals if iv.kind == "io_write")

    right_info_lines = [
        f"run: {args.run}",
        f"window_s: {wall:.3f}",
        f"intervals: {len(intervals)}",
        f"gpu_gpu_ovl_s: {gpu_gpu_overlap:.3f}",
        f"cpu_gpu_ovl_s: {cpu_gpu_overlap:.3f}",
        "-- counts --",
        f"io read/write: {io_read_count}/{io_write_count}",
        f"cpu prep/scat: {cpu_prep_count}/{cpu_scatter_count}",
        f"gpu slot1/2: {gpu_slot1_count}/{gpu_slot2_count}",
    ]

    plot_clean_swimlane(
        intervals,
        out_path,
        title,
        time_axis=args.time_axis,
        right_info_lines=right_info_lines,
    )

    print(f"log: {log_path}")
    print(f"runs found: {len(runs)}")
    print(f"selected run events: {len(run_events)}")
    print(f"intervals plotted: {len(intervals)}")
    print(f"window_s: {wall:.3f}")
    print(f"gpu_gpu_overlap_s: {gpu_gpu_overlap:.3f}")
    print(f"cpu_gpu_overlap_s: {cpu_gpu_overlap:.3f}")
    if phase_rows:
        total_h2d = sum(r.h2d_proxy_s for r in phase_rows)
        total_kernel = sum(r.kernel_s for r in phase_rows)
        total_d2h = sum(r.d2h_proxy_s for r in phase_rows)
        print("phase_breakdown_proxy_seconds:")
        for r in phase_rows:
            print(
                f"  sub={r.sub} h2d_proxy={r.h2d_proxy_s:.3f} "
                f"kernel={r.kernel_s:.3f} d2h_proxy={r.d2h_proxy_s:.3f}"
            )
        print(
            f"phase_totals_proxy_s: h2d_proxy={total_h2d:.3f} "
            f"kernel={total_kernel:.3f} d2h_proxy={total_d2h:.3f}"
        )
    print(f"wrote: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
