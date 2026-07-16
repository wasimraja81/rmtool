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
import textwrap
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple


LINE_RE = re.compile(
    r"^(?P<ts>\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}).*"
    r"\[(?P<cat>[^\]]+)\]\s+\[tid=(?P<tid>\d+)\]\s+(?P<msg>.*)$"
)

ASYNC_RE = re.compile(
    r"async\s+(?P<kind>enqueue|start|done)\s+"
    r"(?P<stage>compute|scatter)\s+slot=(?P<slot>\d+)\s+sub=(?P<sub>\d+)"
)

SUBBLOCK_RE = re.compile(
    r"(?P<label>prep|send|start|done)\s+(?P<sub>\d+)/(?P<tot>\d+)\s+y:\[(?P<yb>\d+),(?P<ye>\d+)\]"
)

STAGE_BOUNDS_RE = re.compile(r"\b(?P<ev>start|done)\b.*x:\[")

THREAD_CPU_RE = re.compile(
    r"thread_timing\s+stage=cpu_extract\s+event=(?P<event>start|done)\s+"
    r"tid=(?P<tid>\d+)\s+rm_block=(?P<rm_block>\d+)\s+nrm_now=(?P<nrm_now>\d+)"
    r"(?:\s+dur_ms=(?P<dur_ms>[0-9]+\.[0-9]+))?"
)


@dataclass(frozen=True)
class Event:
    ts: datetime
    category: str
    tid: int
    message: str
    label: Optional[str] = None
    sub: Optional[int] = None
    kind: Optional[str] = None
    stage: Optional[str] = None
    slot: Optional[int] = None


@dataclass(frozen=True)
class ThreadInterval:
    tid: int
    rm_block: int
    nrm_now: int
    start: datetime
    end: datetime
    dur_ms: Optional[float] = None


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
            tid = int(m.group("tid"))
            msg = m.group("msg")

            if cat == "startup":
                if "run started" in msg:
                    run_starts.append(ts)
                events.append(
                    Event(
                        ts=ts,
                        category=cat,
                        tid=tid,
                        message=msg,
                    )
                )
                continue

            if cat == "tile_async":
                ma = ASYNC_RE.search(msg)
                if ma:
                    events.append(
                        Event(
                            ts=ts,
                            category=cat,
                            tid=tid,
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
                            tid=tid,
                            message=msg,
                            label=ms.group("label"),
                            sub=int(ms.group("sub")),
                        )
                    )
                    continue

            if cat in {
                "tile_read",
                "tile_write",
                "tile_mask",
                "tile_prep",
                "tile_compute",
                "tile_cubestat",
            }:
                mb = STAGE_BOUNDS_RE.search(msg)
                if mb:
                    events.append(
                        Event(
                            ts=ts,
                            category=cat,
                            tid=tid,
                            message=msg,
                            label=mb.group("ev"),
                        )
                    )

            if cat == "tile_thread":
                mt = THREAD_CPU_RE.search(msg)
                if mt:
                    events.append(
                        Event(
                            ts=ts,
                            category=cat,
                            tid=tid,
                            message=msg,
                            label=mt.group("event"),
                            sub=int(mt.group("rm_block")),
                            kind="cpu_extract",
                            slot=int(mt.group("nrm_now")),
                        )
                    )

    events.sort(key=lambda e: e.ts)
    run_starts.sort()
    return events, run_starts


def build_cpu_thread_intervals(events: List[Event]) -> List[ThreadInterval]:
    starts: Dict[Tuple[int, int, int], List[datetime]] = {}
    intervals: List[ThreadInterval] = []

    for ev in events:
        if ev.category != "tile_thread" or ev.kind != "cpu_extract":
            continue
        if ev.sub is None or ev.slot is None:
            continue

        m = THREAD_CPU_RE.search(ev.message)
        dur_ms = float(m.group("dur_ms")) if (m and m.group("dur_ms")) else None
        key = (ev.tid, ev.sub, ev.slot)

        if ev.label == "start":
            starts.setdefault(key, []).append(ev.ts)
        elif ev.label == "done":
            bucket = starts.get(key)
            if bucket:
                st = bucket.pop(0)
                if ev.ts > st:
                    intervals.append(
                        ThreadInterval(
                            tid=ev.tid,
                            rm_block=ev.sub,
                            nrm_now=ev.slot,
                            start=st,
                            end=ev.ts,
                            dur_ms=dur_ms,
                        )
                    )

    intervals.sort(key=lambda x: (x.start, x.tid, x.rm_block))
    return intervals


def plot_cpu_thread_timeline(
    thread_intervals: List[ThreadInterval],
    io_intervals: List[Interval],
    cpu_stage_intervals: List[Interval],
    out_path: Path,
    title: str,
    time_axis: str,
    right_info_lines: Optional[List[str]] = None,
) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    if not thread_intervals:
        raise ValueError("No CPU thread intervals available to plot")

    all_start = [iv.start for iv in thread_intervals] + [iv.start for iv in io_intervals] + [iv.start for iv in cpu_stage_intervals]
    all_end = [iv.end for iv in thread_intervals] + [iv.end for iv in io_intervals] + [iv.end for iv in cpu_stage_intervals]
    t0 = min(all_start)
    t1 = max(all_end)
    wall = max(0.001, (t1 - t0).total_seconds())

    tids = sorted(set(iv.tid for iv in thread_intervals))
    lane_labels = [f"T{tid}" for tid in tids]
    if cpu_stage_intervals:
        lane_labels = ["CPU stage"] + lane_labels
    if io_intervals:
        lane_labels = ["I/O"] + lane_labels
    y_map = {lane: idx for idx, lane in enumerate(lane_labels)}

    fig_h = 3.2 + 0.22 * len(tids)
    fig, ax = plt.subplots(figsize=(14, fig_h))

    for iv in thread_intervals:
        lane = f"T{iv.tid}"
        y = y_map[lane]
        if time_axis == "absolute":
            left = mdates.date2num(iv.start)
            width = max(0.02, (iv.end - iv.start).total_seconds()) / 86400.0
        else:
            left = (iv.start - t0).total_seconds()
            width = max(0.02, (iv.end - iv.start).total_seconds())

        color = "#e97827" if (iv.rm_block % 2 == 1) else "#f0be64"
        hatch = "//" if (iv.rm_block % 2 == 0) else None
        ax.barh(
            y=y,
            left=left,
            width=width,
            height=0.72,
            color=color,
            edgecolor="black",
            linewidth=0.45,
            alpha=0.95,
            hatch=hatch,
        )

    for iv in io_intervals:
        y = y_map["I/O"]
        if time_axis == "absolute":
            left = mdates.date2num(iv.start)
            width = max(0.02, (iv.end - iv.start).total_seconds()) / 86400.0
        else:
            left = (iv.start - t0).total_seconds()
            width = max(0.02, (iv.end - iv.start).total_seconds())

        color = "#2a78d6" if iv.kind == "io_read" else "#1f4f99"
        ax.barh(
            y=y,
            left=left,
            width=width,
            height=0.55,
            color=color,
            edgecolor="black",
            linewidth=0.45,
            alpha=0.9,
        )

    stage_colors = {
        "cpu_stage_mask": "#85b86f",
        "cpu_stage_prep": "#b8b36a",
        "cpu_stage_compute": "#8a8a8a",
        "cpu_stage_cubestat": "#5aa3a5",
    }
    for iv in cpu_stage_intervals:
        y = y_map["CPU stage"]
        if time_axis == "absolute":
            left = mdates.date2num(iv.start)
            width = max(0.02, (iv.end - iv.start).total_seconds()) / 86400.0
        else:
            left = (iv.start - t0).total_seconds()
            width = max(0.02, (iv.end - iv.start).total_seconds())

        ax.barh(
            y=y,
            left=left,
            width=width,
            height=0.55,
            color=stage_colors.get(iv.kind, "#7f7f7f"),
            edgecolor="black",
            linewidth=0.45,
            alpha=0.9,
        )

    ax.set_yticks([y_map[k] for k in lane_labels])
    ax.set_yticklabels(lane_labels)
    ax.set_ylim(-0.8, len(lane_labels) - 0.2)
    if time_axis == "absolute":
        ax.set_xlabel("absolute time")
        ax.set_xlim(mdates.date2num(t0), mdates.date2num(t1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    else:
        ax.set_xlabel("seconds since run start")
        ax.set_xlim(0.0, wall)
    ax.set_ylabel("CPU thread")
    ax.set_title(title)
    ax.grid(axis="x", linestyle="--", alpha=0.35)
    for tick in ax.get_xticklabels():
        tick.set_rotation(30)
        tick.set_ha("right")

    handles = [
        plt.Rectangle((0, 0), 1, 1, color="#e97827", ec="black", lw=0.45),
        plt.Rectangle((0, 0), 1, 1, color="#f0be64", ec="black", lw=0.45, hatch="////"),
    ]
    labels = ["cpu_extract rm_block odd", "cpu_extract rm_block even"]
    if io_intervals:
        handles.extend(
            [
                plt.Rectangle((0, 0), 1, 1, color="#2a78d6", ec="black", lw=0.45),
                plt.Rectangle((0, 0), 1, 1, color="#1f4f99", ec="black", lw=0.45),
            ]
        )
        labels.extend(["I/O read", "I/O write"])
    if cpu_stage_intervals:
        handles.extend(
            [
                plt.Rectangle((0, 0), 1, 1, color="#85b86f", ec="black", lw=0.45),
                plt.Rectangle((0, 0), 1, 1, color="#b8b36a", ec="black", lw=0.45),
                plt.Rectangle((0, 0), 1, 1, color="#8a8a8a", ec="black", lw=0.45),
                plt.Rectangle((0, 0), 1, 1, color="#5aa3a5", ec="black", lw=0.45),
            ]
        )
        labels.extend(["CPU mask", "CPU prep", "CPU compute", "CPU cubestat"])
    fig.tight_layout(rect=(0.0, 0.0, 0.84, 1.0))
    _layout_right_panel(fig, ax, handles, labels, right_info_lines)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=170, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


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


def build_gpu_sync_intervals(events: List[Event]) -> List[Interval]:
    send_ts: Dict[int, List[datetime]] = {}
    done_ts: Dict[int, List[datetime]] = {}
    intervals: List[Interval] = []

    for ev in events:
        if ev.category != "tile_compute" or ev.sub is None:
            continue
        if ev.label == "send":
            send_ts.setdefault(ev.sub, []).append(ev.ts)
        elif ev.label == "done":
            done_ts.setdefault(ev.sub, []).append(ev.ts)

    for sub in sorted(set(send_ts.keys()) | set(done_ts.keys())):
        starts = sorted(send_ts.get(sub, []))
        ends = sorted(done_ts.get(sub, []))
        n = min(len(starts), len(ends))
        for i in range(n):
            if ends[i] > starts[i]:
                intervals.append(
                    Interval(
                        "GPU",
                        "gpu_compute_sync",
                        f"compute {sub} (sync)",
                        starts[i],
                        ends[i],
                        sub=sub,
                        slot=None,
                    )
                )

    intervals.sort(key=lambda x: (x.start, x.lane, x.kind))
    return intervals


def build_cpu_stage_intervals(events: List[Event]) -> List[Interval]:
    starts: Dict[str, List[datetime]] = {
        "tile_mask": [],
        "tile_prep": [],
        "tile_compute": [],
        "tile_cubestat": [],
    }
    ends: Dict[str, List[datetime]] = {
        "tile_mask": [],
        "tile_prep": [],
        "tile_compute": [],
        "tile_cubestat": [],
    }

    for ev in events:
        if ev.category not in starts or ev.label not in {"start", "done"}:
            continue
        if ev.label == "start":
            starts[ev.category].append(ev.ts)
        else:
            ends[ev.category].append(ev.ts)

    kind_map = {
        "tile_mask": "cpu_stage_mask",
        "tile_prep": "cpu_stage_prep",
        "tile_compute": "cpu_stage_compute",
        "tile_cubestat": "cpu_stage_cubestat",
    }
    label_map = {
        "tile_mask": "mask",
        "tile_prep": "prep",
        "tile_compute": "compute",
        "tile_cubestat": "cubestat",
    }

    intervals: List[Interval] = []
    for cat in starts:
        pairs = _pair_intervals(sorted(starts[cat]), sorted(ends[cat]))
        for s, e in pairs:
            intervals.append(
                Interval(
                    lane="CPU",
                    kind=kind_map[cat],
                    label=label_map[cat],
                    start=s,
                    end=e,
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


def _bbox_overlap_area_axes(a, b) -> float:
    if a is None or b is None:
        return 0.0
    w = max(0.0, min(a.x1, b.x1) - max(a.x0, b.x0))
    h = max(0.0, min(a.y1, b.y1) - max(a.y0, b.y0))
    return w * h


def _format_info_block(
    rows: List[Tuple[Optional[str], Optional[str]]],
    key_width: int = 22,
    value_wrap: int = 52,
) -> List[str]:
    lines: List[str] = []
    cont_indent = " " * (key_width + 3)

    for key, value in rows:
        if key is None:
            lines.append(value or "")
            continue

        vtxt = "" if value is None else str(value)
        wrapped = textwrap.wrap(vtxt, width=value_wrap) or [""]
        lines.append(f"{key:<{key_width}} : {wrapped[0]}")
        for part in wrapped[1:]:
            lines.append(f"{cont_indent}{part}")

    return lines


def _layout_right_panel(
    fig,
    ax,
    handles,
    labels,
    right_info_lines: Optional[List[str]],
):
    x_anchor = 1.05
    y_bottom = 0.02
    y_top = 0.98
    min_gap = 0.02

    if not handles and not right_info_lines:
        return

    legend_candidates = [(1, 8), (2, 8), (3, 8), (2, 7), (3, 7)]
    info_font_candidates = [9, 8, 7, 6]

    chosen_legend = None
    chosen_info = None
    chosen_metrics = None

    for ncol, legend_fs in legend_candidates:
        for info_fs in info_font_candidates:
            legend_artist = None
            info_artist = None

            if handles:
                legend_artist = ax.legend(
                    handles,
                    labels,
                    loc="lower left",
                    bbox_to_anchor=(x_anchor, y_bottom),
                    ncol=ncol,
                    fontsize=legend_fs,
                    handlelength=2.0,
                    handleheight=1.1,
                    columnspacing=1.0,
                    borderaxespad=0.0,
                    framealpha=0.9,
                )

            if right_info_lines:
                info_artist = ax.text(
                    x_anchor,
                    y_top,
                    "\n".join(right_info_lines),
                    transform=ax.transAxes,
                    va="top",
                    ha="left",
                    fontsize=info_fs,
                    fontfamily="DejaVu Sans Mono",
                    linespacing=1.05,
                    clip_on=False,
                    bbox={
                        "facecolor": "white",
                        "edgecolor": "#666666",
                        "alpha": 0.9,
                        "boxstyle": "round,pad=0.5",
                    },
                )

            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()

            legend_bbox = None
            info_bbox = None
            if legend_artist is not None:
                legend_bbox = legend_artist.get_window_extent(renderer=renderer).transformed(
                    ax.transAxes.inverted()
                )
            if info_artist is not None:
                info_bbox = info_artist.get_window_extent(renderer=renderer).transformed(
                    ax.transAxes.inverted()
                )

            legend_h = legend_bbox.height if legend_bbox is not None else 0.0
            info_h = info_bbox.height if info_bbox is not None else 0.0
            required_h = legend_h + info_h + (min_gap if legend_bbox is not None and info_bbox is not None else 0.0)
            available_h = y_top - y_bottom
            overflow = max(0.0, required_h - available_h)
            overlap_area = _bbox_overlap_area_axes(legend_bbox, info_bbox)

            fits = overflow <= 1e-6 and overlap_area <= 1e-6
            chosen_legend = legend_artist
            chosen_info = info_artist
            chosen_metrics = {
                "legend_bbox": legend_bbox,
                "info_bbox": info_bbox,
                "legend_area": (legend_bbox.width * legend_bbox.height) if legend_bbox is not None else 0.0,
                "info_area": (info_bbox.width * info_bbox.height) if info_bbox is not None else 0.0,
                "overlap_area": overlap_area,
                "overflow": overflow,
                "ncol": ncol,
                "legend_fs": legend_fs,
                "info_fs": info_fs,
            }

            if fits:
                break

            if legend_artist is not None:
                legend_artist.remove()
            if info_artist is not None:
                info_artist.remove()
            chosen_legend = None
            chosen_info = None

        if chosen_legend is not None or chosen_info is not None:
            break

    if chosen_metrics is not None:
        lb = chosen_metrics["legend_bbox"]
        ib = chosen_metrics["info_bbox"]
        if lb is not None:
            print(
                "layout_legend_bbox_axes: "
                f"x0={lb.x0:.3f} y0={lb.y0:.3f} x1={lb.x1:.3f} y1={lb.y1:.3f} "
                f"area={chosen_metrics['legend_area']:.4f}"
            )
        if ib is not None:
            print(
                "layout_info_bbox_axes: "
                f"x0={ib.x0:.3f} y0={ib.y0:.3f} x1={ib.x1:.3f} y1={ib.y1:.3f} "
                f"area={chosen_metrics['info_area']:.4f}"
            )
        print(
            "layout_overlap_axes_area: "
            f"{chosen_metrics['overlap_area']:.6f} "
            f"overflow={chosen_metrics['overflow']:.6f} "
            f"legend_ncol={chosen_metrics['ncol']} "
            f"legend_fs={chosen_metrics['legend_fs']} info_fs={chosen_metrics['info_fs']}"
        )


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

    lane_order = [
        lane
        for lane in ["GPU", "CPU", "I/O"]
        if any(iv.lane == lane for iv in intervals)
    ]
    if not lane_order:
        raise ValueError("No active lanes available to plot")
    y_map = {lane: idx for idx, lane in enumerate(lane_order)}

    colors = {
        "gpu_compute_slot1": "#4a3aa7",
        "gpu_compute_slot2": "#7d6fd1",
        "gpu_compute_sync": "#5b4bc0",
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
    for tick in ax.get_xticklabels():
        tick.set_rotation(30)
        tick.set_ha("right")

    kinds_present = {iv.kind for iv in intervals}
    legend_items: List[Tuple[str, str, Optional[str]]] = []
    if "io_read" in kinds_present:
        legend_items.append(("I/O read", "#2a78d6", None))
    if "io_write" in kinds_present:
        legend_items.append(("I/O write", "#1f4f99", None))
    if "cpu_prep" in kinds_present:
        legend_items.append(("CPU prep (odd blocks)", "#eda100", None))
        legend_items.append(("CPU prep (even blocks)", "#f4c65a", "////"))
    if "cpu_scatter" in kinds_present:
        legend_items.append(("CPU scatter", "#e34948", None))
    if "gpu_compute_slot1" in kinds_present:
        legend_items.append(("GPU compute async slot 1", "#4a3aa7", None))
    if "gpu_compute_slot2" in kinds_present:
        legend_items.append(("GPU compute async slot 2", "#7d6fd1", "//"))
    if "gpu_compute_sync" in kinds_present:
        legend_items.append(("GPU compute (synchronous fallback)", "#5b4bc0", None))

    handles = []
    for name, color, hatch in legend_items:
        handles.append(plt.Rectangle((0, 0), 1, 1, color=color, ec="black", lw=0.5, hatch=hatch))
    labels = [n for n, _, _ in legend_items]
    # Keep room on the right for the out-of-axes legend.
    fig.tight_layout(rect=(0.0, 0.0, 0.84, 1.0))
    _layout_right_panel(fig, ax, handles, labels, right_info_lines)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=170, bbox_inches="tight", pad_inches=0.12)
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
    cpu_stage_intervals = build_cpu_stage_intervals(run_events)
    thread_intervals = build_cpu_thread_intervals(run_events)
    phase_rows = build_phase_rows(run_events)

    gpu_enabled_hint = any(
        ev.category == "startup" and "GPU requested and enabled" in ev.message
        for ev in run_events
    )

    has_async_gpu_compute = any(iv.kind.startswith("gpu_compute_slot") for iv in intervals)
    if gpu_enabled_hint and not has_async_gpu_compute:
        intervals.extend(build_gpu_sync_intervals(run_events))
        intervals.sort(key=lambda x: (x.start, x.lane, x.kind))

    if not intervals and not thread_intervals:
        raise SystemExit("No intervals produced for selected run")

    all_start = [iv.start for iv in intervals] + [iv.start for iv in thread_intervals]
    all_end = [iv.end for iv in intervals] + [iv.end for iv in thread_intervals]
    t_start = min(all_start)
    t_end = max(all_end)
    plot_window_s = (t_end - t_start).total_seconds()
    total_wall_s = (
        (run_events[-1].ts - run_events[0].ts).total_seconds()
        if len(run_events) > 1
        else 0.0
    )

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
    gpu_sync_count = sum(1 for iv in intervals if iv.kind == "gpu_compute_sync")
    cpu_prep_count = sum(1 for iv in intervals if iv.kind == "cpu_prep")
    cpu_scatter_count = sum(1 for iv in intervals if iv.kind == "cpu_scatter")
    io_read_count = sum(1 for iv in intervals if iv.kind == "io_read")
    io_write_count = sum(1 for iv in intervals if iv.kind == "io_write")

    execution_context = "GPU run inferred" if gpu_enabled_hint else "CPU only run inferred"
    gpu_marker_status = "found" if gpu_enabled_hint else "not found"
    plot_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    cpu_only_thread_mode = (
        bool(thread_intervals)
        and gpu_slot1_count == 0
        and gpu_slot2_count == 0
        and gpu_sync_count == 0
    )
    view_mode = "CPU thread detail" if cpu_only_thread_mode else "Pipeline timeline"

    info_rows: List[Tuple[Optional[str], Optional[str]]] = [
        (None, "Run metadata"),
        ("Run log file", log_path.name),
        ("Run selector", args.run),
        ("Plot date", plot_date),
        ("Total wall time (s)", f"{total_wall_s:.3f}"),
        ("Execution context", execution_context),
        ("GPU startup marker", gpu_marker_status),
        ("View", view_mode),
        (None, ""),
        (None, "Event inventory"),
        ("Pipeline intervals", str(len(intervals))),
        ("CPU thread intervals", str(len(thread_intervals))),
        ("CPU stage intervals", str(len(cpu_stage_intervals))),
        (None, ""),
        (None, "Overlap metrics"),
        ("GPU-GPU overlap (s)", f"{gpu_gpu_overlap:.3f}"),
        ("CPU-GPU overlap (s)", f"{cpu_gpu_overlap:.3f}"),
        (None, ""),
        (None, "Category counts"),
        ("I/O read / write", f"{io_read_count} / {io_write_count}"),
        ("CPU prep / scatter", f"{cpu_prep_count} / {cpu_scatter_count}"),
        (
            "GPU async s1/s2/sync-fb",
            f"{gpu_slot1_count} / {gpu_slot2_count} / {gpu_sync_count}",
        ),
    ]

    if cpu_only_thread_mode:
        io_only = [iv for iv in intervals if iv.kind in {"io_read", "io_write"}]
        cpu_stage_only = [iv for iv in cpu_stage_intervals if iv.kind != "cpu_stage_compute"]
        thread_tids = sorted(set(iv.tid for iv in thread_intervals))
        info_rows.extend(
            [
                (None, ""),
                (None, "Thread layout"),
                ("Threads active", str(len(thread_tids))),
                ("Thread IDs", ", ".join(str(t) for t in thread_tids)),
            ]
        )
        right_info_lines = _format_info_block(info_rows)
        plot_cpu_thread_timeline(
            thread_intervals,
            io_only,
            cpu_stage_only,
            out_path,
            title=title,
            time_axis=args.time_axis,
            right_info_lines=right_info_lines,
        )
    else:
        right_info_lines = _format_info_block(info_rows)
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
    print(f"cpu_thread_intervals: {len(thread_intervals)}")
    print(f"cpu_only_thread_mode: {cpu_only_thread_mode}")
    print(f"view_mode: {view_mode}")
    print(f"execution_context: {execution_context}")
    print(f"total_wall_s: {total_wall_s:.3f}")
    print(f"plot_window_s: {plot_window_s:.3f}")
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
