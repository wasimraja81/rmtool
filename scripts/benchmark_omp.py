#!/home/wasim/venv/rmtool/bin/python3
"""Run OMP thread sweep benchmarks and generate a runtime/speedup plot.

Usage:
    python scripts/benchmark_omp.py \
            --config ../cfg/benchmark.cfg \
            --exe ../bin/rm_synthesis \
            --threads 1,2,3,4,5,6,7,8 \
            --outdir .
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import List, Tuple


def parse_threads(spec: str) -> List[int]:
    vals = []
    for part in spec.split(','):
        part = part.strip()
        if not part:
            continue
        vals.append(int(part))
    vals = sorted(set(vals))
    if not vals:
        raise ValueError("No valid thread counts provided")
    return vals


def run_one(exe: Path, cfg: Path, threads: int, cwd: Path, outbase: str) -> float:
    outputs = [
        cwd / f"{outbase}.AMP.RMCUBE.FITS",
        cwd / f"{outbase}.REAL.RMCUBE.FITS",
        cwd / f"{outbase}.PHA.RMCUBE.FITS",
        cwd / f"{outbase}.POLA.RMCUBE.FITS",
        cwd / f"{outbase}.IMAG.RMCUBE.FITS",
        cwd / f"{outbase}.MASK.CUBE.FITS",
        cwd / f"{outbase}.NVALID.MAP.FITS",
    ]
    for p in outputs:
        if p.exists():
            p.unlink()

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_PROC_BIND"] = env.get("OMP_PROC_BIND", "close")
    env["OMP_PLACES"] = env.get("OMP_PLACES", "cores")

    cmd = ["/usr/bin/time", "-f", "ELAPSED=%e", str(exe), str(cfg)]
    proc = subprocess.run(
        cmd,
        cwd=str(cwd),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    if proc.returncode != 0:
        raise RuntimeError(
            f"Benchmark failed for OMP_NUM_THREADS={threads}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )

    m = re.search(r"ELAPSED=([0-9]+(?:\.[0-9]+)?)", proc.stderr)
    if not m:
        raise RuntimeError(
            f"Could not parse elapsed time for OMP_NUM_THREADS={threads}\n"
            f"stderr:\n{proc.stderr}"
        )
    return float(m.group(1))


def write_csv(rows: List[Tuple[int, float, float]], csv_path: Path) -> None:
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["threads", "elapsed_sec", "speedup_vs_1"])
        for t, elapsed, speedup in rows:
            w.writerow([t, f"{elapsed:.3f}", f"{speedup:.3f}"])


def plot_results(rows: List[Tuple[int, float, float]], png_path: Path) -> None:
    import matplotlib.pyplot as plt

    threads = [r[0] for r in rows]
    elapsed = [r[1] for r in rows]
    speedup = [r[2] for r in rows]

    fig, ax1 = plt.subplots(figsize=(8.5, 5.0))
    ax2 = ax1.twinx()

    ln1 = ax1.plot(threads, elapsed, marker="o", color="#1565c0", label="Runtime (s)")
    ln2 = ax2.plot(threads, speedup, marker="s", color="#2e7d32", label="Speedup vs 1T")

    ax1.set_xlabel("OMP_NUM_THREADS")
    ax1.set_ylabel("Runtime (seconds)", color="#1565c0")
    ax2.set_ylabel("Speedup vs 1 thread", color="#2e7d32")
    ax1.grid(True, linestyle="--", alpha=0.4)

    lines = ln1 + ln2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="best")
    ax1.set_title("RMTool OMP Benchmark")

    fig.tight_layout()
    fig.savefig(png_path, dpi=140)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run OMP benchmark sweep for rm_synthesis")
    parser.add_argument("--config", required=True, help="Path to benchmark cfg (relative to cwd or absolute)")
    parser.add_argument("--exe", required=True, help="Path to rm_synthesis executable")
    parser.add_argument("--threads", default="1,2,3,4,5,6,7,8", help="Comma-separated thread counts")
    parser.add_argument("--outdir", default=".", help="Directory for CSV/plot outputs")
    parser.add_argument("--cwd", default=".", help="Working directory to run benchmark in")
    parser.add_argument("--outbase", default="MY_CASA_BENCHMARK_SUBIM", help="Output FITS basename from cfg outfile=")

    args = parser.parse_args()

    exe = Path(args.exe).expanduser().resolve()
    cfg = Path(args.config).expanduser().resolve()
    cwd = Path(args.cwd).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()

    if not exe.exists():
        raise FileNotFoundError(f"Executable not found: {exe}")
    if not cfg.exists():
        raise FileNotFoundError(f"Config not found: {cfg}")
    if not cwd.exists():
        raise FileNotFoundError(f"Working directory not found: {cwd}")

    outdir.mkdir(parents=True, exist_ok=True)

    threads = parse_threads(args.threads)

    print(f"[bench] exe={exe}")
    print(f"[bench] cfg={cfg}")
    print(f"[bench] cwd={cwd}")
    print(f"[bench] threads={threads}")
    print(f"[bench] OMP_PLACES=cores (physical cores only, no hyperthreads)")

    raw: List[Tuple[int, float]] = []
    for t in threads:
        elapsed = run_one(exe, cfg, t, cwd, args.outbase)
        raw.append((t, elapsed))
        print(f"[bench] threads={t:2d} elapsed={elapsed:.3f}s")

    base = next((e for t, e in raw if t == 1), raw[0][1])
    rows = [(t, e, base / e) for t, e in raw]

    csv_path = outdir / "benchmark_omp_results.csv"
    png_path = outdir / "benchmark_omp_plot.png"
    write_csv(rows, csv_path)

    if shutil.which("python3") is None:
        print("[bench] WARNING: python3 not found for plotting; wrote CSV only")
        return 0

    try:
        plot_results(rows, png_path)
        print(f"[bench] wrote {csv_path}")
        print(f"[bench] wrote {png_path}")
    except Exception as exc:  # pragma: no cover
        print(f"[bench] WARNING: plotting failed ({exc}); CSV is available at {csv_path}")
        print("[bench] Install matplotlib to enable plot output: python3 -m pip install matplotlib")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
