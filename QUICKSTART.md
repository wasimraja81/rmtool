# RM-Synthesis Quick Start Guide

Build the tools and run a complete match → RM synthesis → RM-CLEAN
pipeline end to end, with a real result to look at -- no external data
or accounts required to get there.

## Contents
- [Where things go](#where-things-go)
- [Prerequisites](#prerequisites)
- [1. Build](#1-build)
- [2. Validate your build](#2-validate-your-build)
- [3. Run it end to end](#3-run-it-end-to-end)
- [4. Where to go next](#4-where-to-go-next)

---

## Where things go

Set `$RMTOOL_HOME` first, then clone into it -- one disciplined way of
doing this, not an optional convenience. Every command below, and in
every other doc in this project, assumes you're sitting inside it.

```bash
export RMTOOL_HOME=~/rmtool   # any path you like -- set once per shell session
git clone https://github.com/wasimraja81/rmtool.git "$RMTOOL_HOME"
cd "$RMTOOL_HOME"
```

> [!TIP]
> Opening a new shell later? Export `$RMTOOL_HOME` again before running
> anything here -- it isn't remembered between sessions unless you add
> the `export` line to your own shell profile (`~/.bashrc` or similar).

- `./build.sh` and every `make` command are run from `$RMTOOL_HOME`
  itself.
- Every binary lands in `$RMTOOL_HOME/bin/` (`match_cubes`,
  `reproject_cubes`, `convolve_cubes`, `rmclean_cubes`, and the
  `rm_synthesis_release_*` variants) -- always there, regardless of
  what directory you happened to run `build.sh`/`make` from.
- The commands below run from `$RMTOOL_HOME` using relative paths
  (`bin/...`, `cfg/...`, `tests/...`).

---

## Prerequisites

| Package | Notes |
|---|---|
| Fortran compiler, CFITSIO, GPU compiler, Starlink AST + FFTW3 | See [BUILD.md](BUILD.md#requirements) for exact install commands per platform and which tool needs which dependency (`rmclean_cubes` needs only FFTW3 of that last group, not Starlink AST). |
| `python3`, `python3-venv`, `python3-pip`, plus this project's own Python dependencies (astropy, astroquery, matplotlib, numpy, pyvo, requests) | Needed for the commands below, the test suite, and the plotting tools -- not for the compiled binaries themselves. The first three are system packages (`python3 -m venv` fails without `python3-venv` even when `python3` itself is installed); once those are present, `./build.sh` sets everything else up for you automatically. |
| ffmpeg (diagnostic animation only) | via your package manager -- used only by `scripts/animate_fits_cube.py` |

```bash
# Minimal Ubuntu install
sudo apt-get install gfortran libcfitsio-dev make python3 python3-venv python3-pip
```

---

## 1. Build

```bash
# from $RMTOOL_HOME
./build.sh
```

Checks every dependency above by name, builds all five tools, and sets
up this project's own Python environment automatically (best-effort --
a missing or failed venv never blocks the binaries themselves; see
[Prerequisites](#prerequisites)).

> [!NOTE]
> **Check:** you should see `Build complete` followed by a list of
> binaries under `bin/`, with no `[MISSING]` line above it.

### Diagnosing build failures

`./build.sh` names exactly what's missing and how to install it,
rather than failing deep inside a compile. If a dependency check
passed but the build itself still fails -- a compilation error, an AST
link error, or the Python venv step erroring out -- see
[BUILD.md's Troubleshooting section](BUILD.md#troubleshooting), which
covers the specific failure modes seen so far, each with its exact
fix.

---

## 2. Validate your build

```bash
# from $RMTOOL_HOME
bash tests/run_tests.sh
```

Optional, but strongly recommended before you rely on this build for
real work: rebuilds every variant and runs the full regression suite
(~180 checks: RM peak recovery, OMP/GPU bit-identical comparisons, I/O
parallelism, RM-CLEAN correctness, and more) -- the most thorough
confirmation that this checkout, on this machine, actually works.

> [!IMPORTANT]
> Takes a few minutes, not seconds -- it rebuilds everything from
> scratch first. Worth the wait: this is the same suite that gates
> every release.

> [!NOTE]
> **Check:** the summary ends with `RESULT: ALL PASSED`. `Skip` counts
> vary by platform -- GPU sections skip cleanly when no GPU-capable
> binary could be built. That's expected, not a failure.

---

## 3. Run it end to end

The fastest way to see the whole pipeline work, with no external data
or accounts needed:

```bash
# from $RMTOOL_HOME
~/venv/rmtool/bin/python3 tests/make_test_cubes.py
scripts/run_pipeline.sh cfg/pipeline-e2e-smalltest.cfg
```

The first command writes a small synthetic Stokes Q/U cube pair (two
point sources at known RMs) into `tests/data/`. The second chains
match → `rm_synthesis` → `rmclean_cubes` against it in one call.

> [!NOTE]
> **Check:** the last line printed is `[pipeline] Pipeline finished.`,
> with no `ERROR` line above it. Output lands in
> `tests/output/pipeline_e2e_smalltest/` --
> `smalltest_cleaned.RESTORED.AMP.RMCUBE.FITS` is the final CLEANed
> result.

### Want to try it on real data?

`scripts/casda_fetch.py` checks the CSIRO ASKAP Science Data Archive
(CASDA) for calibrated POSSUM Stokes Q/U cubes around a sky position
and fetches them, generating a ready-to-run pipeline config in the
same style as above.

> [!TIP]
> Needs a free CASDA/OPAL account (self-register at
> `opal.atnf.csiro.au`). Don't have one yet? The synthetic run above
> already gave you a working pipeline -- come back to this once you're
> registered.

**Check what's there** (the default mode; downloads nothing, safe to
run any time):
```bash
# from $RMTOOL_HOME
~/venv/rmtool/bin/python3 scripts/casda_fetch.py --target dancingghosts --username you@example.org
```

**Fetch, and generate a pipeline config:**
```bash
# from $RMTOOL_HOME
~/venv/rmtool/bin/python3 scripts/casda_fetch.py --target dancingghosts --username you@example.org \
    --run-mode=auto --outdir ./dancingghosts_data
```

**Run it**, the same way as the synthetic example above:
```bash
# from $RMTOOL_HOME
scripts/run_pipeline.sh ./dancingghosts_data/dancingghosts_pipeline.cfg
```

`--target` also accepts any name resolvable via SIMBAD/NED (e.g.
`--target "Kes 27"`), or `--ra`/`--dec`/`--radius` for a raw position
not resolvable by name. A starter `rmclean_cubes` cfg is generated
alongside the pipeline one too, but left out of `stages=` by default --
CLEAN's stopping criteria are a choice to make after looking at the
dirty cube, not a safe default (see
[EXAMPLES.md §3](docs/user/EXAMPLES.md#3-choosing-rm-clean-stopping-criteria)
for choosing them, then `--chain-rmclean` at fetch time to opt in).
Full flag reference, every `--run-mode`, and re-run behavior via
`--if-exists`: `scripts/casda_fetch.py --help`.

---

## 4. Where to go next

| Doc | For |
|---|---|
| [docs/user/TUTORIAL.md](docs/user/TUTORIAL.md) | The same synthetic run above, explained step by step -- including inspecting the output pixel by pixel |
| [docs/user/EXAMPLES.md](docs/user/EXAMPLES.md) | Recipes by scenario: multi-band grid/resolution mismatches, RM-CLEAN stopping criteria, memory/IO tuning, GPU vs. CPU, subimage extraction |
| [docs/user/APP_REFERENCE.md](docs/user/APP_REFERENCE.md) | Every parameter, every tool, fully explained -- including `reproject_cubes`/`convolve_cubes`/`match_cubes`/`rmclean_cubes` run standalone, and the OMP/GPU environment variables |
| [BUILD.md](BUILD.md) | Every build variant, troubleshooting, performance tuning |

Swim-lane I/O visualization (`scripts/plot_tile_async_swimlane.py`) and
internal architecture/design documentation live under `docs/dev/` and
in the tools' own `--help`.
