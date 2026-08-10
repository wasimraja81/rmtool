# Tutorial: From a Fresh Checkout to a Cleaned RM Cube

A step-by-step walkthrough for a brand-new user: build the tools,
generate a small synthetic dataset with known answers, run RM synthesis
and RM-CLEAN on it, and inspect the result. Every command below has
been run and verified — this is not a hypothetical walkthrough.

Run everything from the repository root.

## 1. Build

```bash
make OMP=1 GPU=0
make rmclean_cubes
```

This produces `bin/rm_synthesis_release_cpu_omp` (also symlinked as
`bin/rm_synthesis`) and `bin/rmclean_cubes`. See
[QUICKSTART.md](../../QUICKSTART.md) if `make` fails — it covers
prerequisites and every build variant.

## 2. Generate a small synthetic dataset

```bash
python3 tests/make_test_cubes.py
```

This writes a 32×32-pixel, 200-channel synthetic Stokes Q/U cube pair
to `tests/data/TEST.Q.FITSCUBE`/`TEST.U.FITSCUBE` — no real telescope
data needed. Two point sources are injected at known positions and
known rotation measures, recorded in `tests/data/truth.json`:

| Source | Pixel (0-indexed x, y) | RM (rad/m²) |
|---|---|---|
| `src_A` | (12, 10) | −5.0 |
| `src_B` | (22, 20) | +22.0 |

This is the exact same fixture `tests/run_tests.sh` uses throughout the
regression suite, so everything below is also continuously verified by
that suite, not just this document.

## 3. Run RM synthesis

Create a config file (`tutorial-rmsynth.cfg`) with the content below,
or copy `cfg/rmsynth-e2e-smalltest.cfg` and edit its `path=` to a
relative path (`tests/data/`) if you're not running from this exact
checkout location:

```cfg
path                = tests/data/
infileQ             = TEST.Q.FITSCUBE
infileU             = TEST.U.FITSCUBE
outfile             = tutorial_out

remove_badchan      = n
global_badchan_file = /dev/null
subim               = n
rem_mean            = 0
remove_qu_bias      = n
resiQ               = 0.0
slopeQ              = 0.0
resiU               = 0.0
slopeU              = 0.0
ofac                = 1
fac                 = 3.14159265358979
use_auto_rm_range   = 0
beg_rm              = -50.0
end_rm              = 50.0
nrm                 = 201
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = n
```

Run it:

```bash
bin/rm_synthesis_release_cpu_omp tutorial-rmsynth.cfg
```

This writes `tutorial_out.AMP.RMCUBE.FITS`, `tutorial_out.PHA.RMCUBE.FITS`,
`tutorial_out.MASK.CUBE.FITS`, and `tutorial_out.NVALID.MAP.FITS` — the
dirty (pre-CLEAN) RM cube. Every key here is explained in full in
[docs/user/APP_REFERENCE.md](APP_REFERENCE.md#1-rm_synthesis); the same
annotated example lives at `cfg/rmsynth-e2e-smalltest.cfg`.

**Check the known sources were recovered:**

```bash
python3 tests/check_rm_peak.py tutorial_out.AMP.RMCUBE.FITS tests/data/truth.json
```

Expected output:

```
[OK] src_A: expected RM=-5.0, found RM=-5.00 (err=0.00, tol=1.00)
[OK] src_B: expected RM=+22.0, found RM=+22.00 (err=0.00, tol=1.00)
```

## 4. Run RM-CLEAN

```bash
bin/rmclean_cubes ampfile=tutorial_out.AMP.RMCUBE.FITS \
  phafile=tutorial_out.PHA.RMCUBE.FITS \
  maskfile=tutorial_out.MASK.CUBE.FITS \
  outfile=tutorial_cleaned \
  abs_flux_floor=0.05 auto_nsigma=1.0 niter=200 gain=0.1
```

This writes 6 output cubes:
`tutorial_cleaned.CLEAN.AMP/PHA.RMCUBE.FITS` (the CLEAN components),
`tutorial_cleaned.RESID.AMP/PHA.RMCUBE.FITS` (what's left after
subtracting them), and `tutorial_cleaned.RESTORED.AMP/PHA.RMCUBE.FITS`
(components reconvolved with a restoring beam, plus the residual — the
final, publication-style result). CLEAN stops per-pixel on whichever of
`abs_flux_floor`/`auto_nsigma`/`niter` fires first; every key is
explained in full in
[docs/user/APP_REFERENCE.md](APP_REFERENCE.md#5-rmclean_cubes).

For this small, low-noise synthetic cube, expect a stop-reason summary
close to:

```
CLEAN stop-reason summary:
  hit niter cap:        0 (  0.00%)
  stopped at abs_flux:  1022 ( 99.80%)
  stopped at auto_nsigma:2 (  0.20%)
  n_iter_used: mean=1.01 min=1 max=8 (niter cap=200)
```

Nearly every pixel stops in 1-2 iterations because most of the image is
empty sky (no source at all — CLEAN correctly finds nothing to remove
and stops immediately at `abs_flux_floor`). Real data with a genuine
noise floor behaves differently — see `cfg/rmclean-example.cfg` for a
fully worked-through, annotated real-world configuration.

**Confirm both sources are still correctly recovered after CLEAN:**

```bash
python3 tests/check_rm_peak.py tutorial_cleaned.RESTORED.AMP.RMCUBE.FITS tests/data/truth.json
```

## 5. Inspect a single pixel directly

A short Python/astropy snippet to look at one source's own restored
spectrum directly (`src_A`, pixel `(12, 10)` per `truth.json` — note
these are already 0-indexed numpy array positions, not 1-indexed FITS
pixel numbers):

```python
from astropy.io import fits
import numpy as np

amp = fits.getdata('tutorial_cleaned.RESTORED.AMP.RMCUBE.FITS')
hdr = fits.getheader('tutorial_cleaned.RESTORED.AMP.RMCUBE.FITS')

spectrum = amp[:, 10, 12]          # data[:, y, x]
peak_idx = np.argmax(spectrum)
rm = hdr['CRVAL3'] + peak_idx * hdr['CDELT3']
print(f"Peak amplitude {spectrum[peak_idx]:.4f} at RM={rm:.2f} rad/m^2")
# Peak amplitude 0.9995 at RM=-5.00 rad/m^2
```

## 6. Running the full pipeline in one command

Steps 3-4 above can also be run as a single chained call via
`scripts/run_pipeline.sh`, which additionally supports a `match_cubes`
pre-processing stage for data that doesn't already share a common sky
grid and/or resolution — whether that's bands that don't match each
other, or a single band whose own channels don't share one resolution:

```bash
scripts/run_pipeline.sh --help          # full option list, read this first
scripts/run_pipeline.sh cfg/pipeline-example.cfg
```

`cfg/pipeline-example.cfg` is a fully annotated template — edit its
paths for your own data, or see `cfg/pipeline-e2e-smalltest.cfg` for a
complete working example using this same synthetic fixture. See
`scripts/run_pipeline.sh`'s own header comment (or `--help`) for the
full chaining/provenance behaviour.

## Where to go next

- [docs/user/EXAMPLES.md](EXAMPLES.md) — recipes for real scenarios: multi-
  band with matched/mismatched grid or resolution, choosing RM-CLEAN
  stopping criteria, memory/IO tuning, GPU vs. CPU, subimage extraction.
- [docs/user/APP_REFERENCE.md](APP_REFERENCE.md) — every parameter,
  every tool, fully explained.
- [QUICKSTART.md](../../QUICKSTART.md) — build variants, GPU builds, the
  full validation test suite, running on real data.
- [README.md](../../README.md) — top-level overview: what each tool
  does, and which one to use for your data.
- `cfg/rmsynth.cfg`, `cfg/rmclean-example.cfg` — fully annotated
  templates to copy and adapt for real data.
