# RMTool Config Format

This project now uses a single KEY=VALUE config format.

**Scope:** this document covers `rm_synthesis`'s own config format only.
`reproject_cubes` and `convolve_cubes` (see the README's "Multi-Band
Preprocessing Toolchain" section) are separate standalone tools with their
own similar-but-separate KEY=VALUE parsers (deliberately not sharing this
one — see each source file's own top comment for why) — run either with
`--help` for its own full option list, or use `cfg/example_beamLog.txt`/
`.csv` for `convolve_cubes`' ASCII beam-log format specifically.

## Parser Source

The config parser is implemented in src/rm_synthesis_mod.f90:
- subroutine: read_cfg_keyval
- helper: split_key_value
- helper: flag_from_value

The main program src/rm_synthesis.f90 calls this parser.

## Strict Validation Rules

The parser enforces all of the following:
- Unknown keys are rejected.
- Duplicate keys are rejected.
- Numeric parsing errors are rejected.
- Required keys must be present.
- Range checks:
  - use_auto_rm_range must be 0 or 1
  - ofac must be >= 1
  - nrm must be >= 1

Optional key:
- output_mode: ap or ri
  - ap: write amplitude + angle cubes
  - ri: write real + imaginary cubes
- ap_angle_mode: phase or pol
  - phase: write PACUBE as phase angle arg(F)
  - pol: write PACUBE as polarization angle 0.5*arg(F)

Output filenames:
- AP + phase: OUTBASE.AMP.RMCUBE.FITS and OUTBASE.PHA.RMCUBE.FITS
- AP + pol:   OUTBASE.AMP.RMCUBE.FITS and OUTBASE.POLA.RMCUBE.FITS
- RI:         OUTBASE.REAL.RMCUBE.FITS and OUTBASE.IMAG.RMCUBE.FITS

## Required Keys

All keys below are required in strict mode:
- path
- infileQ
- infileU
- outfile
- remove_badchan
- badchan_file
- subim
- rem_mean
- remove_qu_bias
- resiQ
- slopeQ
- resiU
- slopeU
- ofac
- fac
- beg_rm
- end_rm
- nrm
- use_auto_rm_range

Optional keys (may be omitted):
- output_mode
- ap_angle_mode

Subimage configuration:
When subim=true, the following keys control the spatial/channel subset:
- subim_ra_blc: RA first pixel (>= 1)
- subim_ra_trc: RA last pixel (>= subim_ra_blc, or 0 for max)
- subim_ra_inc: RA step (>= 1)
- subim_dec_blc: Dec first pixel (>= 1)
- subim_dec_trc: Dec last pixel (>= subim_dec_blc, or 0 for max)
- subim_dec_inc: Dec step (>= 1)
- subim_chan_blc: Channel first pixel (>= 0, or 0 for first)
- subim_chan_trc: Channel last pixel (>= subim_chan_blc, or 0 for max)
- subim_chan_inc: Channel step (>= 1)

When subim=false, these keys are optional and default to full-range values.
(Legacy subim_parfile is still supported but superseded by KEY=VALUE pairs.)

RM sampling model:
- nrm_out is computed internally as nrm * ofac
- use_auto_rm_range=0:
  - user provides beg_rm, end_rm, nrm, ofac
  - code samples nrm_out points uniformly between beg_rm and end_rm
- use_auto_rm_range=1:
  - beg_rm/end_rm/nrm are derived from channel data
  - ofac controls oversampling, with nrm_out = nrm * ofac

Conditionally required keys:
- path_I and infileI are required only when remove_qu_bias is enabled

## Example

See:
- The full annotated, sectioned cfg reference in the top-level
  [`README.md`](../README.md) ("Configuration" section) -- every key the
  parser accepts, marked required/required-if/optional with its real
  default.
- `cfg/rmsynth-subim.cfg` for a runnable example.

(`cfg/myfits_spec2rm.cfg` and `cfg/example_myfits_spec2rm.cfg`, previously
linked here, no longer exist -- superseded by the above.)

## Adding New Config Variables

If you add a new KEY=VALUE variable in cfg files, you must also update src/rm_synthesis_mod.f90:
1. Add handling in read_cfg_keyval (select case block).
2. Add duplicate tracking for strict mode.
3. Add required-key validation if it is mandatory.
4. Add parse/type checks.

If this is not done, strict validation will fail on unknown keys.

## I-Cube Behavior

- I-cube is opened/read only when remove_qu_bias is enabled.
- If remove_qu_bias is disabled, I-cube keys may be omitted.

## Current Output Convention

For the AP case, the output type is controlled by ap_angle_mode:
- phase -> amplitude + phase angle cubes
- pol -> amplitude + polarization-angle cubes

The amplitude cube is always named with the AMP tag.
