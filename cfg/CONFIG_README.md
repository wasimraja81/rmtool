# RMTool Config Format

This project now uses a single KEY=VALUE config format.

**Scope:** this document covers `rm_synthesis`'s own config *parser
mechanics* only — the strict validation rules the parser enforces, not
a key-by-key list of what each key does (see
[docs/user/APP_REFERENCE.md](../docs/user/APP_REFERENCE.md) for that —
the single source of truth for every key, every default, and every
output file across all 5 tools).

`reproject_cubes`, `convolve_cubes`, and `match_cubes` are separate
standalone tools with their own similar-but-separate KEY=VALUE parsers
(deliberately not sharing this one — see each source file's own top
comment for why) — run any of them with `--help` for its own full
option list, or use `cfg/example_beamLog.txt`/`.csv` for
`convolve_cubes`'/`match_cubes`' ASCII beam-log format specifically.

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

## Example

See:
- [`cfg/rmsynth.cfg`](rmsynth.cfg) for a full annotated, sectioned cfg
  template -- every key the parser accepts, marked required/required-if/
  optional with its real default. Also documented in
  [`docs/user/APP_REFERENCE.md`](../docs/user/APP_REFERENCE.md)'s
  rm_synthesis section.
- `cfg/rmsynth-subim.cfg` for a runnable example.

(`cfg/myfits_spec2rm.cfg` and `cfg/example_myfits_spec2rm.cfg`, previously
linked here, no longer exist -- superseded by the above.)
