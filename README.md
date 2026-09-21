# rmtool

The rmtool package is intended for performing Faraday tomography (or
rotation-measure (RM) synthesis) on radio spectro-polarimetric data.
In addition to processing data from a single contiguous band, rmtool
provides the necessary applications to appropriately combine
heterogeneous multi-band data from different telescopes/bands with
different native projections and angular resolutions. Deconvolution
becomes critical for data sparsely sampled in the spectral domain; the
RM-CLEAN algorithm provided handles this appropriately.

One key motivation for this package is to address Faraday tomography
of the big data cubes that modern observatories are generating at
present. Current packages are unable to handle this data, either
because of the narrow scope they were designed for (i.e. aimed at
small-to-moderate sized data cubes), or because of limited HPC
strategies built into these algorithms. This package addresses both —
every application in it can process large data cubes optimally on
both PCs with small core counts and RAM, as well as HPC-grade nodes
with hundreds of cores and TBs worth of shared memory.

**New here?** → **[QUICKSTART.md](QUICKSTART.md)** builds the tools
and runs a complete match → RM synthesis → RM-CLEAN pipeline end to
end, in a couple of commands.

## Demo

RM-CLEAN's restored output from a WALLABY/EMU multi-band cube reveals line of sight components of the magnetic fields that are in opposite direction at the two edges of this supernova remnant PKS1209-52.

![RM-CLEAN Faraday-depth channel-map animation](docs/user/images/rmclean_channelmap_demo.gif)

**Highlights**

- **Handles cubes larger than RAM.** Data is processed in tiles sized
  to a user-set fraction of system memory (`mem_frac_ram`) — the same
  config runs on a laptop or an HPC node without need for changes.
- **Multi-band pre-processing built in.** Provides tools for matching
  sky-grid projections between multi-band data, and angular resolution
  across multiple bands as well as across the spectral axis of a
  single band. These two stages can be chained together to avoid
  writing intermediate files, saving huge amounts of disk space for
  large cubes.
- **RM-CLEAN included**, sharing the same tiled, memory-budgeted I/O
  engine as RM synthesis. Correctly handles "gaps" in data due to RFI,
  and from multiple widely separated bands.
- **Parallel, overlapped I/O** for HPC filesystems: multiple concurrent
  read/write channels, with a tile's write overlapped against the next
  tile's read and compute.
- **Optional GPU offload** for the core RM synthesis transform
  (work in progress — see [GPU Support](#gpu-support-work-in-progress)).
- **Configuration files.** Flexible run configuration using plain-text
  `KEY=VALUE` pair files to control memory budget, threading, I/O
  parallelism, GPU use, etc., in addition to RM-synthesis-specific
  parameters.

## What's Included

| Tool | Binary | Does |
|---|---|---|
| RM synthesis | `rm_synthesis` | The core λ²→Faraday-depth transform |
| Reproject | `reproject_cubes` | Aligns cubes onto one common sky grid (WCS reprojection) |
| Convolve | `convolve_cubes` | Convolves cubes to one common angular resolution |
| Match | `match_cubes` | Reproject + convolve chained through memory, no intermediate file |
| RM-CLEAN | `rmclean_cubes` | Högbom-style CLEAN deconvolution of RMSF sidelobes |
| Pipeline | `scripts/run_pipeline.sh` | Chains match → rm_synthesis → rmclean from one config |

## Choosing the Right Tool for Your Data

| Your situation | Run this | Example |
|---|---|---|
| Single band, channels already at consistent resolution, and (if multiple bands) same sky grid | `rm_synthesis` | [EXAMPLES.md §1](docs/user/EXAMPLES.md#1-single-band-quickstart) |
| Different sky grid — one band's pointing vs. another's | `reproject_cubes`, then `rm_synthesis` | [EXAMPLES.md §2b](docs/user/EXAMPLES.md#2b-sky-grid-mismatched-resolution-already-matched--reproject_cubes-only) |
| Different angular resolution — across bands, or across channels of a single band (native beam varies with frequency) | `convolve_cubes`, then `rm_synthesis` | [EXAMPLES.md §2c](docs/user/EXAMPLES.md#2c-resolution-mismatched-sky-grid-already-matched--convolve_cubes-only) |
| Multiple bands, both grid and resolution differ | `match_cubes`, then `rm_synthesis` | [EXAMPLES.md §2d](docs/user/EXAMPLES.md#2d-both-mismatched--match_cubes-chained-through-memory) |
| Need CLEAN deconvolution of a dirty RM cube | `rmclean_cubes` | [EXAMPLES.md §3](docs/user/EXAMPLES.md#3-choosing-rm-clean-stopping-criteria) |
| Want one command for match → rmsynth → rmclean | `scripts/run_pipeline.sh` | [TUTORIAL.md §6](docs/user/TUTORIAL.md#6-running-the-full-pipeline-in-one-command) |

Not sure which case applies? [docs/user/EXAMPLES.md §2](docs/user/EXAMPLES.md#2-multi-band-which-preprocessing-do-i-need) walks through how to tell.

## Quick Start

See [QUICKSTART.md](QUICKSTART.md) — build the tools and run a
complete pipeline end to end in a couple of commands, including
fetching your own data from CASDA. For build requirements and every
build variant, see [BUILD.md](BUILD.md); for the same run explained
step by step, see [docs/user/TUTORIAL.md](docs/user/TUTORIAL.md).

## Configuration

Every tool reads a plain-text `KEY=VALUE` config file (`cfg/`), one key
per line, `#` for comments. The parser is strict: unknown keys,
duplicate keys, and unparsable values are all rejected outright, so a
config that loads at all is already validated in that sense.
[cfg/rmsynth.cfg](cfg/rmsynth.cfg) is a fully annotated template; the
complete key-by-key reference for every tool is in
[docs/user/APP_REFERENCE.md](docs/user/APP_REFERENCE.md).

## GPU Support (Work in Progress)

`rm_synthesis` can optionally offload its core transform to a GPU via
OpenMP target offload (`make GPU=1`, `use_gpu=y`). This works, but is
early-stage: validated on one GPU model so far
(`NVIDIA GeForce RTX 3050`), not yet the recommended default, and not
supported by the other four tools. See
[EXAMPLES.md §6](docs/user/EXAMPLES.md#6-gpu-vs-cpu-which-build-should-i-use)
for guidance on choosing CPU vs. GPU.

## Documentation

| Doc | Covers |
|---|---|
| [QUICKSTART.md](QUICKSTART.md) | Build, then run the full pipeline end to end in two commands |
| [BUILD.md](BUILD.md) | Build system details, requirements, release tagging policy |
| [docs/user/TUTORIAL.md](docs/user/TUTORIAL.md) | Step-by-step: build → synthetic data → rm_synthesis → rmclean → inspect output |
| [docs/user/EXAMPLES.md](docs/user/EXAMPLES.md) | Recipes by scenario: single/multi-band, CLEAN stopping criteria, memory/IO tuning, GPU vs. CPU, subimage extraction |
| [docs/user/APP_REFERENCE.md](docs/user/APP_REFERENCE.md) | Full parameter reference for every tool: every key, every default, output files |
| [cfg/CONFIG_README.md](cfg/CONFIG_README.md) | `rm_synthesis` config parser rules |
| [RELEASE_NOTES.md](RELEASE_NOTES.md) | What's in the current release |
| [CHANGELOG.md](CHANGELOG.md) | Version history |

Internal architecture/engineering documentation lives under `docs/dev/`.

## Development

- **Branches:** `main` (stable releases), `develop` (active development).
- **Release tags:** `R`-prefixed `MAJOR.MINOR` on `main` (e.g. `R1.0`).
  See [BUILD.md](BUILD.md#release-tagging-policy) for why the `R`
  prefix.

## License

See [LICENSE](LICENSE). For questions or contributions, please open an issue.

---

*This project was developed with assistance from LLM tools (code
refactorisation, debugging, documentation, and testing support).*
