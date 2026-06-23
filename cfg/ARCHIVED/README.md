# Archived Configuration Files

This directory contains 63 legacy configuration files that were specific to historical research projects using the GMRT (Giant Meterwave Radio Telescope) and CASA (Common Astronomy Software Applications) pipeline.

## Contents

These files are organized by the type of analysis they performed:

### RM-Synthesis Configs (16 files)
- `myfits_spec2rm_*.cfg` - Configurations for Faraday RM tomography
  - Files for calibration sources: 3C147, 3C286, 3C303, 3C345, 3C468.1, 3C48
  - CASA pipeline variants with different processing options
  - Specialized variants (ADHOC, RLCOR, mean-removed versions)

### Image Extraction Configs (3 files)
- `extract_image_from_cube*.cfg` - Extracting 2D images from 3D RM cubes at specific RM values

### Statistics Configs (4 files)
- `rmstat_*.cfg`, `rmpspec_*.cfg` - Statistical analysis of RM spectra

### General Purpose Configs (40 files)
- FITS image operations (align, combine, compare, compute statistics)
- SILO format conversion
- Frequency plane extraction

## Why Archived?

These configurations contain:
- **Hardcoded absolute paths** to `/media/GMRT_CASA/...` which won't work on other systems
- **Project-specific datasets** (source names, observing parameters, proprietary data paths)
- **Historical research parameters** that are specific to individual observations

## For Current Development

See `../example_myfits_spec2rm.cfg` for a **portable, well-documented example** of the correct configuration file format.

## If You Need One of These Files

If you need to adapt one of the archived configs for your own work:

1. Find the appropriate file (e.g., `myfits_spec2rm.cfg` for basic RM-synthesis)
2. Copy it to the parent `cfg/` directory
3. Update paths and filenames to match your data
4. Refer to `example_myfits_spec2rm.cfg` for parameter explanations

Example:
```bash
cp ARCHIVED/myfits_spec2rm.cfg ../my_project.cfg
# Then edit my_project.cfg with your paths and filenames
```

## File Statistics

Total archived files: 63
Total size: ~70 KB (original data files referenced are in various directories)

---

**Last updated:** 2026-06-23  
**Archive reason:** Repository cleanup and modernization
