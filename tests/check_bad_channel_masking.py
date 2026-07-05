#!/usr/bin/env python3
"""
Validation script for bad channel masking and NaN output handling.

Tests:
1. Pixels with one bad channel should have valid RM values (not masked out)
2. Fully-masked pixels should output NaN in RM cube
3. Mask cube should show correct per-channel masking
"""
import sys
import json
from pathlib import Path
import numpy as np
from astropy.io import fits

TESTS_DIR = Path(__file__).resolve().parent
DATA_DIR = TESTS_DIR / "data"
OUTPUT_DIR = TESTS_DIR / "output"

def load_truth():
    truth_path = DATA_DIR / "truth.json"
    with open(truth_path) as f:
        return json.load(f)

def check_bad_pixels_have_valid_rm(rm_cube, bad_pixels, fully_masked_pixels):
    """
    Verify that pixels with one bad channel still have valid (non-NaN) RM values.
    """
    print("Checking bad pixels have valid RM values...")
    errors = []
    
    for x, y in bad_pixels:
        # Check all RM bins for this pixel
        for irm in range(rm_cube.shape[0]):
            val = rm_cube[irm, y, x]
            if np.isnan(val):
                errors.append(f"  Bad pixel ({x},{y}) at RM bin {irm}: got NaN "
                            f"(should have valid RM value from good channels)")
    
    if errors:
        print("  [FAIL] Bad pixels have NaN in RM:")
        for e in errors[:5]:  # Show first 5
            print(e)
        return False
    else:
        print(f"  [OK] All {len(bad_pixels)} bad pixels have valid RM values")
        return True

def check_fully_masked_output_nan(rm_cube, fully_masked_pixels):
    """
    Verify that fully-masked pixels output NaN in RM cube.
    """
    print("Checking fully-masked pixels output NaN...")
    errors = []
    
    for x, y in fully_masked_pixels:
        # Check all RM bins for this pixel
        any_valid = False
        for irm in range(rm_cube.shape[0]):
            val = rm_cube[irm, y, x]
            if not np.isnan(val):
                any_valid = True
                errors.append(f"  Fully-masked pixel ({x},{y}) at RM bin {irm}: "
                            f"got {val} (should be NaN)")
        if not any_valid:
            print(f"  [OK] Fully-masked pixel ({x},{y}) outputs NaN across all RM bins")
    
    if errors:
        print("  [FAIL] Fully-masked pixels have valid values:")
        for e in errors[:5]:
            print(e)
        return False
    else:
        print(f"  [OK] All {len(fully_masked_pixels)} fully-masked pixels are NaN")
        return True

def check_mask_cube_correctness(mask_cube, bad_channel_idx, bad_pixels, fully_masked_pixels):
    """
    Verify that mask cube correctly shows per-channel masking.
    - At bad_channel_idx, pixels in bad_pixels should be masked (0)
    - At fully_masked_pixels, all channels should be masked (0)
    """
    print("Checking mask cube per-channel masking...")
    errors = []
    
    # Check bad pixels at bad channel
    for x, y in bad_pixels:
        val = mask_cube[bad_channel_idx, y, x]
        if val != 0:
            errors.append(f"  Bad pixel ({x},{y}) at chan {bad_channel_idx}: "
                        f"mask={val} (should be 0=masked)")
        # Check good channels at same pixel
        for chan in [0, 1, bad_channel_idx - 1, bad_channel_idx + 1]:
            if chan == bad_channel_idx or chan < 0 or chan >= mask_cube.shape[0]:
                continue
            val = mask_cube[chan, y, x]
            if val != 1:
                errors.append(f"  Bad pixel ({x},{y}) at good chan {chan}: "
                            f"mask={val} (should be 1=good)")
    
    # Check fully masked pixels at all channels
    for x, y in fully_masked_pixels:
        for chan in range(min(10, mask_cube.shape[0])):  # Check first 10 channels
            val = mask_cube[chan, y, x]
            if val != 0:
                errors.append(f"  Fully-masked pixel ({x},{y}) at chan {chan}: "
                            f"mask={val} (should be 0=masked)")
    
    if errors:
        print("  [FAIL] Mask cube has incorrect values:")
        for e in errors[:5]:
            print(e)
        return False
    else:
        print(f"  [OK] Mask cube correctly shows per-channel masking")
        return True

def main():
    # Load test metadata
    truth = load_truth()
    bad_chan_info = truth.get("bad_channel_test", {})
    
    if not bad_chan_info:
        print("[SKIP] No bad channel test data in truth.json")
        return 0
    
    bad_chan_idx = bad_chan_info["bad_channel_idx"]
    bad_pixels = bad_chan_info["bad_pixels"]
    fully_masked = bad_chan_info["fully_masked_pixels"]
    q_file = bad_chan_info["files_badchan"]["Q"]
    
    # Expected RM cube output file (from test suite)
    prefix = OUTPUT_DIR / "badchan"
    rm_cube_path = prefix.parent / f"{prefix.name}.AMP.RMCUBE.FITS"
    mask_cube_path = prefix.parent / f"{prefix.name}.MASK.CUBE.FITS"
    
    print(f"\nBad channel masking test")
    print(f"  Bad channel index: {bad_chan_idx}")
    print(f"  Bad pixels: {bad_pixels}")
    print(f"  Fully masked: {fully_masked}")
    print()
    
    # Check if RM cube exists
    if not rm_cube_path.exists():
        print(f"[SKIP] RM cube not found at {rm_cube_path}")
        print(f"       (Make sure test has been run with BADCHAN config)")
        return 0
    
    # Load RM and mask cubes
    # FITS axis order: (RM, Dec, RA) [no Stokes axis in output]
    with fits.open(rm_cube_path) as hdul:
        rm_data = hdul[0].data
        # rm_data is (RM, Dec, RA) - use directly
        rm_cube = rm_data  # (RM, Dec, RA)
    
    all_pass = True
    
    # Test 1: Bad pixels have valid RM
    if not check_bad_pixels_have_valid_rm(rm_cube, bad_pixels, fully_masked):
        all_pass = False
    print()
    
    # Test 2: Fully masked pixels output NaN
    if not check_fully_masked_output_nan(rm_cube, fully_masked):
        all_pass = False
    print()
    
    # Test 3: Mask cube correctness (if it exists)
    if mask_cube_path.exists():
        with fits.open(mask_cube_path) as hdul:
            # Mask cube is (Freq, Dec, RA)
            mask_data = hdul[0].data
            mask_cube = mask_data  # (Freq, Dec, RA)
        
        if not check_mask_cube_correctness(mask_cube, bad_chan_idx, bad_pixels, fully_masked):
            all_pass = False
    else:
        print("[SKIP] Mask cube not found (write_mask_output may be disabled)")
    print()
    
    if all_pass:
        print("[PASS] Bad channel masking test passed")
        return 0
    else:
        print("[FAIL] Bad channel masking test failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())
