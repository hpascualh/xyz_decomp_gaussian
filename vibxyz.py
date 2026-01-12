"""
VibXYZ
Cartesian decomposition of IR and Raman intensities from Gaussian frequency calculations.

This script provides a clean entry point while preserving the original
legacy implementation for full reproducibility.
"""

import sys
import os

# Ensure legacy module is discoverable
sys.path.append(os.path.join(os.path.dirname(__file__), "legacy"))

from gaus_vib_xyz_raman_hacky_v2 import main


if __name__ == "__main__":
    main()
