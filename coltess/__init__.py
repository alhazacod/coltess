"""
COLTESS: Curves of Light from TESS.

A lightweight Python package for extracting light curves from TESS
(Transiting Exoplanet Survey Satellite) Full Frame Images (FFIs).

Main entry points:

- :func:`create_catalog` — build a Gaia DR3 catalog around a target.
- :func:`process_images_parallel` — download and analyze TESS FFIs in parallel.
- :func:`process_local_images_parallel` — analyze local FITS files in parallel.
- :func:`analyze_image` — analyze a single local FITS image.
- :func:`load_photometry_data` — load the light curve from the output CSV.
- :func:`compute_periodogram` — Lomb-Scargle periodogram with FAP.
"""

from .core import StarData
from .photometry import TessPhotometry, analyze_image
from .catalog import create_catalog, get_star, query_gaia_catalog
from .download import (
    get_tess_sectors,
    download_tess_sector_script,
    download_tess_image,
    download_tess_images,
)
from .analysis import load_photometry_data, compute_periodogram
from .parallel import process_images_parallel, process_local_images_parallel

__all__ = [
    "StarData",
    "TessPhotometry",
    "analyze_image",
    "create_catalog",
    "get_star",
    "query_gaia_catalog",
    "get_tess_sectors",
    "download_tess_sector_script",
    "download_tess_image",
    "download_tess_images",
    "load_photometry_data",
    "compute_periodogram",
    "process_images_parallel",
    "process_local_images_parallel",
]

__version__ = "0.1.0"
