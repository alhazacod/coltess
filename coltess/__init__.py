from .core import StarData
from .photometry import TessPhotometry
from .catalog import create_catalog, get_star, query_gaia_catalog
from .download import get_tess_sectors, download_tess_sector_script, download_tess_image, download_tess_images
from .analysis import load_photometry_data, compute_periodogram
from .parallel import process_images_parallel

__all__ = [
    "StarData",
    "TessPhotometry",
    "create_catalog", "get_star", "query_gaia_catalog",
    "get_tess_sectors", "download_tess_sector_script", "download_tess_image", "download_tess_images",
    "load_photometry_data", "compute_periodogram",
    "process_images_parallel"
]

__version__ = "0.1.0"
