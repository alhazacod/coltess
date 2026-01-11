from pathlib import Path

from coltess.photometry import TessPhotometry
from coltess.catalog import create_catalog, get_star_coordinates
from coltess.download import download_tess_image
from coltess.analysis import load_photometry_data

star_name = "lambda tau"

base_dir = Path("run_single")
base_dir.mkdir(exist_ok=True)

catalog_file = base_dir / "gaia_catalog.csv"

fits_dir = base_dir / "fits"
fits_dir.mkdir(exist_ok=True)

csv_dir = base_dir / "csv"
csv_dir.mkdir(exist_ok=True)

lines = [
        "curl -C - -L -o tess2018338165938-s0005-1-4-0125-s_ffic.fits https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess2018338165938-s0005-1-4-0125-s_ffic.fits"
        ]

for line in lines:
    fits_file = download_tess_image(line, fits_dir)
    #print(f"fits for {line} saved")
    photometry = TessPhotometry()
    ra, dec, gaia_id = get_star_coordinates(star_name)
    success = photometry.process_image(fits_file, catalog_file, ra, dec, output_dir = csv_dir)
    print(success)
    if success:
        print(f"Star found at {line}")

