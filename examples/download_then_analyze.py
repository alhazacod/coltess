from pathlib import Path

from coltess.photometry import TessPhotometry
from coltess.download import download_tess_image
from coltess.catalog import create_catalog

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
    star = create_catalog(star_name, radius_arcmin=10.0, output_file=catalog_file)
    success = photometry.process_image(fits_file, catalog_file, star, output_dir = csv_dir)
    print(success)
    if success:
        print(f"Star found at {line}")

