from coltess.photometry import TessPhotometry
from coltess.download import download_tess_image
from coltess.catalog import create_catalog

star_name = "lambda tau"

catalog_file = "gaia_catalog.csv"
fits_dir = "fits"
csv_dir = "csv"

# ---------------------------------------------------
# Curl commands to download the FITS' from MAST
# ---------------------------------------------------
# These lines can be automatically downloaded
# with coltess.download.download_tess_secto_script.
# Please refer to lambda_tau_light_curve.py to 
# see an example with the full pipeline.
lines = [
        "curl -C - -L -o tess2018338165938-s0005-1-4-0125-s_ffic.fits https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess2018338165938-s0005-1-4-0125-s_ffic.fits"
        ]   # For the example we only download and analyze one image.

# ---------------------------------------------------
# Download and analyze the FITS from the curl commands
# ---------------------------------------------------
for line in lines:
    # 1. Download the image from the curl command.
    fits_file = download_tess_image(line, fits_dir)
    # 2. Initialize the photometry class.
    photometry = TessPhotometry()
    # 3. Create a catalog with the stars within a 10arcmin radius around star_name
    star = create_catalog(star_name, radius_arcmin=10.0, output_file=catalog_file)
    # 4. Check if the star is in the image, do the photometry and save the result inside csv_dir
    success = photometry.process_image(fits_file, catalog_file, star, output_dir = csv_dir)
    print(success)
    if success:
        print(f"Star found at {line}")

