import os

from coltess import create_catalog, load_photometry_data, process_local_images_parallel, compute_periodogram, download_tess_image
import matplotlib.pyplot as plt

star_name = "lambda tau"

catalog_file = "gaia_catalog.csv"
fits_dir = "lambda tau"
output_file = "photometry.csv"

# ---------------------------------------------------
# Curl commands to download the FITS' from MAST
# ---------------------------------------------------
# These lines can be automatically downloaded
# with coltess.download.download_tess_sector_script.
# Please refer to lambda_tau_light_curve.py to
# see an example with the full pipeline.
lines = [
    "curl -C - -L -o tess2018338165938-s0005-1-4-0125-s_ffic.fits https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess2018338165938-s0005-1-4-0125-s_ffic.fits"
]  # For the example we only download and analyze one image.

# ---------------------------------------------------
# 1. Download the FITS files into a local folder
# ---------------------------------------------------
os.makedirs(fits_dir, exist_ok=True)

for line in lines:
    download_tess_image(line, fits_dir)

# ---------------------------------------------------
# 2. Create a catalog with the stars within a 10arcmin radius around star_name
# ---------------------------------------------------
star = create_catalog(star_name, radius_arcmin=10.0, output_file=catalog_file)

# ---------------------------------------------------
# 3. Analyze all local FITS images in parallel.
#    Results are appended to a single CSV file and
#    re-running resumes automatically.
# ---------------------------------------------------
process_local_images_parallel(
    fits_dir=fits_dir,
    catalog_file=catalog_file,
    star=star,
    output_file=output_file,
)

# ---------------------------------------------------
# 4. Load the light curve
# ---------------------------------------------------
times, fluxes, flux_errors = load_photometry_data(output_file, star)

# ---------------------------------------------------
# 5. Load light curve, compute periodogram, plot both
# ---------------------------------------------------
times, fluxes, flux_errors = load_photometry_data(output_file, star)
result = compute_periodogram(times, fluxes, flux_errors)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), tight_layout=True)

ax1.errorbar(times, fluxes, yerr=flux_errors, fmt="o", ms=2, capsize=0, elinewidth=0.5)
ax1.set(xlabel="Julian Date", ylabel="Flux (e⁻/s)", title=f"Light Curve: {star.name}")

ax2.plot(result["periods"], result["power"], lw=0.5)
ax2.axvline(
    result["primary_period"],
    color="r",
    ls="--",
    alpha=0.7,
    label=(
        f"P = {result['primary_period']:.4f} ± "
        f"{result['primary_period_uncertainty']:.4f} d\n"
        f"FAP = {result['primary_fap']:.2e}"
    ),
)
ax2.legend(fontsize=8)
ax2.set(xlabel="Period (d)", ylabel="Power", title="Lomb-Scargle Periodogram")

plt.show()
