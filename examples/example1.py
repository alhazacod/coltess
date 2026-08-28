from coltess import create_catalog, get_tess_sectors, download_tess_sector_script
from coltess import process_images_parallel, load_photometry_data
import matplotlib.pyplot as plt

# 1. Create catalog and get star info
star = create_catalog(
    "lambda tau",
    radius_arcmin=10.0,
    output_file="catalog.csv",
)

# 2. Find available TESS sectors
sectors = get_tess_sectors(star)
sector = int(sectors["sector"][0])

# 3. Download sector script
script_path = download_tess_sector_script(sector)

# 4. Process images in parallel, appending results to a single CSV
process_images_parallel(
    script_file=script_path,
    catalog_file="catalog.csv",
    output_file="photometry_results.csv",
    star=star,
    # keep_images_dir="lambda tau" # Optional field if need to save the images with the star in it
)

# 5. Load and plot light curve
times, fluxes, flux_errors = load_photometry_data("photometry_results.csv", star)

plt.scatter(times, fluxes)
plt.xlabel("Julian Date")
plt.ylabel("Flux (e⁻/s)")
plt.title(f"Light Curve: {star.name}")
plt.show()
