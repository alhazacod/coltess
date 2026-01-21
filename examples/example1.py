from coltess import create_catalog, get_tess_sectors, download_tess_sector_script
from coltess import process_images_parallel, load_photometry_data
import matplotlib.pyplot as plt

# 1. Create catalog and get star info
star = create_catalog("lambda tau", radius_arcmin=10.0, output_file="catalog.csv")

# 2. Find available TESS sectors
sectors = get_tess_sectors(star)
sector = int(sectors["sector"][0])

# 3. Download sector script
script_path = download_tess_sector_script(sector)

# 4. Process images in parallel
process_images_parallel(
    script_file=script_path,
    catalog_file="catalog.csv",
    output_dir="photometry_results",
    star=star
)

# 5. Load and plot light curve
times, fluxes = load_photometry_data("photometry_results", star)

plt.scatter(times, fluxes)
plt.xlabel("Julian Date")
plt.ylabel("Flux (e⁻/s)")
plt.title(f"Light Curve: {star.name}")
plt.show()


