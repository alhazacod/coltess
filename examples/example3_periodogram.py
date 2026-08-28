from coltess import create_catalog, get_tess_sectors, download_tess_sector_script
from coltess import process_images_parallel, load_photometry_data, compute_periodogram
import matplotlib.pyplot as plt

# 1. Create catalog and get star info.
star = create_catalog("HD 2655", radius_arcmin=10.0, output_file="catalog.csv")

# 2. Find available TESS sectors
sectors = get_tess_sectors(star)
print(f"Sectors: {sectors}")
sector = int(sectors["sector"][0])

# 3. Download sector script
script_path = download_tess_sector_script(sector)

# 4. Process images in parallel, appending results to a single CSV
process_images_parallel(
    script_file=script_path,
    catalog_file="catalog.csv",
    output_file="photometry_results.csv",
    star=star,
    start_idx=0,
)

# 5. Load light curve, compute periodogram, plot both
times, fluxes, flux_errors = load_photometry_data("photometry_results.csv", star)
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
