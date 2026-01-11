#!/usr/bin/env python3
"""
Parallel TESS photometry example using Coltess
"""

from pathlib import Path
import matplotlib.pyplot as plt

from coltess.catalog import create_catalog
from coltess.download import get_tess_sectors, download_tess_sector_script
from coltess.parallel import process_images_parallel
from coltess.analysis import load_photometry_data
from coltess.core import StarData


def main():
    star_name = "lambda tau"

    base_dir = Path("run_lambda_tau")
    base_dir.mkdir(exist_ok=True)

    catalog_file = base_dir / "gaia_catalog.csv"
    csv_dir = base_dir / "csv"

    csv_dir.mkdir(exist_ok=True)

    # --------------------------------------------------
    # 1. Create Gaia catalog + StarData
    # --------------------------------------------------
    print(f"Creating catalog for {star_name}")

    star: StarData = create_catalog(
        star_name,
        radius_arcmin=10.0,
        output_file=str(catalog_file),
    )

    print(
        f"Resolved {star.name}: "
        f"RA={star.ra:.6f}, DEC={star.dec:.6f}, Gaia={star.gaia_id}"
    )

    # --------------------------------------------------
    # 2. Find TESS sector
    # --------------------------------------------------
    sectors = get_tess_sectors(star)
    print(sectors)

    sector = int(sectors["sector"][0])
    print(f"Using TESS sector {sector}")

    script_path = download_tess_sector_script(sector)

    # --------------------------------------------------
    # 3. Parallel download + photometry
    # --------------------------------------------------
    process_images_parallel(
        script_file=script_path,
        catalog_file=str(catalog_file),
        output_dir=str(csv_dir),
        star=star,
        start_idx=19000,
        max_workers=6,
    )

    # --------------------------------------------------
    # 4. Load photometry and plot
    # --------------------------------------------------
    jd, flux = load_photometry_data(
        csv_dir=str(csv_dir),
        target_star=star,
        max_sep_arcsec=1.0,
    )

    plt.figure(figsize=(8, 4))
    plt.scatter(jd, flux, s=20)
    plt.xlabel("Julian Date (JD)")
    plt.ylabel("Flux")
    plt.title(f"Light curve: {star.name}")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
