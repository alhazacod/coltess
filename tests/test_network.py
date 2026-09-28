"""Opt-in end-to-end test against ONE real TESS FFI from MAST.

This is the README Quick Start on a single image. It downloads roughly
200 MB, so it is skipped unless explicitly enabled:

    COLTESS_NETWORK_TESTS=1 pytest tests/test_network.py -s

Optional environment variables:
    COLTESS_TEST_STAR    target name resolvable by SIMBAD (default: "lambda tau")
    COLTESS_TEST_RADIUS  Gaia cone-search radius in arcmin (default: 10)

Detection depends on whether the target happens to fall inside that
particular FFI tile, so the test asserts the *pipeline* runs correctly on
real data and, when the target is present, that the light curve round-trips.
"""

from __future__ import annotations

import os

import pandas as pd
import pytest

from coltess import (
    StarData,
    TessPhotometry,
    create_catalog,
    download_tess_images,
    download_tess_sector_script,
    get_tess_sectors,
    load_photometry_data,
    process_images_parallel,
)

NETWORK_ENABLED = os.environ.get("COLTESS_NETWORK_TESTS", "").strip().lower() in {
    "1",
    "true",
    "yes",
}

pytestmark = [
    pytest.mark.network,
    pytest.mark.skipif(
        not NETWORK_ENABLED,
        reason="set COLTESS_NETWORK_TESTS=1 to run real-download tests",
    ),
]


@pytest.fixture(scope="module")
def quick_start(tmp_path_factory):
    """Run the Quick Start up to (but not including) parallel processing."""
    workdir = tmp_path_factory.mktemp("network_quickstart")

    star_name = os.environ.get("COLTESS_TEST_STAR", "lambda tau")
    radius = float(os.environ.get("COLTESS_TEST_RADIUS", "10"))

    catalog_file = str(workdir / "catalog.csv")
    star: StarData = create_catalog(
        star_name, radius_arcmin=radius, output_file=catalog_file
    )

    sectors = get_tess_sectors(star)
    assert len(sectors) > 0, f"no TESS sectors found for {star_name}"
    sector = int(sectors["sector"].to_numpy()[0])

    # The download script is written to the current directory.
    original_cwd = os.getcwd()
    try:
        os.chdir(workdir)
        script_path = download_tess_sector_script(sector)
    finally:
        os.chdir(original_cwd)

    # Keep only the first image so the run takes minutes, not days.
    with open(script_path) as f:
        first_line = f.readline().strip()
    one_image_script = workdir / f"tesscurl_sector_{sector}_ffic_one.sh"
    one_image_script.write_text(first_line + "\n")

    return star, catalog_file, str(one_image_script)


def test_quick_start_on_one_real_image(quick_start, tmp_path):
    star, catalog_file, script_file = quick_start
    output_file = str(tmp_path / "lightcurve.csv")

    process_images_parallel(
        script_file=script_file,
        catalog_file=catalog_file,
        star=star,
        output_file=output_file,
        max_workers=1,
    )

    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
        pytest.skip(
            f"{star.name} is not in the first FFI tile of this sector; "
            "rerun with more images or another target to exercise detection"
        )

    df = pd.read_csv(output_file)
    assert len(df) >= 1
    assert {"flux", "flux_err", "RA", "DEC", "DATE-OBS", "SECTOR"} <= set(df.columns)
    assert bool(df["flux"].notna().all())

    times, fluxes, flux_errors = load_photometry_data(output_file, star)
    assert len(times) == len(fluxes) == len(flux_errors) >= 1
    assert fluxes[0] > 0
    assert flux_errors[0] > 0


def test_real_ffi_pipeline_runs(quick_start, tmp_path):
    """`TessPhotometry` must handle a real FFI header/WCS without crashing."""
    star, catalog_file, script_file = quick_start
    image_dir = str(tmp_path / "images")

    download_tess_images(script_file, start_idx=0, num_images=1, output_dir=image_dir)

    fits_files = [f for f in os.listdir(image_dir) if f.endswith(".fits")]
    assert fits_files, "no FFI was downloaded"

    processor = TessPhotometry()
    catalog = processor.load_catalog(catalog_file)
    result = processor.process_fits(os.path.join(image_dir, fits_files[0]), catalog)

    # Real TESS FFIs are in e-/s, so the gain must be the exposure time.
    assert processor.epadu > 0

    if result is not None:
        assert len(result) > 0
        assert set(result.colnames) >= {"flux", "mag", "mag_err", "flux_err"}
