"""End-to-end photometry pipeline tests on synthetic TESS-like FFIs.

These mirror the README Quick Start (catalog -> per-image photometry ->
light curve -> periodogram) but run against tiny local images, so they are
fast and fully offline.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from astropy.time import Time

from coltess import (
    StarData,
    TessPhotometry,
    analyze_image,
    compute_periodogram,
    load_photometry_data,
    process_local_images_parallel,
)
from conftest import EPOCH_JD, N_IMAGES, target_sky_position, write_ffi


@pytest.fixture
def target_star() -> StarData:
    ra, dec = target_sky_position()
    return StarData(name="synthetic target", ra=ra, dec=dec)


# --- Single image ----------------------------------------------------------


def test_process_fits_measures_target(single_ffi, synthetic_catalog):
    """TessPhotometry.process_fits returns the catalog sources in frame."""
    processor = TessPhotometry()
    catalog = processor.load_catalog(synthetic_catalog)
    result = processor.process_fits(single_ffi, catalog)

    assert result is not None
    assert len(result) == 1 + 3  # target + 3 extra sources
    assert {"flux", "mag", "mag_err", "flux_err"} <= set(result.colnames)
    assert np.all(np.isfinite(result["flux"]))


def test_analyze_image_finds_target(single_ffi, synthetic_catalog, target_star):
    """The Quick Start single-image path returns exactly the target row."""
    row = analyze_image(single_ffi, synthetic_catalog, target_star)

    assert row is not None
    assert len(row) == 1
    assert row["flux"][0] > 0
    assert row["flux_err"][0] > 0
    assert row["ID"][0] == "1"


def test_target_not_detected_returns_none(
    single_ffi, synthetic_catalog_without_target, target_star
):
    """Frames where only distant sources exist must be rejected cleanly."""
    row = analyze_image(single_ffi, synthetic_catalog_without_target, target_star)
    assert row is None


def test_flux_scales_with_amplitude(tmp_path, synthetic_catalog, target_star):
    """A brighter synthetic star yields a proportionally larger flux."""
    date_obs = Time(EPOCH_JD, format="jd", scale="utc").isot
    dim = write_ffi(
        tmp_path / "dim.fits",
        date_obs=date_obs,
        target_amplitude=2500.0,
        seed=7,
    )
    bright = write_ffi(
        tmp_path / "bright.fits",
        date_obs=date_obs,
        target_amplitude=5000.0,
        seed=7,
    )

    dim_flux = analyze_image(dim, synthetic_catalog, target_star)["flux"][0]
    bright_flux = analyze_image(bright, synthetic_catalog, target_star)["flux"][0]

    assert bright_flux > dim_flux
    assert 1.5 < (bright_flux / dim_flux) < 2.5


# --- Effective gain detection ---------------------------------------------


def test_gain_from_exptime_for_electrons_per_second(single_ffi, synthetic_catalog):
    processor = TessPhotometry()
    processor.process_fits(single_ffi, processor.load_catalog(synthetic_catalog))
    assert processor.epadu == pytest.approx(200.0)


def test_gain_is_nominal_for_adu(single_ffi_adu, synthetic_catalog):
    processor = TessPhotometry()
    processor.process_fits(single_ffi_adu, processor.load_catalog(synthetic_catalog))
    assert processor.epadu == pytest.approx(5.22)


# --- Multi-image Quick Start chain ----------------------------------------


def test_local_images_to_periodogram(
    multi_epoch_dir, synthetic_catalog, target_star, tmp_path
):
    """Full offline Quick Start: photometry -> light curve -> periodogram."""
    output_file = str(tmp_path / "lightcurve.csv")

    process_local_images_parallel(
        fits_dir=multi_epoch_dir,
        catalog_file=synthetic_catalog,
        star=target_star,
        output_file=output_file,
        max_workers=2,
    )

    times, fluxes, flux_errors = load_photometry_data(output_file, target_star)

    assert len(times) == N_IMAGES
    assert len(fluxes) == len(flux_errors) == N_IMAGES
    assert np.all(np.isfinite(times))
    assert np.all(np.isfinite(fluxes))
    assert np.all(np.isfinite(flux_errors))
    assert np.all(np.diff(times) > 0)
    assert np.all(flux_errors > 0)

    results = compute_periodogram(times, fluxes, flux_errors)
    assert results["primary_period"] is not None
    # Recover the injected period (generous tolerance: 12 points, 2.2 d baseline).
    assert 0.7 < results["primary_period"] < 1.4
    assert results["primary_period_uncertainty"] is not None
    assert results["primary_fap"] is not None
    assert len(results["periods"]) == len(results["power"])


def test_rerun_does_not_duplicate_rows(
    multi_epoch_dir, synthetic_catalog, target_star, tmp_path
):
    """Re-running over the same folder must resume, not append duplicates."""
    output_file = str(tmp_path / "resume.csv")

    for _ in range(2):
        process_local_images_parallel(
            fits_dir=multi_epoch_dir,
            catalog_file=synthetic_catalog,
            star=target_star,
            output_file=output_file,
            max_workers=2,
        )

    df = pd.read_csv(output_file)
    assert len(df) == N_IMAGES
    assert df["CURL"].nunique() == N_IMAGES


def test_output_metadata_columns(
    multi_epoch_dir, synthetic_catalog, target_star, tmp_path
):
    output_file = str(tmp_path / "metadata.csv")
    process_local_images_parallel(
        fits_dir=multi_epoch_dir,
        catalog_file=synthetic_catalog,
        star=target_star,
        output_file=output_file,
        max_workers=2,
    )

    df = pd.read_csv(output_file)
    assert {"flux", "flux_err", "RA", "DEC", "DATE-OBS", "SECTOR", "CURL"} <= set(
        df.columns
    )
    assert set(df["SECTOR"]) == {57}


# --- Periodogram edge cases ------------------------------------------------


def test_periodogram_without_errors(
    multi_epoch_dir, synthetic_catalog, target_star, tmp_path
):
    output_file = str(tmp_path / "noerr.csv")
    process_local_images_parallel(
        fits_dir=multi_epoch_dir,
        catalog_file=synthetic_catalog,
        star=target_star,
        output_file=output_file,
        max_workers=2,
    )
    times, fluxes, _ = load_photometry_data(output_file, target_star)
    results = compute_periodogram(times, fluxes, flux_errors=None)
    assert results["primary_period"] is not None


def test_periodogram_too_few_points():
    results = compute_periodogram(np.array([1.0, 2.0]), np.array([1.0, 2.0]))
    assert results["primary_period"] is None
    assert results["primary_fap"] is None
    assert len(results["periods"]) == 0


def test_periodogram_ignores_nans():
    times = np.arange(0.0, 10.0, 0.5)
    fluxes = 10.0 + 3.0 * np.sin(2.0 * np.pi * times / 2.0)

    # Inject NaNs into both arrays; they must be masked, not crash the fit.
    times = np.append(times, np.nan)
    fluxes = np.append(fluxes, 5.0)
    fluxes = np.where(np.isclose(times, 4.0), np.nan, fluxes)

    results = compute_periodogram(times, fluxes, min_period=0.5, max_period=5.0)
    assert results["primary_period"] is not None
    assert results["primary_period"] == pytest.approx(2.0, rel=0.1)


def test_load_photometry_missing_file(target_star):
    with pytest.raises(RuntimeError):
        load_photometry_data("does_not_exist.csv", target_star)
