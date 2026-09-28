"""Shared fixtures for the coltess test suite.

The offline suite builds small synthetic TESS-like Full Frame Images so the
whole photometry pipeline can be exercised in seconds without downloading
the ~200 MB real FFIs. The real end-to-end test lives in
``test_network.py`` and is opt-in (``COLTESS_NETWORK_TESTS=1``).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

# --- Synthetic field geometry ---------------------------------------------
# A small 512x512 image with the TESS pixel scale (~21 arcsec/pixel).
IMAGE_SIZE = 512
PIXEL_SCALE_DEG = 21.0 / 3600.0

# Field centre and the exact pixel the target star is placed on.
FIELD_RA = 120.0
FIELD_DEC = -30.0
TARGET_X = 256.0
TARGET_Y = 256.0

# A few extra catalog sources, far enough from the target that they are
# never confused with it (30+ pixels ~= 10 arcmin).
EXTRA_SOURCES = [(100.0, 120.0), (400.0, 130.0), (300.0, 440.0)]

# Multi-epoch dataset used by the periodogram integration test.
N_IMAGES = 12
CADENCE_DAYS = 0.20
PERIOD_DAYS = 1.0
EPOCH_JD = 2457000.0
BACKGROUND = 5.0
NOISE_SIGMA = 2.0
TARGET_AMPLITUDE = 5000.0
GAUSS_SIGMA = 2.0


def build_wcs() -> WCS:
    """Return a simple TAN WCS covering the synthetic field."""
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [IMAGE_SIZE / 2.0, IMAGE_SIZE / 2.0]
    wcs.wcs.crval = [FIELD_RA, FIELD_DEC]
    wcs.wcs.cdelt = [-PIXEL_SCALE_DEG, PIXEL_SCALE_DEG]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return wcs


def target_sky_position() -> tuple[float, float]:
    """Sky coordinates (deg) of the synthetic target star."""
    world = build_wcs().pixel_to_world(TARGET_X, TARGET_Y)
    return float(world.ra.deg), float(world.dec.deg)


def extra_sky_positions() -> list[tuple[float, float]]:
    """Sky coordinates (deg) of the non-target synthetic sources."""
    wcs = build_wcs()
    positions = []
    for x, y in EXTRA_SOURCES:
        world = wcs.pixel_to_world(x, y)
        positions.append((float(world.ra.deg), float(world.dec.deg)))
    return positions


def _add_gaussian(image: np.ndarray, x0: float, y0: float, amplitude: float) -> None:
    yy, xx = np.mgrid[0 : image.shape[0], 0 : image.shape[1]]
    image += amplitude * np.exp(
        -((xx - x0) ** 2 + (yy - y0) ** 2) / (2.0 * GAUSS_SIGMA**2)
    )


def write_ffi(
    path,
    *,
    date_obs: str,
    target_amplitude: float = TARGET_AMPLITUDE,
    bunit: str = "e-/s",
    exptime: float = 200.0,
    sector: int | None = 57,
    seed: int = 12345,
) -> str:
    """Write one synthetic TESS-like FFI (image in HDU 1, as TESS does)."""
    rng = np.random.default_rng(seed)
    image = np.full((IMAGE_SIZE, IMAGE_SIZE), BACKGROUND, dtype=np.float64)
    image += rng.normal(0.0, NOISE_SIGMA, image.shape)

    _add_gaussian(image, TARGET_X, TARGET_Y, target_amplitude)
    for (x, y), scale in zip(EXTRA_SOURCES, (0.4, 0.7, 0.5), strict=True):
        _add_gaussian(image, x, y, TARGET_AMPLITUDE * scale)

    header = build_wcs().to_header()
    header["BUNIT"] = bunit
    header["EXPTIME"] = exptime
    header["DATE-OBS"] = date_obs
    if sector is not None:
        header["SECTOR"] = sector

    hdul = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(data=image, header=header)])
    hdul.writeto(path, overwrite=True)
    return str(path)


def write_catalog(path, *, include_target: bool = True) -> str:
    """Write a Gaia-like catalog CSV (ra, dec, source_id)."""
    rows = []
    if include_target:
        ra, dec = target_sky_position()
        rows.append({"ra": ra, "dec": dec, "source_id": "1"})
    for i, (ra, dec) in enumerate(extra_sky_positions(), start=2):
        rows.append({"ra": ra, "dec": dec, "source_id": str(i)})

    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def target_amplitude_for_epoch(epoch: int) -> float:
    """Sinusoidal target flux used to create a recoverable signal."""
    phase = 2.0 * np.pi * (epoch * CADENCE_DAYS) / PERIOD_DAYS
    return TARGET_AMPLITUDE * (1.0 + 0.25 * np.sin(phase))


@pytest.fixture(scope="session")
def synthetic_catalog(tmp_path_factory):
    """Single Gaia-like catalog shared by the whole session."""
    path = tmp_path_factory.mktemp("catalog") / "catalog.csv"
    return write_catalog(path)


@pytest.fixture(scope="session")
def synthetic_catalog_without_target(tmp_path_factory):
    """Catalog whose nearest source is ~10 arcmin away from the target."""
    path = tmp_path_factory.mktemp("catalog_notarget") / "catalog.csv"
    return write_catalog(path, include_target=False)


@pytest.fixture(scope="session")
def single_ffi(tmp_path_factory):
    """One synthetic FFI with the target present."""
    path = tmp_path_factory.mktemp("single") / "tess-s0057-1-1-0001.fits"
    date_obs = Time(EPOCH_JD, format="jd", scale="utc").isot
    return write_ffi(path, date_obs=date_obs)


@pytest.fixture(scope="session")
def single_ffi_adu(tmp_path_factory):
    """One synthetic FFI in ADU units (tests the fixed-gain branch)."""
    path = tmp_path_factory.mktemp("single_adu") / "tess-s0057-1-1-0002.fits"
    date_obs = Time(EPOCH_JD, format="jd", scale="utc").isot
    return write_ffi(path, date_obs=date_obs, bunit="ADU")


@pytest.fixture(scope="session")
def multi_epoch_dir(tmp_path_factory):
    """A directory of N_IMAGES synthetic FFIs with a periodic signal."""
    directory = tmp_path_factory.mktemp("epochs")
    for epoch in range(N_IMAGES):
        jd = EPOCH_JD + epoch * CADENCE_DAYS
        date_obs = Time(jd, format="jd", scale="utc").isot
        write_ffi(
            directory / f"tess-s0057-1-1-{epoch:04d}.fits",
            date_obs=date_obs,
            target_amplitude=target_amplitude_for_epoch(epoch),
            seed=1000 + epoch,
        )
    return str(directory)
