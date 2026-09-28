"""Smoke tests: the package imports, exposes its public API, ships py.typed."""

from __future__ import annotations

import importlib.resources

import coltess

PUBLIC_API = [
    "StarData",
    "TessPhotometry",
    "analyze_image",
    "create_catalog",
    "get_star",
    "query_gaia_catalog",
    "get_tess_sectors",
    "download_tess_sector_script",
    "download_tess_image",
    "download_tess_images",
    "load_photometry_data",
    "compute_periodogram",
    "process_images_parallel",
    "process_local_images_parallel",
]


def test_version_is_a_nonempty_string():
    assert isinstance(coltess.__version__, str)
    assert coltess.__version__.strip()


def test_public_api_is_importable():
    for name in PUBLIC_API:
        assert hasattr(coltess, name), f"missing public export: {name}"


def test_all_matches_public_api():
    assert set(coltess.__all__) == set(PUBLIC_API)


def test_py_typed_marker_is_shipped():
    """PEP 561: the marker must be inside the installed package."""
    marker = importlib.resources.files("coltess").joinpath("py.typed")
    assert marker.is_file(), "coltess/py.typed is missing from the package"
