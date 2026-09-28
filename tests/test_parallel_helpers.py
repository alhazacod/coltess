"""Unit tests for the parallel/resume bookkeeping in ``coltess.parallel``."""

from __future__ import annotations

import pandas as pd

from coltess.parallel import (
    _append_row,
    _processed_sources,
    _resume_start_idx,
    _sector_from_fits,
    _sector_from_script,
)


def test_sector_from_script():
    assert _sector_from_script("tesscurl_sector_57_ffic.sh") == 57
    assert _sector_from_script("some/nested/tesscurl_sector_14_ffic.sh") == 14
    assert _sector_from_script("random.sh") is None


def test_sector_from_fits_header(tmp_path, single_ffi):
    assert _sector_from_fits(single_ffi) == 57


def test_sector_from_fits_filename(tmp_path):
    path = tmp_path / "tess2019213215929-s0014-3-1-0150-s_ffic.fits"
    path.write_bytes(b"not a real fits file")
    assert _sector_from_fits(str(path)) == 14


def test_sector_from_fits_unknown(tmp_path):
    path = tmp_path / "mystery.fits"
    path.write_bytes(b"not a real fits file")
    assert _sector_from_fits(str(path)) is None


def test_processed_sources_missing_or_empty(tmp_path):
    assert _processed_sources(str(tmp_path / "missing.csv")) is None
    empty = tmp_path / "empty.csv"
    empty.write_text("")
    assert _processed_sources(str(empty)) is None


def test_resume_start_idx(tmp_path):
    script = tmp_path / "tesscurl_sector_57_ffic.sh"
    lines = [
        "curl -o a.fits URL/A\n",
        "curl -o b.fits URL/B\n",
        "curl -o c.fits URL/C\n",
    ]
    script.write_text("".join(lines))

    assert _resume_start_idx(lines, str(tmp_path / "new.csv")) == 0

    output = tmp_path / "out.csv"
    pd.DataFrame({"CURL": ["curl -o a.fits URL/A", "curl -o b.fits URL/B"]}).to_csv(
        output, index=False
    )
    assert _resume_start_idx(lines, str(output)) == 2


def test_resume_falls_back_when_last_image_unknown(tmp_path):
    lines = ["curl -o a.fits URL/A\n"]
    output = tmp_path / "out.csv"
    pd.DataFrame({"CURL": ["curl -o zzz.fits URL/Z"]}).to_csv(output, index=False)
    assert _resume_start_idx(lines, str(output)) == 0


def test_append_row_writes_header_once_and_deduplicates(tmp_path):
    output = tmp_path / "out.csv"
    df = pd.DataFrame(
        {"flux": [1.0], "flux_err": [0.1], "SECTOR": [57], "CURL": ["src-a"]}
    )

    _append_row(df, str(output), "src-a")
    _append_row(df, str(output), "src-a")  # duplicate: skipped
    second = df.copy()
    second["flux"] = [2.0]
    second["CURL"] = ["src-b"]
    _append_row(second, str(output), "src-b")

    result = pd.read_csv(output)
    assert len(result) == 2
    assert list(result["CURL"]) == ["src-a", "src-b"]
    assert list(result["flux"]) == [1.0, 2.0]
