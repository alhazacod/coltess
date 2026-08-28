#!/usr/bin/env python3
"""
Parallel pipeline for TESS photometry.

This module provides functions to download and process large numbers of
TESS Full Frame Images (FFIs) in parallel for a single target star, or to
analyze FITS images that are already stored locally. Each image is
analyzed independently and one photometry row per detection is appended
to a single output CSV file.

Parallelization is implemented using ``multiprocessing.Pool``.
"""

import glob
import os
import re
import sys
import shutil
import tempfile
import multiprocessing as mp
from functools import partial

import pandas as pd
from astropy.io import fits

try:
    import fcntl
except ImportError:
    fcntl = None  # type: ignore[assignment]

from coltess.core import StarData
from coltess.photometry import TessPhotometry
from coltess.download import download_tess_images


def process_images_parallel(
    script_file: str,
    catalog_file: str,
    star: StarData,
    output_file: str | None = None,
    start_idx: int | None = None,
    max_workers: int | None = None,
    keep_images_dir: str | None = None,
):
    """
    Download and process TESS images in parallel for a target star.

    This function coordinates the parallel execution of photometry over
    a list of TESS image URLs or commands contained in a script file.
    Each worker downloads exactly one FITS file, performs photometry,
    appends one row per detection to a single output CSV file, and
    deletes temporary files.

    Parameters
    ----------
    script_file : str
        Path to a TESS download script (one image per line).
    catalog_file : str
        Path to a Gaia catalog CSV.
    star : StarData
        Target star information.
    output_file : str or None, optional
        Path of the single CSV file that photometry rows are appended to.
        Defaults to ``<star name>_<sector>.csv``.
    start_idx : int or None, optional
        Line index in the script file from which to start processing.
        If None (default) and the output file exists with data,
        processing resumes automatically after the last processed image.
        Passing an explicit value overrides the automatic resume.
    max_workers : int or None, optional
        Number of parallel worker processes. Defaults to the number of
        available CPU cores.
    keep_images_dir : str or None, optional
        Directory where FITS images in which the star was detected are
        saved. If None, no images are kept.

    Notes
    -----
    - Each FITS file is handled independently.
    - Rows are appended to the output file under an exclusive file lock,
      so concurrent workers never corrupt it. Duplicate rows (same curl
      command) are skipped.
    - Temporary FITS files are stored in per-worker directories and
      deleted after processing.
    - Pressing ``Ctrl+C`` terminates all workers immediately and exits
      with status code 130.
    """
    if max_workers is None:
        max_workers = mp.cpu_count()

    sector = _sector_from_script(script_file)
    if sector is not None:
        star.sector = sector

    if output_file is None:
        safe_name = star.name.replace(" ", "_")
        if sector is not None:
            output_file = f"{safe_name}_{sector}.csv"
        else:
            output_file = f"{safe_name}.csv"

    output_dir = os.path.dirname(output_file) or "."
    os.makedirs(output_dir, exist_ok=True)

    with open(script_file) as f:
        lines = f.readlines()

    n_images = len(lines)

    if start_idx is None:
        start_idx = _resume_start_idx(lines, output_file)

    print(f"Processing images {start_idx} -> {n_images - 1}")
    print(f"Using {max_workers} workers")

    indices = list(range(start_idx, n_images))

    worker = partial(
        worker_process_fits,
        script_file,
        catalog_file=catalog_file,
        star=star,
        output_file=output_file,
        sector=sector,
        keep_images_dir=keep_images_dir,
    )

    _run_pool(worker, indices, max_workers)


def _run_pool(worker, tasks: list, max_workers: int) -> None:
    """
    Run ``worker`` over ``tasks`` in a multiprocessing pool.

    Parameters
    ----------
    worker : callable
        Function taking one task item and returning a result tuple.
    tasks : list
        Task items to process.
    max_workers : int
        Number of parallel worker processes.

    Notes
    -----
    Pressing ``Ctrl+C`` terminates all workers immediately and exits
    with status code 130.
    """
    try:
        pool = mp.get_context("fork").Pool(processes=max_workers, maxtasksperchild=100)
    except ValueError:
        # Windows without WSL - needs __main__ guard
        print(
            "WARNING: Using 'spawn' method. Scripts should use if __name__ == '__main__' or process images sequentially"
        )
        pool = mp.get_context("spawn").Pool(processes=max_workers, maxtasksperchild=100)

    try:
        for _ in pool.imap_unordered(worker, tasks):
            pass

    except KeyboardInterrupt:
        print("\nCtrl+C detected — terminating workers immediately...")
        pool.terminate()
        pool.join()
        sys.exit(130)

    else:
        pool.close()
        pool.join()


def _sector_from_script(script_file: str) -> int | None:
    """
    Extract the TESS sector number from the download script filename.

    Parameters
    ----------
    script_file : str
        Path to a TESS download script.

    Returns
    -------
    int or None
        Sector number, or None if it cannot be determined.
    """
    match = re.search(r"tesscurl_sector_(\d+)_ffic\.sh", os.path.basename(script_file))
    return int(match.group(1)) if match else None


def _sector_from_fits(fits_file: str) -> int | None:
    """
    Extract the TESS sector number from a FITS file.

    The sector is read from the ``SECTOR`` header keyword, falling back
    to the ``-sNNNN-`` pattern in the TESS filename.

    Parameters
    ----------
    fits_file : str
        Path to a TESS FITS image.

    Returns
    -------
    int or None
        Sector number, or None if it cannot be determined.
    """
    try:
        with fits.open(fits_file) as hdul:
            sector = hdul[1].header.get("SECTOR")
            if sector is not None:
                return int(sector)
    except Exception:
        pass

    match = re.search(r"-s(\d{4})-", os.path.basename(fits_file))
    return int(match.group(1)) if match else None


def _processed_sources(output_file: str) -> list | None:
    """
    Return the stored source values from the output file, in append order.

    Parameters
    ----------
    output_file : str
        Path of the output CSV file.

    Returns
    -------
    list of str or None
        Stored ``CURL`` column values (download commands or local file
        paths), or None if the file is missing or empty.
    """
    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
        return None

    try:
        return list(pd.read_csv(output_file, usecols=["CURL"])["CURL"].astype(str))
    except Exception:
        return None


def _resume_start_idx(lines: list[str], output_file: str) -> int:
    """
    Determine the starting script line index for resuming a run.

    If the output file exists and contains data, the image filename from
    its last row is matched against the script lines and processing
    resumes from the line after it.

    Parameters
    ----------
    lines : list of str
        Lines of the TESS download script.
    output_file : str
        Path of the output CSV file.

    Returns
    -------
    int
        Line index from which to start processing.
    """
    processed = _processed_sources(output_file)
    if processed is None:
        return 0

    last_curl = processed[-1]
    last_filename = os.path.basename(last_curl.split()[-1])

    for i, line in enumerate(lines):
        tokens = line.split()
        if not tokens:
            continue
        if os.path.basename(tokens[-1]) == last_filename:
            resume = i + 1
            print(
                f"Resuming after image {last_filename} "
                f"(script line {i + 1}); starting at line {resume}"
            )
            return resume

    print(
        f"WARNING: last processed image {last_filename} not found in "
        "the script; starting from line 0"
    )
    return 0


def worker_process_fits(
    script_file: str,
    index: int,
    catalog_file: str,
    star: StarData,
    output_file: str,
    sector: int | None = None,
    keep_images_dir: str | None = None,
):
    """
    Process a single TESS FITS image.

    This worker function performs the following steps:
    1. Downloads exactly one FITS file specified by ``index`` in the
       download script.
    2. Runs aperture photometry for the target star.
    3. Appends one row (with sector and curl command) to the output CSV
       file, protected by an exclusive file lock.
    4. Optionally copies the FITS file to ``keep_images_dir`` when the
       star was detected.
    5. Deletes all temporary files and directories.

    Parameters
    ----------
    script_file : str
        Path to the TESS download script.
    index : int
        Line index in the script file corresponding to the image to process.
    catalog_file : str
        Path to a Gaia catalog CSV.
    star : StarData
        Target star information.
    output_file : str
        Path of the single CSV file the photometry row is appended to.
    sector : int or None, optional
        TESS sector number stored in the output row.
    keep_images_dir : str or None, optional
        Directory where the FITS image is saved if the star was detected.
        If None, no image is kept.

    Returns
    -------
    tuple
        (index, success) where ``success`` is True if photometry was
        successfully performed and saved, False otherwise.

    Notes
    -----
    - Each worker runs in its own process and uses a unique temporary
      directory.
    - All temporary files are deleted even if an exception occurs.
    - Appends are serialized with ``fcntl.flock`` and rows with an
      already present curl command are skipped.
    """

    success = False

    print(f"[PID {os.getpid()}] " f"Processing script line {index + 1}")

    with open(script_file) as f:
        lines = f.readlines()
    curl_line = lines[index].strip()

    tmp_dir = tempfile.mkdtemp(prefix="tess_")

    try:
        download_tess_images(
            script_file,
            start_idx=index,
            num_images=1,
            output_dir=tmp_dir,
        )

        fits_files = [
            os.path.join(tmp_dir, f) for f in os.listdir(tmp_dir) if f.endswith(".fits")
        ]

        row = None
        if fits_files:
            processor = TessPhotometry()
            row = processor.process_image(
                fits_files[0],
                catalog_file=catalog_file,
                target_star=star,
            )
            success = row is not None

        if row is not None:
            _append_photometry_row(
                row,
                fits_files[0],
                output_file,
                sector=sector,
                source=curl_line,
                source_label=f"line {index + 1}",
                keep_images_dir=keep_images_dir,
            )

        return index, success

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def process_local_images_parallel(
    fits_dir: str,
    catalog_file: str,
    star: StarData,
    output_file: str | None = None,
    max_workers: int | None = None,
    keep_images_dir: str | None = None,
    pattern: str = "*.fits",
):
    """
    Analyze FITS images already stored locally, in parallel.

    This function performs photometry over all FITS files in a directory
    without downloading anything. One row per detection is appended to a
    single output CSV file, with the FITS header ``SECTOR`` and the local
    file path (in the ``CURL`` column) as metadata.

    Parameters
    ----------
    fits_dir : str
        Directory containing the FITS images.
    catalog_file : str
        Path to a Gaia catalog CSV.
    star : StarData
        Target star information.
    output_file : str or None, optional
        Path of the single CSV file that photometry rows are appended to.
        Defaults to ``<star name>.csv``.
    max_workers : int or None, optional
        Number of parallel worker processes. Defaults to the number of
        available CPU cores.
    keep_images_dir : str or None, optional
        Directory where FITS images in which the star was detected are
        copied. If None, no images are kept.
    pattern : str, optional
        Glob pattern used to select files inside ``fits_dir``. Defaults
        to ``"*.fits"``.

    Notes
    -----
    - Images whose path is already stored in the output file are skipped,
      so re-running resumes automatically. The deduplication keys on the
      exact path string; different path spellings for the same file will
      not be recognized as duplicates.
    - Rows are appended under an exclusive file lock.
    - Pressing ``Ctrl+C`` terminates all workers immediately and exits
      with status code 130.
    """
    if max_workers is None:
        max_workers = mp.cpu_count()

    fits_files = sorted(glob.glob(os.path.join(fits_dir, pattern)))

    if not fits_files:
        raise FileNotFoundError(
            f"No FITS files found in {fits_dir} (pattern {pattern!r})"
        )

    if output_file is None:
        output_file = f"{star.name.replace(' ', '_')}.csv"

    output_path_dir = os.path.dirname(output_file) or "."
    os.makedirs(output_path_dir, exist_ok=True)

    processed = _processed_sources(output_file)
    if processed:
        already = set(processed)
        remaining = [f for f in fits_files if f not in already]
        print(f"Skipping {len(fits_files) - len(remaining)} already-processed images")
        fits_files = remaining

    print(f"Processing {len(fits_files)} local images")
    print(f"Using {max_workers} workers")

    worker = partial(
        worker_process_local_fits,
        catalog_file=catalog_file,
        star=star,
        output_file=output_file,
        keep_images_dir=keep_images_dir,
    )

    _run_pool(worker, fits_files, max_workers)


def worker_process_local_fits(
    fits_file: str,
    catalog_file: str,
    star: StarData,
    output_file: str,
    keep_images_dir: str | None = None,
):
    """
    Process a single local FITS image.

    This worker function performs the following steps:
    1. Runs aperture photometry for the target star.
    2. Appends one row (with sector and file path) to the output CSV
       file, protected by an exclusive file lock.
    3. Optionally copies the FITS file to ``keep_images_dir`` when the
       star was detected.

    Parameters
    ----------
    fits_file : str
        Path to the local FITS image.
    catalog_file : str
        Path to a Gaia catalog CSV.
    star : StarData
        Target star information.
    output_file : str
        Path of the single CSV file the photometry row is appended to.
    keep_images_dir : str or None, optional
        Directory where the FITS image is copied if the star was detected.
        If None, no image is kept.

    Returns
    -------
    tuple
        (fits_file, success) where ``success`` is True if photometry was
        successfully performed and saved, False otherwise.

    Notes
    -----
    - Appends are serialized with ``fcntl.flock`` and rows with an
      already present file path are skipped.
    """

    print(f"[PID {os.getpid()}] " f"Processing {os.path.basename(fits_file)}")

    processor = TessPhotometry()
    row = processor.process_image(
        fits_file,
        catalog_file=catalog_file,
        target_star=star,
    )
    success = row is not None

    if row is not None:
        sector = _sector_from_fits(fits_file)
        _append_photometry_row(
            row,
            fits_file,
            output_file,
            sector=sector,
            source=fits_file,
            source_label=os.path.basename(fits_file),
            keep_images_dir=keep_images_dir,
        )

    return fits_file, success


def _append_photometry_row(
    row,
    fits_file: str,
    output_file: str,
    sector: int | None = None,
    source: str = "",
    source_label: str = "",
    keep_images_dir: str | None = None,
) -> None:
    """
    Append one photometry row to the output file and optionally keep the FITS.

    The row is enriched with the ``SECTOR`` and ``CURL`` columns, appended
    under an exclusive file lock, and the FITS file is copied to
    ``keep_images_dir`` when requested.

    Parameters
    ----------
    row : astropy.table.Table
        Single-row Table with the photometry measurements.
    fits_file : str
        Path of the FITS image the row was measured from.
    output_file : str
        Path of the single CSV file the row is appended to.
    sector : int or None, optional
        TESS sector number stored in the output row.
    source : str, optional
        Download command or local file path stored in the ``CURL``
        column; used for deduplication.
    source_label : str, optional
        Human-readable source description used in the log message.
    keep_images_dir : str or None, optional
        Directory where the FITS image is copied. If None, no copy is made.
    """
    df = row.to_pandas()
    df["SECTOR"] = sector
    df["CURL"] = source
    _append_row(df, output_file, source)

    print(
        f"[PID {os.getpid()}] "
        f"Star detected in {source_label} appended to {output_file}"
    )

    if keep_images_dir is not None:
        os.makedirs(keep_images_dir, exist_ok=True)
        shutil.copy2(
            fits_file,
            os.path.join(keep_images_dir, os.path.basename(fits_file)),
        )


def _append_row(df: pd.DataFrame, output_file: str, source: str) -> None:
    """
    Append one photometry row to the output CSV under an exclusive lock.

    The header is written only when the file is empty and rows whose
    source (curl command or file path) is already present are skipped.

    Parameters
    ----------
    df : pandas.DataFrame
        Single-row DataFrame with the photometry measurements.
    output_file : str
        Path of the output CSV file.
    source : str
        Download command or local file path of the image; used for
        deduplication.
    """
    with open(output_file, "a") as f:
        if fcntl is not None:
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)

        try:
            if os.path.getsize(output_file) == 0:
                header = True
            else:
                existing = set(
                    pd.read_csv(output_file, usecols=["CURL"])["CURL"].astype(str)
                )
                if source in existing:
                    return
                header = False

            df.to_csv(f, header=header, index=False)
            f.flush()

        finally:
            if fcntl is not None:
                fcntl.flock(f.fileno(), fcntl.LOCK_UN)
