#!/usr/bin/env python3
"""
Parallel pipeline for TESS photometry.

This module provides functions to download and process large numbers of
TESS Full Frame Images (FFIs) in parallel for a single target star.
Each image is downloaded and analyzed independently then removed after
photometry results are saved to disk.

Parallelization is implemented using ``multiprocessing.Pool``.
"""

import os
import re
import sys
import shutil
import tempfile
import multiprocessing as mp
from functools import partial

import pandas as pd

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
    start_idx: int = 0,
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
    start_idx : int, optional
        Line index in the script file from which to start processing.
        If the output file exists and contains data, processing resumes
        automatically after the last processed image, overriding this
        value.
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

    start_idx = _resume_start_idx(lines, output_file, start_idx)

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

    try:
        pool = mp.get_context("fork").Pool(processes=max_workers, maxtasksperchild=100)
    except ValueError:
        # Windows without WSL - needs __main__ guard
        print(
            "WARNING: Using 'spawn' method. Scripts should use if __name__ == '__main__' or process images sequentially"
        )
        pool = mp.get_context("spawn").Pool(processes=max_workers, maxtasksperchild=100)

    try:
        for _ in pool.imap_unordered(worker, indices):
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


def _resume_start_idx(lines: list[str], output_file: str, start_idx: int) -> int:
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
    start_idx : int
        User-supplied starting line index.

    Returns
    -------
    int
        Line index from which to start processing.
    """
    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
        return start_idx

    try:
        last_curl = pd.read_csv(output_file, usecols=["CURL"]).iloc[-1]["CURL"]
        last_filename = os.path.basename(str(last_curl).split()[-1])

        for i, line in enumerate(lines):
            tokens = line.split()
            if not tokens:
                continue
            if os.path.basename(tokens[-1]) == last_filename:
                resume = i + 1
                if resume > start_idx:
                    print(
                        f"Resuming after image {last_filename} "
                        f"(script line {i + 1}); starting at line {resume}"
                    )
                    return resume
                return start_idx

        print(
            f"WARNING: last processed image {last_filename} not found in "
            f"the script; starting from line {start_idx}"
        )
    except Exception:
        print(
            f"WARNING: could not resume from {output_file}; starting from line {start_idx}"
        )

    return start_idx


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
            df = row.to_pandas()
            df["SECTOR"] = sector
            df["CURL"] = curl_line
            _append_row(df, output_file, curl_line)

            print(
                f"[PID {os.getpid()}] "
                f"Star detected at line {index + 1} appended to {output_file}"
            )

            if keep_images_dir is not None:
                os.makedirs(keep_images_dir, exist_ok=True)
                shutil.copy2(
                    fits_files[0],
                    os.path.join(keep_images_dir, os.path.basename(fits_files[0])),
                )

        return index, success

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _append_row(df: pd.DataFrame, output_file: str, curl_line: str) -> None:
    """
    Append one photometry row to the output CSV under an exclusive lock.

    The header is written only when the file is empty and rows whose curl
    command is already present are skipped.

    Parameters
    ----------
    df : pandas.DataFrame
        Single-row DataFrame with the photometry measurements.
    output_file : str
        Path of the output CSV file.
    curl_line : str
        Curl command used to download the image; used for deduplication.
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
                if curl_line in existing:
                    return
                header = False

            df.to_csv(f, header=header, index=False)
            f.flush()

        finally:
            if fcntl is not None:
                fcntl.flock(f.fileno(), fcntl.LOCK_UN)
