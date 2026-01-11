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
import sys
import shutil
import tempfile
import multiprocessing as mp
from functools import partial

from coltess.core import StarData
from coltess.photometry import TessPhotometry
from coltess.download import download_tess_images


def process_images_parallel(
    script_file: str,
    catalog_file: str,
    output_dir: str,
    star: StarData,
    start_idx: int = 0,
    max_workers: int | None = None,
):
    """
    Download and process TESS images in parallel for a target star.

    This function coordinates the parallel execution of photometry over
    a list of TESS image URLs or commands contained in a script file.
    Each worker downloads exactly one FITS file, performs photometry,
    writes results to disk, and deletes temporary files.

    Parameters
    ----------
    script_file : str
        Path to a TESS download script (one image per line).
    catalog_file : str
        Path to a Gaia catalog CSV.
    output_dir : str
        Directory where photometry CSV files will be written.
    star : StarData
        Target star information.
    start_idx : int, optional
        Line index in the script file from which to start processing.
        Useful for resuming interrupted runs.
    max_workers : int or None, optional
        Number of parallel worker processes. Defaults to the number of
        available CPU cores.

    Notes
    -----
    - Each FITS file is handled independently.
    - Temporary FITS files are stored in per-worker directories and
      deleted after processing.
    - Pressing ``Ctrl+C`` terminates all workers immediately and exits
      with status code 130.
    """
    if max_workers is None:
        max_workers = mp.cpu_count()

    with open(script_file) as f:
        n_images = sum(1 for _ in f)

    print(f"Processing images {start_idx} -> {n_images - 1}")
    print(f"Using {max_workers} workers")

    indices = list(range(start_idx, n_images))

    worker = partial(
        worker_process_fits,
        script_file,
        catalog_file=catalog_file,
        output_dir=output_dir,
        star=star,
    )

    pool = mp.Pool(processes=max_workers)

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


def worker_process_fits(
    script_file: str,
    index: int,
    catalog_file: str,
    output_dir: str,
    star: StarData,
):
    """
    Process a single TESS FITS image.

    This worker function performs the following steps:
    1. Downloads exactly one FITS file specified by ``index`` in the
       download script.
    2. Runs aperture photometry for the target star.
    3. Writes photometry results to a CSV file.
    4. Deletes all temporary files and directories.

    Parameters
    ----------
    script_file : str
        Path to the TESS download script.
    index : int
        Line index in the script file corresponding to the image to process.
    catalog_file : str
        Path to a Gaia catalog CSV.
    output_dir : str
        Directory where the resulting photometry CSV will be saved.
    star : StarData
        Target star information.

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
    """

    success = False

    print(
        f"[PID {os.getpid()}] "
        f"Processing script line {index + 1}"
    )

    tmp_dir = tempfile.mkdtemp(prefix="tess_")

    try:
        download_tess_images(
            script_file,
            start_idx=index,
            num_images=1,
            output_dir=tmp_dir,
        )

        fits_files = [
            os.path.join(tmp_dir, f)
            for f in os.listdir(tmp_dir)
            if f.endswith(".fits")
        ]

        if fits_files:
            processor = TessPhotometry()
            success = processor.process_image(
                fits_files[0],
                catalog_file=catalog_file,
                target_star = star,
                output_dir=output_dir,
            )

        if success:
            filename = os.path.basename(fits_files[0]).replace(".fits", ".csv")
            output_path = os.path.join(output_dir, filename)
            print(
                f"[PID {os.getpid()}] "
                f"Star detected at line {index + 1} saved at {output_path}"
            )

        return index, success

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
