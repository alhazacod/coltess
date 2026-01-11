#!/usr/bin/env python3
"""
Download sectors, bash script 
for certain sector and imaes.
"""

import os
import pandas as pd

from coltess.core import StarData

import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.mast import Tesscut

import requests
import subprocess

def get_tess_sectors(target_star: StarData) -> pd.DataFrame:
    """
    Return the list of TESS observing sectors covering a sky position.
    
    Parameters
    ----------
    target_star : StarData
        StarData with right acension and declination.
    
    Returns
    -------
    pandas.DataFrame
        Table describing available TESS sectors for the target position.
    """

    ra = target_star.ra 
    dec = target_star.dec 

    coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
    sectors = Tesscut.get_sectors(coordinates=coord)
    return sectors


def download_tess_sector_script(sector: int) -> str:
    """
    Download the official shell script for a TESS sector from:
    https://archive.stsci.edu/missions/tess/download_scripts/sector/
    
    The script contains curl commands for downloading all FFIs
    associated with the specified sector.

    Parameters
    ----------
    sector : int
        TESS sector number.
    
    Returns
    -------
    str
        Path to the downloaded shell script.
    """
    script_name = f"tesscurl_sector_{sector}_ffic.sh"
    
    if not os.path.exists(script_name):
        url = f"https://archive.stsci.edu/missions/tess/download_scripts/sector/{script_name}"
        response = requests.get(url)
        response.raise_for_status()
        
        with open(script_name, "wb") as f:
            f.write(response.content)
    
    return script_name

def download_tess_image(
        shell_command: str, 
        output_dir: str
        ) -> str:
    """
    Download a single TESS FFI using the passed shell command.

    Parameters
    ----------
    shell_command : str 
        Shell command used to download the FITS file.
    output_dir : str 
        Directory where FITS file will be writter

    Returns 
    -------
    str 
        Path to the downloaded FITS file.
    """
    filename = os.path.basename(shell_command.split()[-1])
    output_path = os.path.join(output_dir, filename)

    cmd = f"{shell_command} -o {output_path}"
    #print(cmd)
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, 
                  stderr=subprocess.DEVNULL, cwd=output_dir)
    
    return output_path


def download_tess_images(
        script_path: str, 
        start_idx: int,
        num_images: int, 
        output_dir: str = "./tess_images"
    ) -> None:
    """
    Download a subset of TESS FFIs listed in a sector shell script.
    
    Parameters
    ----------
    script_path : str
        Path to the TESS download shell script.
    start_idx : int
        Index of the first line to process.
    num_images : int
        Number of images to download.
    output_dir : str, optional
        Directory where FITS files will be written.
    
    Notes
    -----
    This function executes shell commands and suppresses stdout/stderr.
    It is safe to call concurrently for non-overlapping index ranges.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    with open(script_path, "r") as f:
        lines = [line.strip() for line in f]
    
    # Get the requested range
    download_lines = lines[start_idx:start_idx + num_images]
    
    for i, line in enumerate(download_lines, 1):
        download_tess_image(line, output_dir)        
