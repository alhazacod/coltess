#!/usr/bin/env python3
"""
Analysis tools.
"""

import pandas as pd
import numpy as np
from scipy.signal import find_peaks

from coltess.core import StarData

from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.timeseries import LombScargle
from astropy import units as u
from typing import Tuple, List

from pathlib import Path

def load_photometry_data(
        csv_dir: str,
        target_star: StarData,
        max_sep_arcsec: float = 0.5 
    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load per-frame photometry CSV files and extract a target light curve.
    
    For each frame, the source below a maximum angular separation is selected.
    
    Parameters
    ----------
    csv_dir : str
        Directory containing per-frame CSV photometry files.
    target_ra : float
        Target right ascension in degrees.
    target_dec : float
        Target declination in degrees.
    max_sep_arcsec : float, optional
        Maximum allowed separation for a valid detection.
    
    Returns
    -------
    times : numpy.ndarray
        Observation times in Julian Date.
    fluxes : numpy.ndarray
        Measured fluxes corresponding to the target.
    
    Raises
    ------
    RuntimeError
        If no CSV files are found or the target is not detected in any frame.
    """

    target_ra = target_star.ra 
    target_dec = target_star.dec
    
    target_coord = SkyCoord(target_ra, target_dec, unit=u.deg)
    
    times = []
    fluxes = []
    
    csv_files = sorted(Path(csv_dir).glob("*.csv"))
    
    if not csv_files:
        raise RuntimeError("No CSV files found.")
    
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        
        # Build coordinates for detected sources in this frame
        coords = SkyCoord(df["RA"].values, df["DEC"].values, unit=u.deg)
        
        seps = target_coord.separation(coords).arcsec
        idx = np.argmin(seps)
        
        if seps[idx] > max_sep_arcsec:
            continue  # target not detected in this frame
        
        flux = df.loc[idx, "flux"]
        date_obs = df.loc[idx, "DATE-OBS"]
        
        jd = Time(date_obs, format="isot", scale="utc").jd
        
        fluxes.append(flux)
        times.append(jd)
        
    if not times:
        raise RuntimeError("Target not found in any CSV file.")
    
    return np.array(times), np.array(fluxes)


def compute_periodogram(times: List[float], fluxes: List[float], 
                       min_period: float = 0.1, max_period: float = 10.0) -> dict:
    """Compute Lomb-Scargle periodogram"""
    times_array = np.array(times)
    fluxes_array = np.array(fluxes)
    
    # Remove NaNs
    mask = ~np.isnan(fluxes_array)
    times_array = times_array[mask]
    fluxes_array = fluxes_array[mask]
    
    # Compute periodogram
    frequency = np.linspace(1/max_period, 1/min_period, 20000)
    ls = LombScargle(times_array, fluxes_array)
    power = ls.power(frequency)
    periods = 1/frequency
    
    # Find peaks
    peaks, _ = find_peaks(power)
    sorted_peaks = sorted(peaks, key=lambda x: power[x], reverse=True)
    
    results = {
        'periods': periods,
        'power': power,
        'peaks': sorted_peaks
    }
    
    if len(sorted_peaks) >= 2:
        results['primary_period'] = periods[sorted_peaks[0]]
        results['secondary_period'] = periods[sorted_peaks[1]]
        
        # Estimate uncertainty
        sigma_p = estimate_period_uncertainty(periods, power, sorted_peaks[1])
        results['period_uncertainty'] = sigma_p
    
    return results


def estimate_period_uncertainty(periods: np.ndarray, power: np.ndarray, peak_idx: int) -> float:
    """Estimate period uncertainty using FWHM"""
    half_max = power[peak_idx] / 2
    above_half = np.where(power >= half_max)[0]
    
    if len(above_half) < 2:
        return periods[peak_idx] * 0.01  # Default 1% uncertainty
    
    fwhm = periods[above_half[0]] - periods[above_half[-1]]
    sigma_p = fwhm / 2.35
    
    return sigma_p
