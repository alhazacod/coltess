#!/usr/bin/env python3
"""
Analysis tools: light-curve loading and Lomb-Scargle periodogram analysis.
"""

import os

import pandas as pd
import numpy as np
from scipy.signal import find_peaks
from typing import Optional, Tuple

from coltess.core import StarData

from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.timeseries import LombScargle
from astropy import units as u

from pathlib import Path


def load_photometry_data(
    csv_path: str, target_star: StarData, max_sep_arcsec: float = 0.5
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load photometry data and extract a target light curve.

    ``csv_path`` may be either a single CSV file containing one row per
    frame (appended photometry output) or a directory of per-frame CSV
    files. Sources below a maximum angular separation are selected.

    Parameters
    ----------
    csv_path : str
        Path to a combined photometry CSV file or a directory of
        per-frame CSV files.
    target_star : StarData
        Target star with RA/Dec coordinates.
    max_sep_arcsec : float, optional
        Maximum allowed separation for a valid detection.

    Returns
    -------
    times : numpy.ndarray
        Observation times in Julian Date, sorted chronologically.
    fluxes : numpy.ndarray
        Measured fluxes corresponding to the target.
    flux_errors : numpy.ndarray
        1-sigma flux uncertainties corresponding to the target.

    Raises
    ------
    RuntimeError
        If no photometry data is found or the target is not detected
        in any frame.
    """

    target_ra = target_star.ra
    target_dec = target_star.dec

    target_coord = SkyCoord(target_ra, target_dec, unit=u.deg)

    times = []
    fluxes = []
    flux_errors = []

    if os.path.isdir(csv_path):
        csv_files = sorted(Path(csv_path).glob("*.csv"))

        if not csv_files:
            raise RuntimeError("No CSV files found.")

        for csv_file in csv_files:
            df = pd.read_csv(csv_file)

            coords = SkyCoord(df["RA"].values, df["DEC"].values, unit=u.deg)

            seps = target_coord.separation(coords).arcsec
            idx = np.argmin(seps)

            if seps[idx] > max_sep_arcsec:
                continue

            flux = df.loc[idx, "flux"]
            flux_err = df.loc[idx, "flux_err"]
            date_obs = df.loc[idx, "DATE-OBS"]

            jd = Time(date_obs, format="isot", scale="utc").jd

            fluxes.append(flux)
            flux_errors.append(flux_err)
            times.append(jd)

    else:
        if not os.path.exists(csv_path):
            raise RuntimeError(f"Photometry file not found: {csv_path}")

        df = pd.read_csv(csv_path)

        coords = SkyCoord(df["RA"].values, df["DEC"].values, unit=u.deg)

        seps = target_coord.separation(coords).arcsec

        for i in range(len(df)):
            if seps[i] > max_sep_arcsec:
                continue

            flux = df.loc[i, "flux"]
            flux_err = df.loc[i, "flux_err"]
            date_obs = df.loc[i, "DATE-OBS"]

            jd = Time(date_obs, format="isot", scale="utc").jd

            fluxes.append(flux)
            flux_errors.append(flux_err)
            times.append(jd)

    if not times:
        raise RuntimeError("Target not found in any frame.")

    times_arr = np.array(times)
    fluxes_arr = np.array(fluxes)
    flux_errors_arr = np.array(flux_errors)

    order = np.argsort(times_arr)

    return times_arr[order], fluxes_arr[order], flux_errors_arr[order]


def compute_periodogram(
    times: np.ndarray,
    fluxes: np.ndarray,
    flux_errors: Optional[np.ndarray] = None,
    min_period: float = 0.1,
    max_period: float = 10.0,
    oversampling: float = 5.0,
) -> dict:
    """Compute Lomb-Scargle periodogram with weighted fitting and FAP.

    Parameters
    ----------
    times : np.ndarray
        Observation times (e.g., Julian Date).
    fluxes : np.ndarray
        Measured fluxes.
    flux_errors : np.ndarray, optional
        1-sigma flux uncertainties. If provided, passed to LombScargle
        as `dy` for weighted periodogram (generalized LS).
    min_period : float, optional
        Minimum period to search (days).
    max_period : float, optional
        Maximum period to search (days).
    oversampling : float, optional
        Oversampling factor relative to the frequency resolution 1/T.
        Default 5 follows VanderPlas (2018) recommendation.

    Returns
    -------
    dict
        Dictionary containing:
        - 'periods': array of periods (days)
        - 'power': array of periodogram power
        - 'frequency': array of frequencies (1/days)
        - 'peak_indices': list of peak indices sorted by power
        - 'primary_period': best period or None
        - 'primary_period_uncertainty': 1-sigma uncertainty or None
        - 'primary_fap': False alarm probability (Baluev 2008) or None
        - 'secondary_period': second best period or None
    """
    times = np.asarray(times)
    fluxes = np.asarray(fluxes)

    mask = ~np.isnan(times) & ~np.isnan(fluxes)
    if flux_errors is not None:
        flux_errors = np.asarray(flux_errors)
        mask &= ~np.isnan(flux_errors)

    times = times[mask]
    fluxes = fluxes[mask]
    if flux_errors is not None:
        flux_errors = flux_errors[mask]

    if len(times) < 3:
        return {
            "periods": np.array([]),
            "power": np.array([]),
            "frequency": np.array([]),
            "peak_indices": [],
            "primary_period": None,
            "primary_period_uncertainty": None,
            "primary_fap": None,
            "secondary_period": None,
        }

    baseline = times.max() - times.min()
    df = 1.0 / baseline / oversampling
    f_min = 1.0 / max_period
    f_max = 1.0 / min_period
    frequency = np.arange(f_min, f_max, df)
    periods = 1.0 / frequency

    ls = LombScargle(times, fluxes, dy=flux_errors, normalization="standard")
    power = ls.power(frequency)

    # Peak finding with prominence and distance thresholds
    power_std = np.std(power)
    min_distance = max(1, int(1.0 / df))  # at least one independent resolution element
    peaks, _ = find_peaks(power, prominence=3.0 * power_std, distance=min_distance)
    sorted_peaks = sorted(peaks, key=lambda x: power[x], reverse=True)

    results = {
        "periods": periods,
        "power": power,
        "frequency": frequency,
        "peak_indices": sorted_peaks,
        "primary_period": None,
        "primary_period_uncertainty": None,
        "primary_fap": None,
        "secondary_period": None,
    }

    if len(sorted_peaks) >= 1:
        primary_idx = sorted_peaks[0]
        results["primary_period"] = float(periods[primary_idx])
        results["primary_period_uncertainty"] = estimate_period_uncertainty(
            frequency, power, primary_idx, baseline
        )
        # Baluev (2008) analytic FAP
        fap = ls.false_alarm_probability(
            power[primary_idx],
            method="baluev",
            minimum_frequency=frequency.min(),
            maximum_frequency=frequency.max(),
        )
        results["primary_fap"] = float(fap)

    if len(sorted_peaks) >= 2:
        results["secondary_period"] = float(periods[sorted_peaks[1]])

    return results


def estimate_period_uncertainty(
    frequency: np.ndarray,
    power: np.ndarray,
    peak_idx: int,
    baseline_T: float,
) -> float:
    """Estimate period uncertainty from FWHM in frequency space.

    The FWHM is measured locally around the peak in frequency space,
    then propagated to period space via sigma_P = P^2 * sigma_f.

    References
    ----------
    VanderPlas 2018, ApJ Suppl., §7.2
    """
    peak_power = power[peak_idx]
    baseline = np.median(power)
    half_max = baseline + (peak_power - baseline) / 2.0

    # Search locally for left/right half-max crossings
    left_slice = power[:peak_idx]
    right_slice = power[peak_idx:]

    left_above = np.where(left_slice >= half_max)[0]
    right_above = np.where(right_slice >= half_max)[0]

    if len(left_above) == 0 or len(right_above) == 0:
        # Fallback: frequency resolution limit sigma_f = 1/T
        # Propagate: sigma_P = P^2 * sigma_f = P^2 / T
        P = 1.0 / frequency[peak_idx]
        return P**2 / baseline_T

    left_idx = left_above[-1]  # nearest to peak on left
    right_idx = peak_idx + right_above[0]  # nearest to peak on right

    fwhm_f = frequency[right_idx] - frequency[left_idx]
    # Lorentzian-like peak: sigma_f ~ FWHM / 2
    sigma_f = fwhm_f / 2.0

    P = 1.0 / frequency[peak_idx]
    sigma_P = P**2 * sigma_f

    return float(sigma_P)
