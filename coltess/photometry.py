#!/usr/bin/env python3
"""
Photometry analysis for TESS data
"""

import os
import numpy as np
from typing import List, Optional
import pandas as pd
import functools

from coltess.core import StarData

from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.stats import sigma_clipped_stats

from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from photutils.centroids import centroid_sources, centroid_com

import warnings
warnings.filterwarnings(
    "ignore",
    category=FITSFixedWarning
)

class TessPhotometry:
    """
    Aperture photometry pipeline for TESS Full Frame Images (FFIs).

    This class provides utilities to load a Gaia-based source catalog,
    identify catalog sources falling within a TESS image footprint,
    and perform aperture photometry with local background subtraction.

    The class is designed to be safe for batch and parallel processing
    of large TESS datasets.
    """
    
    def __init__(
            self,
            aperture_radius: int = 10, 
            annulus_inner: int = 12,
            annulus_outer: int = 14,
            zeropoint: float = 20.4402281476,
        ):
        """
        Initialize the photometry configuration.

        Parameters
        ----------
        aperture_radius : int, optional
            Radius of the circular photometric aperture in pixels.
        annulus_inner : int, optional
            Inner radius of the background annulus in pixels.
        annulus_outer : int, optional
            Outer radius of the background annulus in pixels.
        zeropoint : float, optional
            Photometric zeropoint used to convert fluxes to instrumental
            magnitudes.

        Notes
        -----
        The gain value (e-/ADU) is fixed to the nominal TESS FFI gain.
        """
        self.aperture_radius = aperture_radius
        self.annulus_inner = annulus_inner
        self.annulus_outer = annulus_outer
        self.zeropoint = zeropoint
        self.epadu = 5.22  # TESS gain

    @functools.lru_cache(maxsize=128)
    def load_catalog(self, catalog_file: str):

        """
        Load a Gaia source catalog from a CSV file.

        Parameters
        ----------
        catalog_file : str
            Path to a CSV file containing at least `ra`, `dec`, and `source_id`
            columns in degrees (ICRS).

        Returns
        -------
        list of tuple
            List of `(ra, dec, source_id)` entries for all catalog sources.
        """
        df = pd.read_csv(catalog_file, dtype={"source_id": str})

        coords = SkyCoord(df["ra"], df["dec"], unit=u.deg)

        catalog = [
            (coords[i].ra.deg, coords[i].dec.deg, df["source_id"].iloc[i])
            for i in range(len(df))
        ]

        return catalog       


    def process_fits(self, fits_path: str, catalog: List[tuple]) -> Optional[Table]:
        """
        Process a single TESS FITS image and perform photometry on catalog sources.

        Parameters
        ----------
        fits_path : str
            Path to the TESS FFI FITS file.
        catalog : list of tuple
            Catalog entries as `(ra, dec, source_id)`.

        Returns
        -------
        astropy.table.Table or None
            Table containing photometric measurements and metadata for all
            detected catalog sources in the frame, or None if no sources fall
            within the image footprint.
        """
        with fits.open(fits_path) as hdul:
            image = hdul[1].data
            header = hdul[1].header
            wcs = WCS(header)
        
        # Filter objects in frame
        objects_in_frame = []
        for ra, dec, obj_id in catalog:
            try:
                x, y = SkyCoord(ra, dec, unit=u.deg).to_pixel(wcs)
                if (0 <= x < image.shape[1] and 0 <= y < image.shape[0]):
                    objects_in_frame.append((ra, dec, obj_id, x, y))
            except Exception as e:
                continue
        
        if not objects_in_frame:
            #print("no object in frame")
            return None
        
        # Extract positions
        positions = [(pos[3], pos[4]) for pos in objects_in_frame]
        
        # Perform photometry
        result_table = self._perform_aperture_photometry(image, positions)
        
        # Add metadata
        result_table['RA'] = [pos[0] for pos in objects_in_frame]
        result_table['DEC'] = [pos[1] for pos in objects_in_frame]
        result_table['ID'] = [pos[2] for pos in objects_in_frame]
        result_table['DATE-OBS'] = header.get('DATE-OBS', '')

        return result_table
    
    def _perform_aperture_photometry(self, image: np.ndarray, 
                                    positions: List[tuple]) -> Table:
        """
        Perform aperture photometry at specified pixel positions.

        Parameters
        ----------
        image : numpy.ndarray
            2D image array extracted from the FITS file.
        positions : list of tuple
            Initial `(x, y)` pixel coordinates for photometry.
        header : dict
            FITS header (currently unused but retained for extensibility).

        Returns
        -------
        astropy.table.Table
            Table containing fluxes, magnitudes, and magnitude uncertainties
            for all successfully measured sources.
        """
        # Centroid refinement
        x_init, y_init = zip(*positions)
        x_cent, y_cent = centroid_sources(image, x_init, y_init, 
                                         box_size=3, centroid_func=centroid_com)
        
        # Remove failed centroids
        valid = ~np.isnan(x_cent)
        positions = [(x_cent[i], y_cent[i]) for i in range(len(x_cent)) if valid[i]]
        
        # Aperture definitions
        aperture = CircularAperture(positions, r=self.aperture_radius)
        annulus = CircularAnnulus(positions, r_in=self.annulus_inner, 
                                 r_out=self.annulus_outer)
        
        # Photometry
        phot_table = aperture_photometry(image, [aperture, annulus])
        
        # Background subtraction
        bkg_mean = phot_table['aperture_sum_1'] / annulus.area
        bkg_sum = bkg_mean * aperture.area
        final_flux = phot_table['aperture_sum_0'] - bkg_sum
        
        # Calculate magnitudes
        magnitudes = self.zeropoint - 2.5 * np.log10(np.abs(final_flux))
        
        # Error estimation
        _, _, std = sigma_clipped_stats(image, sigma=3.0)

        # Flux uncertainty (in electrons)
        flux_uncertainty = np.sqrt(
            np.abs(final_flux) / self.epadu +
            aperture.area * std**2 +
            aperture.area * std**2 / annulus.area
        )

        # Magnitude uncertainty
        mag_error = 1.0857 * flux_uncertainty / np.abs(final_flux)
        
        # Create output table
        result = Table()
        result['flux'] = final_flux
        result['mag'] = magnitudes
        result['mag_err'] = mag_error
        result['flux_err'] = flux_uncertainty
        
        return result
    
    def process_image(
            self,
            fits_file: str,
            catalog_file: str,
            target_star: StarData,
            output_dir: str = "./csv_results",
            max_sep_arcsec: float = 0.5
        ) -> bool:
        """
        Process a single FITS file and extract photometry for a target source.

        This method filters all detected sources in the frame and keeps only
        the source below a maximum angular separation threshold.

        Parameters
        ----------
        fits_file : str
            Path to a single TESS FITS image.
        catalog_file : str
            Path to the Gaia catalog CSV file.
        target_ra : float
            Target right ascension in degrees.
        target_dec : float
            Target declination in degrees.
        output_dir : str, optional
            Directory where the output CSV file will be written.
        max_sep_arcsec : float, optional 
            Maximum angular separation threshold.

        Returns
        -------
        bool
            True if the target source was detected and saved,
            False otherwise.

        Notes
        -----
        This method is safe for parallel execution and is intended to be used
        inside multiprocessing workers.
        """
        
        target_ra = target_star.ra 
        target_dec = target_star.dec

        os.makedirs(output_dir, exist_ok=True)

        catalog = self.load_catalog(catalog_file)


        target_coord = SkyCoord(
            target_ra,
            target_dec,
            unit=u.deg
        )

        try:
            result = self.process_fits(fits_file, catalog)
            if result is None or len(result) == 0:
                return False

            coords = SkyCoord(
                result["RA"],
                result["DEC"],
                unit=u.deg
            )

            separations = target_coord.separation(coords)
            idx = np.argmin(separations)

            # Reject frame if the star is not detected
            if separations[idx].arcsec > max_sep_arcsec:
                print(
                    f"[PID {os.getpid()}] "
                    f"Target not found (sep={separations[idx].arcsec:.2f}\") "
                    f"in {os.path.basename(fits_file)}",
                    flush=True
                )
                return False

            # Keep only if it contains the star
            lambda_tau_row = result[idx:idx+1]

            filename = os.path.basename(fits_file).replace(".fits", ".csv")
            output_path = os.path.join(output_dir, filename)
            lambda_tau_row.write(output_path, overwrite=True)

            print(
                f"[PID {os.getpid()}] "
                f"Saved photometry at {output_path}",
                flush=True
            )

            return True
        except Exception as e:
            print(f"[PID {os.getpid()}] Error processing {fits_file}: {e}")
            return False
