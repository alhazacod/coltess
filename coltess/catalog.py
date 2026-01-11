#!/usr/bin/env python3
"""
Create and query gaiadr3 catalog 
and get star coordinates.
"""

from coltess.core import StarData

import pandas as pd

from astroquery.gaia import Gaia
from astroquery.simbad import Simbad

# Suppress warnings
Gaia.ROW_LIMIT = -1


def get_star(name: str) -> StarData:
    """
    Resolve an object name using SIMBAD and return a 
    StarData with its sky coordinates and id.
    
    Parameters
    ----------
    name : str
        Object identifier resolvable by SIMBAD.
    
    Returns
    -------
    star : StarData 
        Returns the found right ascension and declination, can
        be accessed as star.ra and star.dec respectively.
    
    """
    Simbad.add_votable_fields("ids")
    result = Simbad.query_object(name)
    
    if result is None:
        raise ValueError(f"Object '{name}' not found in SIMBAD")
    
    ra = result['ra'][0]
    dec = result['dec'][0]
    
    # Get Gaia ID
    gaia_id = None
    for id_str in result["ids"][0].split("|"):
        if "Gaia DR3" in id_str:
            gaia_id = id_str.replace("Gaia DR3 ", "").strip()
            break

    star = StarData(
        name=name,
        ra=ra,
        dec=dec,
        gaia_id=gaia_id,
    )
    
    return star


def query_gaia_catalog(star: StarData, radius_arcmin: float = 1.0) -> pd.DataFrame:
    """
    Query the Gaia DR3 catalog around a given sky position.
    
    Parameters
    ----------
    star : StarData 
        Target star.
    radius_arcmin : float, optional
        Search radius in arcminutes.
    
    Returns
    -------
    pandas.DataFrame
        Table containing Gaia sources within the search radius.
    """
    ra = star.ra 
    dec = star.dec

    radius_deg = radius_arcmin / 60.0
    
    query = f"""
    SELECT ra, dec, source_id
    FROM gaiadr3.gaia_source
    WHERE CONTAINS(POINT('ICRS', ra, dec), 
                   CIRCLE('ICRS', {ra}, {dec}, {radius_deg})) = 1
    """
    
    job = Gaia.launch_job_async(query)
    results = job.get_results()
    
    return results.to_pandas()


def create_catalog(star_name: str, radius_arcmin: float = 10.0, 
                  output_file: str = "gaia_catalog.csv") -> StarData:
    """
    Create and save a Gaia source catalog centered on a target star.
    
    Parameters
    ----------
    star_name : str
        Object name resolvable by SIMBAD.
    radius_arcmin : float, optional
        Search radius in arcminutes.
    output_file : str, optional
        Output CSV filename.
    
    Returns
    -------
    star : StarData 
        StarData with name, right acension, declination and gaiadr3 id.
    """

    star = get_star(star_name)
    catalog = query_gaia_catalog(star, radius_arcmin)
    
    catalog.to_csv(output_file, index=False)
    print(f"Catalog saved to {output_file} with {len(catalog)} stars")
    
    return star
