# Coltess - TESS FFI Photometry Pipeline

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue)](https://www.python.org/downloads/)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL_v3.0-blue)](https://opensource.org/licenses/MIT)

**Coltess** is a lightweight Python package for extracting light curves from TESS (Transiting Exoplanet Survey Satellite) Full Frame Images (FFIs). 

## Features

- 🌟 **Automated catalog generation** from Gaia DR3
- 📥 **Direct FFI downloads** from MAST archive
- 🔭 **Aperture photometry** with local background subtraction
- ⚡ **Parallel processing** for analyzing thousands of images efficiently
- 📊 **Periodogram analysis** using Lomb-Scargle with error weighting and false alarm probability
- 🎯 **Simple API** designed for both interactive and scripted workflows


## Comparison with Lightkurve

| Feature | Coltess | Lightkurve |
|---------|---------|------------|
| Target | Raw FFI photometry | Pre-processed light curves + TPFs |
| Use case | Custom apertures, faint targets | Quick analysis of cataloged targets |
| Data products | DIY light curves | Official SPOC/QLP products |
| Flexibility | Full control | Standardized pipeline |
| Speed | Slower (raw processing) | Faster (pre-computed) |

**When to use Coltess:**
- You need custom aperture sizes
- Your target isn't in the TESS Input Catalog
- You want complete control over the photometry
- You're analyzing very faint sources

**When to use Lightkurve:**
- Your target has official light curves
- You want quick exploratory analysis
- You need TESS pipeline systematics corrections

## Installation

### From PyPI (not published yet)
```bash
pip install coltess
```

### From source
```bash
git clone https://github.com/yourusername/coltess.git
cd coltess
pip install -e .
```

### Dependencies

Core requirements:
- numpy
- pandas
- astropy >= 5.0
- astroquery
- photutils >= 1.5
- scipy
- matplotlib
- requests

## Quick Start

Extract a light curve for Lambda Tau in just a few lines:

**Note for Windows users:** Parallel processing requires `if __name__ == '__main__'` guard in scripts. Windows + Jupyter users should use WSL or process images sequentially.

```python
from coltess import create_catalog, get_tess_sectors, download_tess_sector_script
from coltess import process_images_parallel, load_photometry_data
import matplotlib.pyplot as plt

# 1. Create catalog and get star info
star = create_catalog("lambda tau", radius_arcmin=10.0, output_file="catalog.csv")

# 2. Find available TESS sectors
sectors = get_tess_sectors(star)
sector = int(sectors["sector"][0])

# 3. Download sector script
script_path = download_tess_sector_script(sector)

# 4. Process images in parallel
process_images_parallel(
    script_file=script_path,
    catalog_file="catalog.csv",
    output_dir="photometry_results",
    star=star
)

# 5. Load and plot light curve with error bars
times, fluxes, flux_errors = load_photometry_data("photometry_results", star)

plt.errorbar(times, fluxes, yerr=flux_errors, fmt="o", ms=2, capsize=0, elinewidth=0.5)
plt.xlabel("Julian Date")
plt.ylabel("Flux (e⁻/s)")
plt.title(f"Light Curve: {star.name}")
plt.show()
```

## Usage Examples

Please refer to the examples folder.


## API Reference

### Core Classes

#### `StarData`
Container for star information and photometry data.

**Attributes:**
- `name` (str): Star identifier
- `ra` (float): Right ascension in degrees
- `dec` (float): Declination in degrees
- `gaia_id` (str, optional): Gaia DR3 source ID
- `times` (np.ndarray, optional): Observation times (JD)
- `fluxes` (np.ndarray, optional): Measured fluxes

#### `TessPhotometry`
Main photometry processing class.

**Parameters:**
- `aperture_radius` (int): Aperture radius in pixels (default: 10)
- `annulus_inner` (int): Inner annulus radius (default: 12)
- `annulus_outer` (int): Outer annulus radius (default: 14)
- `zeropoint` (float): Magnitude zeropoint (default: 20.44)

**Methods:**
- `process_image(fits_file, catalog_file, target_star, output_dir)`: Process single FITS image
- `load_catalog(catalog_file)`: Load Gaia catalog from CSV

### Catalog Functions

#### `create_catalog(star_name, radius_arcmin, output_file)`
Create Gaia DR3 catalog centered on target star.

**Returns:** `StarData` object with resolved coordinates

#### `get_star(name)`
Resolve star name via SIMBAD.

**Returns:** `StarData` object

#### `query_gaia_catalog(star, radius_arcmin)`
Query Gaia DR3 around sky position.

**Returns:** pandas DataFrame with sources

### Download Functions

#### `get_tess_sectors(star)`
Find TESS sectors covering target.

**Returns:** pandas DataFrame with sector information

#### `download_tess_sector_script(sector)`
Download official MAST download script for sector.

**Returns:** Path to shell script

#### `download_tess_image(shell_command, output_dir)`
Download single TESS FFI using curl command.

**Returns:** Path to downloaded FITS file

### Analysis Functions

#### `load_photometry_data(csv_dir, target_star, max_sep_arcsec)`
Load light curve from photometry CSV files.

**Returns:** Tuple of (times, fluxes, flux_errors) as numpy arrays

#### `compute_periodogram(times, fluxes, flux_errors=None, min_period=0.1, max_period=10.0, oversampling=5.0)`
Compute Lomb-Scargle periodogram with optional error weighting and false alarm probability.

**Parameters:**
- `times`: Observation times (numpy array)
- `fluxes`: Measured fluxes (numpy array)  
- `flux_errors`: 1-sigma flux uncertainties (optional, numpy array)
- `min_period`: Minimum period to search (days, default: 0.1)
- `max_period`: Maximum period to search (days, default: 10.0)
- `oversampling`: Oversampling factor for frequency grid (default: 5.0)

**Returns:** Dictionary containing:
- `periods`: Array of periods (days)
- `power`: Array of periodogram power
- `frequency`: Array of frequencies (1/days)
- `peak_indices`: List of peak indices sorted by power
- `primary_period`: Best period or None
- `primary_period_uncertainty`: 1-sigma uncertainty or None  
- `primary_fap`: False alarm probability (Baluev 2008) or None
- `secondary_period`: Second best period or None

### Parallel Processing

#### `process_images_parallel(script_file, catalog_file, output_dir, star, start_idx, max_workers)`
Process TESS images in parallel.

Automatically downloads, analyzes, and cleans up temporary files for each image.

## How It Works

### Photometry Pipeline

1. **Catalog Creation**: Queries Gaia DR3 for all sources within a radius around your target
2. **Image Download**: Retrieves TESS FFI files from MAST archive
3. **Source Matching**: Identifies catalog sources within each image's field of view using WCS
4. **Centroid Refinement**: Refines positions using center-of-mass centroiding
5. **Aperture Photometry**: 
   - Measures flux in circular aperture around each source
   - Estimates local background from surrounding annulus
   - Subtracts background and calculates uncertainties
6. **Target Selection**: Matches photometry to target star by position

### Error Propagation

Flux uncertainties are calculated including:
- Poisson noise from the source
- Sky background noise
- Background estimation uncertainty

The formula used is:
```
σ_flux = √(flux/gain + A_ap×σ_sky² + A_ap×σ_sky²/A_ann)
```

Where:
- `gain = 5.22 e⁻/ADU` (TESS nominal)
- `A_ap` = aperture area
- `A_ann` = annulus area
- `σ_sky` = background standard deviation


### Periodogram Analysis

The Lomb-Scargle periodogram is used to detect periodic signals in unevenly sampled time series data. Coltess implements a generalized Lomb-Scargle periodogram that:

- Accepts measurement uncertainties for weighted fitting
- Uses the Baluev (2008) method to compute false alarm probabilities (FAP)
- Estimates period uncertainties from the peak width in frequency space
- Automatically selects an appropriate frequency grid based on the time baseline

The periodogram power is computed using Astropy's LombScargle implementation with options for:
- Weighted fitting when flux errors are provided
- Various normalization schemes (default: 'standard')
- Adaptive frequency grid spacing based on observational baseline

Significant peaks are identified using prominence and distance thresholds to avoid noise fluctuations, and the false alarm probability provides a statistical measure of peak significance.


## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Citation

If you use Coltess in your research, please cite:

```bibtex
@software{coltess2025,
  author = {Manuel Garcia},
  title = {Curves of Light from TESS (COLTESS): photometry tool for TESS FFIs},
  year = {2026},
  url = {https://github.com/alhazacod/coltess}
}
```

Please also cite the relevant TESS papers:
- Ricker et al. 2015 (TESS Mission): [2015JATIS...1a4003R](https://ui.adsabs.harvard.edu/abs/2015JATIS...1a4003R)

And the data sources:
- Gaia DR3: [2023A&A...674A...1G](https://ui.adsabs.harvard.edu/abs/2023A%26A...674A...1G)

## Acknowledgments

This package uses:
- [Astropy](https://www.astropy.org/) for astronomical calculations
- [Photutils](https://photutils.readthedocs.io/) for aperture photometry
- [Astroquery](https://astroquery.readthedocs.io/) for catalog access
- Data from the [TESS mission](https://tess.mit.edu/) and [Gaia DR3](https://www.cosmos.esa.int/gaia)

## License

GPL-3.0-or-later - see LICENSE file for details

## Support

- 📧 Email: mangarciama@unal.edu.co
- 🐛 Issues: [GitHub Issues](https://github.com/alhazacod/coltess/issues)
- 📖 Documentation: (coming soon)

## Roadmap

- [ ] Write documentation.
- [ ] Use logging instead of prints.

---

**Note**: This is scientific research software. Always verify results and report any issues you encounter!
