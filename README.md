# Coltess - TESS FFI Photometry Pipeline

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue)](https://www.python.org/downloads/)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL_v3.0-blue)](https://www.gnu.org/licenses/gpl-3.0.en.html)

**Coltess** is a lightweight Python package for extracting light curves from TESS (Transiting Exoplanet Survey Satellite) Full Frame Images (FFIs). 

## Features

- 🌟 **Automated catalog generation** from Gaia DR3
- 📥 **Direct FFI downloads** from MAST archive
- 🔭 **Aperture photometry** with local background subtraction and centroid refinement
- ⚡ **Parallel processing** for analyzing thousands of images efficiently and **automatic resume** after interruption
- 🖼️ **Local-image analysis** of your own FITS folders or single files (`process_local_images_parallel`, `analyze_image`)
- 📊 **Periodogram analysis** using Lomb-Scargle with error weighting, false alarm probability and period uncertainty
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

### From PyPI
```bash
pip install coltess
```

### From source
```bash
git clone https://github.com/alhazacod/coltess.git
cd coltess
pip install -e .
```

### Dependencies

Core requirements (installed automatically):
- numpy
- pandas
- astropy >= 5.0
- astroquery
- photutils >= 1.5
- scipy
- requests

Development and examples: matplotlib, ipython, black, ruff, mypy (`pip install -e ".[dev]"`)

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

# 4. Process images in parallel, appending results to a single CSV
process_images_parallel(
    script_file=script_path,
    catalog_file="catalog.csv",
    output_file="lambda_tau_photometry.csv",
    star=star
)

# 5. Load and plot light curve with error bars
times, fluxes, flux_errors = load_photometry_data("lambda_tau_photometry.csv", star)

plt.errorbar(times, fluxes, yerr=flux_errors, fmt="o", ms=2, capsize=0, elinewidth=0.5)
plt.xlabel("Julian Date")
plt.ylabel("Flux (e⁻/s)")
plt.title(f"Light Curve: {star.name}")
plt.show()
```

## Usage Examples

Please refer to the examples folder.

### Analyzing local FITS images

For images you already have on disk, use `process_local_images_parallel` (whole folder, in
parallel) or `analyze_image` (single file):

```python
from coltess import process_local_images_parallel, analyze_image, load_photometry_data

# Analyze every FITS file in a folder, appending results to one CSV:
process_local_images_parallel(
    fits_dir="my_fits",
    catalog_file="catalog.csv",
    star=star,
    output_file="photometry.csv",
)

# Analyze a single image:
row = analyze_image("my_fits/image.fits", "catalog.csv", star)

times, fluxes, flux_errors = load_photometry_data("photometry.csv", star)
```


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
σ_flux = √(F/g + n_ap × (1 + n_ap/n_ann) × σ_sky²)
```

Where:
- `F` = background-subtracted flux of the source
- `g` = effective gain: 5.22 e⁻/ADU for ADU images, or the exposure time when the FFI is in e⁻/s units (detected automatically from the header)
- `n_ap` = aperture area (pixels)
- `n_ann` = annulus area (pixels)
- `σ_sky` = background standard deviation (σ-clipped, σ=3)


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
@software{coltess2026,
  author = {Manuel Garcia},
  title = {Curves of Light from TESS (COLTESS): photometry tool for TESS FFIs},
  year = {2026},
  version = {0.1.0},
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
- 📖 Documentation: [API reference](https://alhazacod.github.io/coltess) (generated with pdoc)

---

**Note**: This is scientific research software. Always verify results and report any issues you encounter!
