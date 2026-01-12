# Ultra-minimal Sphinx configuration
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# Project info
project = 'Coltess'
copyright = '2025, Manuel Garcia'
author = 'Manuel Garcia'
release = '0.1.0'

# Only essential extensions
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
]

# Napoleon for NumPy docstrings
napoleon_numpy_docstring = True
napoleon_google_docstring = False

# Theme
html_theme = 'sphinx_rtd_theme'
