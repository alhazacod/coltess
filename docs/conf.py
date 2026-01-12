# Minimal Sphinx configuration for Coltess

import os
import sys

# Add package to path
sys.path.insert(0, os.path.abspath('../..'))

# Project information
project = 'Coltess'
copyright = '2026, Manuel Garcia'
author = 'Manuel Garcia'
release = '0.1.0'

# Extensions
extensions = [
    'sphinx.ext.autodoc',      # Auto-generate API docs from docstrings
    'sphinx.ext.napoleon',     # Support NumPy docstring style
    'sphinx.ext.viewcode',     # Add source code links
    'myst_parser',             # Support Markdown files
]

# MyST (Markdown) settings
myst_enable_extensions = [
    "colon_fence",
]

# Source files
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': True,
}

# Napoleon settings (NumPy docstrings)
napoleon_numpy_docstring = True
napoleon_google_docstring = False

# HTML theme
html_theme = 'sphinx_rtd_theme'

# Don't show source in output
html_show_sourcelink = False
html_copy_source = False
