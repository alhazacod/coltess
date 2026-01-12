import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
project = 'Coltess'
copyright = '2026, Manuel Garcia'
author = 'Manuel Garcia'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
]
templates_path = []
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
html_theme = 'classic'  
html_static_path = []  
