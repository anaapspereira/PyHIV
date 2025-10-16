# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))

project = 'PyHIV'
author = 'Ana Santos-Pereira; Joao Correia'
release = '0.0.3'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',       # for Google/NumPy docstrings
    'sphinx.ext.viewcode',       # links to source code
    'sphinx.ext.autosummary',    # auto-generate module summaries
    'myst_parser',               # Markdown support if needed
]

templates_path = ['_templates']
exclude_patterns = []

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}
autodoc_typehints = "description"

# Autosummary settings
autosummary_generate = True