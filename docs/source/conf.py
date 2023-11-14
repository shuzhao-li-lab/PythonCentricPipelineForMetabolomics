# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys, os

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath('../..'))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PythonCentricPipelineForMetabolomics'
copyright = '2023, Joshua M. Mitchell, Yuanye Chi, Maheshwor Thapa, Shujian Zheng, Shuzhao Li'
author = 'Joshua M. Mitchell, Yuanye Chi, Maheshwor Thapa, Shujian Zheng, Shuzhao Li'
release = '0.0.9'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.imgconverter',
    'myst_parser',
]

napoleon_custom_sections = [('Updates','params_style'),('Input','params_style'), ('Outputs','rubric')]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

source_suffix = ['.rst', '.md']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
# if on_rtd:
#     html_theme = 'default'
# else:
#     # html_theme = 'default'
#     html_theme = 'sphinx_rtd_theme'

html_theme = 'sphinx_rtd_theme'