# Configuration file for the Sphinx documentation builder.

# -- Project information -----------------------------------------------------
project = "Nexus-cat"
author = "Julien Perradin"
release = "1.0.1"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",  # Include documentation from docstrings
    "sphinx.ext.napoleon",  # Support for NumPy and Google style docstrings
    "sphinx.ext.viewcode",  # Add links to the source code
    "sphinx.ext.mathjax",  # Render math via JavaScript
    "myst_parser",  # Support for markdown
    "sphinx_rtd_theme",  # Read the Docs theme
    "sphinx.ext.autosummary",  # Generate autodoc summaries
    "sphinx.ext.autosectionlabel",  # Generate unique section labels
    "sphinx.ext.todo",  # Support for todo items
    "sphinx.ext.intersphinx",  # Link to other project's documentation
    "sphinx.ext.extlinks",  # Add external links
    "sphinx.ext.githubpages",  # Publish the documentation on GitHub pages
    "sphinx.ext.doctest",  # Test snippets in the documentation
    "sphinx.ext.inheritance_diagram",  # Generate inheritance diagrams
    "sphinx.ext.graphviz",  # Support for graphviz
    "sphinx.ext.ifconfig",  # Include content based on configuration
    "sphinx_togglebutton",  # Add toggle buttons
    "sphinx_copybutton",  # Add copy buttons
]

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
