# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

project = "ProLint"
copyright = "2026, Daniel P. Ramirez-Echemendia and Besian I. Sejdiu"
author = "Daniel P. Ramirez-Echemendia and Besian I. Sejdiu"
version = "2.0.0"
release = "2.0.0"

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
    "sphinx_design",
    "autoapi.extension",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

html_theme = "furo"
html_static_path = ["_static"]
html_title = "ProLint Documentation"
html_favicon = "_static/logo.svg"

# Furo theme options
html_theme_options = {
    # Logo - Furo light/dark variants
    "light_logo": "logo.svg",
    "dark_logo": "logo.svg",

    # Sidebar
    "sidebar_hide_name": False,

    # Navigation
    "navigation_with_keys": True,

    # Footer
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/ProLint/prolint",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],

    # Brand colors - Modern purple palette
    "light_css_variables": {
        "color-brand-primary": "#5C3D8C",
        "color-brand-content": "#5C3D8C",
        "color-admonition-background": "#f8f9fa",
    },
    "dark_css_variables": {
        "color-brand-primary": "#9D7EC9",
        "color-brand-content": "#9D7EC9",
        "color-admonition-background": "#1e1e1e",
    },
}

# Code syntax highlighting
pygments_style = "sphinx"
pygments_dark_style = "monokai"

# Custom CSS (minimal overrides)
html_css_files = ["custom.css"]

# -- MyST Parser configuration -----------------------------------------------

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "tasklist",
]

# Suppress warnings for included files that start with H2 and ambiguous refs
suppress_warnings = ["myst.header", "ref.python"]

# Allow .md files
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# -- Napoleon settings (for NumPy docstrings) --------------------------------

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True

# -- Intersphinx mapping -----------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "MDAnalysis": ("https://docs.mdanalysis.org/stable/", None),
}

# -- AutoAPI configuration ---------------------------------------------------

autoapi_dirs = ["../prolint"]
autoapi_ignore = ["*/web/*", "*/__pycache__/*", "*/tests/*"]
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "imported-members",
]
autoapi_python_class_content = "class"
autoapi_member_order = "bysource"
autoapi_keep_files = True
autoapi_add_toctree_entry = True
