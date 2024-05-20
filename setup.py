"""
ProLint2: Lipid-Protein Interaction Analysis
============================================

:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License

ProLint2 calculates distance-based lipid-protein interactions from molecular dynamics trajectories.
"""

import sys
from setuptools import setup, find_packages
import versioneer

short_description = "ProLint2: Lipid-Protein Interaction Analysis.".split("\n")[0]

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest-runner"] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = None

setup(
    # Self-descriptive entries which should always be present
    name="prolint2",
    author="Daniel P. Ramirez & Besian I. Sejdiu",
    author_email="daniel.ramirezecheme@ucalgary.ca",
    description=short_description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="MIT",
    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),
    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,
    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,
    scripts=[
        "bin/prolint2",
    ],
    # Additional entries you may want simply uncomment the lines you want and fill in the data
    url="https://github.com/ProLint/prolint2",  # Website
    install_requires=[
        "pandas",
        "mdanalysis",
        "bottle",
        "jsonargparse[signatures]",
        "requests",
        "networkx",
        "logomaker",
        "seaborn",
    ],  # Required packages, pulls from pip if needed; do not use for Conda deployment
    platforms=[
        "Linux",
        "Mac OS-X",
        "Windows",
    ],  # Valid platforms your code works on, adjust to your flavor
    python_requires=">=3.10",  # Python version restrictions
    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,
)
