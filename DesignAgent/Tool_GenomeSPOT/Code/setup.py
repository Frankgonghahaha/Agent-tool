#!/usr/bin/env python

import sys

from genome_spot._version import __version__
from setuptools import (
    find_packages,
    setup,
)


if sys.version_info < (3, 8, 16) or sys.version_info >= (3, 12):
    sys.exit("Python version must be >=3.8.16 and <3.12")

setup(
    name="genome_spot",
    version=__version__,
    description="Predict oxygen, temp, salinity, and pH preferences of bacteria and archaea from a genome",
    url="https://github.com/cultivarium/GenomeSPOT",
    author="Tyler Barnum",
    author_email="tyler@cultivarium.org",
    license="MIT License",
    package_data={"genome_spot": ["bioinformatics/hmm/hmm_signal_peptide.joblib"]},
    packages=find_packages(exclude=["tests"]),
    scripts=["genome_spot/genome_spot.py"],
    install_requires=[
        "biopython>=1.83",
        "hmmlearn==0.3.0",
        "scikit-learn==1.2.2",
        "bacdive>=0.2",
        "numpy>=1.23.5",
        "pandas>=1.5.3",
        "pytest>=7.4.3",
    ],
    zip_safe=False,
)
