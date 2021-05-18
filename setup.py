"""Setup script for bio-bookshelf"""

import os.path
from setuptools import setup, find_packages

# The directory containing this file
HERE = os.path.abspath(os.path.dirname(__file__))

# The text of the README file
with open(os.path.join(HERE, "README.md")) as fid:
    README = fid.read()

setup(
    name='biobookshelf',
    version='0.1.15',
    author="Hyunsu An",
    author_email="ahs2202@gm.gist.ac.kr",
    description="a collection of python scripts and functions for exploratory analysis of bioinformatic data in Python",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/ahs2202/biobookshelf",
    license="GPLv3",
    packages=find_packages(include=['biobookshelf', 'biobookshelf.*']),
    include_package_data=True,
    install_requires=[
        'ipython>=7.19.0',
        'plotly>=4.11.0',
        'plotnine>=0.7.1',
        'pandas>=1.1.4',
        'numpy>=1.18.5',
        'pysam>=0.16.0.1',
        'scanpy>=1.6.0',
        'intervaltree>=3.1.0',
        'umap-learn>=0.4.6',
        'regex>=2020.10.15',
        'scipy>=1.5.2',
        'matplotlib>=3.3.2',
        'leidenalg>=0.8.3',
        'numba>=0.52.0',
        'scikit-learn>=0.24.1',
        'bokeh>=2.2.3',
        'seaborn>=0.11.1',
        'statsmodels>=0.12.1',
        'bitarray>=1.6.1',
        'xmltodict>=0.12.0',
        'beautifulsoup4>=4.9.3',
        "parse>=1.18.0",
        "UpSetPlot>=0.4.1"
    ],
    entry_points={
        "console_scripts": [
            "biobook-check_plasmid_with_nanopore_sequencing=biobookshelf.ONT.check_plasmid_with_nanopore_sequencing:main",
            "biobook-run_guppy_and_combine_output=biobookshelf.ONT.run_guppy_and_combine_output:main",
        ]
    },
)
