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
    version='0.1.42',
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
        'plotnine>=0.10.1',
        'pandas>=1.1.4',
        'numpy>=1.18.5',
        'pysam>=0.16.0.1',
        'scanpy>=1.6.0',
        'intervaltree>=3.1.0',
        'umap-learn>=0.4.6',
        'regex>=2020.10.15',
        'scipy>=1.4.1',
        'matplotlib>=3.3.2',
        'leidenalg>=0.8.3',
        'numba>=0.52.0',
        'scikit-learn>=1.0.2',
        'bokeh>=2.2.3',
        'seaborn>=0.11.1',
        'statsmodels>=0.12.1',
        'bitarray>=1.6.1',
        'xmltodict>=0.12.0',
        'beautifulsoup4>=4.9.3',
        "parse>=1.18.0",
        "UpSetPlot>=0.4.1",
        "seqfold>=0.7.7",
        'mappy>=2.21',
        'primer3-py>=0.6.1',
        'biopython>=1.79',
        'edlib>=1.3.9',
        'h5py>=2.1.0',
        'psutil>=5.7.0',
        'zarr>=2.11.0',
    ],
    entry_points={
        "console_scripts": [
            "biobook-check_plasmid_with_nanopore_sequencing=biobookshelf.ONT.nanopore_functions:Check_plasmid_with_nanopore_sequencing",
            "biobook-run_guppy_and_combine_output=biobookshelf.ONT.nanopore_functions:Guppy_Run_and_Combine_Output",
            'biobook-Server_Status=biobookshelf.CLI.unclassified_programs:Server_Status',
            'biobook-Stop_a_job=biobookshelf.CLI.unclassified_programs:Stop_a_job',
        ]
    },
)
