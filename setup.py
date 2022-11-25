from distutils.core import setup
from setuptools import find_packages
from _version import __version__  # noqa

setup(
    name='analytics',
    version=__version__,
    description='Analysis Package for Molecular Dynamics on Conjugated Polymers',
    author='Nicholas Siemons',
    author_email='nicholas9182@gmail.com',
    package_dir={
        'MDAnalysis': '/Users/nicholassiemons/opt/miniconda3/envs/md_analysis/lib/python3.10/site-packages/'
    },
    packages=find_packages(),
    scripts=['cli_tools/plot_hills.py']
)
