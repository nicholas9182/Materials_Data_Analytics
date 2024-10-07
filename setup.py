from distutils.core import setup
from setuptools import find_packages

setup(
    name='Materials_Data_Analytics',
    version="6.1.10", 
    description='Data analysis package for materials characterisation at Stanford University',
    author='Material Science Stanford',
    author_email='nsiemons@stanford.edu',
    url="https://github.com/nicholas9182/Materials_Data_Analytics/",
    packages=find_packages(),
    install_requires=[
        "scipy >= 1.11.2",
	    "pandas >= 2.1.1",
        "plotly >= 5.17.0",
        "matplotlib >= 3.8.0",
        "typer >= 0.9.0",
        "click >= 8.1.7",
        "numpy >= 1.26.0",
        "torch >= 2.2.0",
        "dash >= 2.17.1",
        "networkx >= 3.1",
        "MDAnalysis >= 2.6.1"
    ],
    scripts=[
        'cli_tools/plot_hills.py',
	    'cli_tools/colvar_plotter.py',
	    'cli_tools/get_cv_sample.py',
	    'cli_tools/get_polymer_contacts.py'
    ],
    classifiers=[ 
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.4"
)
