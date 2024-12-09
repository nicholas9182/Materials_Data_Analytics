from distutils.core import setup
from setuptools import find_packages

setup(
    name='Materials_Data_Analytics',
    version="6.2.0", 
    description='Data analysis package for materials characterisation at Stanford University',
    author='Material Science Stanford',
    author_email='nsiemons@stanford.edu',
    url="https://github.com/nicholas9182/Materials_Data_Analytics/",
    packages=find_packages(),
    install_requires=[
        "scipy",
        "pandas",
        "plotly",
        "matplotlib",
        "typer",
        "click",
        "numpy",
        "torch",
        "networkx",
        "MDAnalysis",
        "dash",
        "kaleido",
        "pyFAI",
        "pygix",
        "holoviews",
        "Datetime",
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
    python_requires=">=3.10",
    entry_points={
        'console_scripts': [
            'launch_dash_app=dash_app.app:run',
        ],
    }
)
