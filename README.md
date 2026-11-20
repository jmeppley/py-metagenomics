py-metagenomics
===============

- Author: John Eppley <jmeppley@hawaii.edu>
- License: MIT

About
-----
This is a collection of python (3.10+) modules and scripts for metagenomic analyses for use on the command-line or in [http://ipython.org/](iPython).

The primary scripts are useful for filtering sequence search results from common bioinformatics tools (blast, last, hmmer, etc) by hit metadata.

Additionally, there are tools to turn hits against common reference databases (EG RefSeq or KEGF) into taxonomic or functional annotations.

For use in Galaxy, look at the python2 (py2) branch. I hope to have a separate galaxy repository set up eventually, but I don't have time for that at the moment. Please feel free to ask for help (jmeppley@hawaii.edu) if you would like to try it for yourself.

Differences from the Python2 version
------------------------------------
A number of scripts and tools have been dropped in the move to python3. I've tried to limit the py3 version to the most generally useful tools, so I dropped many things were too specialized including:

 * The galaxy tools
 * plotting tools
 * QC tools specific to our lab's inhouse workflow
 * lastWrapper (because lastal added most of the needed features after v693)

I also consolodated and simplified the assignXXX scripts, so some features will be missing.

Installation
------------
### Conda ###
The easiest way to get py-metagenomics (with the scripts on your PATH and the external tools lastal, tantan and hmmer) is the conda package from the jmeppley channel:

    $ conda install -c jmeppley -c conda-forge -c bioconda py-metagenomics

The database download/formatting workflows (the `databases` folder) are installed to `$CONDA_PREFIX/share/py-metagenomics/databases`.

### Dependencies ###
You must have a number of python modules installed to run all the included scripts. These are listed in requirements.yaml (a conda environment file). 

The most critical are biopython and numpy. Pandas and matplotlib are only used in some obscure corners of the modules, but they are invaluable if you are planning to use pymg within ipython.

### Comand-line ###
This version of py-metagenomics requires python 3.10 or greater and has been tested with python 3.10 through 3.14, numpy 2.2-2.5, pandas 2.3-3.1 and biopython 1.88. There is a python 2 branch, with more features, but less active development. 

The scripts will be usable from anywhere. As long as the "edl" folder is in the same location as the scripts, they will be able to find the necessary modules.

It is recommended that you put the root directory of py-metagenomics in your PATH. This way, you won't have to type the full path each time.

### iPython ###
You may want to use the modules yourself from within iPython. If you are working in a [https://virtualenv.pypa.io/en/latest/](virtual environment) or wish to install the modules to the central library, you can use pip to install the modules and scripts into your python environment:

    $ pip install .

This will also install the core dependencies (numpy, pandas, and biopython). Add matplotlib (`pip install .[plots]`) if you want the plotting functions in edl.assembly.

### Releasing ###
The version number is set in one place: `__version__` in `edl/__init__.py` (pyproject.toml reads it from there). To make a release, update that value, commit, and tag the commit with the same version prefixed by "v" (e.g. `v1.4.2`).

The conda recipe is in the `conda` folder and reads the version from the same place. To build and upload the conda package for a release:

    $ git checkout v1.4.2
    $ conda build -c conda-forge -c bioconda conda
    $ anaconda upload --user jmeppley $(conda build -c conda-forge -c bioconda conda --output)

Only bump `build: number` in conda/meta.yaml when rebuilding a version that has already been uploaded.

New command-line scripts must also be added to the `script-files` list in pyproject.toml (a test checks this).

### Testing ###
The unit tests use pytest and the command-line tests use bats. From the root of the repository:

    $ python -m pytest test
    $ bats test/bats

Overview
--------



Known Issues
------------
