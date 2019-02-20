py-metagenomics
===============

- Author: John Eppley <jmeppley@hawaii.edu>
- License: GPL v2

About
-----
This is a collection of python (3.5) modules and scripts for metagenomic analyses for use on the command-line or in [http://ipython.org/](iPython).

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
### Dependencies ###
You must have a number of python modules installed to run all the included scripts. These are listed in requirements.txt. 

The most critical are biopython and numpy. Pandas and matplotlib are only used in some obscure corners of the modules, but they are invaluable if you are planning to use pymg within ipython.

### Comand-line ###
This version of py-metagenomics is compatible with python 3 or greater. There is a python 2 branch, with more features, but less active development. 

The scripts will be usable from anywhere. As long as the "edl" folder is in the same location as the scripts, they will be able to find the necessary modules.

It is recommended that you put the root directory of py-metagenomics in your PATH. This way, you won't have to type the full path each time.

### iPython ###
You may want to use the modules yourself from within iPython. If you are working in a [https://virtualenv.pypa.io/en/latest/](virtual environment) or wish to install the modules to the central library, you can use the setup.py script to install the modules into your python environment:

    $ python setup.py install

Don't forget to install the dependencies listed in requirements.txt!

Overview
--------



Known Issues
------------
