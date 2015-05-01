py-metagenomics
===============

- Author: John Eppley <jmeppley@hawaii.edu>
- License: GPL v2

About
-----
This is a collection of python (2.7) modules and scripts for metagenomic analyses for use on the command-line, in [http://galaxyproject.org](Galaxy), or in [http://ipython.org/](iPython).

Installation
------------
### Dependencies ###
You must have a number of python modules installed to run all the included scripts. These are listed in requirements.txt. 

### Comand-line ###
py-metagenomics is compatible with python 2.7 only. The scripts will be usable from anywhere. As long as the "edl" folder is in the same location as the scripts, they will be able to find the necessary modules.

It is recommended that you put the root directory of py-metagenomics in your PATH. This way, you won't have to type the full path each time.

### iPython ###
You may want to use the modules yourself from within iPython. If you are working in a [https://virtualenv.pypa.io/en/latest/](virtual environment) or wish to install the modules to the central library, you can use the setup.py script to install the modules:

    $ python setup.py install

Don't forget to install the dependencies listed in requirements.txt!

### Galaxy ###
## Environment ##
To use the galaxy tools, the root directory must be in the execution PATH of jobs spawned by Galaxy. In simple instalations, just add to the PATH of the process that launches Galaxy. 

More complex Galaxy installations have a shell script that runs before every job. This is configured by the "environment_setup_file" parameter in the main galaxy config file. Find this shell script and add the py-metagenomics location to the PATH environemnt variable. NOTE: the path to the environment_setup_file should be absolute.

If you installed the python requirements/dependencies in a virtual environment, you must invoke this in the environment_setup_file, too.

To sum up, the following lines must be in the job execution environment (with the paths edited appropriately):

	export PATH=/data/galaxy/py-metagenomics:$PATH
	export PYMG_PATH=/data/galaxy/py-metagenomics

## Tools ##
The tools are not yet ready to be released as toolshed repositories. In the meantime, they must be installed manually. Many files have paths that assume the galaxy root diretory (usually galaxy-dist) and this repository (py-metagenomics) are in the same parent folder. If this is not the case, you'll have to edit a few extra files.

# Tool definitions #
Add the file galaxy/pymg_tool_conf.xml (using the full path or path 
relative to the galaxy root) to the list of tool config files 
(tool_config_file=) in the main galaxy configuration file config/galaxy.ini. 
Also, if the galaxy root is not in the same parent directory as this 
repository, you'll need to edit the tool_path value set in 
galaxy/pymg_tool_conf.xml.

# Data tables #
First, copy the contents of galaxy/tool-data/tool_data_table_conf.xml (from this repository) into the galaxy data table config file, usually config/tool_data_table_conf.xml (If it doesn't exist yet, copy the .sample version into place and edit from that. The lines to copy are also here:

    <!-- Locations of last(al) databases -->
    <table name="lastdb" comment_char="#">
        <columns>value, name, type, path, taxmapfile, taxkey, kofile, descfile, taxdump</columns>
        <file path="tool-data/lastdb.loc" />
    </table>
    <!-- Information about the chemistries available to MiSeq -->
    <table name="illuminaChemistries" comment_char="#">
        <columns>value, name, primerTemplate</columns>
        <file path="tool-data/illuminaChemistries.loc" />
    </table>

Next, link or copy the .loc.sample files into the tool-data folder in galaxy and drop the .sample extension from the names. You'll need to edit the paths in illuminaChemistries.loc to point to the absolute paths of the primer files.

To use any database searching tools, you'll need to build your own last databases from some or all of the public repositories REfSeq, Silva, or KEGG. See the instructions for downloading and building these in the databases subfolder. Once built, add the locations into lastdb.loc.

# from the tool shed #

To run the included workflow, you'll need to install the following tools from the Galaxy toolshed:

 * sortmerna
 * fastq_to_fasta

Overview
--------


