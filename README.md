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

More complex Galaxy installations have a shell script that runs before every job. This is configured by the "environment_setup_file" parameter in the main galaxy config file. Find this shell script and add the py-metagenomics location to the PATH environemnt variable. 

If you installed the python requirements/dependencies in a virtual environment, you must invoke this in the environment_setup_file, too.

## Tools ##
The tools are not yet ready to be released as toolshed repositories. In the meantime, they must be installed manually. 

First, copy the galaxy/pymg folder into the tools subfolder of your Galaxy installation. 

Second, add the following XML snippet to the tool_conf.xml file:

```
  <section id="pymg" name="py-metagenomics">
    <tool file="pymg/prinseq_lite.xml" />
    <tool file="pymg/readLengths.xml" />
    <tool file="pymg/screen_list.xml" />
    <tool file="pymg/createPrimerFile.xml" />
    <tool file="pymg/illuminaPrep.xml" />
    <tool file="pymg/last_wrapper.xml" />
    <tool file="pymg/filter_blast.xml" />
    <tool file="pymg/translate_column.xml" />
    <tool file="pymg/assignTaxa.xml" />
    <tool file="pymg/assignPaths.xml" />
    <tool file="pymg/countHits.xml" />
    <tool file="pymg/countTaxa.xml" />
    <tool file="pymg/countTaxaSimple.xml" />
    <tool file="pymg/countPaths.xml" />
    <tool file="pymg/sunburstMeganCSV.xml" />
  </section>
```

Thrid, set up the lastdb filetype. This is a little involved (sorry, the toolshed version will fix this)

Normally you would install datatypes via a ToolShed, which would move
the provided lastdb.py file into a suitable location and process the
datatypes_conf.xml entry to be combined with your local configuration. If you have a private toolshed running, you can upload the galaxy/last_dataypes folder into a toolshed repository and use it from there.

However, if you don't have a private toolshed already set up, it's much easier to install the datatypes manually. (1) Add
the following lines from galaxy/last_datatypes/dataypes_conf.xml in this package to the datatypes_conf.xml file in the Galaxy main folder::

        <datatype extension="lastdbn" type="galaxy.datatypes.lastdb:LastalNucDb" mimetype="text/html" display_in_upload="false"/>
        <datatype extension="lastdbp" type="galaxy.datatypes.lastdb:LastalProtDb" mimetype="text/html" display_in_upload="false"/>

(2) Also create the file lib/galaxy/datatypes/lastdb.py by moving, 
copying or linking the lastdb.py file in galaxy/last_datatypes. 
(3) Finally add 'import lastdb' near the start of file 
lib/galaxy/datatypes/registry.py (after the other import
lines).

Overview
--------


