py-metagenomics: databases
==========================

Instructions for downloading and building local copies of public databases. Downloads are driven by makefiles written in either GNU make or Snakemake.

Overview
--------
In this folder is a makefile for each database. You should be able to simply `make` each database with the commands:

    make -f download_XXX.make SEQDB_ROOT=/path/in/which/to/create/seqdbs

or

    snakemake -s download_XXX.snake --config seqdb_root=/path/in/which/to/create/seqdbs

Two databases (GTDB and PhyloDB) require a second step:

    snakemake -s format_XXX.snake
	
If you are impatient and are working on a multi-core computer, you can add "-j N" to the command to use N threads when possible.

### snake vs make

The following DBs are downloaded using make:

 * RefSeq
 * Silva
 * KEGG
 * HMMs: COG, TIGR, pFAM

These use snakemake:

 * GTDB
 * PhyloDB
 * EGGNOG

## Requirements ##

You need to have a number of things installed for this to work properly, including:
 
 * gnu make or snakemake
 * perl
 * lastal/lastdb
 * tantan
 * python with modules: biopython, numpy

You can install all of these with conda. The included conda/yaml file includes
everything but make and perl which are standard on most unix-like systems. To
set up an environemnt run:

    conda env create -p ./conda.env -f conda.yaml

This will create a folder named conda.env that contains all the dependencies
you need.

Next, run:

    conda activate ./conda.env

This will "activate" the environment and make all the tools accessible to your
command shell. If your system does not have make and/or perl installed, you can
also install those with conda. Once the environment is activated, run:

    conda install -c conda-forge perl make

(You can omit "-c conda-forge" if you have the channel already configured.)

# DB specific options #
## KEGG ##
To download the KEGG database, you'll need a paid account. Sorry. To tell the makefile what your account is, you'll need to set the KEGG_USER and KEGG_PASSWORD variables. You can set these in the file, in your environment, or in the make call:

    make -f download_kegg.make KEGG_USER=myname@myschool.edu KEGG_PASSWORD=mybadpassword

## RefSeq and KEGG ##
You can build a mapping file from RefSeq accessions to KEGG ko ids if you have download the KEGG database link files. These are automatically downloaded with the KEGG db. If you don't want the whold KEGG db, but do want to build the RefSeq to KO map, you can download just the link files with:

	make -f download_kegg.make links

Remember to define KEGG_USER and KEGG_PASSWORD either in the command or in the makefile!

To generate the RefSeq mapping, set BUILD_KO_MAP=True (either in the makefile, the environment, or command line) when making the RefSeq database. EG: 

	make -f download_refseq.make BUILD_KO_MAP=True

One of the nice things about the `make` program is that you can re-run the database build command with the new setting, and only the newly requested files will be processed. It won'e try to download and build the RefSeq sequences again.


## HMMs ##
The hmm makefile will download PFAM HMMs by default and can also download TIGR
and COG with, EG:

    make -f download_hmms.make GET_TIGR=True

# GTDB and PhyloDB #
You can generate a merged GTDB (prokaryotes) and PhyloDB (everything else) db using format_mergedb.snake.
