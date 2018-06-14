py-metagenomics: databases
==========================

Instructions for downloading and building local copies of RefSeq, Silva, or KEGG sequence databases and metadata

Overview
--------
In this folder is a makefile for each database. You should be able to simply `make` each database with the command:

    make -f download_XXX.make

If you are impatient and are working on a multi-core computer, you can add "-j N" to the command to use N threads when possible.

You need to have a number of things installed for this to work properly, including:
 
 * gnu make
 * perl
 * lastal/lastdb
 * tantan
 * python with modules: biopython, numpy

The above (except for make and perl, which are standard on most systems) will be installed already if you installe py-metagenomics using conda.

# LASTDBCHUNK #
Unless you have a lot of RAM available (>100GB), you should limit the size of the generated last databases so they fit into memory. To do so, set the make variable: LASTDBCHUNK=50G (or whatever is appropriate for your system).

Larger chunks will be processed faster, but there is a architectural size limit
of about 100GB. The lastal8 and lastdb8 binaries can work around this limit, but
that's not yet in the makefiles.

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

## RefSeq Additions ##
The add_to_refseq.make makefile will create a version of the RefSeq database
with added genomes. Added genomes must be in two files:

  - additions.protein.fasta: fasta file of protein sequences
  - acc.to.taxid.protein.additions: map from sequence ids to taxid

Run with

    make -f add_to_refseq.make ADDITIONS_SOURCE=${ADDITIONS_FOLDER}

...where ADDITIONS_FOLDER is the folder containing the two files mentioned
above.

## HMMs ##
The hmm makefile will download PFAM HMMs by default and can also download TIGR
and COG with, EG:

    make -f download_hmms.make GET_TIGR=True

