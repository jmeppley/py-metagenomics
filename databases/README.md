py-metagenomics: databases
==========================

Instructions for downloading and building local copies of RefSeq, Silva, or KEGG sequence databases and metadata

Overview
--------
In this folder is a makefile for each database. You should be able to simply `make` each database with the command:

    make -f download_XXX.make

If you are impatient and are working on a multi-core computer, you can add "-j N" to the command to use N threads when possible.

KEGG
----
To download the KEGG database, you'll need a paid account. Sorry.


