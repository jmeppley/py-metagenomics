#!/usr/bin/python

#################
# This script takes any file that can be divided into records and breaks it up
#	 into the given number of chunks or into chunks of the given size
#	 size can be measuered in number of records or by record size
#
# Records can be fasta, fastq, genbank, or something described by a simple RegExp
#       user defined records will be sized on the total contained text
#       pre-defined types will use sequence length
################

from optparse import OptionParser
from os import path
from edl.util import *
from edl.batch import *

import re, sys

def main():
## set up CLI
    description = """
    This script takes any file that can be divided into records and breaks it up
     into the given number of chunks or into chunks of the given size
     size can be measuered in number of records or by record size

     Records can be fasta, fastq, genbank, or something described by a simple RegExp
      user defined records will be sized on the total contained text
      pre-defined types will use sequence length
    """

    parser = OptionParser(description=description)
    addFragmentingOptions(parser)

    parser.add_option("-i", "--infile", dest="infile",
                  metavar="FILE", help="Read data from FILE"),
    parser.add_option("-o", "--outName", dest="outname",
                      metavar="OUTNAME", help="Create fragments with this name, except 'fragNNN' is inserted before any extension. Defaults to input file name (or './stdin')")

    addUniversalOptions(parser)

    (options, cmdargs) = parser.parse_args()

    setupLogging(options, description)

    logging.info("Reading sequences from: " + options.infile)

    # allow input from STDIN only if chunk given
    if options.infile is None and options.chunk is None:
        parser.error("Cannot fragment STDIN unless a chunk size is specified. Please give an input file (-i) or chunk size (-C)")

    # if no output prefix given, use infile name in current directory
    if options.outname is None:
        outdir='.'
        if options.infile is None:
            outpref='stdin'
            outsuff=''
        else:
            (inpath,inname)=path.split(options.infile)
            (outpref,outsuff)=path.splitext(inname)
    else:
        if path.isdir(options.outname):
            outdir=options.outname
            outname=''
            if options.infile is None:
                outname='stdin'
            else:
                outname=path.split(options.infile)[1]
        else:
            (outdir,outname)=path.split(options.outname)
        (outpref,outsuff)=path.splitext(outname)

    if outdir.strip()=="":
        outdir='.'
    logging.debug("outdir: '%s' :: outpref: '%s' :: outsuff: '%s'" % (outdir, outpref, outsuff))

    # fragment!
    num=fragmentInput(options.infile, options, outdir, outpref, suffix=outsuff)
    logging.info("Created %s fragments" % num)

if __name__ == '__main__':
    main()
