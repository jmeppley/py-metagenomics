#!/usr/bin/env python
"""
This script takes any file that can be divided into records and breaks it up
 into the given number of chunks or into chunks of the given size
 size can be measuered in number of records or by record size

 Records can be fasta, fastq, genbank, or something described by a simple RegExp
  user defined records will be sized on the total contained text
  pre-defined types will use sequence length
"""

from os import path
from edl.util import *
from edl.batch import *

import re, sys, argparse

def main():
## set up CLI
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    add_fragmenting_arguments(parser, defaults={"splits":None})

    parser.add_argument("-i", "--infile", dest="infile",
                  metavar="FILE", help="Read data from FILE"),
    parser.add_argument("-o", "--outName", dest="outname",
                      metavar="OUTNAME", help="Create fragments with this name, except 'fragNNN' is inserted before any extension. Defaults to input file name (or './stdin')")
    parser.add_argument("-Z", "--zero-padding", type=int, default=None, dest="padding",
                    help="Zero pad fragment number in output file names. Leave blank to calculat from requested splits. Falls back to 5 if neither option used.")

    add_universal_arguments(parser)
    arguments = parser.parse_args()
    setup_logging(arguments)

    logging.info("Reading sequences from: " + arguments.infile)

    # allow input from STDIN only if chunk given
    if arguments.infile is None and arguments.chunk is None:
        parser.error("Cannot fragment STDIN unless a chunk size is specified. Please give an input file (-i) or chunk size (-C)")

    # if no output prefix given, use infile name in current directory
    if arguments.outname is None:
        outdir='.'
        if arguments.infile is None:
            outpref='stdin'
            outsuff=''
        else:
            (inpath,inname)=path.split(arguments.infile)
            (outpref,outsuff)=path.splitext(inname)
    else:
        if path.isdir(arguments.outname):
            outdir=arguments.outname
            outname=''
            if arguments.infile is None:
                outname='stdin'
            else:
                outname=path.split(arguments.infile)[1]
        else:
            (outdir,outname)=path.split(arguments.outname)
        (outpref,outsuff)=path.splitext(outname)

    if outdir.strip()=="":
        outdir='.'
    logging.debug("outdir: '%s' :: outpref: '%s' :: outsuff: '%s'" % (outdir, outpref, outsuff))

    # fragment!
    num=fragmentInput(arguments.infile, arguments, outdir, outpref, suffix=outsuff)
    logging.info("Created %s fragments" % num)

if __name__ == '__main__':
    main()
