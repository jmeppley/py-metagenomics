#!/usr/bin/env python
"""
"""

import sys, logging, argparse
from edl.blastm8 import filterM8, add_hit_table_arguments, FilterParams
from edl.util import add_universal_arguments, setup_logging

def main():
    description = """
    Take a blast result table and output a subset of hits based on the chosen filtering options. If more than one blast file given, use -O to get multiple output files, otherwise all output data will be concatenated into one output.
    """

# command line arguments
    parser = argparse.ArgumentParser(description=description, conflict_handler='resolve')
    add_hit_table_arguments(parser, flags='all')
    parser.add_argument("-o", "--outfilenome", dest="outfilename", default=None,
                      metavar="OUTFILENAME", help="Write masked fasta output to OUTFILENAME.")
    parser.add_argument('-O', '--autoOutName', default=False,
                      action='store_true',
                      help="Automatically generate output file name from input name and options. Overridden by -o, cannot be used with data from STDIN.")
    parser.add_argument('hit_table', nargs='*', 
            type=argparse.FileType('rU'), default=[sys.stdin,],
            help="Table of search results to be filtered. If absent, \
                    data will be read from STDIN")

    add_universal_arguments(parser)

    arguments = parser.parse_args()

    setup_logging(arguments)

    # check that we have blast file as argument

    # if we're not doing auto file names, wriate all outputs to same file
    if not arguments.autoOutName:
        if arguments.outfilename is not None:
            logging.info("Writing data to %s" % (arguments.outfilename))
            outfile_handle=open(arguments.outfilename,'w')
        else:
            logging.info("writing data to STDOUT")
            outfile_handle=sys.stdout

    # loop over inputs
    for infile_handle in arguments.hit_table:
        logging.info("reading data from %s" % (infile_handle.name))
        if arguments.autoOutName:
            outfile_handle=open(getOutputFile(infile_handle.name,arguments),'w')

        # filter
        params=FilterParams.create_from_arguments(arguments)
        filterM8(infile_handle,outfile_handle,params)

        if arguments.autoOutName:
            outfile_handle.close()
        infile_handle.close()

#############
# Functions #
#############
def getOutputFile(infile, arguments):
    """
    Use the requested arguments to name the output file
    """
    outfile = infile
    if arguments.filterPctid > 0:
        outfile += ".i%g" % arguments.filterPctid
    if arguments.filterLength > 0:
        outfile += ".l%d" % arguments.filterLength
    if arguments.filterBits > 0:
        outfile += ".b%g" % arguments.filterBits
    if arguments.filterEvalue is not None:
        outfile += ".e%g" % arguments.filterEvalue
    if arguments.filterAln is not None and arguments.filterAln>0:
        outfile += ".a%g" % arguments.filterAln
    if arguments.filterHspsPerHit!=1:
        outfile += ".h%d" % arguments.filterHspsPerHit
    if arguments.filterTopPct >= 0:
        outfile += ".p%g" % arguments.filterTopPct
    if arguments.filterNonoverlapping:
        outfile += ".u"
    if arguments.filterHitsPerRead > 0:
        outfile += ".n%d" % arguments.filterHitsPerRead

    if outfile == infile:
        sys.exit("outfile and infile are the same!!\n%s" % infile)

    return outfile

if __name__ == '__main__':
    main()
