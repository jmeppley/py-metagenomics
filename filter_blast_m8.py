#!/usr/bin/env python
"""
"""

import sys
import logging
import argparse
from edl.blastm8 import filterM8, add_hit_table_arguments, FilterParams
from edl.util import add_universal_arguments, setup_logging


def main():
    description = """
    Take a blast result table and output a subset of hits based on the
    chosen filtering options. If more than one blast file given, use -O
    to get multiple output files, otherwise all output data will be
    concatenated into one output.
    """

# command line arguments
    parser = argparse.ArgumentParser(
        description=description,
        conflict_handler='resolve')
    add_hit_table_arguments(parser, flags='all')
    parser.add_argument(
        "-o",
        "--outfilenome",
        dest="outfilename",
        default=None,
        metavar="OUTFILENAME",
        help="Write masked fasta output to OUTFILENAME.")
    parser.add_argument(
        '-O',
        '--autoOutName',
        default=False,
        action='store_true',
        help="Automatically generate output file name from input name "
             "and options. Overridden by -o, cannot be used with data "
             "from STDIN.")
    parser.add_argument('-G', '--gff', default=False, action='store_true',
                        help="output GFF format instead of input format")
    parser.add_argument('hit_table', nargs='*',
                        type=argparse.FileType('r'), default=[sys.stdin, ],
                        help="Table of search results to be filtered. "
                             "If absent, data will be read from STDIN")

    add_universal_arguments(parser)

    arguments = parser.parse_args()

    setup_logging(arguments)

    # check that we have blast file as argument

    # if we're not doing auto file names, wriate all outputs to same file
    if not arguments.autoOutName:
        if arguments.outfilename is not None:
            logging.info("Writing data to %s" % (arguments.outfilename))
            outfile_handle = open(arguments.outfilename, 'w')
        else:
            logging.info("writing data to STDOUT")
            outfile_handle = sys.stdout

    if arguments.gff:
        logging.info("Converting to GFF")

    # loop over inputs
    for infile_handle in arguments.hit_table:
        logging.info("reading data from %s" % (infile_handle.name))
        if arguments.autoOutName:
            outfile_handle = open(
                getOutputFile(
                    infile_handle.name,
                    arguments),
                'w')

        # filter
        params = FilterParams.create_from_arguments(arguments)
        filterM8(infile_handle, outfile_handle, params, to_gff=arguments.gff)

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
    if arguments.filter_pctid > 0:
        outfile += ".I%g" % arguments.filter_pctid
    if arguments.filter_length > 0:
        outfile += ".L%d" % arguments.filter_length
    if arguments.filter_bits > 0:
        outfile += ".B%g" % arguments.filter_bits
    if arguments.filter_evalue is not None:
        outfile += ".E%g" % arguments.filter_evalue
    if arguments.filter_aln is not None and arguments.filter_aln > 0:
        outfile += ".N%g" % arguments.filter_aln
    if arguments.filter_hsps_per_hit != 0:
        outfile += ".P%d" % arguments.filter_hsps_per_hit
    if arguments.filter_top_pct >= 0:
        outfile += ".F%g" % arguments.filter_top_pct
    if arguments.filter_nonoverlapping >= 0:
        outfile += ".U"
        if arguments.filter_nonoverlapping >= 1:
            outfile += str(arguments.filter_nonoverlapping)
    if arguments.filter_hits_per_read > 0:
        outfile += ".H%d" % arguments.filter_hits_per_read

    if outfile == infile:
        sys.exit("outfile and infile are the same!!\n%s" % infile)

    logging.info("Writing output to: %s", outfile)
    return outfile


if __name__ == '__main__':
    main()
