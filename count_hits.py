#! /usr/bin/env python
"""
Count hits in a tabular blast output. By default, first hit for each read
is used.
"""

import sys
import re
import logging
import argparse
from edl.hits import add_count_arguments, getAllMethod, add_weight_arguments, applyFractionalCutoff
from edl.util import add_universal_arguments, setup_logging

# a habit that stuck in Perl
die = sys.exit


def main():
    # set up CLI
    description = __doc__

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--infile", dest="infile",
                        metavar="FILE", help="Read raw table from INFILE")
    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        metavar="OUTFILE",
        help="Write collapsed table to OUTFILE")
    parser.add_argument("-d", "--delim", dest="delim", default="\t",
                        help="Input table delimiter", metavar="DELIM")
    parser.add_argument("-D", "--delimOut", dest="delimOut", default="\t",
                        help="Output table delimiter", metavar="DELIM")
    parser.add_argument(
        '-F',
        '--countFirst',
        action='store_true',
        default=False,
        help="Don't skip the first line, it's NOT a header")
    parser.add_argument(
        "-R",
        "--readColumn",
        dest="readCol",
        type=int,
        default=0,
        help="Index (starting at 0) of column with read name, 0 is default",
        metavar="READCOL")
    parser.add_argument(
        "-H",
        "--hitColumn",
        dest="hitCol",
        type=int,
        default=2,
        help="Index (starting at 0) of column with hit name (for counting), "
             "2 is default, if less than zero, all (non-read) columns will "
             "be used as multiple hits",
        metavar="HITCOL")
    parser.add_argument(
        '-s',
        '--hitSep',
        default=None,
        help="Use this string to split multiple values in single hit cell. "
             "Default is 'None' to leave hits as is, use 'eval' to parse "
             "as python repr strings")
    add_weight_arguments(parser, multiple=False)
    parser.add_argument("-T", "--total", default=False, action="store_true",
                        help="Report 'Total' in the first row")

    # cutoff options
    add_count_arguments(parser, {'cutoff': 0})

    # log level and help
    add_universal_arguments(parser)
    arguments = parser.parse_args()
    setup_logging(arguments)

    # make sure we have something to do
    if (arguments.infile is None):
        logging.info("Reading table from: STDIN")
    else:
        logging.info("Reading table from: " + arguments.infile)

    if (arguments.outfile is None):
        logging.info("Writing counts to: STDOUT")
    else:
        logging.info("Writing counts to: " + arguments.outfile)

    # process arguments
    takeFirst = (arguments.allMethod == 'first')
    splitHits = (arguments.hitSep is not None and arguments.hitSep != 'None')
    uncluster = (arguments.weights is not None)

    if arguments.hitSep == 'eval':
        parser.error("Sorry, parsing with eval is not yet supported!")

    # inform the curious user
    logging.info("Delimiter: '" + arguments.delim)
    logging.info("Read names in col: '" + str(arguments.readCol))
    logging.info("Hit names in col: '" + str(arguments.hitCol))
    if splitHits:
        logging.info("Splitting hits with: %s" % (arguments.hitSep))
        logging.warn(
            "Splitting hits has not been tested yet! Let me know how it goes.")
    if takeFirst:
        logging.info("Taking first hit for each read.")
    else:
        if arguments.allMethod == 'portion':
            logging.info("Dividing count among all hits for each read.")
        else:
            logging.info("Adding 1 to every hit for each read")
    if uncluster:
        logging.info(
            "Getting read cluster sizes from: %s" %
            (arguments.weights))
    if arguments.countFirst:
        logging.info("First line is data")
    else:
        logging.info("Skipping first line")

    # Do the counting!
    counts = {}
    countHitsForRead = getAllMethod(arguments.allMethod)

    clusteredReadCounts = {}
    if uncluster:
        clusteredReadCounts = parseMapFile(
            arguments.clusterFile, valueType=int)

    currentRead = ''
    readCount = 1
    hits = []

    if arguments.infile is None:
        infile = sys.stdin
    else:
        infile = open(arguments.infile)

    # loop over lines
    if not arguments.countFirst:
        # skip first line
        try:
            next(infile)
        except StopIteration:
            raise Exception("No lines in %s" % str(infile))

    for line in infile:
        line = line.rstrip('\r\n')
        rowcells = line.split(arguments.delim)
        # get read
        read = rowcells[arguments.readCol]

        # if it's a new read, process previous read
        if currentRead == '':
            currentRead = read
        elif read != currentRead and currentRead != '':
            readCount += 1
            logging.info("Checking hits for %s" % currentRead)

            # was it part of a cluster?
            multiplier = 1
            if uncluster:
                multiplier = clusteredReadCounts[currentRead]

            # where does the count for this read go
            countHitsForRead(hits, counts, multiplier=multiplier)

            hits = []
            currentRead = read

        # get hit from this line
        if arguments.hitCol >= 0:
            hit = rowcells[arguments.hitCol]
            if splitHits:
                hits.extend(hit.split(arguments.hitSep))
            else:
                hits.append(hit)
        else:
            rowcells.pop(arguments.readCol)
            hits.extend(rowcells)

    # check last read!
    logging.info("Checking hits for %s" % currentRead)
    # was it part of a cluster?
    multiplier = 1
    if uncluster:
        multiplier = clusteredReadCounts[currentRead]
    # where does the count for this read go
    countHitsForRead(hits, counts, multiplier=multiplier)

    # apply cutoff
    if arguments.cutoff > 0:
        applyFractionalCutoff(counts, threshold=arguments.cutoff * readCount)

    # print output
    if arguments.outfile is None:
        outhandle = sys.stdout
    else:
        outhandle = open(arguments.outfile, 'w')

    if arguments.total:
        outhandle.write("Total%s%d\n" % (arguments.delimOut, readCount))

    if arguments.allMethod == 'portion':
        outFmtString = "%s%s%f\n"
    else:
        outFmtString = "%s%s%d\n"

    delimRE = re.compile(arguments.delimOut)
    for hit in sorted(counts.keys()):
        count = counts[hit]
        hit = delimRE.sub('_', hit)
        outhandle.write(outFmtString % (hit, arguments.delimOut, count))

if __name__ == '__main__':
    main()
