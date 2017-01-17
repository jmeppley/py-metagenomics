#!/usr/bin/env python
"""
A simplified taxon assignment script.
"""

import argparse
import sys
import re
import logging
from edl.taxon import ranks, TaxNode
from edl.hits import *
from edl.util import *
from edl.expressions import accessionRE, nrOrgRE
from count_taxa import cleanRanks, formatTaxon

def main():
    description = """
Takes a hit table (reads searched against a database) and assigns
each read to a taxon. Hit table may be specified with -i or piped to STDIN.

    Notes:

     * Specifying a top score precent (-F) will force hits to be sorted
       by score within each read. However, it is assumed that the hits in
       the input table(s) are already grouped by read. This program does
       not attempt to sort the entire input.
    """
    parser = argparse.ArgumentParser(description)
    add_IO_arguments(parser)
    parser.add_argument("-T", "--taxids", default=False, action="store_true",
                        help="Output taxids instead of names")
    add_taxon_arguments(parser)
    parser.add_argument(
        "-r",
        "--rank",
        dest="rank",
        default=None,
        metavar="RANK",
        help=" Rank to collect counts on.  Defaults to None (whatever "
             "the annotation was). Corresponds to rank names in nodes.dmp. "
             "To see list run: 'cut -f5 nodes.dmp | uniq | sort | uniq' in "
             "ncbi tax dir")
    parser.add_argument(
        "-R",
        "--printRank",
        dest="printRanks",
        action="append",
        help="Include indeicated rank(s) in lineage of printed taxa. "
             "Will be ignored if beyond the rank of the taxa "
             "(IE We can't include species if the taxon being counted "
             "is genus)")
    add_universal_arguments(parser)
    arguments = parser.parse_args()
    setup_logging(arguments)

    logging.debug("Parsing style is: %s" % (arguments.parseStyle))

    # Handle the case where Galaxy tries to set None as a string
    arguments.printRanks = checkNoneOption(arguments.printRanks)

    # check arguments
    if arguments.taxids and arguments.taxdir is None:
        parser.error("Only use -T when a taxonomy is specified")
    if arguments.rank is not None and arguments.taxdir is None:
        parser.error(
            "Please supply NCBI phylogeny(-n) if specifying a rank(-r).")
    if arguments.printRanks is not None and arguments.taxdir is None:
        parser.error(
            "Please supply NCBI phylogeny(-n) if specifying a rank(-R).")
    if arguments.rank is not None:
        if arguments.rank == 'domain':
            logging.warn('translating domain to superkingdom')
            arguments.rank = 'superkingdom'
        if arguments.rank not in ranks:
            parser.error("Unknown rank: %s" % (arguments.rank))

    try:
        # Make sure the rank lists make sense
        if arguments.printRanks is not None:
            arguments.printRanks = cleanRanks(arguments.printRanks)
    except Exception as e:
        parser.error(str(e))

    # load necessary maps
    (taxonomy, valueMap) = readMaps(arguments)

    # loop over inputs
    for (inhandle, outhandle) in inputIterator(arguments):
        logging.debug(
            "Reading from %s and writing to %s" %
            (inhandle, outhandle))
        hitIter = parseM8FileIter(
            inhandle,
            valueMap,
            arguments.hitTableFormat,
            arguments.filterTopPct,
            arguments.parseStyle,
            arguments.countMethod,
            taxonomy=taxonomy,
            rank=arguments.rank)

        ##
        # print output
        # choose output method
        if arguments.taxids:
            printer = taxidPrinter
        else:
            if arguments.printRanks is None:
                printer = defaultPrinter
            else:
                def printer(read, hits):
                    return tax_table_printer(read,
                                             hits,
                                             arguments.rank,
                                             arguments.printRanks)

        # loop over reads
        outhandle.write("Read\tHit\n")
        for (read, hits) in hitIter:
            outhandle.write(printer(read, hits))

#############
# Functions #
#############


def taxidPrinter(read, hits):
    """
    Convert hits from TaxNodes to taxid ints before printing with default
    """
    if hits is not None:
        newHits = []
        for h in hits:
            if isinstance(h, TaxNode):
                newHits.append(h.id)
            else:
                if isinstance(h, str):
                    logging.debug("Failed to translate: %s" % (h))
                newHits.append(h)
        hits = newHits
    return defaultPrinter(read, hits)


def defaultPrinter(read, hit):
    """
    return the read and all the hit strings separated by tabs
    """
    if isinstance(hit, list):
        hitString = ""
        for h in sorted(hit):
            hitString = "%s\t%s" % (hitString, h)
    else:
        hitString = "\t%s" % (hit)
    return "%s%s\n" % (read, hitString)


def tax_table_printer(read, hit, leaf_rank, displayed_ranks):
    """
    return a line for each hit and the lineage separated by tabs
    """
    if not isinstance(hit, list):
        hit = [hit,]
    line_string = ""
    for h in sorted(hit):
        line_string += read + "\t"
        line_string += formatTaxon(hit, displayed_ranks, leaf_rank,
                                   delim="\t")
        line_string += "\n"

    return line_string


if __name__ == '__main__':
    main()
