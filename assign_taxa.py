#! /usr/bin/python
"""
A simplified taxon assignment script.
"""

import argparse
import sys, re, logging
from edl.taxon import ranks, TaxNode
from edl.hits import *
from edl.util import *
from edl.expressions import accessionRE, nrOrgRE

def main():
    description = """
Takes a hit table (reads searched against a database) and assigns each read to a taxon. Hit table may be specified with -i or piped to STDIN.

    Notes:

     * Specifying a top score precent (-F) will force hits to be sorted by score within each read. However, it is assumed that the hits in the input table(s) are already grouped by read. This program does not attempt to sort the entire input.
    """
    parser = argparse.ArgumentParser(description)
    parser.add_argument("-i", "--inputfile", dest="infile",
                      metavar="INFILE", help="Read data table from INFILE. (this is here for legacy, you can also just put list file names without this options)"),
    add_IO_arguments(parser)
    parser.add_argument("-T", "--taxids", default=False, action="store_true",
                      help="Output taxids instead of names")
    add_taxon_arguments(parser)
    parser.add_argument("-r", "--rank", dest="rank", default=None,
                      metavar="RANK",
                      help=" Rank to collect counts on.  Defaults to None (whatever the annotation was). Corresponds to rank names in nodes.dmp. To see list run: 'cut -f5 nodes.dmp | uniq | sort | uniq' in ncbi tax dir")
    add_universal_arguments(parser)
    arguments = parser.parse_args()
    setup_logging(arguments)

    logging.debug("Parsing style is: %s" % (arguments.parseStyle))

    # check arguments
    if arguments.taxids and arguments.taxdir is None:
        parser.error("Only use -T when a taxonomy is specified")
    if arguments.rank is not None and arguments.taxdir is None:
        parser.error("Please supply NCBI phylogeny(-n) if specifying a rank(-r).")
    if arguments.rank is not None:
        if arguments.rank == 'domain':
            logging.warn('translating domain to superkingdom')
            arguments.rank='superkingdom'
        if arguments.rank not in ranks:
            parser.error("Unknown rank: %s" % (arguments.rank))

    # load necessary maps
    (taxonomy,valueMap)=readMaps(arguments)

    # loop over inputs
    for (inhandle,outhandle) in inputIterator(arguments):
        logging.debug("Reading from %s and writing to %s" % (inhandle, outhandle))
        hitIter = parseM8FileIter(inhandle, valueMap, arguments.hitTableFormat, arguments.filterTopPct, arguments.parseStyle, arguments.countMethod,taxonomy=taxonomy,rank=arguments.rank)

        ##
        # print output
        # choose output method
        if arguments.taxids:
            printer = taxidPrinter
        else:
            printer = defaultPrinter

        # loop over reads
        outhandle.write("Read\tHit\n")
        for (read,hits) in hitIter:
            outhandle.write(printer(read,hits))

#############
# Functions #
#############
def taxidPrinter(read,hits):
    """
    Convert hits from TaxNodes to taxid ints before printing with default
    """
    if hits is not None:
        newHits = []
        for h in hits:
            if isinstance(h,TaxNode):
                newHits.append(h.id)
            else:
                if isinstance(h,str):
                    logging.debug("Failed to translate: %s" % (h))
                newHits.append(h)
        hits=newHits
    return defaultPrinter(read,hits)

def defaultPrinter(read,hit):
    """
    return the read and all the hit strings separated by tabs
    """
    if isinstance(hit,list):
        hitString = ""
        for h in sorted(hit):
            hitString="%s\t%s" % (hitString,h)
    else:
        hitString="\t%s" % (hit)
    return "%s%s\n" % (read,hitString)

if __name__ == '__main__':
    main()

