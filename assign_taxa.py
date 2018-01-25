#!/usr/bin/env python
"""
A simplified taxon assignment script.
"""

import argparse
import logging
from edl.taxon import ranks, TaxNode
from edl import hits as edlhits, util
from count_taxa import cleanRanks, formatTaxon


def main():
    """ The CLI """
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
    util.add_IO_arguments(parser)
    parser.add_argument("-T", "--taxids", default=False, action="store_true",
                        help="Output taxids instead of names")
    edlhits.add_taxon_arguments(parser)
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
    parser.add_argument(
        "--no-header",
        dest="no_header",
        default=False,
        action='store_true',
        help="do not write header line")

    util.add_universal_arguments(parser)
    arguments = parser.parse_args()
    util.setup_logging(arguments)

    logging.debug("Parsing style is: %s", arguments.parseStyle)

    # Handle the case where Galaxy tries to set None as a string
    arguments.printRanks = util.checkNoneOption(arguments.printRanks)

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
            logging.warning('translating domain to superkingdom')
            arguments.rank = 'superkingdom'
        if arguments.rank not in ranks:
            parser.error("Unknown rank: %s" % (arguments.rank))

    try:
        # Make sure the rank lists make sense
        if arguments.printRanks is not None:
            arguments.printRanks = cleanRanks(arguments.printRanks)
    except Exception as exc:
        parser.error(str(exc))

    # load necessary maps
    (taxonomy, value_map) = edlhits.readMaps(arguments)

    # loop over inputs
    for (inhandle, outhandle) in util.inputIterator(arguments):
        logging.debug(
            "Reading from %s and writing to %s",
            inhandle, outhandle)
        hit_iter = edlhits.parseM8FileIter(
            inhandle,
            value_map,
            arguments.hitTableFormat,
            arguments.filter_top_pct,
            arguments.parseStyle,
            arguments.countMethod,
            taxonomy=taxonomy,
            rank=arguments.rank)

        ##
        # print output
        # choose output method
        if arguments.taxids:
            hit_header = 'taxid'
            printer = taxid_printer
        else:
            if arguments.printRanks is None:
                hit_header = 'Hit(s)'
                printer = default_printer
            else:
                hit_header = '\t'.join(arguments.printRanks)

                def printer(read, hits):
                    return tax_table_printer(read,
                                             hits,
                                             arguments.rank,
                                             arguments.printRanks)

        # loop over reads
        if not arguments.no_header:
            outhandle.write("Read\t{}\n".format(hit_header))
        for (read, hits) in hit_iter:
            outhandle.write(printer(read, hits))

#############
# Functions #
#############


def taxid_printer(read, hits):
    """
    Convert hits from TaxNodes to taxid ints before printing with default
    """
    if hits is not None:
        new_hits = []
        for h in hits:
            if isinstance(h, TaxNode):
                new_hits.append(h.id)
            else:
                if isinstance(h, str):
                    logging.debug("Failed to translate: %s", h)
                new_hits.append(h)
        hits = new_hits
    return default_printer(read, hits)


def default_printer(read, hit):
    """
    return the read and all the hit strings separated by tabs
    """
    if isinstance(hit, list):
        hit_string = ""
        for h in sorted(hit):
            hit_string = "%s\t%s" % (hit_string, h)
    else:
        hit_string = "\t%s" % (hit)
    return "%s%s\n" % (read, hit_string)


def tax_table_printer(read, hit, leaf_rank, displayed_ranks):
    """
    return a line for each hit and the lineage separated by tabs
    """
    if not isinstance(hit, list):
        hit = [hit, ]
    line_string = ""
    for h in sorted(hit):
        line_string += read + "\t"
        line_string += formatTaxon(h, displayed_ranks, leaf_rank,
                                   delim="\t")
        line_string += '\t'*(1 + len(displayed_ranks) - line_string.count('\t'))
        line_string += "\n"

    return line_string


if __name__ == '__main__':
    main()
