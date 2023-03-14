#! /usr/bin/env python
"""
Takes an m8 blast and picks the best hit for each.

There are two count methods available, and these are not listed in the
options below: 'tophit' and 'toporg'. The other count methods (first,
most, etc) are unavailable.

First, only the best scores are used, but if there is a tie (aka ambiguous
hit), than a winner is assigned using the abundances found from unambiguous
hits. When 'tophit' is selected, the hit with the highest overal abundance
(from unambigous reads) is used. When 'toporg' is used, the hit which belongs
to the most abundant taxon is used.

If the -P or --proportinal flag is given, then ambiguous hits are resolved
so that the overall proportion of hit (when using tophit) or taxon (for
toporg) abundance is changed the least.

filter_top_pct defaults to 0, but can be altered, but I don't recommend it.
"""

import argparse
import logging
import os
import sys
from urllib.parse import unquote_plus
from edl import redistribute
from edl.hits import ACCS, FilterParams, add_taxon_arguments, \
        readMaps
from edl.util import add_IO_arguments, add_universal_arguments, \
        inputIterator, setup_logging


def main():
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    add_IO_arguments(parser)
    add_taxon_arguments(
        parser,
        defaults={
            'filter_top_pct': 0,
            'parseStyle': ACCS,
            'countMethod': 'tophit'},
        choices={
            'countMethod': (
                'tophit',
                'toporg')})
    parser.add_argument(
        "-P",
        "--proportional",
        dest="proportional",
        default=False,
        action="store_true",
        help="Assign reads that have multiple equal top hits to taxa such "
             "that the overal proportion of taxa is consistent with the "
             "unambiguious hits. This is meant for use with the 'toporg' "
             "count method.")
    parser.add_argument(
        "-i",
        "--individualFiles",
        dest="individual",
        default=False,
        action="store_true",
        help="Use this flag to process files independently. Normally, "
             "counts from all files are pooled for making choices.")

    add_universal_arguments(parser)
    arguments = parser.parse_args()
    setup_logging(arguments)

    # load necessary maps
    params = FilterParams.create_from_arguments(arguments)
    if arguments.countMethod == 'toporg':
        (taxonomy, hitStringMap) = readMaps(arguments)

    wta = not (arguments.proportional)

    if len(arguments.input_files) <= 1 or arguments.individual:
        # loop over input
        for (inhandle, outhandle) in inputIterator(arguments):
            logging.debug(
                "Reading from %s and writing to %s" %
                (inhandle, outhandle))

            if arguments.countMethod == 'tophit':
                # don't give any taxonomy, just map to accessions for
                # redistribution
                readHits = redistribute.pickBestHitByAbundance(
                    inhandle,
                    filterParams=params,
                    return_lines=True,
                    winnerTakeAll=wta,
                    parseStyle=arguments.parseStyle)
            else:
                # translate to organism before finding most abundant
                readHits = redistribute.pickBestHitByAbundance(
                    inhandle,
                    filterParams=params,
                    return_lines=True,
                    winnerTakeAll=wta,
                    taxonomy=taxonomy,
                    hitStringMap=hitStringMap,
                    parseStyle=arguments.parseStyle)

            for line in readHits:
                outhandle.write(line)

    else:
        # process all files at once
        multifile = redistribute.multipleFileWrapper(arguments.input_files)

        # Build a map from input file name to output handle
        outputMap = {}
        for infile_name in arguments.input_files:
            if arguments.output_file is None:
                outputMap[infile_name] = sys.stdout
            elif len(arguments.input_files) <= 1:
                outputMap[infile_name] = open(arguments.output_file, 'w')
            else:
                # use outfileName as suffix
                if arguments.cwd:
                    # strip path info first
                    (infilePath, infileFile) = os.path.split(infile_name)
                    outfile = "./" + infileFile + arguments.output_file
                else:
                    outfile = infile_name + arguments.output_file
                outputMap[infile_name] = open(outfile, 'w')

        if arguments.countMethod == 'tophit':
            # don't give any taxonomy, just map to accessions for
            # redistribution
            readHits = redistribute.pickBestHitByAbundance(
                multifile,
                filterParams=params,
                return_lines=False,
                winnerTakeAll=wta,
                parseStyle=arguments.parseStyle)
        else:
            # translate to organism before finding most abundant
            readHits = redistribute.pickBestHitByAbundance(
                multifile,
                filterParams=params,
                return_lines=False,
                winnerTakeAll=wta,
                taxonomy=taxonomy,
                hitStringMap=hitStringMap,
                parseStyle=arguments.parseStyle)

        for (read, hit) in readHits:
            infile_name, read = read.split("/", 1)
            outhandle = outputMap[unquote_plus(infile_name)]
            outhandle.write(hit.line.split("/", 1)[1])

        if arguments.output_file is not None:
            for outhandle in outputMap.values():
                outhandle.close()


if __name__ == '__main__':
    main()
