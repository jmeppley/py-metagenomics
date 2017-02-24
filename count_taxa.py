#! /usr/bin/env python
"""
Takes m8 blast files and generates a table of taxon hit counts for the
given rank. Columns are input files and rows are taxa. If multiple ranks
given (the default), multiple output files are produced, each with the
rank name appended to the output file name.
"""
import sys
import re
import argparse
import logging
from urllib.parse import unquote_plus
from edl.taxon import *
from edl.hits import *
from edl.util import add_universal_arguments, setup_logging, parseMapFile, \
        checkNoneOption, parseMapFile, passThrough
from edl.expressions import accessionRE, nrOrgRE

ORG_RANK = 'organism'


def main():
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("input_files", nargs="+",
                        default=[],
                        metavar="INFILE",
                        help="List of hit tables to process")
    parser.add_argument("-o", "--outfile", dest="outfile",
                        metavar="OUTFILE",
                        help="Write count table to OUTFILE")
    parser.add_argument("-r", "--rank", dest="ranks", default=None,
                        metavar="RANK", action="append",
                        help=""" Rank(s) to collect counts on. Use flag
                        multiple
                        times to specify multiple ranks. If multiple values
                        given, one table produced for each with rank name
                        appended to file name. Defaults to all major ranks
                        between phylum and species. Corresponds to rank names
                        in nodes.dmp. To see list run:
                        'cut -f5 nodes.dmp | uniq | sort | uniq'
                        in ncbi tax dir. Will also accept 'organism' to mean
                        no rank (ie, just the organism name).""")
    parser.add_argument(
        "-s",
        "--collapseToDomain",
        default=False,
        action="store_true",
        help="Collapse all taxa below given rank down to "
             "superkingdom/domain. EG: in the genus output, anything "
             "assigned to Cyanobactia, will be lumped in with all "
             "other bacteria")
    parser.add_argument(
        "-R",
        "--printRank",
        dest="printRanks",
        action="append",
        help="Include indeicated rank(s) in lineage of printed taxa. "
             "Will be ignored if beyond the rank of the taxa "
             "(IE We can't include species if the taxon being counted "
             "is genus)")

    # option for deconvoluting clusters or assemblies
    add_weight_arguments(parser, multiple=True)

    # cutoff options
    add_count_arguments(parser)

    # format, tax dir, and more
    add_taxon_arguments(
        parser,
        choices={
            'countMethod': (
                'LCA',
                'all',
                'first',
                'most',
                'tophit',
                'toporg',
                'consensus')})

    # log level and help
    add_universal_arguments(parser)
    arguments = parser.parse_args()
    setup_logging(arguments)

    if len(arguments.input_files) == 0:
        parser.error("Must supply at least one m8 file to parse")

    # Handle the case where Galaxy tries to set None as a string
    arguments.ranks = checkNoneOption(arguments.ranks)
    arguments.printRanks = checkNoneOption(arguments.printRanks)

    logging.info("Printing out ranks: {}".format(arguments.ranks))

    # Set defaults and check for some conflicts
    if arguments.ranks is None and arguments.taxdir is None:
        # using hit names only
        arguments.ranks = [ORG_RANK]
        if arguments.printRanks is not None:
            parser.error("Display ranks are not used without taxonomic info")
    else:
        if arguments.taxdir is None:
            parser.error("Cannot select ranks without a taxonomy")
        if arguments.ranks is None:
            # set a default
            arguments.ranks = [
                'phylum',
                'class',
                'order',
                'family',
                'genus',
                'species']

        try:
            # Make sure the rank lists make sense
            arguments.ranks = cleanRanks(arguments.ranks)
            if arguments.printRanks is not None:
                arguments.printRanks = cleanRanks(arguments.printRanks)
        except Exception as e:
            parser.error(str(e))

    # load weights file
    sequenceWeights = loadSequenceWeights(arguments.weights)

    # only print to stdout if there is a single rank
    if len(arguments.ranks) > 1 and arguments.outfile is None:
        parser.error("STDOUT only works if a single rank is chosen!")

    cutoff = arguments.cutoff

    # Because rank is used in parsing hits, we can only do multiple ranks for
    # certain kinds of count methods
    if len(arguments.ranks) > 1:
        rank = None
        if arguments.countMethod in ['consensus', 'most']:
            parser.error(
                "Using multiple ranks does not work with the 'consensus' "
                "or 'most' counting methods. LCA should give the same "
                "results as consensus. If you really want to do this, "
                "use a bash loop:'for rank in phylum order genus; do "
                "COMMAND -r ${rank}; done'")
    else:
        rank = arguments.ranks[0]

    # load necessary maps
    (taxonomy, hitStringMap) = readMaps(arguments)

    # parse input files
    fileCounts = {}
    totals = {}
    fileLabels = {}
    sortedLabels = []

    # Allow for file names to be preceded with TAG=
    for filename in arguments.input_files:
        bits = filename.split("=", 1)
        if len(bits) > 1:
            (filetag, filename) = bits
        else:
            filetag = filename
        fileLabels[filename] = filetag
        # keep order so that column order matches arguments
        sortedLabels.append(filetag)
        fileCounts[filetag] = {}
        totals[filetag] = 0

    if arguments.countMethod == 'tophit' or arguments.countMethod == 'toporg':
        # Process all files at once and use overall abundance to pick best hits
        from edl import redistribute
        params = FilterParams.create_from_arguments(arguments)
        multifile = redistribute.multipleFileWrapper(fileLabels.keys())

        if arguments.countMethod == 'tophit':
            # don't give any taxonomy, just map to accessions for
            # redistribution
            readHits = redistribute.pickBestHitByAbundance(
                multifile,
                filterParams=params,
                returnLines=False,
                winnerTakeAll=True,
                parseStyle=arguments.parseStyle,
                sequenceWeights=sequenceWeights)
            # define method to turn Hits into orgnaisms
            hitTranslator = getHitTranslator(parseStyle=arguments.parseStyle,
                                             taxonomy=taxonomy,
                                             hitStringMap=hitStringMap)

            translateHit = lambda hit: hitTranslator.translateHit(hit=hit)[0]

        else:
            # translate to organism before finding most abundant
            readHits = redistribute.pickBestHitByAbundance(
                multifile,
                filterParams=params,
                returnLines=False,
                returnTranslations=True,
                winnerTakeAll=True,
                taxonomy=taxonomy,
                hitStringMap=hitStringMap,
                parseStyle=hits.ACCS)

            # Organisms will be returned, make translator trivial:
            translateHit = passThrough

        # use read->file mapping and hit translator to get file based counts
        #  from returned (read,Hit) pairs
        increment = 1
        for (read_name, hit) in readHits:
            file_name, read_name = read_name.split("/", 1)
            file_tag = fileLabels[unquote_plus(file_name)]
            taxon = translateHit(hit)
            taxcount = fileCounts[file_tag].setdefault(taxon, 0)
            if sequenceWeights is not None:
                increment = sequenceWeights.get(read_name, 1)
            fileCounts[file_tag][taxon] = taxcount + increment
            totals[file_tag] += increment
        logging.debug(str(totals))

    else:
        # Original way, just process each file separately
        for (filename, filetag) in fileLabels.items():
            infile = open(filename, 'rU')

            hitIter = parseM8FileIter(infile,
                                      hitStringMap,
                                      arguments.hitTableFormat,
                                      arguments.filterTopPct,
                                      arguments.parseStyle,
                                      arguments.countMethod,
                                      taxonomy=taxonomy,
                                      rank=rank)

            (total, counts, hitMap) = \
                countIterHits(hitIter,
                              allMethod=arguments.allMethod,
                              weights=sequenceWeights)
            fileCounts[filetag] = counts
            totals[filetag] = total

            logging.info(
                "parsed %d hits (%d unique) for %d reads from %s" %
                (total, len(counts), len(hitMap), filename))

            infile.close()

    printCountTablesByRank(fileCounts, totals, sortedLabels, arguments)


def cleanRanks(rankList):
    if ORG_RANK not in ranks:
        ranks.insert(0, ORG_RANK)

    # don't allow duplicates
    rankList = list(set(rankList))

    # translate domain to superkingdom
    if 'domain' in rankList:
        rankList.remove('domain')
        rankList.append('superkingdom')

    # make sure the ranks are real
    badRanks = []
    for rank in rankList:
        if rank not in ranks:
            badRanks.append(rank)
    if len(badRanks) > 0:
        parser.error("Unknown rank(s): %s" % (badRanks))

    # return ranks in proper order
    return sorted(rankList, key=ranks.index, reverse=True)


def printCountTablesByRank(fileCounts, totals, fileNames, options):
    """
    Create a new file for each rank witha tab separated table of counts
    """
    cutoff = options.cutoff

    # create an output table for each requested rank
    for rank in options.ranks:

        # For each rank, try to force all counts to be at that rank
        fileRankTotals = {}
        rankCounts = {}
        rankTaxa = {}
        thresholds = {}
        for (filename, counts) in fileCounts.items():
            fileRankTotals[filename] = 0
            thresholds[filename] = totals[filename] * cutoff
            fileRankCounts = rankCounts.setdefault(filename, {})

            fileTotal = 0
            for taxon in counts.keys():
                # get the counts from this node
                taxonCount = counts[taxon]
                fileTotal += taxonCount

                # get parent taxon at the given rank
                if taxon is None:
                    ranked = None
                elif rank is None or rank == ORG_RANK:
                    ranked = taxon
                else:
                    if options.collapseToDomain:
                        # If we are beyond this rank already, fall back to SK
                        fallback = taxon.getAncestorAtRank('superkingdom')
                        if fallback is None:
                            fallback = taxon.getRootNode()
                    else:
                        fallback = taxon
                    ranked = getAncestorClosestToRank(
                        taxon,
                        rank,
                        default=fallback,
                        useChildOfFirstRankedAncestor=not(
                            options.collapseToDomain))
                    if ranked is None:
                        # This shouldn't happen...
                        logging.warn("getAncestorClosestRoRank return None!")
                        # ...but if it doesn, leave unchanged
                        ranked = taxon

                # update counts
                fileRankCounts[ranked] = fileRankCounts.get(
                    ranked, 0) + taxonCount
                rankTaxa[ranked] = True

            logging.debug(
                "File %s has %d hits (had %d)" %
                (filename, fileTotal, totals[filename]))

        # logging.debug(repr(rankTaxa))
        # logging.debug(repr(rankCounts))
        if logging.getLogger().level <= logging.DEBUG:
            for (filename, counts) in rankCounts.items():
                logging.debug("File %s hs %d ranked counts" %
                              (filename, sum(counts.values())))

        # apply cutoff
        for taxon in list(rankTaxa.keys()):
            # check to see if taxon is over cutoff in any file
            over = False
            for (filename, fileRankCount) in rankCounts.items():
                frTaxonCount = fileRankCount.get(taxon, 0)
                fileRankTotals[filename] += frTaxonCount
                if frTaxonCount > thresholds[filename]:
                    over = True
            if not over:
                # this taxon is not over the cutoff for any file
                rankTaxa.pop(taxon)
                if taxon is not None:
                    if options.taxdir is None:
                        other = 'Other'
                    else:
                        other = taxon.getAncestorAtRank('superkingdom')
                        if other is None:
                            other = taxon.getRootNode()
                else:
                    other = None

                rankTaxa[other] = True
                for (filename, fileRankCount) in rankCounts.items():
                    fileRankCount[other] = fileRankCount.get(
                        other, 0) + fileRankCount.pop(taxon, 0)

        if logging.getLogger().level <= logging.DEBUG:
            for (filename, counts) in rankCounts.items():
                logging.debug("File %s hs %d ranked counts" %
                              (filename, sum(counts.values())))
                missed = False
                for taxa in counts.keys():
                    if taxa not in rankTaxa:
                        missed = True
                        logging.debug(
                            "Missing taxon %s has %d counts for %s" %
                            (taxa, counts[taxa], filename))
                if not missed:
                    logging.debug(
                        "There are no missing taxa from %s" %
                        (filename))

        logging.debug("Final file counts: %r" % (fileRankTotals))

        # output file
        if options.outfile is None:
            outs = sys.stdout
        else:
            if len(options.ranks) > 1:
                outfile = "%s.%s" % (options.outfile, rank)
            else:
                outfile = options.outfile
            outs = open(outfile, 'w')

        # write to file(s?)
        # header
        outs.write("Taxon\t%s\n" % ('\t'.join(fileNames)))

        taxonFormatter = getTaxonFormatter(options.printRanks, rank)
        for taxon in sorted(rankTaxa.keys(), key=taxonFormatter):
            outs.write(taxonFormatter(taxon))
            for filename in fileNames:
                outs.write("\t")
                outs.write(str(rankCounts[filename].get(taxon, 0)))
            outs.write("\n")

        # close out stream
        if options.outfile is not None:
            outs.close()


def getTaxonFormatter(displayedRanks, leafRank):
    if displayedRanks is None:
        return str
    else:
        return lambda t: formatTaxon(t, displayedRanks, leafRank)


def formatTaxon(taxon, displayedRanks, leafRank, delim=';'):
    """
    Generates lineage using all display ranks that are less than the 
    leaf rank. This is probably ineffecient, as we have to figure out
    which ranks to display for every item. 

    This method is also used by assign_taxa.py!
    """
    if isinstance(taxon,list):
        if len(taxon)==0:
            taxon = None
        elif len(taxon)==1:
            taxon = taxon[0]
        else:
            raise Exception("taxon should not be a list:\n{}"\
                             .format(repr(taxon)))
    if taxon is None:
        logging.debug("Taxon is None")
        return 'None'
    if taxon is taxon.getRootNode():
        return str(taxon)

    lineage = ""
    logging.debug(
        "Creating lineage with: %s, %s, %s" %
        (taxon, displayedRanks, leafRank))
    for rank in displayedRanks:
        if ranks.index(rank) <= ranks.index(leafRank):
            logging.debug(
                "Rank of %s (%d) is less than %s (%d)" %
                (rank, ranks.index(rank), leafRank, ranks.index(leafRank)))
            break
        ancestor = taxon.getAncestorAtRank(rank)
        if ancestor is taxon:
            logging.debug(
                "ancestor at %s of %s is %s" %
                (taxon, rank, ancestor))
            break
        if ancestor is None:
            ancestor = ""
        lineage += str(ancestor) + delim
        logging.debug("Lineage: %s" % lineage)
    lineage += str(taxon)
    logging.debug("Lineage: %s" % lineage)
    return lineage

if __name__ == '__main__':
    main()
