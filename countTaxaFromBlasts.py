#! /usr/bin/python
"""
"""
from optparse import OptionParser
import sys, re
from edl.taxon import *
from edl.hits import *
from edl.util import addUniversalOptions, setupLogging, checkNoneOption, parseMapFile
from edl.expressions import accessionRE, nrOrgRE

ORG_RANK='organism'

def main():
    usage = "usage: %prog [OPTIONS] BLAST_M8_FILE[S]"
    description = """
Takes m8 blast files and generates a table of taxon hit counts for the given rank. Columns are input files and rows are taxa. If multiple ranks given (the default), multiple output files are produced, each with the rank name appended to the output file name.
    """
    parser = OptionParser(usage, description=description)
    parser.add_option("-o", "--outfile", dest="outfile",
                      metavar="OUTFILE", help="Write count table to OUTFILE")
    parser.add_option("-r", "--rank", dest="ranks", default=None,
                      metavar="RANK", action="append",
                      help=""" Rank(s) to collect counts on. Use flag multiple
                      times to specify multiple ranks. If multiple values
                      given, one table produced for each with rank name
                      appended to file name. Defaults to all major ranks
                      between phylum and species. Corresponds to rank names 
                      in nodes.dmp. To see list run: 
                      'cut -f5 nodes.dmp | uniq | sort | uniq' 
                      in ncbi tax dir. Will also accept 'organism' to mean 
                      no rank (ie, just the organism name).""")
    parser.add_option("-s","--collapseToDomain", default=False, action="store_true",
                      help="Collapse all taxa below given rank down to superkingdom/domain. EG: in the genus output, anything assigned to Cyanobactia, will be lumped in with all other bacteria")
    parser.add_option("-R","--printRank",dest="printRanks",action="append",
                      help="Include indeicated rank(s) in lineage of printed taxa. Will be ignored if beyond the rank of the taxa (IE We can't include species if the taxon being counted is genus)")

    # option for deconvoluting clusters or assemblies
    addWeightOption(parser, multiple=True)

    # cutoff options
    addCountOptions(parser)

    # format, tax dir, and more
    addTaxonOptions(parser,choices={'countMethod':('LCA','all','first','most','tophit','toporg','consensus')})

    # log level and help
    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options, description)

    if len(args)==0:
        parser.error("Must supply at least one m8 file to parse")

    # Handle the case where Galaxy tries to set None as a string
    options.ranks=checkNoneOption(options.ranks)
    options.printRanks=checkNoneOption(options.printRanks)

    # Set defaults and check for some conflicts
    if options.ranks is None and options.taxdir is None:
        # using hit names only
        options.ranks=[ORG_RANK]
        if options.printRanks is not None:
            parser.error("Display ranks are not used without taxonomic info")
    else:
        if options.taxdir is None:
            parser.error("Cannot select ranks without a taxonomy")
        if options.ranks is None:
            # set a default
            options.ranks=['phylum','class','order','family','genus','species']

        try:
            # Make sure the rank lists make sense
            options.ranks=cleanRanks(options.ranks)
            if options.printRanks is not None:
                options.printRanks=cleanRanks(options.printRanks)
        except Exception as e:
            parser.error(str(e))

    # load weights file
    sequenceWeights = loadSequenceWeights(options.weights)

    # only print to stdout if there is a single rank
    if len(options.ranks)>1 and options.outfile is None:
        parser.error("STDOUT only works if a single rank is chosen!")

    cutoff=options.cutoff

    # Because rank is used in parsing hits, we can only do multiple ranks for
    # certain kinds of count methods
    if len(options.ranks)>1:
        rank=None
        if options.countMethod in ['consensus','most']:
            parser.error("Using multiple ranks does not work with the 'consensus' or 'most' counting methods. LCA should give the same results as consensus. If you really want to do this, us a bash loop:'for rank in phylum order genus; do COMMAND -r ${rank}; done'")
    else:
        rank=options.ranks[0]

    # load necessary maps
    (taxonomy,hitStringMap)=readMaps(options)

    # parse input files
    fileCounts={}
    totals={}
    fileLabels={}
    sortedLabels=[]

    # Allow for file names to be preceded with TAG=
    for filename in args:
        bits=filename.split("=",1)
        if len(bits) > 1:
            (filetag,filename)=bits
        else:
            filetag=filename
        fileLabels[filename]=filetag
        # keep order so that column order matches arguments
        sortedLabels.append(filetag)
        fileCounts[filetag]={}
        totals[filetag]=0

    if options.countMethod == 'tophit' or options.countMethod == 'toporg':
        # Process all files at once and use overall abundance to pick best hits
        from edl import redistribute
        params = FilterParams.createFromOptions(options)
        (multifile,readFileDict) = redistribute.multipleFileWrapper(fileLabels.keys(), params, returnLines=True)

        if options.countMethod == 'tophit':
            # don't give any taxonomy, just map to accessions for redistribution
            readHits = redistribute.pickBestHitByAbundance(multifile,
                                              filterParams=params,
                                              returnLines=False,
                                              winnerTakeAll=True,
                                              parseStyle=options.parseStyle,
                                              sequenceWeights=sequenceWeights)
            # define method to turn Hits into orgnaisms
            hitTranslator=getHitTranslator(parseStyle=options.parseStyle,
                                           taxonomy=taxonomy,
                                           hitStringMap=hitStringMap)
            translateHit = lambda hit: hitTranslator.translateHit(hit)[0]
        else:
            # translate to organism before finding most abundant
            readHits = redistribute.pickBestHitByAbundance(multifile,
                                                      filterParams=params,
                                                      returnLines=False,
                                                      returnTranslations=True,
                                                      winnerTakeAll=True,
                                                      taxonomy=taxonomy,
                                                      hitStringMap=hitStringMap,
                                                      parseStyle=hits.ACCS)
            # Organisms will be returned, make translator trivial:
            translateHit = lambda hit: hit


        # use read->file mapping and hit translator to get file based counts
        #  from returned (read,Hit) pairs
        increment=1
        for (read, hit) in readHits:
            filename = readFileDict[read]
            filetag = fileLabels[filename]
            taxon = translateHit(hit)
            taxcount = fileCounts[filetag].setdefault(taxon,0)
            if sequenceWeights is not None:
                increment = sequenceWeights.get(read,1)
            fileCounts[filetag][taxon]=taxcount+increment
            totals[filetag]+=increment
        logging.debug(str(totals))

    else:
        # Original way, just process each file separately
        for (filename, filetag) in fileLabels.iteritems():
            infile = open(filename,'rU')

            hitIter = parseM8FileIter(infile, hitStringMap, options.hitTableFormat, options.filterTopPct, options.parseStyle, options.countMethod, taxonomy=taxonomy, rank=rank, sortReads=options.hitTableSortReads)

            (total,counts,hitMap) = countIterHits(hitIter,allMethod=options.allMethod, weights=sequenceWeights)
            fileCounts[filetag] = counts
            totals[filetag]=total

            logging.info("parsed %d hits (%d unique) for %d reads from %s" % (total, len(counts), len(hitMap),filename))

            infile.close()

    printCountTablesByRank(fileCounts, totals, sortedLabels, options)

def cleanRanks(rankList):
    if ORG_RANK not in ranks:
        ranks.insert(0,ORG_RANK)

    # don't allow duplicates
    rankList=list(set(rankList))

    # translate domain to superkingdom
    if 'domain' in rankList:
        rankList.remove('domain')
        rankList.append('superkingdom')

    # make sure the ranks are real
    badRanks=[]
    for rank in rankList:
        if rank not in ranks:
            badRanks.append(rank)
    if len(badRanks)>0:
        parser.error("Unknown rank(s): %s" % (badRanks))

    # return ranks in proper order
    return sorted(rankList, key=ranks.index, reverse=True)

def printCountTablesByRank(fileCounts, totals, fileNames, options):
    """
    Create a new file for each rank witha tab separated table of counts
    """
    cutoff=options.cutoff

    # create an output table for each requested rank
    for rank in options.ranks:

        # For each rank, try to force all counts to be at that rank
        fileRankTotals={}
        rankCounts = {}
        rankTaxa = {}
        thresholds = {}
        for (filename, counts) in fileCounts.iteritems():
            fileRankTotals[filename]=0
            thresholds[filename] = totals[filename]*cutoff
            fileRankCounts = rankCounts.setdefault(filename, {})

            fileTotal=0
            for taxon in counts.keys():
                # get the counts from this node
                taxonCount = counts[taxon]
                fileTotal+=taxonCount

                # get parent taxon at the given rank
                if taxon is None:
                    ranked=None
                elif rank is None or rank == ORG_RANK:
                    ranked=taxon
                else:
                    if options.collapseToDomain:
                        # If we are beyond this rank already, fall back to SK
                        fallback = taxon.getAncestorAtRank('superkingdom')
                        if fallback is None:
                            fallback = taxon.getRootNode()
                    else:
                        fallback = taxon
                    ranked = getAncestorClosestToRank(taxon,rank,default=fallback, useChildOfFirstRankedAncestor=not(options.collapseToDomain))
                    if ranked is None:
                        # This shouldn't happen...
                        logging.warn("getAncestorClosestRoRank return None!")
                        # ...but if it doesn, leave unchanged
                        ranked=taxon

                # update counts
                fileRankCounts[ranked] = fileRankCounts.get(ranked,0) + taxonCount
                rankTaxa[ranked]=True

            logging.debug("File %s has %d hits (had %d)" % (filename, fileTotal, totals[filename]))

        #logging.debug(repr(rankTaxa))
        #logging.debug(repr(rankCounts))
        if logging.getLogger().level <= logging.DEBUG:
            for (filename,counts) in rankCounts.iteritems():
                logging.debug("File %s hs %d ranked counts" % (filename,sum(counts.values())))

        # apply cutoff
        for taxon in rankTaxa.keys():
            # check to see if taxon is over cutoff in any file
            over=False
            for (filename,fileRankCount) in rankCounts.iteritems():
                frTaxonCount=fileRankCount.get(taxon,0)
                fileRankTotals[filename]+=frTaxonCount
                if frTaxonCount > thresholds[filename]:
                    over=True
            if not over:
                # this taxon is not over the cutoff for any file
                rankTaxa.pop(taxon)
                if taxon is not None:
                    if options.taxdir is None:
                        other='Other'
                    else:
                        other=taxon.getAncestorAtRank('superkingdom')
                        if other is None:
                            other = taxon.getRootNode()
                else:
                    other=None

                rankTaxa[other]=True
                for (filename, fileRankCount) in rankCounts.iteritems():
                    fileRankCount[other]=fileRankCount.get(other,0) + fileRankCount.pop(taxon,0)

        if logging.getLogger().level <= logging.DEBUG:
            for (filename,counts) in rankCounts.iteritems():
                logging.debug("File %s hs %d ranked counts" % (filename,sum(counts.values())))
                missed=False
                for taxa in counts.iterkeys():
                    if taxa not in rankTaxa:
                        missed=True
                        logging.debug("Missing taxon %s has %d counts for %s" % (taxa, counts[taxa], filename))
                if not missed:
                    logging.debug("There are no missing taxa from %s" % (filename))

        logging.debug("Final file counts: %r" % (fileRankTotals))

        # output file
        if options.outfile is None:
            outs = sys.stdout
        else:
            if len(options.ranks)>1:
                outfile = "%s.%s" % (options.outfile, rank)
            else:
                outfile = options.outfile
            outs = open(outfile,'w')

        # write to file(s?)
        # header
        outs.write("Taxon\t%s\n" % ('\t'.join(fileNames)))

        taxonFormatter = getTaxonFormatter(options.printRanks, rank)
        for taxon in sorted(rankTaxa.keys(),key=taxonFormatter):
            outs.write(taxonFormatter(taxon))
            for filename in fileNames:
                outs.write("\t")
                outs.write(str(rankCounts[filename].get(taxon,0)))
            outs.write("\n")

        # close out stream
        if options.outfile is not None:
            outs.close()

def getTaxonFormatter(displayedRanks, leafRank):
    if displayedRanks is None:
        return str
    else:
        return lambda t: formatTaxon(t, displayedRanks, leafRank)

def formatTaxon(taxon, displayedRanks, leafRank):
    if taxon is None:
        logging.debug("Taxon is None")
        return 'None'
    if taxon is taxon.getRootNode():
        return str(taxon)

    lineage = ""
    logging.debug("Creating lineage with: %s, %s, %s" % (taxon, displayedRanks, leafRank))
    for rank in displayedRanks:
        if ranks.index(rank) <= ranks.index(leafRank):
            logging.debug("Rank of %s (%d) is less than %s (%d)" % (rank, ranks.index(rank), leafRank, ranks.index(leafRank)))
            break
        ancestor=taxon.getAncestorAtRank(rank)
        if ancestor is taxon:
            logging.debug("ancestor at %s of %s is %s" % (taxon, rank, ancestor))
            break
        if ancestor is None:
            ancestor=""
        lineage += str(ancestor) + ";"
        logging.debug("Lineage: %s" % lineage)
    lineage += str(taxon)
    logging.debug("Lineage: %s" % lineage)
    return lineage

if __name__ == '__main__':
    main()

