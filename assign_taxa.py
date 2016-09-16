#! /usr/bin/python
"""
"""

from optparse import OptionParser
import sys, re, logging
from edl.taxon import ranks, TaxNode
from edl.redistribute import redistributeHits
from edl.hits import *
from edl.util import *
from edl.expressions import accessionRE, nrOrgRE

def main():
    usage = "usage: %prog [OPTIONS] HIT_TABLE(S)"
    description = """
Takes an m8 blast and assigns each read to a taxon. Blast may be specified with -i or piped to STDIN.

Attempts to process one record at a time (without filling up RAM) when possible. This is NOT possible if any of the following is selected: a non-zero count cutoff, redistribted LCA (rLCA) counting method, read sorting.
    """
    parser = OptionParser(usage, description=description)
    parser.add_option("-i", "--inputfile", dest="infile",
                      metavar="INFILE", help="Read data table from INFILE. (this is here for legacy, you can also just put list file names without this options)"),
    addIOOptions(parser)
    parser.add_option('-O', "--pythonText", default=False,
                      action='store_true',
                      help="Output is formatted as python strings (hits in quotes and multiple hits in brackets). By default, there are no quotes and multiple hits are in multiple columns")
    parser.add_option("-T", "--taxids", default=False, action="store_true",
                      help="Output taxids instead of names")
    addTaxonOptions(parser)
    parser.add_option("-r", "--rank", dest="rank", default=None,
                      metavar="RANK",
                      help=" Rank to collect counts on.  Defaults to None (whatever the annotation was). Corresponds to rank names in nodes.dmp. To see list run: 'cut -f5 nodes.dmp | uniq | sort | uniq' in ncbi tax dir")
    parser.add_option("-R", "--fallbackRank", default=None, metavar="RANK",
                      help="If applying cutoff with rank, use this second rank to collcet counts under cutoff (instead of just moving up phylogeny step by step). 'superkingdom' is recommended")
    parser.add_option("-c", "--cutoff", dest="cutoff", type="float", default=0.,
                      help="Cutoff for showing taxa. May be given as a fraction or absolute number. If a (fractional) count for a taxa is below this value, it will be folded up into its parent. Rank option will be interpreted as lowest level to split counts to. Default is 0 (no cutoff).",
                  metavar="CUTOFF")
    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options,description)

    logging.debug("PArsing style is: %s" % (options.parseStyle))

    # check arguments
    if options.taxids and options.taxdir is None:
        parser.error("Only use -T when a taxonomy is specified")
    if options.rank is not None and options.taxdir is None:
        parser.error("Please supply NCBI phylogeny(-n) if specifying a rank(-r).")
    if options.rank is not None:
        if options.rank == 'domain':
            logging.warn('translating domain to superkingdom')
            options.rank='superkingdom'
        if options.rank not in ranks:
            parser.error("Unknown rank: %s" % (options.rank))
        if options.fallbackRank == 'domain':
            logging.warn('translating domain to superkingdom')
            options.fallbackRank='superkingdom'

    # load necessary maps
    (taxonomy,valueMap)=readMaps(options)

    # loop over inputs
    for (inhandle,outhandle) in inputIterator(args, options):
        logging.debug("Reading from %s and writing to %s" % (inhandle, outhandle))
        hitIter = parseM8FileIter(inhandle, valueMap, options.hitTableFormat, options.filterTopPct, options.parseStyle, options.countMethod,taxonomy=taxonomy,rank=options.rank,sortReads=options.hitTableSortReads)
        hitMap=None

        # manipulate mappings based on options
        if options.cutoff is not None and options.cutoff>0:
            if taxonomy is None:
                hitMap = applySimpleCutoff(hitIter,
                                           options.cutoff)
            else:
                hitMap = applyCutoff(hitIter,
                                 options.cutoff, options.fallbackRank)

            logging.info("Cutoffs complete for %d reads" % (len(hitMap)))

        # apply rLCA
        if options.countMethod == 'rLCA':
            if hitMap is None:
                # redistribute method can take either map or iterator
                hitMap = hitIter
            hitMap = redistributeHits(hitMap, options.rank)
            logging.info("Redistribution complete for %d reads" % (len(hitMap)))

        # get sorted iterator over map, so counts are reproducible
        if hitMap is not None:
            hitIter = sortedHitIterator(hitMap)

        ##
        # print output
        # choose output method
        if options.taxids:
            printer = taxidPrinter
        elif options.pythonText:
            printer = pythonTextPrinter
        else:
            printer = defaultPrinter

        # loop over reads
        outhandle.write("Read\tHit\n")
        for (read,hits) in hitIter:
            outhandle.write(printer(read,hits))

#############
# Functions #
#############
def applySimpleCutoff(hitMap, cutoff):
    """
    collect hits and remove items under cutoff
    """

    logging.info("Starting applySimpleCutoff (%d)" % (cutoff))
    # bin counts by hit and find total
    (total, counts, hitMap) = countIterHits(hitMap)

    logging.debug("%d different hits from %d total" % (len(counts),total))

    # identify hits under cutoff
    if cutoff<=1.0:
        minCount = float(cutoff)*total
    else:
        minCount = cutoff
    hitTranslation={}
    for (hit,count) in counts.iteritems():
        if count<minCount:
            hitTranslation[hit]='other'


    # update hits and counts from translation
    logging.debug("applySimpleCutoff: done translating hit names")
    translateHits(hitMap, hitTranslation)
    logging.debug("apllyCutoff: done translating hits")
    if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        (total,counts)=countHits(hitMap)
        logging.debug("recounted: %d different hits from %d total" % (len(counts),total))

    return hitMap

def applyCutoff(hitMap, cutoff, fallbackRank):
    """
    apply cutoff using phylogeny to combine hits
    """
    if cutoff<=0:
        logging.debug("no cutoff to apply!")
        return None

    if fallbackRank is not None:
        return applyFallbackCutoff(hitMap, cutoff, fallbackRank)

    logging.info("Starting applyCutoff (%d)" % (cutoff))
    (total, counts, hitMap) = countIterHits(hitMap)

    logging.debug("%d different hits from %d total" % (len(counts),total))

    # apply cutoff
    hitTranslations = {}
    if cutoff <= 1.0:
        threshold = float(cutoff)*total
    else:
        threshold = cutoff
    logging.info("Collapsing counts under threshold: %s" % (threshold))

    # first apply cutoff to any hits that aren't TaxNode objects
    for hit in counts:
        if not isinstance(hit,TaxNode):
            logging.debug("untranslated hit: %s: %s" % (hit, counts[hit]))
            if float(counts[hit])<threshold:
                hitTranslations[hit]="other"

    #get the root of the tax tree
    for hit in counts.iterkeys():
        if isinstance(hit,TaxNode):
            root=hit.getRootNode()
            break
    else:
        sys.exit("No hits were translated into TaxNode objects!")
    logging.debug("Root node is %s %s" % (root.id, root.name))

    # now walk up tax tree from root, collecting counts
    (leftoverCount,leftoverNodes) = root.getCollapsedCounts(counts,threshold,hitTranslations)

    # deal with anything that fell out to the root node
    if len(leftoverNodes)>0:
        for node in leftoverNodes:
            hitTranslations[node]='root'
        logging.warn("%d hits from %d nodes assigned to root node!" % (leftoverCount,len(leftoverNodes)))

    # apply collapse via hit translations
    translateHits(hitMap, hitTranslations)
    logging.debug("apllyCutoff: done applying Cutoff")

    if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        (total, counts) = countHits(hitMap)
        logging.debug("%d different hits from %d total" % (len(counts),total))

    return hitMap

def applyFallbackCutoff(hitMap, cutoff, fallbackRank):
    logging.info("Starting applyCutoffFallback: (%d,%s)" % (cutoff,fallbackRank))
    (total, counts, hitMap) = countIterHits(hitMap)

    logging.debug("%d different hits from %d total" % (len(counts),total))

    # set up
    hitTranslations = {}
    if cutoff <= 1.0:
        threshold = float(cutoff)*total
    else:
        threshold = cutoff
    logging.info("Collapsing counts under threshold: %s" % (threshold))

    #get the root of the tax tree
    for hit in counts.iterkeys():
        if isinstance(hit,TaxNode):
            root=hit.getRootNode()
            break
    else:
        sys.exit("No hits were translated into TaxNode objects!")
    logging.debug("Root node is %s %s" % (root.id, root.name))

    # apply cutoff
    for hit in counts:
        if float(counts[hit])<threshold:
            if isinstance(hit,TaxNode):
                other=hit.getAncestorAtRank('superkingdom')
                if other is None:
                    other = root
            else:
                other='other'
            hitTranslations[hit]=other

    # apply collapse via hit translations
    translateHits(hitMap, hitTranslations)
    logging.debug("applyCutoff: done applying Cutoff")

    if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        (total, counts) = countHits(hitMap)
        logging.debug("%d different hits from %d total" % (len(counts),total))

    return hitMap

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

def pythonTextPrinter(read,hit):
    """
    Simply print the repr of each hit list
    """
    return("%s\t%r\n" % (read,hit))

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

