#! /usr/bin/python
"""
"""

from optparse import OptionParser
import sys, logging
from edl.util import parseMapFile, inputIterator, addUniversalOptions, setupLogging, addIOOptions
from edl.kegg import *
from edl.hits import *
from edl.expressions import *
from edl.blastm8 import addHitTableOptions

def main():
    usage = "usage: %prog [OPTIONS] INPUT_FILE(S)"
    description = """
Takes an m8 blast and assigns each read to a pathway or gene family. Blast may be specified with -i or piped to STDIN.
    """
    parser = OptionParser(usage, description=description)
    parser.add_option("-i", "--inputfile", dest="infile",
                      metavar="INFILE", help="Read data table from INFILE"),
    addIOOptions(parser)
    parser.add_option('-O', "--pythonText", default=False,
                      action='store_true',
                      help="Output is formatted as python strings (hits in quotes and multiple hits in brackets). By default, there are no quotes and multiple hits are in multiple columns")
    parser.add_option("-m", "--mapFile", dest="mapFile",
                      metavar="MAPFILE", help="Location of file containing table of with db hit name as first column and geneIDs (Knumber) in second column.")
    parser.add_option("-M", "--mapStyle", default='auto', choices=['auto','kegg','tab'],
                      help="What type of mapping file are you using: simple tab separated list of IDs and kos, or the genes_ko.list file from KEGG (which adds ko: to the K numbers and can have multiple records for each gene id). By default, this script will inspect the file name and guess, but you can force either 'kegg' or 'tab' with this option.")
    parser.add_option("-p", "--parseStyle",
                      default=KEGG,
                      choices=[ACCS,GIS,KEGG,HITID,HITDESC],
                      help="What should be parsed from the hit table: accessions('accs'), 'gis', K numbers in description ('kegg'), the full hit name('hitid'), or the full hit description('hitdesc'). (defaults to '%default')")
    parser.add_option("-c", "--cutoff", dest="cutoff", type="float", default=0.01,
            help="Cutoff for showing paths or genes. If a fractional count for a path/gene is below this value, it will be labelled None.",
                  metavar="CUTOFF")

    # format and filterPct
    addHitTableOptions(parser)

    parser.add_option("-C", "--countMethod", dest="countMethod", default="all", choices=('first','most','all','consensus'),
                      help="How to deal with counts from multiple hits. (first, most: can return multiple hits, all (default): return every hit, consensus: return None unless all the same)",
                    metavar="COUNTMETHOD")
    parser.add_option("-r","--filterForKO",action="store_true", dest="koHitsOnly", default=False, help="ignore hits with no KO assignment. This means reads with no hits to KO tagged sequences will not be in the output.")
    parser.add_option("-l","--level", dest="level", default="ko", choices=('ko','NAME','DEFINITION','EC','PATHWAY','1','2','3'), help="Either 'ko'; a string to look for in ko file ('PATHWAY','NAME', 'DEFINITION', or 'EC'); or level in kegg class heirarchy (1, 2, or 3 (should be same as PATHWAY))")
    parser.add_option("-k", "--koFile", dest="ko", metavar="KOFILE", default=None,
                      help="File containing kegg heirarchy (either ko or ko00001.keg)")
    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options, description)

    if options.infile is None:
        infile = sys.stdin
    else:
        infile = open(options.infile)

    if options.parseStyle==KEGG:
        if options.mapFile is not None:
            logging.warn("Do you REALLY want to apply a mapping to KOs?")

    if options.level != 'ko':
        if options.ko is None:
            options.error("Please supply KEGG file if sepcifying a level other than 'ko' ")

        # read KEGG file
        koTranslation = readKEGGFile(options.ko, options.level)
    else:
        koTranslation = None

    # map reads to hits
    if options.mapFile is not None:
        if options.mapStyle=='kegg' or ( options.mapStyle=='auto' and len(options.mapFile)>=13 and options.mapFile[-13:]=='genes_ko.list'):
            valueMap=parseLinkFile(options.mapFile)
        else:
            if options.parseStyle == GIS:
                keyType=int
            else:
                keyType=None
            valueMap = parseMapFile(options.mapFile,valueType=None,keyType=keyType)
    else:
        valueMap=None

    for (inhandle,outhandle) in inputIterator(args, options):
        logging.debug("Reading from %s and writing to %s" % (inhandle, outhandle))
        hitMap = parseM8File(inhandle, valueMap, options.hitTableFormat, options.filterTopPct, options.parseStyle, options.countMethod, ignoreEmptyHits=options.koHitsOnly,sortReads=options.hitTableSortReads)

        # manipulate mappings
        hitMap = applySimpleCutoff(hitMap, options.cutoff, koTranslation)

        log("maps complete for %d reads" % (len(hitMap)))

        # print out hit table
        outhandle.write("Read\tHit\n")
        if options.pythonText:
            for read in sorted(hitMap.keys()):
                hit=hitMap[read]
                outhandle.write(str(read))
                outhandle.write("\t")
                outhandle.write(repr(hit))
                outhandle.write("\n")
        else:
             for read in sorted(hitMap.keys()):
                hit=hitMap[read]
                outhandle.write(str(read))
                if type(hit) is type([]):
                    for h in sorted(hit):
                        outhandle.write("\t")
                        outhandle.write(str(h))
                else:
                    outhandle.write("\t")
                    outhandle.write(str(hit))
                outhandle.write("\n")

#############
# Functions #
#############
log = logging.info
warn = logging.warn
debug = logging.debug

def parseM8FileOrig(infile, mapFile, hitCol, style, scoreCol, countMethod, filterHits):
    """
    return a map from read names to hits. If mapFile given, translate hit names. Use
    scoreCol, etc to pick best hit among many
    """

    translateHit=False
    if mapFile is not None:
        hitMap = parseMapFile(mapFile)
        translateHit = True
    justFirst = countMethod=="first"
    readMap={}
    lastRead=None
    hits={}
    for line in infile:
        cells = line.split('\t')
        read = cells[0].strip()
        hit = cells[hitCol].strip()
        if style is ACCS:
            m=accessionRE.search(hit)
            if (m):
                hit=m.group(1)
        elif style is KEGG:
            m = koRE.search(hit)
            if m:
                hit=m.group(0)
            else:
                hit=None
        if translateHit:
            # return hit if it's not in the map
            hit = hitMap.get(hit,hit)
        if filterHits:
            if hit is None or hit == "None" or hit == '':
                continue

        if read != lastRead:
            if justFirst:
                readMap[read]=hit
                lastRead=read
                continue

            if lastRead is not None:
                readMap[lastRead]=pickBestHit(hits, countMethod)
                hits={}
            lastRead=read

        if not justFirst:
            score = cells[scoreCol].strip()
            scores = hits.setdefault(hit,[])
            scores.append(score)

    if not justFirst:
        readMap[read]=pickBestHit(hits, countMethod)

    return readMap

def pickBestHit(hits, countMethod):
    if countMethod == 'plurality':
        bestHit = pickMostCommonHit(hits)
    elif countMethod == 'best-score':
        bestHit = findBestScore(hits)
    else:
        sys.exit("Unrecognized count method: %s!" % (countMethod))

    # make sure there is not a None or '' in the list
    if type(bestHit) is type([]):
        try:
            bestHit.remove(None)
        except ValueError:
            pass
        try:
            bestHit.remove('')
        except ValueError:
            pass
        if len(bestHit)==1:
            bestHit=bestHit[0]
        elif len(bestHit)==0:
            bestHit=None

    return bestHit

def findBestScore(hits):
    bestHit=None
    bestScore=0
    for hit in hits:
        maxScore = max(hits[hit])
        if bestHit is None or maxScore>bestScore:
            bestHit = hit
            bestScore = maxScore
        elif maxScore==bestScore:
            if type(bestHit)==type([]):
                bestHit.append(hit)
            else:
                bestHit=[bestHit,hit]

    return bestHit

def pickMostCommonHit(hits):
    """
    for all hits, return hit with most scores.
    Use sum of scores as a tiebreaker.
    """

    """
    # find utoff score
    cutoff=0
    if topPct>=0:
        # start with the highest score
        for (hit,scores) in hits.iteritems():
            cutoff=max(max(scores),cutoff)

        # subtract topPct from it
        if topPct>0:
            cutoff = cutoff*(100.0 - float(topPct))/100.0
    """

    # count hits over cutoff
    bestHit=None
    bestCount=0
    bestScores=0
    for (hit,scores) in hits.iteritems():
        count=len(scores)
        scoreSum=0
        for score in scores:
            scoreSum+=float(score)
        if bestHit is None or count>bestCount:
            bestHit=hit
            bestCount=count
            bestScores=scoreSum
        elif count==bestCount:
            if scoreSum>bestScores:
                bestHit=hit
                bestCount=count
                bestScores=scoreSum
            elif scoreSum==bestScores:
                if type(bestHit)==type([]):
                    bestHit.append(hit)
                else:
                    besetHit=[bestHit,hit]

    return bestHit

def applySimpleCutoff(hitMap, cutoff, koTranslation):
    """
    collect hits and remove items under cutoff
    """
    if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        debug("Starting applySimpleCutoff (%d) with hits for %d reads" % (cutoff,len(hitMap)))
        # bin counts by hit and find total
        (total, counts) = countHits(hitMap)
        debug("%d different hits from %d total" % (len(counts),total))

    # rename hits to Gene or Pathway
    if koTranslation is not None:
        translateHits(hitMap, koTranslation)

    # bin counts by hit and find total
    (total, counts) = countHits(hitMap)

    # identify hits under cutoff
    minCount = float(cutoff)*total
    hitTranslation={}
    for (hit,count) in counts.iteritems():
        if count<minCount:
            hitTranslation[hit]='other'

    # alter hit map accordingly
    translateHits(hitMap, hitTranslation)
    log("apllyCutoff: done translating hits")
    if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
        (total,counts)=countHits(hitMap)
        log("recounted: %d different hits from %d total" % (len(counts),total))

    return hitMap

if __name__ == '__main__':
    main()

