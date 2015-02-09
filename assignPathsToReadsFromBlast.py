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
    parser.add_option('-O', "--outputStyle", default="cols",
                      choices=['cols','lines','python'],
                      help="How are multiple assignments displayed in output. By default ('cols'), multiple hits show up in multiple columns. The 'lines' option prints out a new line for each assignment. The 'python' option prints each assignment as a python string (in quotes) or a list of strings (in quotes, separted by commas, surrounded bya  pair of sqaure brackets).")
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
                      help="How to deal with assignments from multiple hits. (first, most: can return multiple hits, all (default): return every hit, consensus: return None unless all the same)",
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
        if options.outputStyle=='python':
            for read in sorted(hitMap.keys()):
                hit=hitMap[read]
                outhandle.write(str(read))
                outhandle.write("\t")
                outhandle.write(repr(hit))
                outhandle.write("\n")
        if options.outputStyle=='lines':
            for read in sorted(hitMap.keys()):
                hit=hitMap[read]
                if type(hit) is type([]):
                    for h in sorted(hit):
                        outhandle.write(str(read))
                        outhandle.write("\t")
                        outhandle.write(str(h))
                        outhandle.write("\n")
                else:
                    outhandle.write(str(read))
                    outhandle.write("\t")
                    outhandle.write(str(hit))
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

