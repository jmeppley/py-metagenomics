#! /usr/bin/python
"""
"""

from optparse import OptionParser
import sys, re, logging
from edl import redistribute
from edl.hits import *
from edl.util import *
from edl.blastm8 import M8Stream

def main():
    usage = "usage: %prog [OPTIONS] HIT_TABLE(S)"
    description = """
Takes an m8 blast and picks the best hit for each. First, only the best scores are used, but if there is a tie (aka ambiguous hit), than a winner is assigned so that the proportion reads assigned to each organism matches the proportion of unambiguos hits.

    FilterPct defaults to 0, but can be altered, but I don't recommend it.
    ParseStyle and countMethod are ignored.
    """
    parser = OptionParser(usage, description=description)
    addIOOptions(parser)
    addTaxonOptions(parser,defaults={'filterPct':0,'parseStyle':ACCS,'countMethod':'tophit'},choices={'countMethod':('tophit','toporg')})
    addUniversalOptions(parser)
    parser.add_option("-i","--individualFiles", dest="individual", default=False, action="store_true", help="Use this flag to process files independently. Normally, counts from all files are pooled for making choices.")

    (options, args) = parser.parse_args()

    setupLogging(options,description)

    # load necessary maps
    params = FilterParams.createFromOptions(options)
    if options.countMethod=='toporg':
        (taxonomy,hitStringMap)=readMaps(options)

    if len(args)<=1 or options.individual:
        # loop over input
        for (inhandle,outhandle) in inputIterator(args, options):
            logging.debug("Reading from %s and writing to %s" % (inhandle, outhandle))

            m8stream=M8Stream(inhandle)
            if options.countMethod == 'tophit':
                # don't give any taxonomy, just map to accessions for redistribution
                readHits = redistribute.pickBestHitByAbundance(m8stream,
                                                      filterParams=params,
                                                      returnLines=True,
                                                      winnerTakeAll=True,
                                                      parseStyle=options.parseStyle)
            else:
                # translate to organism before finding most abundant
                readHits = redistribute.pickBestHitByAbundance(m8stream,
                                                          filterParams=params,
                                                          returnLines=True,
                                                          winnerTakeAll=True,
                                                          taxonomy=taxonomy,
                                                          hitStringMap=hitStringMap,
                                                          parseStyle=options.parseStyle)

            for line in readHits:
                outhandle.write(line)

    else:
        # process all files at once
        (multifile,readFileDict) = redistribute.multipleFileWrapper(args, params, returnLines=True)

        # Build a map from input file name to output handle
        outputMap={}
        for infileName in args:
            if options.outfile is None:
                outputMap[infileName]=sys.stdout
            elif len(args)<=1:
                outputMap[infileName]=open(options.outfile,'w')
            else:
                # use outfileName as suffix
                if options.cwd:
                    # strip path info first
                    (infilePath,infileFile)=os.path.split(infileName)
                    outfile="./"+infileFile+options.outfile
                else:
                    outfile=infileName+options.outfile
                outputMap[infileName]=open(outfile,'w')


        if options.countMethod == 'tophit':
            # don't give any taxonomy, just map to accessions for redistribution
            readHits = redistribute.pickBestHitByAbundance(multifile,
                                                  filterParams=params,
                                                  returnLines=False,
                                                  winnerTakeAll=True,
                                                  parseStyle=options.parseStyle)
        else:
            # translate to organism before finding most abundant
            readHits = redistribute.pickBestHitByAbundance(multifile,
                                                      filterParams=params,
                                                      returnLines=False,
                                                      winnerTakeAll=True,
                                                      taxonomy=taxonomy,
                                                      hitStringMap=hitStringMap,
                                                      parseStyle=options.parseStyle)


        for (read, hit) in readHits:
            outhandle = outputMap[readFileDict[read]]
            outhandle.write(hit.line)

        if options.outfile is not None:
            for outhandle in outputMap.itervalues():
                outhandle.close()

if __name__ == '__main__':
    main()

