#!/usr/bin/python
"""
Take list of reads (or any names) from file. Remove these (or all but these) from each of a list of text tables (e.g. m8 files).

Table to screen can also be piped through STDIN.

To identify reads hitting specific sequences, use screenHitTable.
"""

import sys, logging
from edl.util import *

def main():
    from optparse import OptionParser

    ## set up CLI
    usage = "usage: %prog -l LIST [OPTIONS] TABLE(S)"
    description = __doc__
    parser = OptionParser(description=description)
    addScreenOptions(parser, accs=True)
    addIOOptions(parser)
    parser.add_option("-d", "--delim", dest="delim", default="\t",
                      help="Input table delimiter (tab is default). If set to 'None', split on any whitespace.", metavar="DELIM")
    parser.add_option("-c", "--col", dest="col", type='int', default=0,
                      help="Column to screen (0 is default)", metavar="INDEX")
    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options,description)

    # allow funky characters in delim arguments
    if options.delim == 'None':
        options.delim=None
    elif options.delim != '\t':
        options.delim=options.delim.decode('string-escape')
    if options.listDelim is not None:
        options.listDelim = options.listDelim.decode('string-escape')

    # get read list
    logging.debug("List file: '%s'\nList delim: '%s'" % (options.listFile, options.listDelim))
    readDict = getScreenList(options, accs=options.accs)
    logging.debug("Got list of %d reads" % (len(readDict)))
    if len(readDict)>0:
        logging.debug("For example: %s" % (readDict.iterkeys().next()))

    for (inhandle,outhandle) in inputIterator(args, options):
        scanFileForReads(readDict, inhandle, options.keep, outhandle, options.delim, options.col, options.accs)

################
# Functions
################
def die( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit(1)

def scanFileForReads(reads,inhandle,keep,outhandle,delim,col,accs):
    lineCount=0
    matchCount=0
    if keep:
        logging.info("Keeping matched reads")
    else:
        logging.info("Discarding matched reads")
    for line in inhandle:
        lineCount+=1
        read = line.rstrip('\r\n').split(delim)[col]
        if accs:
            read=parseAcc(read)
        logging.debug("looking for %s in %s" % (read,line))

        match = read in reads
        if match==keep:
            # write line if either
            #  - read matches list AND keep is True
            #  - read not in list AND keep is False
            outhandle.write(line)
            matchCount+=1

    logging.info("Kept %d of %d lines" % (matchCount,lineCount))

if __name__ == "__main__":
    main()