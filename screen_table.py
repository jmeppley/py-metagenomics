#!/usr/bin/python
"""
Take list of reads (or any names) from file. Remove these (or all but these) from each of a list of text tables (e.g. m8 files).

Table to screen can also be piped through STDIN.

To identify reads hitting specific sequences, use screenHitTable.
"""

import sys, logging
from edl.util import *

def main():
    import argparse

    ## set up CLI
    usage = "usage: %prog -l LIST [OPTIONS] TABLE(S)"
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    add_screen_arguments(parser, accs=True)
    add_IO_arguments(parser)
    parser.add_argument("-d", "--delim", dest="delim", default="\t",
                      help="Input table delimiter (tab is default). If set to 'None', split on any whitespace.", metavar="DELIM")
    parser.add_argument("-c", "--col", dest="col", type=int, default=0,
                      help="Column to screen (0 is default)", metavar="INDEX")

    add_universal_arguments(parser)

    arguments = parser.parse_args()

    setup_logging(arguments)

    # allow funky characters in delim arguments
    if arguments.delim == 'None':
        arguments.delim=None
    elif arguments.delim != '\t':
        arguments.delim=arguments.delim.decode('string-escape')
    if arguments.listDelim is not None:
        arguments.listDelim = arguments.listDelim.decode('string-escape')

    # get read list
    logging.debug("List file: '%s'\nList delim: '%s'" % (arguments.listFile, arguments.listDelim))
    screen_set = get_screen_list(arguments, accs=arguments.accs)
    logging.debug("Got list of %d reads" % (len(screen_set)))
    if len(screen_set)>0:
        logging.debug("For example: %s" % (next(iter(screen_set))))

    for (inhandle,outhandle) in inputIterator(arguments):
        scanFileForReads(screen_set, inhandle, arguments.keep, outhandle, arguments.delim, arguments.col, arguments.accs)

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
