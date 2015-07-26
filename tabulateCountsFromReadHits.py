#! /usr/bin/python
"""
"""

from optparse import OptionParser
import sys, re, urllib2, httplib
from edl.hits import parseHits
from edl.util import tupleIteratorToMap

def main():
    usage = "usage: %prog TABLE_1 TABLE_2"
    description = """
Takes two tables of read hits. Produces a table of hit counts with hit types from the first file as rows and the hit types from the second file as columns.

    """
    parser = OptionParser(usage, description=description)
    parser.add_option("-o", "--outfile", dest="outfile",
                      metavar="OUTFILE", help="Write count table to OUTFILE")
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default=False,
                      help="Print status messages to stderr")
    parser.add_option("-H", "--hitCol1", dest="hitCol1", type="int", default=-1,
            help="Index (starting at 0) of column in file 1 with hit name, -1 is default meaning all columns that are not the read name are hit names.",
                  metavar="HITCOL")
    parser.add_option("-I", "--hitCol2", dest="hitCol2", type="int", default=-1,
            help="Index (starting at 0) of column in file 2 with hit name, -1 is default meaning all columns that are not the read name are hit names.",
                  metavar="HITCOL")
    parser.add_option("-S", "--skipFirstRow", action="store_true", default=False, help="hit tables have a header row which needs to be skipped")
    parser.add_option("-A", "--about",
                      action="store_true", dest="about", default=False,
                      help="Print description")

    (options, args) = parser.parse_args()

    # output for scripts page
    if options.about:
        print description
        exit(0)

    # check arguments
    if options.verbose:
        global verbose
        verbose=True

    if len(args) != 2:
        options.error("Please supply two input files")

    (file1,file2)=args
    log("reading hits from %s" % (file1))
    fhandle = open(file1)
    hits1=parseHits(fhandle,0,options.hitCol1,options.skipFirstRow, None)
    log("reading hits from %s" % (file2))
    fhandle = open(file2)
    hits2=parseHits(fhandle,0,options.hitCol2,options.skipFirstRow, None)

    hits1=tupleIteratorToMap(hits1)
    hits2=tupleIteratorToMap(hits2)

    fhandle.close()
    fhandle.close()

    log("counting hits")
    (table,cols) = combineCounts(hits1,hits2)

    # print out hit table
    if options.outfile is None:
        log("printing table to stdout")
        outhandle = sys.stdout
    else:
        log("printing table to options.outfile")
        outhandle = open(options.outfile,'w')

    printTable(outhandle,table,cols)


##############
# Classes    #
##############

################
# compiled REs #
################

#############
# Functions #
#############
verbose=False
def log(msg):
    if verbose:
        sys.stderr.write(msg)
        sys.stderr.write("\n")

def die( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def warn(msg):
    sys.stderr.write("WARNING: %s\n" % (msg))

def combineCounts(hits1,hits2):
    counts={}
    types2={}
    for (read,hit1) in hits1.iteritems():
        try:
            hit2=hits2.pop(read)
        except KeyError:
            warn("%s not in file 2" % (read))
            continue

        if type(hit1) is not type([]):
            hit1=[hit1,]
            logging.warn('DEBUG: changing hit1 for %s to a list')

        if type(hit2) is not type([]):
            hit2=[hit2,]
            logging.warn('DEBUG: changing hit2 for %s to a list')

        for h2 in hit2:
            for h1 in hit1:
                h1counts=counts.setdefault(h1,{})
                h1counts[h2]=h1counts.setdefault(h2,0)+1
            types2[h2]=1

    return (counts,types2.keys())

def printTable(output, table, cols=None):
    if cols==None:
        cols=[]
        for row in table.values():
            cols.extend(row.keys())
        cols=set(cols)

    cols=sorted(cols)

    # print header row
    output.write("\t")
    output.write("\t".join(cols))
    output.write("\n")

    # print data
    for key in sorted(table.keys()):
        row = table[key]
        output.write(key)
        for col in cols:
            output.write("\t")
            output.write(str(row.get(col,"0")))
        output.write("\n")

if __name__ == '__main__':
    main()

