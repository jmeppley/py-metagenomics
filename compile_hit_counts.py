#! /usr/bin/python
"""
Takes two tables of read hits. Produces a table of hit counts with hit types from the first file as rows and the hit types from the second file as columns.
"""

import argparse
import sys, re, urllib2, httplib
from edl.hits import parseHits
from edl.util import tupleIteratorToMap, add_universal_arguments, setup_logging

def main():
    description = __doc__
    parser = argparse.ArgumentParser(description)
    parser.add_argument("input_files", nargs=2, 
                        type=argparse.FileType('r'),
                        metavar=("INPUT_TABLE_1","INPUT_TABLE_2"),
                        help="Input files. Please supply two tables")
    parser.add_argument("-o", "--outfile", dest="outfile",
                      metavar="OUTFILE", help="Write count table to OUTFILE")
    parser.add_argument("-L", "--long_output", default=False,
            action="store_true",
            help="Print one number per row (prefixed by two keys) instead of a table with one seet of keys as column names and one set as row names.")
    parser.add_argument("-H", "--hitCol1", dest="hitCol1", type=int, default=-1,
            help="Index (starting at 0) of column in file 1 with hit name, -1 is default meaning all columns that are not the read name are hit names.",
                  metavar="HITCOL")
    parser.add_argument("-I", "--hitCol2", dest="hitCol2", type=int, default=-1,
            help="Index (starting at 0) of column in file 2 with hit name, -1 is default meaning all columns that are not the read name are hit names.",
                  metavar="HITCOL")
    parser.add_argument("-S", "--skipFirstRow", action="store_true", default=False, help="hit tables have a header row which needs to be skipped")

    add_universal_arguments(parser)
    arguments = parser.parse_args()
    setup_logging(arguments)

    if len(arguments.input_files) != 2:
        parser.error("Please supply two input files")

    (file1,file2)=arguments.input_files
    log("reading hits from %s" % (file1))
    hits1=parseHits(file1,0,arguments.hitCol1,arguments.skipFirstRow, None)
    log("reading hits from %s" % (file2))
    hits2=parseHits(file2,0,arguments.hitCol2,arguments.skipFirstRow, None)

    hits1=tupleIteratorToMap(hits1)
    hits2=tupleIteratorToMap(hits2)

    log("counting hits")
    (table,cols) = combineCounts(hits1,hits2)

    # print out hit table
    if arguments.outfile is None:
        log("printing table to stdout")
        outhandle = sys.stdout
    else:
        log("printing table to arguments.outfile")
        outhandle = open(arguments.outfile,'w')

    printTable(outhandle,table,cols, long_output=arguments.long_output)


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

def combineCounts(hits1,hits2, unmatched_1="Unknown", unmatched_2="Unknown"):
    # compile counts into nested dicts
    counts={}

    # keep track of all unique hit ids from set 2
    types2=set()

    # Start by using reads from hits1 as index
    for (read,hit_list1) in hits1.items():
        hit_list2=hits2.pop(read, [unmatched_2,])

        for hit2 in hit_list2:
            for hit1 in hit_list1:
                h1counts=counts.setdefault(hit1,{})
                h1counts[hit2]=h1counts.get(hit2,0)+1
            types2.add(hit2)

    # remaining reads in hits2 counted as Unknown
    # we know these don't exist in hits1
    h1counts = counts.setdefault(unmatched_1,{})
    for read, hit_list2 in hits2.items():
        for hit2 in hit_list2:
            h1counts[hit2] = h1counts.get(hit2,0)+1
            types2.add(hit2)

    return (counts,types2)

def printTable(output, table, cols=None, long_output=False):
    if long_output:
        # print one row per value
        for row_name, row in table.items():
            for col_name, value in row.items():
                output.write("\t".join([row_name, col_name, str(value)]) + "\n")

    else:
        # We need to know the column names up front
        if cols==None:
            cols = set.union(*[r.keys() for r in table.values()])

        # Lets sort them
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

