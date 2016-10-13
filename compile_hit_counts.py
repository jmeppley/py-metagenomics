#! /usr/bin/env python
"""
Takes two tables of read hits. Produces a table of hit counts with hit types from the first file as rows and the hit types from the second file as columns.
"""

import argparse, logging
import sys, re
from edl.hits import parseHits
from edl.util import tupleIteratorToMap, add_universal_arguments, setup_logging

def main():
    description = __doc__
    parser = argparse.ArgumentParser(description)
    parser.add_argument("-1","--input_file_1",
                        default=None,
                        type=argparse.FileType('r'),
                        metavar=("INPUT_TABLE_1"),
                        help="Input table 1")
    parser.add_argument("-2","--input_file_2",
                        default=None,
                        type=argparse.FileType('r'),
                        metavar=("INPUT_TABLE_2"),
                        help="Input table 2")
    parser.add_argument("-o", "--outfile", dest="outfile",
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        metavar="OUTFILE", 
                        help="Write count table to OUTFILE. (Defaults to STDOUT")
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

    if arguments.input_file_1 is None or arguments.input_file_2 is None:
        parser.error("Please supply two input files")

    logging.info("reading hits from %s" % (arguments.input_file_1.name))
    hits1=parseHits(arguments.input_file_1,
                    0,
                    arguments.hitCol1,
                    arguments.skipFirstRow, 
                    None)
    logging.info("reading hits from %s" % (arguments.input_file_2.name))
    hits2=parseHits(arguments.input_file_2,
                    0,
                    arguments.hitCol2,
                    arguments.skipFirstRow, 
                    None)

    hits1=tupleIteratorToMap(hits1)
    hits2=tupleIteratorToMap(hits2)

    logging.info("counting hits")
    (table,cols) = combineCounts(hits1,hits2)

    # print out hit table
    logging.info("printing table to " + arguments.outfile.name)
    printTable(arguments.outfile,table,cols, long_output=arguments.long_output)

##############
# Classes    #
##############

################
# compiled REs #
################

#############
# Functions #
#############
def die( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

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
        for row_name in sorted(table.keys()):
            row = table[row_name]
            for col_name in sorted(row.keys()):
                value = row[col_name]
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

