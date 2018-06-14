#! /usr/bin/env python
"""
Takes two tables of read hits. Produces a table of hit counts with hit
types from the first file as rows and the hit types from the second file
as columns.

There are two output modes. by default, a full text table of hits is printed.
The long-form output can be invoked with -L. This prints out a three column
(double-index) table with the hit IDs from the two input files as the first
two columns and the combined counts as the third column. EG:

    org_a    gene_1   10
    org_a    gene_2   20
    org_b    gene_2    3

The user can supply a multiplier file that supplies a mulitplier
value for each read or gene in the input tables. This is useful
for contig based counts where you have coverage or abundance data.
The multiplier file is expected to be a two column tab-separated file.

In long form, if a mlutiplier file is given, four columns are printed. The two
data columns are the raw counts and multiplied counts. EG:

    org_a    gene_1   10   35.0
    org_a    gene_2   20   70.0
    org_b    gene_2    3    6.0

"""

import argparse
import logging
import sys
from edl.hits import parseHits
from edl.util import tupleIteratorToMap,\
                     add_universal_arguments,\
                     setup_logging,\
                     parseMapFile


def main():
    """ set up the command line interface """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-1", "--input_file_1",
                        default=None,
                        type=argparse.FileType('r'),
                        metavar=("INPUT_TABLE_1"),
                        help="Input table 1")
    parser.add_argument("-2", "--input_file_2",
                        default=None,
                        type=argparse.FileType('r'),
                        metavar=("INPUT_TABLE_2"),
                        help="Input table 2")
    parser.add_argument("-m", "--multiplier",
                        default=None,
                        metavar=("MULTIPLIER_TABLE"),
                        help=("Table of values to multiply each sequence. "
                              "EG assembly coverages."))
    parser.add_argument("-T", "--total_reads",
                        default=0,
                        metavar="TOTAL_READS",
                        type=int,
                        help="Total number of reads to expect. (This allows "
                             "the reporting of unknown read count)")
    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        metavar="OUTFILE",
        help="Write count table to OUTFILE. (Defaults to STDOUT")
    parser.add_argument(
        "-L",
        "--long_output",
        default=False,
        action="store_true",
        help="Print one number per row (prefixed by two keys) instead "
             "of a table with one seet of keys as column names and one "
             "set as row names.")
    parser.add_argument(
        "-H",
        "--hitCol1",
        dest="hitCol1",
        type=int,
        default=-1,
        help="Index (starting at 0) of column in file 1 with hit name, -1 "
             "is default meaning all columns that are not the read name are "
             "hit names.",
        metavar="HITCOL")
    parser.add_argument(
        "-I",
        "--hitCol2",
        dest="hitCol2",
        type=int,
        default=-
        1,
        help="Index (starting at 0) of column in file 2 with hit name, -1 "
             "is default meaning all columns that are not the read name "
             "are hit names.",
        metavar="HITCOL")
    parser.add_argument(
        "-S",
        "--skipFirstRow",
        action="store_true",
        default=False,
        help="hit tables have a header row which needs to be skipped")

    add_universal_arguments(parser)
    arguments = parser.parse_args()
    setup_logging(arguments)

    if arguments.input_file_1 is None or arguments.input_file_2 is None:
        parser.error("Please supply two input files")

    logging.info("reading hits from %s", arguments.input_file_1.name)
    hits1 = parseHits(arguments.input_file_1,
                      0,
                      arguments.hitCol1,
                      arguments.skipFirstRow,
                      None)
    logging.info("reading hits from %s", arguments.input_file_2.name)
    hits2 = parseHits(arguments.input_file_2,
                      0,
                      arguments.hitCol2,
                      arguments.skipFirstRow,
                      None)

    hits1 = tupleIteratorToMap(hits1)
    hits2 = tupleIteratorToMap(hits2)

    if arguments.multiplier is not None:
        multipliers = parseMapFile(arguments.multiplier, valueType=float)
    else:
        multipliers = None

    logging.info("counting hits")
    (table, cols) = combine_counts(hits1, hits2, multipliers,
                                   total_reads=arguments.total_reads)

    # print out hit table
    logging.info("printing table to " + arguments.outfile.name)
    print_table(arguments.outfile, table, cols,
                is_multiplied=multipliers is not None,
                long_output=arguments.long_output)

#############
# Functions #
#############


def combine_counts(hits1,
                   hits2,
                   multipliers=None,
                   total_reads=0,
                   unmatched_1="Unknown",
                   unmatched_2="Unknown",
                  ):
    """ compile counts into nested dicts """
    total_counted = 0
    counts = {}

    # keep track of all unique hit ids from set 2
    types2 = set()

    if multipliers is None:
        def _get_increment(read):
            return 1

        def _update_hits(increment, counts, hit1, hit2):
            """ Just add 1 to get raw numbers """
            h1counts = counts.setdefault(hit1, {})
            h1counts[hit2] = h1counts.get(hit2, 0) + increment
    else:
        # if a mult table was given, use it to get total count
        if total_reads == 0:
            total_reads = len(multipliers)

        def _get_increment(read):
            """ get multiplier. Use pop to see leftovers """
            return multipliers.pop(read, 1)

        def _update_hits(increment, counts, hit1, hit2):
            """ count both raw numbers and multiplied """
            h1counts = counts.setdefault(hit1, {})
            count_tuple = h1counts.get(hit2, (0, 0))
            count_tuple = (count_tuple[0] + 1,
                           count_tuple[1] + increment)
            h1counts[hit2] = count_tuple

    # Start by using reads from hits1 as index
    for (read, hit_list1) in hits1.items():
        # remove hits from hits2 as we go, so we know what didn't match hits1
        # default to umatched_2
        total_counted += 1
        increment = _get_increment(read)
        hit_list2 = hits2.pop(read, [unmatched_2, ])

        for hit2 in hit_list2:
            for hit1 in hit_list1:
                _update_hits(increment, counts, hit1, hit2)
            types2.add(hit2)

    # remaining reads in hits2 counted as Unknown
    # we know these don't exist in hits1
    hit1 = unmatched_1
    for read, hit_list2 in hits2.items():
        total_counted += 1
        increment = _get_increment(read)
        for hit2 in hit_list2:
            _update_hits(increment, counts, hit1, hit2)
            types2.add(hit2)

    # if a total was given
    if total_reads > 0:
        unknown_counts = counts.setdefault(unmatched_1, {})
        if multipliers is None:
            unknown_counts[unmatched_2] = total_reads - total_counted
        else:
            unknown_counts[unmatched_2] = (total_reads - total_counted,
                                           sum(multipliers.values()))

    return (counts, types2)


def print_table(output,
                table,
                cols=None,
                is_multiplied=False,
                long_output=False):
    """ Render the output """
    logging.debug([k for k in table.keys()])
    if long_output:
        # print one row per value
        for row_name in sorted(table.keys()):
            row = table[row_name]
            for col_name in sorted(row.keys()):
                value = row[col_name]
                if is_multiplied:
                    # in long form we can report raw and multiplied counts
                    value = '\t'.join([str(v) for v in value])
                output.write(
                    "\t".join([row_name, col_name, str(value)]) + "\n")

    else:
        # We need to know the column names up front
        if cols is None:
            cols = set.union(*[r.keys() for r in table.values()])

        # Lets sort them
        cols = sorted(cols)

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
                if is_multiplied:
                    # for the table, just print the multiplied value
                    value = row.get(col, (0, 0))
                    output.write(str(value[1]))
                else:
                    output.write(str(row.get(col, "0")))
            output.write("\n")


if __name__ == '__main__':
    main()
