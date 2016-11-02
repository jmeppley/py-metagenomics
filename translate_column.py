#!/usr/bin/env python
"""
This script takes a tab delimited text table and creates a new column from
 an existing one using an external translation table.

based on stats/columnMaker.py
author: John Eppley
"""
import sys
import re
import os.path
import traceback
import logging
import argparse
from edl.util import *


def main():
    # set up CLI
    description = """
    This script takes a tab delimited text table and creates a new column
    from an existing one using an external translation table.
    """

    parser = argparse.ArgumentParser(description=description)
    add_IO_arguments(parser)
    parser.add_argument(
        "-f",
        "--fillMissing",
        dest="fill",
        metavar="FILL",
        help="Put FILL in column when value not in map. If not used, "
             "entire line is skipped. If set to 'KEY', value in key "
             "column is used."),
    parser.add_argument("-m", "--mapFile", dest="mapFile",
                        metavar="MAPFILE", help="Location of mapping table.")
    parser.add_argument(
        "-c",
        "--column",
        dest="col",
        type=int,
        default=1,
        help="Column number (first column is 1)",
        metavar="COLUMN")
    parser.add_argument(
        "-C",
        "--newColumn",
        dest="newcol",
        type=int,
        default=None,
        help="Column number to insert new column after. Default is the "
             "after the source column. 0=>make it the first column. "
             "-1=>make it the last column.",
        metavar="COLUMN")
    parser.add_argument(
        "-D",
        "--deleteColumn",
        dest="delcols",
        default=[],
        action='append',
        metavar='COLUMN',
        help="Delete this column (starting at 1, after new column "
             "inserted). May be used multiple times for multiple columns")

    # log level and help
    add_universal_arguments(parser)
    arguments = parser.parse_args()
    setup_logging(arguments)

    logging.info("Value map from: " + arguments.mapFile)
    logging.debug("Fill: '%s'" % (arguments.fill))

    translation = parseMapFile(arguments.mapFile)

    for (inhandle, outhandle) in inputIterator(arguments):
        # setup some counters
        ncols = 0
        total_lines = 0
        skipped_lines = 0
        lines_kept = 0
        first_invalid_line = 0
        invalid_line = None

        # loop over lines
        for i, line in enumerate(inhandle):
            total_lines += 1
            line = line.rstrip('\r\n')
            if not line or line.startswith('#'):
                skipped_lines += 1
                continue
            try:
                cells = line.split('\t')
                if ncols == 0:
                    # count columns and check requested column number
                    ncols = len(cells)
                    if arguments.col > ncols:
                        sys.exit("first line has fewer columns (%d) "
                                 "than requested column number(%d)!" %
                                 (ncols, arguments.col))

                # get value from column
                value = cells[arguments.col - 1]
                if value in translation:
                    newCol = translation[value]
                else:
                    if arguments.fill is not None:
                        if arguments.fill == 'KEY':
                            newCol = value
                        else:
                            newCol = arguments.fill
                    else:
                        logging.debug(
                            "skipping value not in translation: %s" %
                            (value))
                        skipped_lines += 1
                        continue

                # insert new value
                if arguments.newcol is None:
                    cells.insert(arguments.col, newCol)
                elif arguments.newcol < 0 or arguments.newcol >= ncols:
                    cells.append(newCol)
                else:
                    cells.insert(arguments.newcol, newCol)

                # perform any requested column deletions
                for delcol in sorted(
                        [int(c) for c in arguments.delcols], reverse=True):
                    cells.pop(delcol - 1)

                new_line = '\t'.join(cells)
                # print >> outhandle, new_line
                print(new_line, file=outhandle)
                lines_kept += 1
            except:
                logging.warn(
                    "Unexpected error (%s): %s" %
                    (sys.exc_info()[0], sys.exc_info()[1]))
                logging.warn(traceback.format_tb(sys.exc_info()[2]))
                logging.warn('Skipping %s' % (line))
                skipped_lines += 1
                if not invalid_line:
                    first_invalid_line = i + 1
                    invalid_line = line

        # set insertion column for logging
        if arguments.newcol is None:
            inserted = arguments.col + 1
        elif arguments.newcol < 0 or arguments.newcol >= ncols:
            inserted = ncols
        else:
            inserted = arguments.newcol + 1

        valid_lines = total_lines - skipped_lines
        message = "Processed: %s\nCreated column %d with mapfile %s \
applied to column %d" % (inhandle,
                         inserted,
                         arguments.mapFile,
                         arguments.col)
        if valid_lines > 0:
            message += '\nkept %d of %d lines.' % (lines_kept, total_lines)
        else:
            message += '\nNo valid lines found'
        if skipped_lines > 0:
            message += '\nSkipped %d lines' % (skipped_lines)
        logging.info(message)
        if invalid_line:
            logging.warn(
                'Invalid lines found! EG line #%d: "%s"' %
                (first_invalid_line, invalid_line))

if (__name__ == '__main__'):
    main()
