#! /usr/bin/env python
"""
add_kegg_pathway.py

Given a tabular file with a KO column, add a new column with KEGG pathways
listed.
"""

import sys
import re
import os.path
import traceback
import logging
import argparse
from edl import kegg, util


def main():
    # set up CLI
    description = """
    Takes a tabular text file and translate a column of KO values into a new
    column of KEGG pathways. KO column can have multiple entries per row.
    Output column will have multiple entries per pathway cell.
    """

    parser = argparse.ArgumentParser(description=description)
    util.add_IO_arguments(parser)
    parser.add_argument("-l", "--level", dest="level",
                        default="PATHWAY",
                        metavar="LEVEL",
                        help=""" Level to collect counts on. Level can
                        be one of:
                          NAME, PATHWAY, EC, DEFINITION,
                          or a level in the CLASS heirachy: 1, 2, or 3
                        """)
    parser.add_argument(
        "-f",
        "--fill_missing",
        dest="fill",
        metavar="FILL",
        default="No Pathway",
        help="Put FILL in column when KO has no pathway assigned. "
             "Defaults to 'No Pathway'")
    # format, ortholog heirarchy, and more
    parser.add_argument("-k", "--ko_file", dest="ko_file",
                        metavar="MAPFILE", help="Location of kegg ko file")
    parser.add_argument(
        "-c",
        "--ko_column",
        type=int,
        default=1,
        help="Column number (first column is 1)",
        metavar="COLUMN")
    parser.add_argument(
        "-C",
        "--new_column",
        type=int,
        default=None,
        help="Column number to insert new column after. Default is the "
             "after the source column. 0=>make it the first column. "
             "-1=>make it the last column.",
        metavar="COLUMN")
    parser.add_argument(
        "-L",
        "--long_output",
        default=False,
        action='store_true',
        help="Insert new duplicate rows if KO maps to multiple values")
    parser.add_argument(
            "-H", "--header",
            default=None,
            metavar='HEADER',
            help="Put HEADER in first row instead of trying to translate")
    parser.add_argument(
            "-Q", "--quotes",
            default=False,
            action="store_true",
            help="Encase translated values in double quotes")
    parser.add_argument(
            "-s", "--sep",
            dest='sep',
            default='\t',
            help="""Character separating table cells. Default is tab""")
    parser.add_argument(
            "-S", "--ko_sep",
            dest='ko_sep',
            default=';',
            help="""Character separating multiple KO values in iput table
            and used to separate multiple values in output column.
            Default is ";". Ignored for output if --longOutput
            requested""")

    # log level and help
    util.add_universal_arguments(parser)
    arguments = parser.parse_args()
    util.setup_logging(arguments)

    logging.info("KO mapping from: " + arguments.ko_file)
    logging.debug("Fill: '%s'" % (arguments.fill))

    translation = kegg.parse_KEGG_file(arguments.ko_file,
                                       arguments.level)

    # switch to zero indexing
    if arguments.new_column:
        arguments.new_column -= 1
    arguments.ko_column -= 1

    for (inhandle, outhandle) in util.inputIterator(arguments):
        for new_line in translate_ko_column(
                inhandle,
                sep=arguments.sep,
                ko_sep=arguments.ko_sep,
                ko_column=arguments.ko_column,
                new_column=arguments.new_column,
                translation=translation,
                default=arguments.fill,
                quotes=arguments.quotes,
                header=arguments.header,
                long_out=arguments.long_output, ):
            outhandle.write(new_line)


def translate_ko_column(line_iter,
                        sep='\t',
                        ko_sep=';',
                        ko_column=1,
                        new_column=None,
                        translation={},
                        default='No Pathway',
                        long_out=False,
                        drop_empty=False,
                        header=None,
                        quotes=False,
                       ):
    """
    Given a table with a KO column
    add a column of translated KOs. EG: A list of pathway names

    Parameters:

 * translation({}): dictionary mapping KOs to new values (EG Pathways).
                    New value entries should be lists.
 * sep('\t'): The value delimiter in the table
 * ko_sep(';'): The separtor delimiting multiple KOs in one cell
 * ko_column(1): column with KO values
 * new_column(None): where to insert the translated values.
                     A column number or None(=>after KO column).
 * default('No Pathway'): The value to be inserted
                     if the KO(s) is/are not in the translation table.
 * drop_empty(False): If true drop lines with no KO
    """
    comment_lines = 0
    dropped_lines = 0
    ncols = 0

    # loop over table lines
    first_line = True
    for i, line in enumerate(line_iter):
        line = line.rstrip('\r\n')
        if not line or line.startswith('#'):
            comment_lines += 1
            continue
        cells = line.split(sep)
        if ncols == 0:
            # count columns and check requested column number
            ncols = len(cells)

        # don't translate header line (if asked)
        if first_line and header:
            added_values = [header, ]
        else:
            # get value from column
            values = [c.strip('" ') for c in cells[ko_column].split(ko_sep)]
            if len(values) == 0 and drop_empty:
                # skip this line
                dropped_lines += 1
                continue

            # translate values, account for mappings to multiple things
            added_values = set()
            for value in values:
                translated_values = translation.get(value, [default, ])
                if isinstance(translated_values, str):
                    translated_values = (translated_values,)
                for translated_value in translated_values:
                    added_values.add(translated_value)

            # sort added values so tests are reproducible
            added_values = sorted(added_values)
        first_line = False

        # insert new value(s)
        if not long_out:
            # make it a one element list with all values concatenated
            added_values = [ko_sep.join(added_values), ]

        # generate a line for each element in list
        for added_value in added_values:
            # slow, but safe. Revisit if we have performance issues
            if quotes:
                added_value = '"{}"'.format(added_value)
            new_cells = list(cells)
            if new_column is None:
                new_cells.insert(ko_column + 1, added_value)
            elif new_column < 0 or new_column >= ncols:
                new_cells.append(added_value)
            else:
                new_cells.insert(new_column, added_value)

            yield sep.join(new_cells) + '\n'

    logging.info("Skipped %d comments and %d empty lines out of %d",
                 comment_lines, dropped_lines, i + 1)


if __name__ == '__main__':
    main()
