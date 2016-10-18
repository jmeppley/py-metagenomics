#!/usr/bin/env python
# This script takes a tab delimited text table and creates a new column from
#  an existing one using an external translation table.
##
# based on stats/columnMaker.py
# author: John Eppley
import sys, re, os.path
import traceback, logging
from edl.util import *

from optparse import OptionParser

def main():
## set up CLI
    description = """
    This script takes a tab delimited text table and creates a new column from an existing one using an external translation table.
    """

    parser = OptionParser(description=description)
    addIOOptions(parser)
    parser.add_option("-f", "--fillMissing", dest="fill",
                      metavar="FILL", help="Put FILL in column when value not in map. If not used, entire line is skipped. If set to 'KEY', value in key column is used."),
    parser.add_option("-m", "--mapFile", dest="mapFile",
                      metavar="MAPFILE", help="Location of mapping table.")
    parser.add_option("-c", "--column", dest="col", type="int", default=1,
                      help="Column number (first column is 1)", metavar="COLUMN")
    parser.add_option("-C", "--newColumn", dest="newcol", type="int", default=None,
                      help="Column number to insert new column after. Default is the after the source column. 0=>make it the first column. -1=>make it the last column.", metavar="COLUMN")
    parser.add_option("-D", "--deleteColumn", dest="delcols", default=[], action='append',
                      metavar='COLUMN',
                      help="Delete this column (starting at 1, after new column inserted). May be used multiple times for multiple columns")
    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options,description)

    logging.info ("Value map from: " + options.mapFile)
    logging.debug("Fill: '%s'" % (options.fill))

    translation = parseMapFile(options.mapFile)

    for (inhandle,outhandle) in inputIterator(args,options):
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
            line = line.rstrip( '\r\n' )
            if not line or line.startswith( '#' ):
                skipped_lines += 1
                continue
            try:
                cells = line.split( '\t' )
                if ncols==0:
                    # count columns and check requested column number
                    ncols = len(cells)
                    if options.col>ncols:
                        sys.exit( "first line has fewer columns (%d) than requested column number(%d)!" % (ncols,options.col))

                # get value from column
                value = cells[options.col-1]
                if value in translation:
                    newCol = translation[value]
                else:
                    if options.fill is not None:
                        if options.fill == 'KEY':
                            newCol = value
                        else:
                            newCol = options.fill
                    else:
                        logging.debug( "skipping value not in translation: %s"% (value))
                        skipped_lines += 1
                        continue

                # insert new value
                if options.newcol is None:
                    cells.insert(options.col,newCol)
                elif options.newcol<0 or options.newcol>=ncols:
                    cells.append(newCol)
                else:
                    cells.insert(options.newcol,newCol)

                # perform any requested column deletions
                for delcol in sorted([int(c) for c in options.delcols], reverse=True):
                    cells.pop(delcol-1)

                new_line = '\t'.join(cells)
                print >> outhandle, new_line
                lines_kept += 1
            except:
                logging.warn( "Unexpected error (%s): %s"% (sys.exc_info()[0],sys.exc_info()[1]))
                logging.warn( traceback.format_tb(sys.exc_info()[2]))
                logging.warn( 'Skipping %s' %(line))
                skipped_lines += 1
                if not invalid_line:
                    first_invalid_line = i + 1
                    invalid_line = line

        # set insertion column for logging
        if options.newcol is None:
            inserted = options.col+1
        elif options.newcol<0 or options.newcol>=ncols:
            inserted = ncols
        else:
            inserted = options.newcol+1

        valid_lines = total_lines - skipped_lines
        message="Processed: %s\nCreated column %d with mapfile %s applied to column %d" % ( inhandle, inserted, options.mapFile, options.col )
        if valid_lines > 0:
            message+= '\nkept %d of %d lines.' % ( lines_kept, total_lines )
        else:
            message+= '\nNo valid lines found'
        if skipped_lines > 0:
            message+= '\nSkipped %d lines' % (skipped_lines)
        logging.info(message)
        if invalid_line:
            logging.warn('Invalid lines found! EG line #%d: "%s"' % ( first_invalid_line, invalid_line ))

if (__name__=='__main__'):
    main()
