#!/bin/env python
"""
Generates an ascii histogram from list of values
"""
import sys
import logging
import argparse
from numpy import histogram
from edl.util import add_universal_arguments, setup_logging, ascii_histogram

def main():
    # set up CLI
    description = __doc__

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'),
            default=sys.stdin,
            help="Input file (defaults to STDIN) containing a value " + \
                 "on each line")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
            default=sys.stdout,
            help="File to which to write histogram")
    parser.add_argument('-b', '--bins', type=int, default=50)
    parser.add_argument('-l', '--label', default="value")
    parser.add_argument('-w', '--width', default=75, type=int)
    parser.add_argument('-L', '--log', action='store_true')
    parser.add_argument('-W', '--max-label-width', type=int, default=10)

    # log level and help
    add_universal_arguments(parser)
    arguments = parser.parse_args()
    setup_logging(arguments)

    values = []
    for line in arguments.infile:
        try:
            values.append(float(line.strip()))
        except ValueError:
            if len(line.strip())>0:
                logging.warn("Skipping bad line:\n{}".format(line.strip()))
    logging.info("Read in {} values".format(len(values)))
    arguments.outfile.write(ascii_histogram(histogram(values,
                                                      bins=arguments.bins),
                                            label=arguments.label,
                                            width=arguments.width,
                                            log=arguments.log,
                                            maxLabelWidth=\
                                                    arguments.max_label_width
                                           ))

if __name__ == "__main__":
    main()
