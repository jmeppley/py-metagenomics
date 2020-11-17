#!/usr/bin/env python
"""
Take list of sequence ids (eg read names) from file. Remove these (or all
but these) from each of a list of sequences files. Sequence files may
be fasta or fastq (or theoretically any biopython format with an 'id'
field).
"""
from Bio import SeqIO
import re
import sys
import logging
from edl.util import *


def main():
    import argparse

    # set up CLI
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    add_IO_arguments(parser)
    add_screen_arguments(parser, accs=True)
    parser.add_argument(
        "-f",
        "--formatIn",
        default='fasta',
        help="Any format string supported by biopython SeqIO. Default: fasta")
    parser.add_argument(
        "-F",
        "--formatOut",
        default=None,
        help="Any format string supported by biopython SeqIO. Default: "
             "same as input")

    add_universal_arguments(parser)
    arguments = parser.parse_args()
    setup_logging(arguments)

    # check formats
    if arguments.formatOut is None:
        arguments.formatOut = arguments.formatIn
    logging.info("Reading %s and writing %s" %
                 (arguments.formatIn, arguments.formatOut))
    if arguments.accs:
        logging.info("Stripping everything out except accessions")

    # get read list
    read_set = get_screen_list(arguments, accs=arguments.accs)
    numReads = len(read_set)
    logging.info("Loaded list of %d reads" % (numReads))

    # seqio wants a fully functional file object, my wrapper doesn't work
    for (fastaIn, fastaOut) in inputIterator(arguments,
                                             use_file_wrapper=False):
        scanFileForReads(
            read_set,
            fastaIn,
            fastaOut,
            arguments.keep,
            arguments.accs,
            arguments.formatIn,
            arguments.formatOut)

    if len(read_set) > 0:
        logging.warn(
            "%d of %d items in list were not found in fasta data" %
            (len(read_set), numReads))
        logging.debug(read_set)


verbose = False


def log(msg):
    if verbose:
        sys.stderr.write(msg)
        sys.stderr.write("\n")


def die(msg):
    sys.stderr.write("INFO: %s\n" % msg)
    sys.exit(1)


def warn(msg):
    sys.stderr.write("WARNING: %s\n" % (msg))


readRE = re.compile(r'^(\S+)')


def scanFileForReads(reads, fastaIn, fastaOut, keep, acc, formatIn, formatOut):
    """
    Given:
        a set of reads
        an input file and format
        keep: True->keep listed reads, False->drop them
        acc: True->Parse accesions from read names to match names in dict
        output file and format
    Print records from input file to output file
    """
    records = SeqIO.parse(fastaIn, formatIn)
    count = 0
    for record in records:
        read = record.id

        if acc:
            # parse accession out of read name
            read = parseAcc(read)

        # is read in list?
        try:
            reads.remove(read)
            match = True
        except KeyError:
            match = False

        if match == keep:
            count += 1
            # write line if either
            #  - read matches list AND keep is True
            #  - read not in list AND keep is False
            SeqIO.write((record,), fastaOut, formatOut)

    logging.info("Printed %d records" % count)


if __name__ == "__main__":
    main()
