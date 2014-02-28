#!/usr/bin/python
"""
Take list of sequence ids (eg read names) from file. Remove these (or all but these) from each of a list of sequences files. Sequence files may be fasta or fastq (or theoretically any biopython format with an 'id' field).
"""
from Bio import SeqIO
import re, sys, logging
from edl.util import *

def main():
    from optparse import OptionParser

    ## set up CLI
    description = __doc__
    parser = OptionParser(description=description)
    addIOOptions(parser)
    addScreenOptions(parser, accs=True)
    parser.add_option("-f", "--formatIn", default='fasta', help="Any format string supported by biopython SeqIO. Default: fasta")
    parser.add_option("-F", "--formatOut", default=None, help="Any format string supported by biopython SeqIO. Default: same as input")
    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options,description)

    # check formats
    if options.formatOut is None:
        options.formatOut=options.formatIn
    logging.info("Reading %s and writing %s" % (options.formatIn, options.formatOut))
    if options.accs:
        logging.info("Stripping everything out except accessions")

    # get read list
    readDict = getScreenList(options, accs=options.accs)
    numReads=len(readDict)
    logging.info("Loaded list of %d reads" % (numReads))

    for (fastaIn,fastaOut) in inputIterator(args, options):
        scanFileForReads(readDict, fastaIn, fastaOut, options.keep, options.accs, options.formatIn, options.formatOut)

    if len(readDict)>0:
        #logging.info("(%s)" % (','.join(readDict.keys())))
        logging.warn("%d of %d items in list were not found in fasta data" % (len(readDict),numReads))
        logging.debug(readDict)

verbose=False
def log(msg):
    if verbose:
        sys.stderr.write(msg)
        sys.stderr.write("\n")

def die( msg ):
    sys.stderr.write( "INFO: %s\n" % msg )
    sys.exit(1)

def warn(msg):
    sys.stderr.write("WARNING: %s\n" % (msg))

readRE=re.compile(r'^(\S+)')

def scanFileForReads(reads,fastaIn,fastaOut,keep,acc,formatIn,formatOut):
    """
    Given:
        a dict of reads
        an input file and format
        keep: True->keep listed reads, False->drop them
        acc: True->Parse accesions from read names to match names in dict
        output file and format
    Print records from input file to output file
    """
    records = SeqIO.parse(fastaIn,formatIn)
    count=0
    for record in records:
        read = record.id

        if acc:
            # parse accession out of read name
            read = parseAcc(read)

        # is read in list?
        try:
            reads.pop(read)
            match=True
        except KeyError:
            match=False

        if match==keep:
            count+=1
            # write line if either
            #  - read matches list AND keep is True
            #  - read not in list AND keep is False
            SeqIO.write((record,),fastaOut, formatOut)

    logging.info("Printed %d records" % count)

if __name__ == "__main__":
    main()
