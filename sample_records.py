#!/usr/bin/python

#################
# This script takes any file that can be divided into records and 
#  returns N randomly selected records
#
# Records can be fasta, fastq, genbank, or something described by a simple RegExp
################

from optparse import OptionParser
from os import path
from edl.util import *
from edl.batch import *

import re, sys

def main():
## set up CLI
    description = """
 This script takes any file that can be divided into records and 
  returns N randomly selected records. 

  NOTE: 
  All sampled records are stored in memory. This requires a good amount of RAM (depending on record size and sample size). 

     Records can be fasta, fastq, genbank, or something described by a simple RegExp
    """

    parser = OptionParser(description=description)
    addIOOptions(parser)
    addRecordParsingOptions(parser)

    parser.add_option("-s", "--sample_size", default=1000000, type=int,
            help="Number of records to pull out. Defaults to 1 million.")

    addUniversalOptions(parser)

    (options, cmdargs) = parser.parse_args()

    setupLogging(options, description)

    for inhandle, outhandle in inputIterator(cmdargs, options):
        # We need the file name to ge the type, get from handle (if not stdin)
        infilename = inhandle.name if len(cmdargs)>0 else None
        fileType = getFileType(options, infilename)
        record_iterator = fileType.recordStreamer(inhandle)

        # read file once with reservoir sampling
        logging.debug("Looking for %d records in %s" % (options.sample_size,
                                                        infilename))
        sampled_records, record_count = reservoirSample(record_iterator,
                                                        N=options.sample_size,
                                                        returnCount=True)
        logging.debug("Sampled %d of %d records" % (len(sampled_records),
                                                    record_count))
        for record in sampled_records:
            outhandle.writelines(record)

if __name__ == '__main__':
    main()
