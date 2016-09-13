#!/usr/bin/env python

#################
# This script takes any file that can be divided into records and 
#  returns N randomly selected records
#
# Records can be fasta, fastq, genbank, or something described by a simple RegExp
################

from os import path
from edl.util import *
from edl.batch import *

import re, sys, argparse

def main():
## set up CLI
    description = """
 This script takes any file that can be divided into records and 
  returns N randomly selected records. 

  NOTE: 
  By default, all sampled records are stored in memory. This requires a good amount of RAM (depending on record size and sample size). To avoid this, specify the number of records or request a count using the -n (population_size) option.

     Records can be fasta, fastq, genbank, or something described by a simple RegExp
    """

    parser = argparse.ArgumentParser(description=description)
    add_IO_arguments(parser)
    add_record_parsing_arguments(parser)

    parser.add_argument("-s", "--sample_size", default=1000000, type=int,
            metavar="SAMPLE_SIZE",
            help="Number of records to pull out. Defaults to 1 million.")
    parser.add_argument("-n", "--population_size", type=int, default=0,
            metavar="POPULATION_SIZE",
            help="Number of records in file. An integer, should be greater \
                    than the SAMPLE_SIZE, except: \
                    0 => do a separate pass to count records first \
                    -1 (default) => reservoir sample to RAM on the fly")

    add_universal_arguments(parser)

    arguments = parser.parse_args()

    setup_logging(arguments)

    # check arguments
    if arguments.input_files == [sys.stdin,] and arguments.population_size == 0:
        parser.error("We cannot count records from STDIN, please specify a \
                positive population size or use reservoir sampling (-n -1)")
    if arguments.population_size>0 and arguments.population_size<arguments.sample_size:
        parser.error("We cannot sample more records then there are in the file!")

    for inhandle, outhandle in inputIterator(arguments):
        # We need the file name to ge the type, get from handle (if not stdin)
        infilename = inhandle.name
        fileType = getFileType(arguments, infilename)
        record_iterator = fileType.recordStreamer(inhandle)

        logging.debug("Looking for %d records in %s" % (arguments.sample_size,
                                                        infilename))
        #if arguments.population_size<0:
        #    # read file once with reservoir sampling
        #    sampled_records, record_count = reservoirSample(record_iterator,
        #                                                N=arguments.sample_size,
        #                                                returnCount=True)
        #else:

        # count records if asked to
        if arguments.population_size == 0:
            arguments.population_size = get_total_size(inhandle, fileType)

        # get sampled record generator (will use reservoir if P is <0)
        sampled_records = indexed_sample_generator(record_iterator,
                N = arguments.sample_size,
                P = arguments.population_size)

        # print out sampled records
        count=0
        for record in sampled_records:
            outhandle.writelines(record)
            count+=1
        logging.debug("Sampled %d records" % (count))

if __name__ == '__main__':
    main()
