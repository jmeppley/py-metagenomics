#!/usr/bin/env python
"""
Merge a hit table into unique query/reference pairs with aggregate values for all non-overlapping hits between each pair: total length, total score, pctid

The output table will have the format:

  query  hit  total_hit_length  total_hit_score  aggregate_pctid

If the input table format does not include pctid info, the pctid column will be all zeros
"""

import sys
import logging
import argparse
from collections import defaultdict
from edl.blastm8 import filterM8Stream, add_hit_table_arguments, FilterParams
from edl.util import add_universal_arguments, setup_logging


def main():
    # command line arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve')

    # default to non-overlapping=0
    add_hit_table_arguments(parser, flags='all', defaults={'nonoverlapping':0})
    parser.add_argument(
        "-o",
        "--outfilenome",
        dest="outfilename",
        default=None,
        metavar="OUTFILENAME",
        help="Write masked fasta output to OUTFILENAME.")
    parser.add_argument('hit_table', nargs='?',
                        type=argparse.FileType('rU'), default=sys.stdin,
                        help="Table of search results to be filtered. "
                             "If absent, data will be read from STDIN")

    add_universal_arguments(parser)

    arguments = parser.parse_args()

    setup_logging(arguments)

    # output file or STDOUT
    if arguments.outfilename is not None:
        logging.info("Writing data to %s" % (arguments.outfilename))
        outfile_handle = open(arguments.outfilename, 'w')
    else:
        logging.info("writing data to STDOUT")
        outfile_handle = sys.stdout

    # input file or STDIN (handled by argparse)
    infile_handle = arguments.hit_table
    logging.info("reading data from %s" % (infile_handle.name))

    # filter, but don't apply nonoverlapping yet
    # non-overlapping should be applied per-reference only
    params = FilterParams.create_from_arguments(arguments)
    # save user supplied value for later
    overlap_buffer = params.nonoverlapping
    # turn off for now
    params.set_nonoverlapping(-1)
    
    # merge
    hit_iter = filterM8Stream(infile_handle, params, return_lines=False)
    for query, query_hits in hit_iter:
        # group by reference hit
        hits_by_ref = defaultdict(list)
        for hit in query_hits:
            hits_by_ref[hit.hit].append(hit)
        
        # one output for query/reference pair
        for ref, ref_hits in hits_by_ref.items():
        
            # remove overlaps unless the buffer has been set to <0
            if overlap_buffer >= 0:
                ref_hits = remove_overlapping_hits(ref_hits, on_hit=True, buffer=params.nonoverlapping)
                ref_hits = remove_overlapping_hits(ref_hits, on_hit=False, buffer=params.nonoverlapping)
            
            # aggregate values
            length, score, identities = 0, 0, 0
            for hit in ref_hits:
                length += hit.mlen
                score += hit.score
                try:
                    # this will be off by 100x
                    identities += hit.pctid * hit.mlen
                except:
                    # just report pctid=0 if no pctid column in input
                    pass
            
            outfile_handle.write("%s\t%s\t%d\t%d\t%0.2f\n" % (query, ref, length, score, 
                                                            identities/length))
                    
    outfile_handle.close()
    infile_handle.close()

#############
# Functions #
#############

def remove_overlapping_hits(hits, on_hit=True, buffer=0):
    """ remove hits similar to the way the non-overlapping flag works but either on query 
    or on refernce """
    regions = []
    for hit in hits:
        start, end = (hit.hstart, hit.hend) \
                     if on_hit else \
                     sorted((hit.qstart, hit.qend))
        
        # if hit is shorter than buffer, then keep it
        if end - start <= buffer:
            regions.append((start, end))
            yield hit
            continue
        
        # adjust for buffer
        buf_start, buf_end = start + buffer, end - buffer

        # check overlap for all ranges
        for i, already_occ_range in enumerate(regions):
            if (buf_start >= already_occ_range[1] 
                or buf_end <= already_occ_range[0]):                    
                # does not overlap this range (try next range) 
                continue
            else:
                # overlaps a previous hit...move on
                break
        else:
            # there was no overlap, save range and yield this hit
            regions.append((start, end))
            yield hit



if __name__ == '__main__':
    main()
