#!/usr/bin/python
"""
Uses biopython and matplotlib to create a histogram of read lengths
"""

import matplotlib, logging
from Bio import SeqIO


def main():
    from optparse import OptionParser

    ## set up CLI
    description = __doc__
    parser = OptionParser(description=description)
    parser.disable_interspersed_args()
    parser.add_option("-o", "--outfile", help="Where to write fugure")
    parser.add_option("-f", "--format", choices=['pdf','png'], help="Figure format")
    parser.add_option("-r", "--reads", help="directory to put output and support files into")
    parser.add_option("-b", "--bins", default=50, type="int", help="Number of histogram bins")

    (options, args) = parser.parse_args()

    # configure matplotlib
    backend = options.format
    if backend=='png':
        backend='agg'
    matplotlib.use(backend)
    import matplotlib.pyplot as plt

    # get list of read lengths
    readLengths=[]
    readsFormat = fastaOrFastq(options.reads)
    for record in SeqIO.parse(options.reads, readsFormat):
        readLengths.append(len(record))

    # plot
    plt.hist(readLengths,bins=options.bins)
    plt.savefig(options.outfile,format=options.format)


def fastaOrFastq(reads):
    with open(reads) as f:
        if f.next()[0]=='>':
            return 'fasta'
        else:
            return 'fastq'

if __name__ == "__main__":
    main()
