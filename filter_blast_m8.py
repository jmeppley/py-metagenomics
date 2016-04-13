#! /usr/bin/python
"""
"""

from optparse import OptionParser
import sys, logging
from edl.blastm8 import filterM8, addHitTableOptions, FilterParams
from edl.util import addUniversalOptions, setupLogging

def main():
    usage = "usage: %prog [OPTIONS] BLAST_FILE"
    description = """
    Take a blast result table and output a subset of hits based on the chosen filtering options. If more than one blast file given, use -O to get multiple output files, otherwise all output data will be concatenated into one output.
    """

# command line options
    parser = OptionParser(usage, description=description, conflict_handler='resolve')
    addHitTableOptions(parser, flags='all')
    parser.add_option("-o", "--outfilenome", dest="outfilename", default=None,
                      metavar="OUTFILENAME", help="Write masked fasta output to OUTFILENAME.")
    parser.add_option('-O', '--autoOutName', default=False,
                      action='store_true',
                      help="Automatically generate output file name from input name and options. Overridden by -o, cannot be used with data from STDIN.")

    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options,description)

    #if options.hitTableFormat=='last':
    #    if options.hitTableSort=='evalue':
    #        parser.error("The last format has no evalue to sort by, sorry")

    # check that we have blast file as argument
    if len(args) <= 1:
        # input
        if len(args) == 1:
            infile = args[0]
            logging.info("reading data from %s" % (infile))
            instream = open(infile,'rU')
        else:
            infile = './stdin'
            logging.info("reading data from STDIN")
            instream=sys.stdin

        # output
        if options.outfilename is not None:
            logging.info("Writing data to %s" % (options.outfilename))
            outstream=open(options.outfilename,'w')
        elif options.autoOutName:
            outfile=getOutputFile(infile,options)
            logging.info("Writing data to %s" % (outfile))
            outstream=open(outfile,'w')
        else:
            logging.info("writing data to STDOUT")
            outstream=sys.stdout

        # filter
        params=FilterParams.createFromOptions(options)
        filterM8(instream,outstream,params)
    else:
        if not options.autoOutName:
            if options.outfilename is not None:
                logging.info("Writing data to %s" % (options.outfilename))
                outstream=open(options.outfilename,'w')
            else:
                logging.info("writing data to STDOUT")
                outstream=sys.stdout
        for infilename in args:
            logging.info("reading data from %s" % (infilename))
            instream=open(infilename,'rU')
            if options.autoOutName:
                outstream=open(getOutputFile(infilename,options),'w')

            # filter
            params=FilterParams.createFromOptions(options)
            filterM8(instream,outstream,params)

            if options.autoOutName:
                outstream.close()
            instream.close()

#############
# Functions #
#############
def getOutputFile(infile, options):
    """
    Use the requested options to name the output file
    """
    outfile = infile
    if options.filterPctid > 0:
        outfile += ".i%g" % options.filterPctid
    if options.filterLength > 0:
        outfile += ".l%d" % options.filterLength
    if options.filterBits > 0:
        outfile += ".b%g" % options.filterBits
    if options.filterEvalue is not None:
        outfile += ".e%g" % options.filterEvalue
    if options.filterAln is not None and options.filterAln>0:
        outfile += ".a%g" % options.filterAln
    if options.filterHspsPerHit!=1:
        outfile += ".h%d" % options.filterHspsPerHit
    if options.filterTopPct >= 0:
        outfile += ".p%g" % options.filterTopPct
    if options.filterNonoverlapping:
        outfile += ".u"
    if options.filterHitsPerRead > 0:
        outfile += ".n%d" % options.filterHitsPerRead

    if outfile == infile:
        sys.exit("outfile and infile are the same!!\n%s" % infile)

    return outfile

if __name__ == '__main__':
    main()
