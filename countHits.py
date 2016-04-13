#! /usr/bin/python
"""
Count hits in a tabular blast output. By default, first hit for each read is used.

usage: %prog [options]
  -h, --help            show this help message and exit
  -i FILE, --infile=FILE
                        Read raw table from INFILE
  -o OUTFILE, --outfile=OUTFILE
                        Write collapsed table to OUTFILE
  -d DELIM, --delim=DELIM
                        Input table
  -D DELIM, --delimOut=DELIM
                        Output table delimiter
  -F, --countFirst      Don't skip the first line, it's NOT a header
  -R READCOL, --readColumn=READCOL
                        Index (starting at 0) of column with read name, 0 is
                        default
  -H HITCOL, --hitColumn=HITCOL
                        Index (starting at 0) of column with hit name (for
                        counting), 2 is default, if less than zero, all (non-
                        read) columns will be used as multiple hits
  -s HITSEP, --hitSep=HITSEP
                        Use this string to split multiple values in single hit
                        cell. Default is 'None' to leave hits as is, use
                        'eval' to parse as python repr strings
  -C CLUSTERFILE, --clusterFile=CLUSTERFILE
                        Location of cd-hit output. Use this to multiply hits
                        to recreate counts of unclustered data.
  -c CUTOFF, --cutoff=CUTOFF
                        Cutoff for showing taxa. If a fractional count for a
                        taxa is below this value, it will be folded up into
                        its parent domain. Defaults to: 0
  -a ALLMETHOD, --allMethod=ALLMETHOD
                        'first' means +1 for every hit found for each read. 'all' means 1 to the first hit for each read. 'portion' means 1/(nhits) for all hits of each read. Defaults to 'all'
"""

from optparse import OptionParser
import sys, re, logging
from edl.hits import addCountOptions, getAllMethod
from edl.util import addUniversalOptions, setupLogging

# a habit that stuck in Perl
die=sys.exit

def main():
    ## set up CLI
    usage = "usage: %prog [options]"
    description = """
    Count hits in a table with read and hit names.
    """

    parser = OptionParser(usage, description=description)
    parser.add_option("-i", "--infile", dest="infile",
                      metavar="FILE", help="Read raw table from INFILE")
    parser.add_option("-o", "--outfile", dest="outfile",
                      metavar="OUTFILE", help="Write collapsed table to OUTFILE")
    parser.add_option("-d", "--delim", dest="delim", default="\t",
                      help="Input table delimiter", metavar="DELIM")
    parser.add_option("-D", "--delimOut", dest="delimOut", default="\t",
                      help="Output table delimiter", metavar="DELIM")
    parser.add_option('-F', '--countFirst', action='store_true', default=False,
                       help="Don't skip the first line, it's NOT a header")
    parser.add_option("-R", "--readColumn", dest="readCol", type="int", default=0,
                      help="Index (starting at 0) of column with read name, 0 is default",
                      metavar="READCOL")
    parser.add_option("-H", "--hitColumn", dest="hitCol", type="int", default=2,
                      help="Index (starting at 0) of column with hit name (for counting), 2 is default, if less than zero, all (non-read) columns will be used as multiple hits",
                      metavar="HITCOL")
    parser.add_option('-s', '--hitSep', default=None,
                      help="Use this string to split multiple values in single hit cell. Default is 'None' to leave hits as is, use 'eval' to parse as python repr strings")
    addWeightOption(parser, multiple=False)
    parser.add_option("-T", "--total", default=False, action="store_true",
                      help="Report 'Total' in the first row")

    # cutoff options
    addCountOptions(parser,{'cutoff':0})

    # logging and help
    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options, description)

    # make sure we have something to do
    if (options.infile==None):
        logging.info("Reading table from: STDIN")
    else:
        logging.info ("Reading table from: " + options.infile )

    if (options.outfile==None):
        logging.info("Writing counts to: STDOUT")
    else:
        logging.info ("Writing counts to: " + options.outfile )

    # process arguments
    takeFirst = (options.allMethod == 'first')
    splitHits = (options.hitSep is not None and options.hitSep != 'None')
    uncluster = (options.weights != None)

    if options.hitSep=='eval':
        parser.error("Sorry, parsing with eval is not yet supported!")

    ## inform the curious user
    logging.info ("Delimiter: '" + options.delim )
    logging.info ("Read names in col: '" + str(options.readCol) )
    logging.info ("Hit names in col: '" + str(options.hitCol) )
    if splitHits:
        logging.info("Splitting hits with: %s" % (options.hitSep))
        logging.warn("Splitting hits has not been tested yet! Let me know how it goes.")
    if takeFirst:
        logging.info("Taking first hit for each read.");
    else:
        if options.allMethod == 'portion':
            logging.info ("Dividing count among all hits for each read.")
        else:
            logging.info ("Adding 1 to every hit for each read")
    if uncluster:
        logging.info("Getting read cluster sizes from: %s" % (options.weights));
    if options.countFirst:
        logging.info("First line is data")
    else:
        logging.info("Skipping first line")

    # Do the counting!
    counts = {}
    countHitsForRead=getAllMethod(options.allMethod)

    clusteredReadCounts={}
    if uncluster:
        clusteredReadCounts = parseMapFile(options.clusterFile, valueType=int)

    currentRead=''
    readCount=1
    hits=[]

    if options.infile is None:
        infile = sys.stdin
    else:
        infile = open(options.infile)

    # loop over lines
    if not options.countFirst:
        # skip first line
        try:
            infile.next()
        except StopIteration:
            raise Exception("No lines in %s" % str(infile))

    for line in infile:
        line=line.rstrip('\r\n')
        rowcells = line.split(options.delim)
        # get read
        read = rowcells[options.readCol]

        # if it's a new read, process previous read
        if currentRead=='':
            currentRead=read
        elif read != currentRead and currentRead != '':
            readCount+=1
            logging.info( "Checking hits for %s" % currentRead)

            # was it part of a cluster?
            multiplier = 1
            if uncluster:
                multiplier = clusteredReadCounts[currentRead]

            # where does the count for this read go
            countHitsForRead(hits, counts, multiplier=multiplier)

            hits=[]
            currentRead=read

        # get hit from this line
        if options.hitCol>=0:
            hit=rowcells[options.hitCol]
            if splitHits:
                hits.extend(hit.split(options.hitSep))
            else:
                hits.append(hit)
        else:
            rowcells.pop(options.readCol)
            hits.extend(rowcells)

    # check last read!
    logging.info( "Checking hits for %s" % currentRead)
    # was it part of a cluster?
    multiplier = 1
    if uncluster:
        multiplier = clusteredReadCounts[currentRead]
    # where does the count for this read go
    countHitsForRead(hits,counts,multiplier=multiplier)

    # apply cutoff
    if options.cutoff>0:
        applyFractionalCutoff(counts, threshold=options.cutoff*readCount)

    # print output
    if options.outfile is None:
        outhandle = sys.stdout
    else:
        outhandle = open(options.outfile,'w')

    if options.total:
        outhandle.write("Total%s%d\n" % (options.delimOut, readCount))

    if options.allMethod=='portion':
        outFmtString = "%s%s%f\n"
    else:
        outFmtString = "%s%s%d\n"

    delimRE = re.compile(options.delimOut)
    for hit, count in counts.iteritems():
        hit=delimRE.sub('_',hit)
        outhandle.write(outFmtString % (hit,options.delimOut,count))

if __name__ == '__main__':
    main()
