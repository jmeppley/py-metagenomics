#! /usr/bin/env python
"""
Generates a simple bar plot from a table of data. By default first two columns are taken as labels and data. Can alternatively transpose and take selected row/column.

usage: %prog
 -v, --verbose: extra output
 -t, --taxDump: directory with ncbi files
 -r, --rank: collaps counts on this rank
"""

from optparse import OptionParser
import sys, re, logging
import matplotlib
import numpy as np
import pandas                         # gives us an R-like data frame
from edl.util import addUniversalOptions, setupLogging

def main():
    usage = "usage: %prog OPTIONS"
    description = """
Generates a simple bar plot from a table of data. By default first two columns are taken as labels and data. Can alternatively transpose and take selected row/column.
    """
    parser = OptionParser(usage, description=description)
    addUniversalOptions(parser)
    parser.add_option("-i", "--inputfile", dest="infile",
                      metavar="INFILE", help="Read data table from INFILE"),
    parser.add_option("-o", "--others", dest="others", default=False, action="store_true", help="include sum of values under cutoff")
    parser.add_option("-p", "--pct", dest="pct", default=False, action="store_true", help="plot percents (before cutoff) instead of counts")
    parser.add_option("-t", "--transpose",
                  action="store_true", dest="transpose", default=False,
                  help="Transpose data")
    parser.add_option("-c", "--cutoff", dest="cutoff", type="float", default=0.025,
		  help="data cutoff value (default .025)",
                  metavar="CUTOFF")
    parser.add_option("-d", "--dataColumn", dest="dataCol", default='0',
            help="Index (starting at 0) or name of column with data to be plotted",
                  metavar="HITCOL")
    parser.add_option("-f", "--format", dest="format", default='pdf', choices=['png','ps','pdf','svg'],
		  help="Format for output image", metavar="FORMAT")
    parser.add_option("-O","--outputFile", default=None, help="Manually set output file name")
    parser.add_option("-N","--outputNameOnly", default=False, action='store_true', help="Only print name of file to create. Don't do anything else")

    (options, args) = parser.parse_args()

    # check arguments
    setupLogging(options, description)

    if options.infile is None:
        options.error("Please supply a blast file name!")

    if options.outputFile is None:
        ofbits = [options.infile,'1dBarPlot']
        if options.transpose:
            ofbits.append('T')
        if options.pct:
            ofbits.append('P')
        if options.others:
            ofbits.append('O')
        ofbits.extend([options.dataCol,str(options.cutoff),options.format])
        outfile = '.'.join(ofbits)
    else:
        outfile=options.outputFile

    backend = options.format
    if backend=='png':
        backend='agg'
    matplotlib.use(backend)
    import matplotlib.pyplot as plt

    print outfile
    if options.outputNameOnly:
        sys.exit(0)

    log("Data in: %s\nImage out: %s" % (options.infile,outfile))
    log("Output format: %s" % options.format)

    series = makeSeriesFromFile(options.infile,options.dataCol, options.transpose, options.cutoff, options.others, options.pct)

    # create list of values over cutoff and sum of others
    labels = series.index

    # plot
    fig = plt.figure()
    ax =fig.add_axes([.1,.4,.8,.5])
    x = range(len(series))
    ax.bar(x,series)
    plt.xticks([i+.5 for i in x],labels,size=7,rotation=-90)

    plt.savefig(outfile,format=options.format)

#############
# Functions #
#############
def log(msg):
    logging.info(msg)

def die( msg ):
    logging.error(msg)
    sys.exit()

def warn(msg):
    logging.warn(msg)

def makeSeriesFromFile(infile, dataCol, transpose, cutoff, others, pct):
    log("loading %s cutoff:%s" % (infile,str(cutoff)))
    df = pandas.io.parsers.read_table(infile,index_col=0)
    if transpose:
        df = df.T

    col = dataCol
    try:
        col = int(col)
        col = df.columns[col]
    except ValueError:
        # column is not an int!
        pass
    log("using colum: %s" % col)
    series = df[col]

    series = cleanSeries(series,['Total','None'])
    series = findRowsOverCutoff(series, cutoff, others, pct)

    return series

def loadMultipleFiles(fileList, dataCol, transpose, cutoff, others, pct):
    data = []
    for f in fileList:
        data.append(makeSeriesFromFile(f,dataCol, transpose, cutoff, others, pct))
    return data

def printMultipleFilesT(fileList, dataCol, transpose, cutoff, others, pct, output):
    data = loadMultipleFiles(fileList, dataCol, transpose, cutoff, others, pct)
    names=[]
    for d in data:
        names.extend(d.index)
    names=sorted(set(names))

    output.write("Counts\t%s\n" % ("\t".join(fileList)))
    for n in names:
        output.write(n)
        for d in data:
            output.write("\t")
            if n in d:
                output.write(str(d[n]))
            else:
                output.write('0')
        output.write("\n")

def printMultipleFiles(fileList, dataCol, transpose, cutoff, others, pct, output):
    data = loadMultipleFiles(fileList, dataCol, transpose, cutoff, others, pct)
    names=[]
    for d in data:
        names.extend(d.index)
    names=sorted(set(names))

    output.write("Counts\t%s\n" % ("\t".join(names)))
    for i in range(len(data)):
        d=data[i]
        f=fileList[i]
        output.write(f)
        for n in names:
            output.write("\t")
            if n in d:
                output.write(str(d[n]))
            else:
                output.write('0')
        output.write("\n")

def cleanSeries(series, removals):
    """
    Given a data series and list of row and column names to ignore,
    returns list of indices to use.
    Also, looks for rows columns to combine if names differ by only whitespace and capitalization
    """
    N=len(series)
    slice = range(N)
    for key in removals:
        try:
            slice.remove(series.index.get_loc(key))
        except KeyError:
            pass

    nameMap={}
    keys = series.keys()
    for name in keys:
        baseName = name.strip().lower()
        if baseName in nameMap:
            # add values to first row with this name
            trueName = nameMap[baseName]
            series[trueName] += series[name]
            # remove this version of the name
            slice.remove(series.index.get_loc(name))
        else:
            # store this as the 'true' name
            nameMap[baseName] = name

    return series[slice]

def findRowsOverCutoff(series, cutoff, others, retpct):
    """
    REturn a modified series with only values over given cutoff
    """
    leftoverSum = 0
    slice = []
    total = series.sum()

    for i in xrange(len(series)):
        pct = float(series[i])/total
        if pct>=cutoff:
            slice.append(i)
        else:
            leftoverSum+=series[i]

    newSeries = series[slice]
    if others:
        newSeries = newSeries.append(pandas.Series([leftoverSum],['Other']))
    if retpct:
        newSeries = newSeries/float(series.sum())

    return newSeries

if __name__ == '__main__':
    main()
