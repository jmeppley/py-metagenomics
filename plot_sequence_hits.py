#! /usr/bin/env python
"""
Generates a pdf file of hit plots, one plot per query sequence. Each plot lays out the HSPs for the top N organisms.

usage: %prog
 -v, --verbose: extra output
 -t, --taxDump: directory with ncbi files
 -r, --rank: collaps counts on this rank
"""

import argparse
import logging
import numpy
import re
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages 
from Bio.Blast import NCBIXML
from edl import blastm8
from edl.util import add_universal_arguments, add_IO_arguments, \
                     inputIterator, setup_logging

def main():
    description = """
Generates a pdf file of hit plots, one plot per query sequence. Each plot lays out the HSPs for the top N organisms.
    """
    parser = argparse.OptionParser(description=description)
    add_universal_arguments(parser)
    add_IO_arguments(parser)
    add_hit_table_arguments(parser, flags=['format',
                                           'sort',
                                           'pctid',
                                           'length',],
                                    defaults={'pctid':75,
                                              'length':750},)

    parser.add_argument('-p','--pctid_cutoff',type=int, default=95,
            help="Only include alignments with at least on HSP with better than this (95%) pctid")
    parser.add_argument('-l','--length_cutoff',type=int, default=1000,
            help="Only include alignents with at least one HSP over this (1000bp) legnth")
    parser.add_argument('-N','--number_lines_per_plot',type=int, default=10,
            help="The number of lines (10) in each plot")
    parser.add_argument('-n','--number_plots_per_page', default=3, type=int, 
            help="The number of queries (4) plotted on each page")
    parser.add_argument('-e','--extension_cutoff',default=0, type=int,
            help="If this is set to a positive integer, hits that are within this distance of the end of the query will be annoated with the remaining length of the hit")
    parser.add_argument('-r','--reference_plot',default=False,action='store_true',
            help="Use reference sequences (aka targets) as the frame of reference for plots, displaying the locations of the query sequences below.")
    parser.add_argument('-C','--color_map_name',default='Reds',
            help="Name of the matplotlib colormap to use for HSP colors. Defaults to 'Reds'")
    parser.add_argument('-R','--label_pattern',default=None,
            help="Regexp pattern to pull label from hit descriptions or query IDs.")
    (arguments, args) = parser.parse_args()

    # check arguments
    setup_logging(arguments, description)

    for (inhandle, outhandle) in inputIterator(arguments):

        logging.info("Data in: %s\nImage out: %s",
                     inhandle,
                     outhandle)

        recordMap = parse_hits_for_plotting_from_arguments(inhandle, arguments)

        logging.debug("Pulling color map %s from matplotlib",
                      arguments.color_map_name)
        colorMap=plt.get_cmap(arguments.color_map_name)

        alignment_sort_key = sortArguments.get(arguments.sort_hits_by,None)
        if arguments.labelPattern is None:
            if arguments.referencePlot:
                labelPattern = r'^.+_([^_]+)\s*$'
            else:
                labelPattern = r'^\S+\s+(\S.*)$'
        else:
            labelPattern = arguments.labelPattern
        labelRE=re.compile(labelPattern)
        
        if arguments.referencePlot:
            drawHitsToReferenceToPDF(recordMap, 
                                     outputPdfFile,
                                     alignmentSortKey=alignment_sort_key,
                                     hspPct=.1, alnPct=.5, overlapPadPct=5.,
                                     color=getPctIDColorFunction(colorMap=colorMap,
                                            pctIDRange=[arguments.pctid_cutoff_hard,
                                                        100]),
                                     label=lambda r: "%s (%s)" % \
                                            (getLabelFromTitle(r.query,
                                                                 labelRE),
                                             getShortLengthString(r.query_length)),
                                     #label=lambda r: "%s (%s)" % (r.query.split("_")[-1],getShortLengthString(r.query_length)),
                                     topAlignmentColor='c',
                                     )
            pass
        else:
            drawQueryHits(recordMap,
                          outputPdfFile,
                          alignment_sort_key=alignment_sort_key,
                          N=arguments.number_orgs_per_plot,
                          label=lambda a: "%s (%s)" % \
                                  (getLabelFromTitle(a.hit, labelRE),
                                   getShortLengthString(a.length)),
                          color=getPctIDColorFunction(colorMap=colorMap,
                                        pctIDRange=[arguments.pctid_cutoff_hard,
                                                    100]),
                          maxDistanceToEnd=arguments.extensionCutoff,
                          queriesPerPage=arguments.number_plots_per_page,
                          )


def die( msg ):
    logging.error(msg)
    sys.exit()

def warn(msg):
    logging.warn(msg)

def parse_hits_for_plotting_from_arguments(inhandle, arguments):
    """
    Parse a hit_table into records
    """
    # build base filter from arguments
    params = blastm8.FilterParams.create_from_arguments(arguments)
    # keep all HSPs
    params.hspsPerHit=0
    record_map = parse_hits_for_plotting(inhandle,
                                    params,
                                    at_least_1_pctid=arguments.pctid_cutoff,
                                    at_least_1_length=arguments.length_cutoff)
    return record_map

def parse_hits_for_plotting(m8file, params, at_least_1_pctid, at_least_1_length):
    """
    Parse a hit table that has one line per hsp into a nested 
    data structure organized by qeury and hit ids.

    m8file: the hit table (as file name or file like object)
    params: a edl.blastm8.FilterParams object that will discard low quality HSPs
    at_least_1_pctid: all hits must have at least one HSP with %ID > this
    at_least_1_length: all hits must have at least one HSP with length > this
    """
    record_map = {}

    if not isinstance(m8file, blastm8.M8Stream):
        m8file = blastm8.M8Stream(m8file)

    query_count = 0
    kept_query_count = 0
    hit_count = 0
    kept_hit_count = 0
    hsp_count = 0
    good_hsp_count = 0
    min_hit_lengths = {}
    # blastm8 will filter out short, bad HSPs
    for query, hits in blastm8.filterM8Stream(m8file, params,
                                              returnLines=False):
        # we don't know the query length , so build up min est. as we go
        min_query_length = 0
        query_count += 1
        # each hit is an HSP, but not grouped by hit id
        hsps_by_hit = {}
        hits_with_good_hsps = set()
        for hit in hits:
            hsp_count += 1
            min_query_length = max(min_query_length, max(hit.qstart, hit.qend))
            min_hit_lengths[hit.hit] = max(min_hit_lengths.get(hit.hit, 0),
                                           max(hit.hstart, hit.hend))
            hsps_by_hit.setdefault(hit.hit, []).append(hit)
            if hit.pctid >= at_least_1_pctid and hit.mlen >= at_least_1_length:
                hits_with_good_hsps.add(hit.hit)
                good_hsp_count += 1
            
        hit_count += len(hsps_by_hit)
        if len(hits_with_good_hsps) == 0:
            # skip this query if no hits were good enough
            continue

        kept_query_count += 1
        kept_hit_count += len(hits_with_good_hsps)
        alignments = [Alignment(h, l, min_hit_lengths[h]) \
                      for h, l in hsps_by_hit.items() \
                      if h in hits_with_good_hsps]
        record_map[query] = Query(query, alignments, min_query_length)
            
    logging.info("Found %d records with %d hits out of %d records with %d hits",
                 kept_query_count, kept_hit_count, query_count, hit_count)
    logging.debug("Processed %d lines and found %d HSPs. %d of these were "
                  "also better than the higher filters",
                  m8file.lines, hsp_count, good_hsp_count)
    return record_map


class Query():
    def __init__(self, query_id, alignments, query_length=0):
        self.query = query_id
        self.alignments = alignments
        try:
            query_length = alignments[0].hsps[0].qlen
        except:
            if query_length == 0:
                for a in alignments:
                    for h in a.hsps:
                        query_length = max(query_length, max(h.qstart, h.qend))
        self.query_length = query_length

    def apply_cutoffs(self,
                      min_pctid=80,
                      min_length=100,
                      at_least_1_pctid=90,
                      at_least_1_length=1000,
                      inplace=False,
                      ):
        """ drop hsps under these thresholds 
        min_*: all hsps under this dropped
        at_least_1_*: all alignments must have at least 1 hsp better 
        inplace: don't create a new query object
        """
        new_alignments = []
        for alignment in self.alignments:
            new_hsps = [h for h in alignment.hsps \
                          if h.pctid >= min_pctid \
                          and h.mlen >= min_length]
            good_hsp_count = sum(1 for h in alignment.hsps \
                                 if h.pctid >= at_least_1_pctid \
                                 and h.mlen >= at_least_1_length) \
                             if at_least_1_pctid > min_pctid \
                             and at_least_1_length > min_length \
                             else len(alignment.hsps)
            if good_hsp_count > 0:
                if inplace:
                    alignment.hsps = new_hsps
                else:
                    alignment = Alignment(alignment.hit, 
                                          new_hsps,
                                          alignment.hit_length)
                new_alignments.append(alignment)
                
        if inplace:
            self.alignments = new_alignments
            return self
        else:
            return Query(self.query, new_alignments, self.query_length)

    def clone(self):
        return Query(self.query,
                     [a.clone() for a in self.alignemnts],
                     self.query_length)

class Alignment():
    def __init__(self, hit, hsps, hit_length=0):
        self.hit = hit
        self.hsps = hsps
        self.hit_length = hit_length

    def clone(self):
        return Alignment(self.hit,
                         list(self.hsps),
                         self.hit_length)

LONGEST_PCTID='pctid_of_longest'
def longestPctid(alignment):
    """
    Returns the pctid of the longest HSP
    """
    longestLength=0
    longestPctid=0
    for hsp in alignment.hsps:
        if hsp.mlen > longestLength:
            longestLength = hsp.mlen
            longestPctid = hsp.pctid
    return longestPctid

BEST_PCTID='best_pctid'
bestPctid = lambda a: max([h.pctid for h in a.hsps])

AVG_PCTID = 'total_pctid'
def combinedPctid(alignment):
    total_length = 0
    total_identities = 0
    for h in alignment.hsps:
        total_length += h.mlen
        total_identities += h.mlen * h.pctid / 100
    return total_identities / total_length

TOTAL_LENGTH='total_length'
totalLength = lambda a: sum([h.mlen for h in a.hsps])

LONGEST_HSP='max_length'
longestHsp = lambda a: max([h.mlen for h in a.hsps])

BEST_SCORE='best_score'
bestScore = lambda a: max([h.score for h in a.hsps])

TOTAL_SCORE='total_score'
totalScore = lambda a: sum([h.score for h in a.hsps])

sortArguments={LONGEST_PCTID:longestPctid,
             BEST_PCTID:bestPctid,
             AVG_PCTID:combinedPctid,
             TOTAL_LENGTH:totalLength,
             LONGEST_HSP:longestHsp,
             BEST_SCORE:bestScore,
             TOTAL_SCORE:totalScore}

def getLabelFromTitle(title, titleRE):
    m = titleRE.search(title)
    logging.debug(titleRE.pattern)
    if m:
        if len(m.groups())>=1:
            logging.debug("Label for %s is %s" % (title, m.group(1)))
            return m.group(1)
        else:
            logging.debug("Label for %s is %s" % (title, m.group()))
            return m.group()
    else:
        logging.debug("Label for %s is itself" % (title))
        return title

def drawQueryHits(queryHitMap,fileName,queriesPerPage=4,alignment_sort_key=None,figureSize=[11,8.5],subplotParams={'hspace':0.25},saveFigArgs={},**kwargs):
    """
    This splits up the plots into a multi-page PDF file.
    """
    
    numQueries=len(queryHitMap)
    queries=sorted(queryHitMap.keys(), 
                   key=lambda c: queryHitMap[c].query_length, 
                   reverse=True)
    firstQuery=0
    #if True:
    with PdfPages(fileName) as pdf:
        #if True:
        while firstQuery<numQueries:
            lastQuery=min(firstQuery+queriesPerPage,numQueries)
            fig,axes = plt.subplots(nrows=queriesPerPage,ncols=1,squeeze=False)
            fig.set_figwidth(figureSize[0])
            fig.set_figheight(figureSize[1])
            for query, ax in zip(queries[firstQuery:lastQuery],axes.flat):
                record = queryHitMap[query]
                plotWidth = record.query_length
                padding = plotWidth/25
                if alignment_sort_key is not None:
                    record.alignments.sort(reverse=True, \
                        key=alignment_sort_key)
                bars = plotTopHSPsWithExtensions(ax, 
                                                 record, 
                                                 extensionPadding=padding/5,
                                                 extensionLength=padding,
                                                 **kwargs)
                xlim=ax.set_xlim([0-padding, plotWidth+1+padding])
            fig.subplots_adjust(**subplotParams) 
            #plt.savefig(fileName)
            pdf.savefig(**saveFigArgs) 
            plt.close()
            firstQuery=lastQuery

def plotTopHSPsWithExtensions(ax, record, N=10, sort=None, label=None, 
                              color=None,
                              maxDistanceToEnd=0, extensionPadding=None, 
                              extensionLength=None, extensionColor='blue'):
    """
    Plots the hits for the first 10 orgnaisms in the given record's alignments

    Parameters:
    N: number of orgnaisms
    sort: a function operating on an alignment object that will return a key    for sorting
    label: a function operating on an alignment opject that returns an organism name for labelling and grouping
    color: either a single color for all HSPs or a function that takes an HSP as a single argument and returns a color.
    maxDistanceToEnd: how close does a hit have to be to the end to be considered for exentsions
    extensionPadding: how much space to leave between ends of plot and extension annotations (defaults to 1% of width)
    extensionLength: how long extension arrows should be (defaults to 5% of width)

    """

    # setup
    rows=[]
    starts=[]
    ends=[]
    if label is None:
        label = lambda a: "%s (%d bp)" % (a.hit, a.length)
    collectColors=not(color is None and isinstance(color,str) and isinstance(color,tuple))
    if collectColors:
        colorVals=[]
    labels={}
    if extensionPadding is None:
        extensionPadding=record.query_length/100
    if extensionLength is None:
        extensionLength=record.query_length/20
    extensions={}
    
    # translate alignments/hsps into coordinates
    for alignment in record.alignments:
        if len(alignment.hsps)==0:
            break
            
        alignmentLabel = label(alignment)
        if alignmentLabel in labels:
            row = labels[alignmentLabel]
        else:
            if len(labels)>=N:
                continue
            else:
                row=len(labels)
                labels[alignmentLabel]=row
                
        for hsp in sorted(alignment.hsps,key=lambda h: h.pctid):
            rows.append(row)
            (start,end) = [hsp.qstart,hsp.qend]
            starts.append(start)
            ends.append(end)
            if collectColors:
                colorVals.append(color(hsp))
                
            # does this hit extend off the end?
            if maxDistanceToEnd<=0:
                continue
            rowExtensions=extensions.setdefault(row,{})
            if start<=maxDistanceToEnd:
                if hsp.hstart>hsp.hend:
                    overhangLength = alignment.length - hsp.hstart
                else:
                    overhangLength = hsp.hstart
                rowExtensions['left']=max(rowExtensions.get('left',0),overhangLength)
                #print "Overhang left: %d(%d) = (%d,%d) of %d to (%d,%d) of %d" % (overhangLength,row,hsp.qstart,hsp.qend,record.query_length,hsp.hstart,hsp.hend,alignment.length)
            elif record.query_length - end < maxDistanceToEnd:
                if hsp.hstart>hsp.hend:
                    overhangLength = hsp.hend
                else:
                    overhangLength = alignment.length - hsp.hend
                rowExtensions['right']=max(rowExtensions.get('right',0),overhangLength)                
                #print "Overhang right: %d(%d) = (%d,%d) of %d to (%d,%d) of %d" % (overhangLength,row,hsp.qstart,hsp.qend,record.query_length,hsp.hstart,hsp.hend,alignment.length)
            
    # generate plot
    rows=numpy.array(rows)
    starts=numpy.array(starts)
    ends=numpy.array(ends)
    if collectColors:
        color=colorVals
    bars=ax.barh(0-rows, ends-starts, left=starts, height=.3, color=color)

    # write labels in the plot, just under each set of HSP bars
    for label, row in labels.items():
        ax.text(int(record.query_length/2.0),-.5-row, label, horizontalalignment='center')
    
    # Draw extensions
    for row, extensionDict in extensions.items():
        if 'left' in extensionDict:
            eString = getShortLengthString(extensionDict['left'])
            ax.arrow(0-extensionPadding, .1-row, 0-extensionLength-extensionPadding, 0, head_width=0.3, head_length=extensionPadding, clip_on=False, fc=extensionColor, ec=extensionColor)
            ax.text(0-(extensionLength/2), -.5-row, eString, horizontalalignment='right', clip_on=False, color=extensionColor)
        if 'right' in extensionDict:
            eString = getShortLengthString(extensionDict['right'])
            ax.arrow(record.query_length+extensionPadding,.1-row, extensionLength-extensionPadding, 0, head_width=0.3, head_length=extensionPadding, clip_on=False, fc=extensionColor, ec=extensionColor)
            ax.text(record.query_length+extensionPadding, -.5-row, eString, horizontalalignment='left', clip_on=False, color=extensionColor)
    
    # tidy and decorate plot
    ax.set_yticks([])
    ax.set_ylim(-N,1)
    ax.set_title(record.query)

    return bars

steps=['','K','M','G','T','P']
def getShortLengthString(length, unit="bp"):
    """
    Given number of kilobytes, return string like 34K or .4M or 4T.
    Will always be 3 characters lone
    """

    # first, figure out ourder of magnitude
    size=float(length)
    step=0
    while size>99.5:
        size/=1024.
        step+=1

    # next, format
    if size>=.95:
        size='%.0f' % size
    else:
        size='%.1f' % size
        size=size[-2:]

    # add suffix
    size+=" " + steps[step] + unit
    return size

def plotTopHSPs(ax, record, N=10, sort=None, label=None, colorMap=plt.get_cmap('Blues'), colorKey=None, colorRange=None):
    """
    Plots the hits for the first 10 orgnaisms in the given record's alignments

    Parameters:
    N: number of orgnaisms
    sort: a function operating on an alignment object that will return a key for sorting
    label: a function operating on an alignment opject that returns an organism name for labelling and grouping
    colorMap: a matplotlib color map for coloring bars
    colorKey: a function operating on a HSP object that returns an scalar value.
    colorRange: the range of values for color mapping
    """
    rows=[]
    starts=[]
    ends=[]
    if label is None:
        label = lambda a: a.hit
    if colorKey is not None:
        colorVals=[]
    labels={}
    
    iterator = record.alignments if sort is None \
                     else sorted(record.alignments, key=sort, reverse=True)
    for alignment in iterator:
        if len(alignment.hsps)==0:
            break
            
        alignmentLabel = label(alignment)
        if alignmentLabel in labels:
            row = labels[alignmentLabel]
        else:
            if len(labels)>=N:
                continue
            else:
                row=len(labels)
                labels[alignmentLabel]=row
                
        for hsp in alignment.hsps:
            rows.append(row)
            starts.append(hsp.qstart)
            ends.append(hsp.qend)
            if colorKey is not None:
                colorVals.append(colorKey(hsp))

    if colorKey is None:
        color=None
    else:
        cVals = numpy.array(colorVals)
        if colorRange is None:
            colorRange=(min(cVals),max(cVals))
        scale = 256/float(colorRange[1]-colorRange[0])
        color=[colorMap(int(scale*(val-colorRange[0]))) for val in cVals]

    rows=numpy.array(rows)
    starts=numpy.array(starts)
    ends=numpy.array(ends)
    bars=ax.barh(0-rows, ends-starts, left=starts, color=color)

    ax.set_yticks(numpy.array(range(0,0-N,-1))+.4)
    ax.set_ylim(-N,1)
    
    labelList = sorted(labels.keys(), key=lambda l: labels[l])
    while (len(labelList)<N):
        labelList.append("")
    ax.set_yticklabels(labelList)
    ax.set_title(record.query)

    return bars
###############
# METHODS for plotting hits against Reference Sequences
def getPctIDColorFunction(colorMap=plt.get_cmap('Reds'),
                          pctIDRange=[75,100],
                          colorMapMax=256,
                          ):
    scale = colorMapMax/float(pctIDRange[1]-pctIDRange[0])
    return lambda hsp: colorMap(int(scale*(pctIDFunction(hsp)-pctIDRange[0])))

def drawHitsToReferenceToPDF(recordMap, fileName,
                              refSeqsPerPage=3,
                              rowsPerPlot=8,
                              hit_sort_key=None,
                              figureSize=[11,8.5],
                              subplotParams={'hspace':0.3},
                              saveFigArgs={},
                              alignmentSortKey=None,
                              topAlignmentColor='k',
                              **kwargs):
    """
    This splits up the plots into a multi-page PDF file.
    """
    alignmentsByRef={}
    topAlignmentsByRef={}
    for contig,record in recordMap.items():
        #filterAlignments(record, keepAllHSPs=False, pctid=75)
        if len(record.alignments)==0:
            continue
        if alignmentSortKey is None:
            alignments=record.alignments
        else:
            alignments=sorted(record.alignments,key=alignmentSortKey,reverse=True)
        topAlignmentsByRef.setdefault(alignments[0].hit_id,[]).append((alignments[0],record))
        for alignment in alignments:
            alignmentsByRef.setdefault(alignment.hit_id,[]).append((alignment,record))
    rankedRefIds = sorted(topAlignmentsByRef.keys(), key=lambda i: len(alignmentsByRef[i]), reverse=True)
    
    refNum=0
    pageCount=0
    with PdfPages(fileName) as pdf:
        while refNum<len(rankedRefIds):
            fig,axes = plt.subplots(nrows=refSeqsPerPage,ncols=1,squeeze=False)
            fig.set_figwidth(figureSize[0])
            fig.set_figheight(figureSize[1])

            for i,ax in enumerate(axes.flat):
                while True:
                    if refNum+i>=len(rankedRefIds):
                        break
                    refId = rankedRefIds[refNum+i]
                    refDesc=alignmentsByRef[refId][0][0].hit_def
                    ax.set_title(refDesc)
                    topAlignments = set([a for a,r in topAlignmentsByRef[refId]])

                    bars=plotReferenceHSPs(ax,alignmentsByRef[refId],
                                           N=rowsPerPlot,
                                           sort=alignmentSortKey,
                                           lineColor=lambda a: topAlignmentColor if a in topAlignments else 'k',
                                           **kwargs)
                    if bars is not None:
                        logging.debug("Plot %s at %d" % (refId, i))
                        break
                    else:
                        logging.debug( "No plot for %s at %d" % ( refId, i))
                        refNum+=1

            fig.subplots_adjust(**subplotParams) 
            pdf.savefig(**saveFigArgs)  # saves the current figure into a pdf page
            plt.close()
           
            refNum+=refSeqsPerPage
            pageCount+=1
            logging.debug("Reference number:" + str( refNum))
    return pageCount

def plotReferenceHSPs(ax, alignmentRecordPairs, hspPct=3., alnPct=5., overlapPadPct=5.,
                      N=8, sort=None, label=None, 
                      color="b", lineColor="k"):
    """
    Reverse hit plot. On the given axis, given the alignments for which this subject was a hitPlots the hits for the first 10 orgnaisms in the given record's alignments

    Parameters:
    N: number of rows in plot
    hspPct: minimum percent of the full span that an HSP must cover to be plotted.
    sort: a function operating on an alignment object that will return a key for sorting
    label: a function operating on an record object that returns a string for labelling hits
    color: a function operating on a HSP object that returns a color OR a single color for all HSPs.
    lineColor: a function operating on Alignments that returns a color OR a single color for all lines.

    """

    # setup
    minHSPlength = hspPct * alignmentRecordPairs[0][0].length / 100.
    minAlnLength = alnPct * alignmentRecordPairs[0][0].length / 100.
    overlapPadding = overlapPadPct * alignmentRecordPairs[0][0].length / 100.
    #print "Cutoffs: HSP: %d Aln: %d" % (minHSPlength, minAlnLength)
    alignments=[]
    alignmentRecordMap={}
    for (a,r) in alignmentRecordPairs:
        alignments.append(a)
        alignmentRecordMap[a]=r
    if sort is not None:
        alignments.sort(key=sort, reverse=True)
    if label is None:
        label = lambda record: record.query
    generateLineColors=not (isinstance(lineColor,str) or isinstance(lineColor,list) or isinstance(lineColor,tuple) or lineColor is None)
    alignmentRowMap=layoutAlignments(alignments, maxDepth=N, minHSPlength=minHSPlength, minAlnLength=minAlnLength, overlapPadding=overlapPadding)
    
    if len(alignmentRowMap)==0:
        return None
    
    # plot line for full extent of alignment
    for alignment, row in alignmentRowMap.items():
        (start,end)=_getRange(alignment, minHSPlength=minHSPlength, minAlnLength=minAlnLength)
        if generateLineColors:
            alignmentLineColor=lineColor(alignment)
        else:
            alignmentLineColor=lineColor
        line = ax.plot((start,end),(.15-row,.15-row),color=alignmentLineColor,linestyle='solid')
        ax.text(int(start+(end-start)/2.0),-.5-row, label(alignmentRecordMap[alignment]), horizontalalignment='center')

    # plot HSPs
    rows=[]
    starts=[]
    ends=[]
    generateColors=not (isinstance(color,str) or isinstance(color,list) or isinstance(color,tuple) or color is None)
    if generateColors:
        colors=[]
    else:
        colors=color
    if generateLineColors:
        edgeColors = []
    else:
        edgeColors=lineColor   
    
    # translate alignments/hsps into coordinates
    for alignment, record in alignmentRecordPairs:
        if alignment not in alignmentRowMap:
            continue
           
        row = alignmentRowMap[alignment]
        for hsp in sorted(alignment.hsps,key=lambda hsp: 100.0*hsp.identities/float(hsp.mlen)):
            (start,end) = (numpy.min((hsp.hstart,hsp.hend)),
                           numpy.max((hsp.hstart,hsp.hend)))
            if end-start+1<minHSPlength:
                continue
            rows.append(row)
            starts.append(start)
            ends.append(end)
            if generateColors:
                colors.append(color(hsp))
            if generateLineColors:
                edgeColors.append(lineColor(alignment))
                
    # generate plot
    rows=numpy.array(rows)
    starts=numpy.array(starts)
    ends=numpy.array(ends)
    bars=ax.barh(0-rows, ends-starts, left=starts, height=.3, color=colors, edgecolor=edgeColors)

    # tidy and decorate plot
    ax.set_ylim(-N,1)

    return bars

# <codecell>

def layoutAlignments(alignments, maxDepth=10, truncate=False, minHSPlength=0, minAlnLength=0):
    """
    Return a map from aligment to row that puts one alignment per row
    """
    row=0
    rowMap={}
    for alignment in alignments:
        if len(alignment.hsps)==0:
            continue
        alignmentRange=_getRange(alignment, minHSPlength=minHSPlength, minAlnLength=minAlnLength)
        if alignmentRange is None:
            continue
        rowMap[alignment]=row
        row+=1
        if row==maxDepth:
            break

    return rowMap

# <codecell>

def layoutAlignments(alignments, maxDepth=10, truncate=False, minHSPlength=0, minAlnLength=0, overlapPadding=0):
    """
    Return a map from aligment to row that will maximize the number of alignments that fit in the given maxDepth
    """
    rows={}
    rangesByRow={}
    for alignment in alignments:
        if len(alignment.hsps)==0:
            continue
        alignmentRange=_getRange(alignment, minHSPlength=minHSPlength, minAlnLength=minAlnLength)
        if alignmentRange is None:
            continue
        for row in range(maxDepth):
            for hitRange in rangesByRow.get(row,[]):
                if _intersects(alignmentRange,hitRange,overlapPadding=overlapPadding):
                    # skip to next row
                    break
            else:
                # no overlaps, exit loop
                hitRow=row
                break
        else:
            # overlap in every row
            if truncate:
                break
            else:
                continue
        rows[alignment]=hitRow
        rangesByRow.setdefault(hitRow,[]).append(alignmentRange)
    
    return rows

def _intersects(range0, range1,overlapPadding=0):
    return (range0[0]<=(range1[1]+overlapPadding)) and (range0[1]>=(range1[0]-overlapPadding))

def _getRange(alignment, minHSPlength=0, minAlnLength=0):
    """
    Return the extent of the alignment
    """
    
    # get list of end points for all HSPs over cutoff
    positions=[]
    for hsp in alignment.hsps:
        alnLength=numpy.absolute(hsp.hend-hsp.hstart)+1
        if alnLength<minHSPlength:
            continue
        positions.append(hsp.hend)
        positions.append(hsp.hstart)
    
    # Did we find any HSPs over the minimum length?
    if len(positions)==0:
        return None
    
    # Does the whole alignment surpass the minimum length?
    (start,end)=(min(positions),max(positions))
    if 1+numpy.absolute(end-start)<minAlnLength:
        return None

    # return the range
    return (start,end)

if __name__ == '__main__':
    main()
