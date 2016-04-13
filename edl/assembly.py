from numpy import histogram
import os, sys, numpy, argparse, logging, re, pandas
import edl.blastm8
from Bio import SeqIO
from edl.util import asciiHistogram
logger=logging.getLogger(__name__)

## 
# the calc_stat, plot_assembly, and getN50 methods originally come from the
# assemlbly_quality_stats.py script by:
# Author: Travis Poulsen
# Date: 09 Feb. 2013
# http://travispoulsen.com/blog/2013/07/basic-assembly-statistics/
# https://gist.github.com/tpoulsen/422b1a19cbd8c0f514fe/raw/assembly_quality_stats.py
def calc_stats(file_in, return_data=False, txt_width=0, log=False, backend=None, format='fasta', minLength=0, **kwargs):
    """
    Given contigs in fastsa format:
     * calculate length stats (including N50)
     * plot histogram (use txt_width and backend to select format)
    """
    with open(file_in, 'r') as seq:
        sizes = [len(record) for record in SeqIO.parse(seq, format) if len(record)>=minLength]

    sizes = numpy.array(sizes)
    data = get_contig_length_stats(sizes)

    if not return_data:
        report = get_contig_length_report(data)

    if backend is not None:
        h = plot_assembly(sizes, file_in, data, log=log, backend=backend, **kwargs)
    if txt_width>0:
        if backend is None:
            h=histogram(sizes, **kwargs)
        histogramText = asciiHistogram(h, log=log, width=txt_width)
        if not return_data:
            if log:
                report += "\n\nContig length histogram (log):\n"
            else:
                report += "\n\nContig length histogram:\n"
            report += histogramText
        else:
            data['histogram']=histogramText

    if return_data:
        return data
    else:
        return report

def get_contig_length_stats(sizes):
    """
    return a dict of useful contig length stats
    """
    return {'min':numpy.min(sizes),
            'max':numpy.max(sizes),
            'mean':numpy.mean(sizes),
            'median':numpy.median(sizes),
            'N50':int(getN50(sizes)),
            'N75':int(getN50(sizes,N=75)),
            'N90':int(getN50(sizes,N=90)),
            'count':len(sizes),
            }

def get_contig_length_report(data):
    """
    return a formatted string summarizing contig length data
    """
    report ='Number of contigs:\t%i' % data['count']
    report += '\nN50:\t%i' % data['N50']
    report += '\nN75:\t%i' % data['N75']
    report += '\nN90:\t%i' % data['N90']
    report += '\nMean contig length:\t%.2f' % data['mean']
    report += '\nMedian contig length:\t%.2f' % data['median']
    report += '\nMinimum contig length:\t%i' % data['min']
    report += '\nMaximum contig length:\t%i' % data['max']
    return report

def mira_stats(contigStatsFile,minLength=0,bins=20,**kwargs):
    """
    Get length, coverage, and GC stats from mira info file
    Returns text with N50 and histograms
    """
    contigStats=pandas.read_csv(contigStatsFile,index_col=0,sep='\t')
    if minLength>0:
        contigStats=contigStats[contigStats.length>=minLength]
    sizes = contigStats['length']
    data = get_contig_length_stats(sizes)
    report = get_contig_length_report(data)

    # add histograms to report
    report += '\nHistograms:\n'
    for key in ['length','GC%','av.cov','mx.cov.','av.qual']:
        report += '\n'
        report += edl.util.asciiHistogram(histogram(contigStats[key],bins=bins),label=key,**kwargs)
    
    return report

def plot_assembly(sizes, file_in, length_data, backend=None,**kwargs):
    min_contig = length_data['min']
    max_contig = length_data['max']
    avg_contig = length_data['mean']
    num_contig = length_data['count']

    if backend:
        import matplotlib
        matplotlib.use(backend)
    from matplotlib import pyplot as plt
    #import pylab
    h=plt.hist(sizes, **kwargs)
    plt.title('%i %s sequences\nLengths %i to %i, Average contig length: %.2f' % (num_contig, file_in, min_contig, max_contig, avg_contig))
    plt.xlabel('Sequence length (bp)')
    plt.ylabel('Count')
    return h

# Get the N50 of the contigs. This is the sequence length at which point
# half of the bases in the entire assembly are contained in contigs of a
# smaller size.
def getN50_old(sizes):
    bases = []
    for read in sorted(sizes):
        for i in range(read):
            bases.append(read)
    return numpy.median(bases)
def getN50(sizes, N=50):
    if N==50:
        compare=getN50_old(sizes)

    totalLength = sum(sizes)
    targetLength = float(totalLength) * N/100.
    totalLength = 0
    for size in sorted(sizes,reverse=True):
        totalLength+=size
        if totalLength>=targetLength:
            if N==50:
                if size != compare:
                    raise Exception("N50s don't match: %s != %s" % (size,compare))
            return size
    else:
        raise Exception("Target length never reached!\nN=%d, target=%d, total=%d" % (N, targetLength, totalLength))

######
# A set of methods for plotting the quality of reseq hits to a set of contigs
def plotHitStats(axes, sequenceFile, hitsFile, 
                 referenceLengths=None,
                 sequenceFormat='fasta', 
                 bins=20, hlog=False, lengthRange=None,
                 barcolor='b', baredgecolor='k', hcolor='r', params=None, 
                 **kwargs):
    """
    Given two or three matplotlib.axes.AxesSubplot objects create plot in each binned by sequence length:
     * overlay a histogram of sequence lengths on the fraction of sequences in each bin that have a hit
     * same bins as above, but use total sequence bases on top of fraction of bases covered by hits
     * if fasta or lengths of reference hits given, plot (using same bins) fraction of reference bases used in hits
     
    Positional Arguments:
     * axes: length 2 list or tuple of ax objects
     * sequenceFile: fasta or similar file of sequence data
     * hitsFile: text hit table
    
    Parameters:
     * hit parsing
      * params=None edl.blatm8.FilterParams object to filter hits
      * **kwargs used to create FilterParams object if params object not given
     * sequence parsing
      * sequenceFormat='fasta'. Can be anything supported by BioPython
      * referenceLengths=None: if give, create 3rd plot using given dictionary of hits. It can also just be the fasta of the reference sequences and the code will look up the lengths.
     * plotting:
      * bins=20 Number of length bins to divide sequence data into
      * barcolor='b' Color of data bars
      * baredgecolor='k' Color of data bar edges
      * hcolor='r' Color of histogram line and axis labels
      * lengthRange=None Can be used to force the x axis to span a specific range 
      * hlog=False If set to True, histogram data plotted in log scale
    """
    
    # get sequence lengths
    lengths=getSequenceLengths(sequenceFile,format=sequenceFormat)
    
    # parse hit file
    if params is None:
        params = edl.blastm8.FilterParams(**kwargs)
    hits=getSequenceHits(hitsFile, params)
    
    # plot data
    plotTranscriptHitRateByLengthBins(axes[0], lengths, hits,
                                      bins=bins, lengthRange=lengthRange,
                                      barcolor=barcolor, 
                                      baredgecolor=baredgecolor,
                                      hcolor=hcolor, hlog=hlog)
    plotTranscriptCoverageByLengthBins(axes[1], lengths, hits,
                                      bins=bins, lengthRange=lengthRange,
                                      barcolor=barcolor, 
                                      baredgecolor=baredgecolor,
                                      hcolor=hcolor, hlog=hlog)
    if referenceLengths is not None:
        plotHitCoverageByLengthBins(axes[2], lengths, hits, referenceLengths,
                                    bins=bins, lengthRange=lengthRange,
                                    barcolor=barcolor, 
                                    baredgecolor=baredgecolor,
                                    hcolor=hcolor, hlog=hlog)

def getSequenceHits(hitsFile, params):
    """
    build a map from sequences to their hits
    """
    sequenceHits={}
    hitCount=0
    m8stream=edl.blastm8.M8Stream(hitsFile)
    for seqid, hits in edl.blastm8.filterM8Stream(m8stream, params, returnLines=False):
        hits=list(hits)
        if len(hits)==0:
            continue
        hitCount+=len(hits)
        sequenceHits[seqid]=hits
    logging.debug("Parsed %d hits for %d sequences fromm %d lines" % (hitCount, len(sequenceHits), m8stream.lines))
    return sequenceHits

def getSequenceLengths(sequenceFile,format='fasta'):
    """
    Get the sequence sizes from a fasta or other file
    """
    if format=='CAF':
        return getContigLengthsFromCAF(sequenceFile)
    sequenceLengths={}
    for record in SeqIO.parse(sequenceFile,format=format):
        sequenceLengths[record.id]=len(record)
    logging.debug("Parsed lengths for %d sequences" % (len(sequenceLengths)))
    return sequenceLengths
    
def plotTranscriptHitRateByLengthBins(ax, lengths, hits, bins=20, lengthRange=None, barcolor='b',baredgecolor='k',hlog=False,hcolor='r'):
    """
    Given a dictionary of transcript lengths and a dictionary of hits,
    Produce a plot of hit rate by length bin.
    """
    
    # Don't try to plot empty data
    if len(lengths)==0:
        raise Exception("Lengths cannot be empty!")
        
    # Draw counts as steppted histogram
    ax2=ax.twinx()
    #transcriptCounts,boundaries=numpy.histogram(lengths.values(), bins=bins, range=lengthRange)
    transcriptCounts,boundaries=ax2.hist(lengths.values(),bins=bins,range=lengthRange, histtype='step',log=hlog,color=hcolor)[:2]
    ax2.set_ylabel('counts',color=hcolor)
    for tl in ax2.get_yticklabels():
        tl.set_color(hcolor)

    # count hits by bin
    hitCounts=numpy.zeros(transcriptCounts.shape)
    for transcript in hits:
        try:
            index=getBin(lengths[transcript],boundaries)
        except ValueError:
            # length was outside range
            continue
        hitCounts[index]+=1

    # normalize hit counts by transcript counts
    hitRate = hitCounts/transcriptCounts
    # remove infinities
    hitRate[transcriptCounts==0]=0
    
    # Draw histogram bars
    lefts=boundaries[:-1]
    widths=[boundaries[i+1]-boundaries[i] for i in range(len(boundaries)-1)]
    ax.bar(lefts,hitRate,width=widths,color=barcolor,edgecolor=baredgecolor)
    ax.set_ylim([0,1])
    ax.set_ylabel('hit rate')
    ax.set_xlabel('transcript length')    
    
def plotTranscriptCoverageByLengthBins(ax, lengths, hits, bins=20, lengthRange=None, barcolor='b',baredgecolor='k',hlog=False,hcolor='r', includeMissed=True):
    """
    Given a dictionary of transcript lengths and a dictionary of hits,
    Produce a plot of coverate rate by length bin. IE: What fracton of total transcript bases were matched.
    """
    
    # Don't try to plot empty data
    if len(lengths)==0:
        raise Exception("Lengths cannot be empty!")
        
    transcriptCounts,boundaries=numpy.histogram(lengths.values(), bins=bins, range=lengthRange)

    # count bases by bin
    hitBaseCounts=numpy.zeros(transcriptCounts.shape)
    totalBaseCounts=numpy.zeros(transcriptCounts.shape)
    for transcript,hitList in hits.iteritems():
        try:
            index=getBin(lengths[transcript],boundaries)
        except ValueError:
            # length was outside range
            continue
        totalBaseCounts[index]+=lengths[transcript]
        hitBaseCounts[index]+=longestHit(hitList)

    if includeMissed:
        for transcript,length in lengths.iteritems():
            if transcript not in hits:
                totalBaseCounts[index]+=length
                
    # Simulate stepped histogram of total bases
    ax2=ax.twinx()
    x,y = getSteppedBars(totalBaseCounts, boundaries)
    if hlog:
        ax2.set_yscale("log",nonposy='clip')
    ax2.plot(x,y,color=hcolor)
    ax2.set_ylabel('total bases',color=hcolor)
    for tl in ax2.get_yticklabels():
        tl.set_color(hcolor)
                
    # normalize hit counts by transcript counts
    hitRate = hitBaseCounts/totalBaseCounts
    # remove infinities
    hitRate[totalBaseCounts==0]=0
    
    # Draw histogram bars
    lefts=boundaries[:-1]
    widths=[boundaries[i+1]-boundaries[i] for i in range(len(boundaries)-1)]
    ax.bar(lefts,hitRate,width=widths,color=barcolor,edgecolor=baredgecolor)
    ax.set_ylim([0,1])
    ax.set_ylabel('bases matched')
    ax.set_xlabel('transcript length')

def plotHitCoverageByLengthBins(ax, lengths, hits, referenceLengths, bins=20, lengthRange=None, barcolor='b',baredgecolor='k',hlog=False,hcolor='r', includeMissed=False):
    """
    Given a dictionary of transcript lengths, a dictionary of hits, and a dict of reference sequence lengths...
    Produce a plot of reference coverate rate by length bin. IE: What fracton of total residues in the reference sequences were matched.

    The param referenceLengths can be a dictionary from hit names to lengths or a fasta file of sequences. The names in both should match the hit names in the "hits" dictionary.
    """
    import screed
    
    # Don't try to plot empty data
    if len(lengths)==0:
        raise Exception("Lengths cannot be empty!")
        
    transcriptCounts,boundaries=numpy.histogram(lengths.values(), bins=bins, range=lengthRange)

    if isinstance(referenceLengths, str):
        # assume we have the path to a fasta file
        # has it been parsed by screed?
        if not os.path.exists("%s_screed" % (referenceLengths)):
            screed.read_fasta_sequences(referenceLengths)
        refScreed=screed.ScreedDB(referenceLengths)
        getHitLength=lambda h: len(refScreed[h]['sequence'])
    else:
        getHitLength=lambda h: referenceLengths[h]

    # count bases by bin
    hitBaseCounts=numpy.zeros(transcriptCounts.shape)
    referenceBaseCounts=numpy.zeros(transcriptCounts.shape)
    totalBaseCounts=numpy.zeros(transcriptCounts.shape)
    for transcript,hitList in hits.iteritems():
        try:
            index=getBin(lengths[transcript],boundaries)
        except ValueError:
            # length was outside range
            continue
        totalBaseCounts[index]+=lengths[transcript]
        firstHit=hitList[0]
        hitLength=getHitLength(firstHit.hit)
        logger.debug("Hit of length %d goes from %d to %d" % (hitLength,
            firstHit.hstart, firstHit.hend))
        referenceBaseCounts[index]+=hitLength
        hitBaseCounts[index]+=numpy.abs(firstHit.hend-firstHit.hstart)+1

    if includeMissed:
        for transcript,length in lengths.iteritems():
            if transcript not in hits:
                totalBaseCounts[index]+=length
                
    # Simulate stepped histogram of total bases
    ax2=ax.twinx()
    x,y = getSteppedBars(totalBaseCounts, boundaries)
    if hlog:
        ax2.set_yscale("log",nonposy='clip')
    ax2.plot(x,y,color=hcolor)
    ax2.set_ylabel('total bases',color=hcolor)
    for tl in ax2.get_yticklabels():
        tl.set_color(hcolor)
                
    # normalize hit counts by transcript counts
    hitRate = hitBaseCounts/referenceBaseCounts
    # remove infinities
    hitRate[totalBaseCounts==0]=0
    
    # Draw histogram bars
    lefts=boundaries[:-1]
    widths=[boundaries[i+1]-boundaries[i] for i in range(len(boundaries)-1)]
    ax.bar(lefts,hitRate,width=widths,color=barcolor,edgecolor=baredgecolor)
    ax.set_ylim([0,1])
    ax.set_ylabel('% reference matched')
    ax.set_xlabel('transcript length')
    
def plotSortedContigLengths(ax, lengths, linecolor='b', log=False, minContigLength=500):
    """
    Given a list of contig lengths, create a stepped plot of ordered lengths.
    """
    
    # Don't try to plot empty data
    if len(lengths)==0:
        raise Exception("Lengths cannot be empty!")
        
    x,y = getSteppedBars(sorted([l for l in lengths if l>minContigLength],
                                reverse=True))
        
    if log:
        ax.set_yscale("log",nonposy='clip')
    objects=ax.plot(x,y,color=linecolor)
    ax.set_ylabel('contig length')
    ax.set_xlabel('contig number')
    return objects

def longestHit(hits):
    return max([numpy.abs(hit.qend-hit.qstart)+1 for hit in hits])
    
def getSteppedBars(values, boundaries=None):
    """
    return lists of x and y coordinates for a stepped line/bar plot given N values (bar heights) and N+1 boundaries (bar edges)
    """
    x=[]
    y=[]
    translate=lambda b: b if numpy.isfinite(b) else 0 

    if boundaries is None:
        boundaries = range(len(values)+1)

    x.append(boundaries[0])
    y.append(0)
    for i in range(len(values)):
        y.append(translate(values[i]))
        x.append(boundaries[i])
        y.append(translate(values[i]))
        x.append(boundaries[i+1])
    x.append(boundaries[-1])
    y.append(0)

    return x,y

def getBin(value, binboundaries):
    """
    Given a value and a set of N+1 ordered values defining N segments or bins,
    return the index of the bin containing value.
    
    Note, this currently uses brute force and could be sped up with a simple dividing by two search
    """
    if binboundaries[0]>value:
        raise ValueError("Value too low")
    for i in range(len(binboundaries)-1):
        if binboundaries[i]<=value and value<=binboundaries[i+1]:
            return i
    raise ValueError("Value too high")

#

sequenceRE=re.compile(r'^Sequence\s+:\s+(\S+)')
dnaRE=re.compile(r'^DNA\s+:\s+(\S+)')
qualityRE=re.compile(r'^BaseQuality\s+:\s+(\S+)')
def getContigLengthsFromCAF(cafFile):
    """
    Return a dictionary of contig names and lengths from a CAF file
    """
    if (isinstance(cafFile,str)):
        cafHandle = open(cafFile)
    else:
        cafHandle = cafFile

    lengths={}
    try:
        while True:
            line=cafHandle.next()
            m = sequenceRE.match(line)
            if m is None:
                continue
            
            # found a sequence header
            seqName=m.group()
            
            if cafHandle.next().strip() != 'Is_contig':
                continue
                
            # it is a contig
            contigName = seqName
            # jump to sequence
            while True:
                line=cafHandle.next()
                if dnaRE.match(line):
                    break

            # parse out sequence
            sequenceLength=0
            while True:
                line=cafHandle.next()
                sequenceLength+=len(line.strip())
                if qualityRE.match(line):
                    break

            lengths[seqName]=sequenceLength
        
    except StopIteration:
        pass
    
    if (isinstance(cafFile,str)):
        cafHandle.close()

    return lengths

if __name__ == '__main__':
    main()
