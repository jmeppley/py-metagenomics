from numpy import histogram
import os
import sys
import numpy
import argparse
import logging
import re
import pandas
from Bio import SeqIO, SeqUtils
logger = logging.getLogger(__name__)


def main():
    """
    Simple hook for running some of the functions below as a script.
    Only works with positional arguments that are strings.

    Examples:
     python assembly.py contig_read_counts file_1
       will run countig_read_counts("file_1")
       """

    # This is a little hack to make this module runnable as a script
    sys.path[0] += "/.."
    from edl.util import asciiHistogram
    import edl.blastm8

    log_level = logging.WARN
    while sys.argv[1] == '-v':
        sys.argv.pop(1)
        log_level -= 10
    logger.setLevel(log_level)
    logger.debug("Log level: {}".format(log_level))

    function = eval(sys.argv[1])
    args = []
    kwargs = {}
    for arg in sys.argv[2:]:
        try:
            param, value = arg.split("=", 1)
            try:
                value = eval(value)
            except NameError:
                pass
            kwargs[param] = value
        except ValueError:
            args.append(arg)

    logger.debug("Function: {}\nArguments: {}\nKWArgs: {}".format(function,
                                                                  args,
                                                                  kwargs))
    function(*args, **kwargs)


def get_contig_stats(contigs_fasta,
                     contig_depth_file=None,
                     contig_read_counts_file=None,
                     contig_stats_file=None,
                     contig_histogram_file=None,
                     **kwargs):
    """
    Extracts GC and length from contigs fasta

    Can optionally merge with read counts and mapped coverage if
    samtools output files given.

    provide a contig_stats_file location to write data to disk instead
    of just returning a pandas DataFrame.

    provide contig_histogram_file to produce a file with summary stats
    and histgrams for each metric. See contig_length_stats() and numpy.
    histogram() for additional kwargs that can be passed when using this
    option.
    """
    # parse contigs fasta
    logger.info("Parsing contig fasta file: {}".format(contigs_fasta))
    contig_stats = get_stats_from_contigs(contigs_fasta)

    # add other files if requested
    if contig_read_counts_file is not None:
        # read counts
        logger.info("Parsing read count file: {}"
                    .format(contig_read_counts_file))
        read_count_table = pandas.read_table(contig_read_counts_file,
                                             delim_whitespace=True,
                                             names=['read count',
                                                    'contig'])\
                                 .set_index('contig')
        contig_stats = contig_stats.join(read_count_table, how='left')

    if contig_depth_file is not None:
        # convert base by base depth data into coverage
        logger.info("Parsing read depth file: {}"
                    .format(contig_depth_file))
        mapping_depth_table = get_samtool_depth_table(contig_depth_file)
        contig_stats = contig_stats.join(mapping_depth_table, how='left')

    # sort and get cumulative length
    contig_stats.fillna(0, inplace=True)
    contig_stats.sort_values(by='length', ascending=False, inplace=True)
    contig_stats['cumul length'] = contig_stats.length.cumsum()
    for col in ['length', 'read count', 'mx cov', 'mn cov', 'cumul length']:
        if col in contig_stats.columns:
            contig_stats[col] = contig_stats[col].astype(int)

    if contig_stats_file is not None:
        logger.info("Writing stats table to: {}".format(contig_stats_file))
        contig_stats.to_csv(contig_stats_file, sep='\t', float_format="%0.2f")

    if contig_histogram_file is not None:
        with open(contig_histogram_file, 'w') as OUTF:
            if 'min_lengths' in kwargs:
                min_lengths = kwargs.pop('min_lengths')
            elif 'min_length' in kwargs:
                min_lengths = [kwargs.pop('min_length'), ]
            else:
                min_lengths = [0, 500, 2000]

            for i, min_length in enumerate(min_lengths):
                logger.info("Making report for contigs >= {}"
                            .format(min_length))
                if i > 0:
                    OUTF.write("\n===============================\
                                ================\n\n")
                OUTF.write("CONTIGS longer or equal to {}bp:\n\n"
                           .format(min_length))
                OUTF.write(contig_length_stats(contig_stats,
                                               return_type='report',
                                               min_length=min_length,
                                               **kwargs
                                               ))

    return contig_stats


def get_stats_from_contigs(contigs_fasta):
    """
    Use BioPython parser and GC calculator to get contig lengths and
    GCs from contigs fasta
    """

    # initialize lists
    contigs = []
    lengths = []
    gcs = []

    # loop over fasta records (this is 2-3 times faster than SeqIO.parse)
    # (and only marginally slower than my custom built parser.)
    with open(contigs_fasta, 'r') as CF:
        for title, sequence in SeqIO.FastaIO.SimpleFastaParser(CF):
            # parse title with RegEx
            contig = title.split(None, 1)[0]
            length = len(sequence)
            contigs.append(contig)
            lengths.append(length)
            gcs.append(SeqUtils.GC(sequence))

    # convert to DataFrame and return
    return pandas.DataFrame({'contig': contigs,
                             'length': lengths,
                             'GC': gcs}).set_index('contig')


def get_samtool_depth_table(depth_file):
    """
    Calculate coverage stats for each contig in an assembly

    Params:
     depth_file: output file from the command:
                    `samtools depth reads.v.contigs.bam`
                 this is a 3 column file with one line per base.
                 columns are:
                     'contig_id base_index base_depth'

    Returns:
     pandas.DataFrame with one row per contig and the three following columns:
            contig  av cov  mx cov
            """
    with open(depth_file, 'r') as DEPTHS:
        return get_samtool_depth_table_from_handle(DEPTHS)


def get_samtool_depth_table_from_handle(depth_stream):
    """
    Calculate coverage stats for each contig in an assembly

    Params:
     depth_stream: output file from the command:
                    `samtools depth reads.v.contigs.bam`

                    passed as an open file-like object (aka a file handle)
                 this is a 3 column file with one line per base.
                 columns are:
                     'contig_id base_index base_depth'

    Returns:
     pandas.DataFrame with one row per contig and the three following columns:
            contig  av cov  mx cov
            """

    # reading into lists is a fast way to build a big DataFrame
    contigs, av_covs, mn_covs, mx_covs = [], [], [], []

    # loop over contig bases
    current_contig = None
    for line in depth_stream:
        contig, base, depth = line.split()
        depth = int(depth)
        if contig != current_contig:
            if current_contig is not None:
                # end of contig, save numbers
                contigs.append(current_contig)
                av_covs.append(depths / bases)
                mn_covs.append(min_depth)
                mx_covs.append(max_depth)
            bases = 0
            depths = 0
            max_depth = depth
            min_depth = depth
            current_contig = contig

        # update contig numbers with current base
        bases += 1
        depths += depth
        min_depth = min(depth, min_depth)
        max_depth = max(depth, max_depth)

    # end of final contig, save numbers
    contigs.append(current_contig)
    av_covs.append(depths / bases)
    mn_covs.append(min_depth)
    mx_covs.append(max_depth)

    return pandas.DataFrame({'contig': contigs,
                             'av cov': av_covs,
                             'mx cov': mx_covs,
                             'mn cov': mn_covs},
                            columns=['contig',
                                     'av cov',
                                     'mn cov',
                                     'mx cov']).set_index('contig')

##
# this evolved from (but now bears little resemblance to) the
# assemlbly_quality_stats.py script by:
# Author: Travis Poulsen
# Date: 09 Feb. 2013
# http://travispoulsen.com/blog/2013/07/basic-assembly-statistics/
# https://gist.github.com/tpoulsen/422b1a19cbd8c0f514fe/raw/assembly_quality_stats.py


def contig_length_stats(contig_stats, return_type=None,
                        txt_width=0,
                        log=False,
                        min_length=0,
                        **kwargs):
    """
    Given contig stats table
     * calculate length stats (including N50)
     * optionally plot histogram (use txt_width and backend to select format)
       (if txt_width is greater than 0 (should be at least 40 for a good plot))
     * return_types:
       None: just print text to STDOUT
       'report': return text
       'data': return dictionary of data
    """
    report_data = {"min_length": min_length}
    contig_stats = contig_stats.loc[contig_stats.length >= min_length]

    if contig_stats.shape[0] == 0:
        report_data['Assembly'] = {'count': 0}
    else:

        report_data['Assembly'] = get_N_stats(contig_stats)
        for column, label in {'length': 'Contig Lengths',
                              'read count': 'Reads per Contig',
                              'av cov': 'Mean Mapped Depth',
                              'mx cov': 'Maximum Mapped Depth',
                              'mn cov': 'Minimum Mapped Depth',
                              'GC': 'GC Content'}.items():
            if column not in contig_stats.columns:
                continue
            report_data[label] = get_column_stats(contig_stats[column])
            if txt_width > 0:
                report_data[label]['log'] = "(log)" if log else ""
                report_data[label]['histogram'] = \
                    asciiHistogram(numpy.histogram(contig_stats[column],
                                                   **kwargs),
                                   log=log,
                                   width=txt_width)

    if return_type == 'data':
        return report_data

    report = get_contig_stats_report(report_data)
    if return_type is None:
        print(report)
    else:
        return report


def get_contig_stats_report(report_data):
    """
    return a formatted string summarizing contig length data
    """
    N_stats = report_data.pop("Assembly")
    if N_stats['count'] == 0:
        return """\
Assembly Summary Stats:
    Contigs: 0
"""

    report = """\
Assembly Summary Stats:
    Contigs: {count}
    N50:     {N50}
    N75:     {N75}
    N90:     {N90}
""".format(**N_stats)

    for column in ['Contig Lengths',
                   'Reads per Contig',
                   'GC Content',
                   'Mean Mapped Depth',
                   'Maximum Mapped Depth',
                   'Minimum Mapped Depth',
                   ]:
        if column not in report_data:
            continue
        report += """
Summary of {column}:
    Min:    {min}
    Max:    {max}
    Mean:   {mean}
    Median: {median}
    StdDev: {std}
""".format(column=column, **report_data[column])
        if 'histogram' in report_data[column]:
            report += """
Histogram of {column} {log}:
{histogram}
""".format(column=column, **report_data[column])
    return report


def get_N_stats(contig_stats, N_levels=[50, 75, 90]):
    """
    uses "length" and "cumulength" columns in contig_stats table
    to quickly get N50 and others for given N_levels (def: 50,75,90)
    """
    N_stats = {'count': contig_stats.shape[0]}

    # tota length is last cumulative length value
    total_length = contig_stats['cumul length'].iloc[-1]
    # set up iterator over just these columns
    cumulen_iter = iter(contig_stats[['length', 'cumul length']].iterrows())

    # Loop over N's. Since they are sorted, we don't need to restart
    #  the length/cumul_length iterator
    cumulength = 0
    for N in sorted(N_levels):
        # looking for N% of the total length
        target = total_length * N / 100

        # stop when we get there
        while cumulength < target:
            contig, (length, cumulength) = next(cumulen_iter)

        # Save the contig length that got use here
        N_key = "N{0:02d}".format(N)
        N_stats[N_key] = length

    return N_stats


def get_column_stats(data):
    """
    return a dict of useful stats
    """
    return {'min': data.min(),
            'max': data.max(),
            'mean': data.mean(),
            'median': data.median(),
            'std': data.std()
            }

##
# the calc_stat, plot_assembly, and getN50 methods originally come from the
# assemlbly_quality_stats.py script by:
# Author: Travis Poulsen
# Date: 09 Feb. 2013
# http://travispoulsen.com/blog/2013/07/basic-assembly-statistics/
# https://gist.github.com/tpoulsen/422b1a19cbd8c0f514fe/raw/assembly_quality_stats.py


def calc_stats(file_in,
               return_type=None,
               txt_width=0,
               log=False,
               backend=None,
               format='fasta',
               minLength=0,
               **kwargs):
    """
    Given contigs in fastsa format:
     * calculate length stats (including N50)
     * plot histogram (use txt_width and backend to select format)
     * return_types:
       None: just print text to STDOUT
       'report': return text
       'data': return dictionary of data
    """
    with open(file_in, 'r') as seq:
        sizes = [len(record) for record in SeqIO.parse(
            seq, format) if len(record) >= minLength]

    sizes = numpy.array(sizes)
    data = get_contig_length_stats(sizes)

    if return_type != 'data':
        report = get_contig_length_report(data)

    if backend is not None:
        h = plot_assembly(sizes, file_in, data, log=log,
                          backend=backend, **kwargs)
    if txt_width > 0:
        if backend is None:
            h = histogram(sizes, **kwargs)
        histogramText = asciiHistogram(h, log=log, width=txt_width)
        if return_type != 'data':
            if log:
                report += "\n\nContig length histogram (log):\n"
            else:
                report += "\n\nContig length histogram:\n"
            report += histogramText
        else:
            data['histogram'] = histogramText

    if return_type == 'data':
        return data
    elif return_type is None:
        print(report)
    else:
        return report


def get_contig_length_stats(sizes):
    """
    return a dict of useful contig length stats
    """
    return {'min': numpy.min(sizes),
            'max': numpy.max(sizes),
            'mean': numpy.mean(sizes),
            'median': numpy.median(sizes),
            'N50': int(getN50(sizes)),
            'N75': int(getN50(sizes, N=75)),
            'N90': int(getN50(sizes, N=90)),
            'count': len(sizes),
            }


def get_contig_length_report(data):
    """
    return a formatted string summarizing contig length data
    """
    report = 'Number of contigs:\t%i' % data['count']
    report += '\nN50:\t%i' % data['N50']
    report += '\nN75:\t%i' % data['N75']
    report += '\nN90:\t%i' % data['N90']
    report += '\nMean contig length:\t%.2f' % data['mean']
    report += '\nMedian contig length:\t%.2f' % data['median']
    report += '\nMinimum contig length:\t%i' % data['min']
    report += '\nMaximum contig length:\t%i' % data['max']
    return report


def mira_stats(contigStatsFile, minLength=0, bins=20, **kwargs):
    """
    Get length, coverage, and GC stats from mira info file
    Returns text with N50 and histograms
    """
    contigStats = pandas.read_csv(contigStatsFile, index_col=0, sep='\t')
    if minLength > 0:
        contigStats = contigStats[contigStats.length >= minLength]
    sizes = contigStats['length']
    data = get_contig_length_stats(sizes)
    report = get_contig_length_report(data)

    # add histograms to report
    report += '\nHistograms:\n'
    for key in ['length', 'GC%', 'av.cov', 'mx.cov.', 'av.qual']:
        report += '\n'
        report += edl.util.asciiHistogram(
            histogram(contigStats[key], bins=bins), label=key, **kwargs)

    return report


def plot_assembly(sizes, file_in, length_data, backend=None, **kwargs):
    min_contig = length_data['min']
    max_contig = length_data['max']
    avg_contig = length_data['mean']
    num_contig = length_data['count']

    if backend:
        import matplotlib
        matplotlib.use(backend)
    from matplotlib import pyplot as plt
    h = plt.hist(sizes, **kwargs)
    plt.title('%i %s sequences\n\
               Lengths %i to %i, \
               Average contig length: %.2f' % (num_contig, file_in,
                                               min_contig, max_contig,
                                               avg_contig))
    plt.xlabel('Sequence length (bp)')
    plt.ylabel('Count')
    return h


def getN50(sizes, N=50):
    """
    Get the N50 of the contigs. This is the sequence length at which point
    half of the bases in the entire assembly are contained in contigs of a
    smaller size.
    """
    totalLength = sum(sizes)
    targetLength = float(totalLength) * N / 100.
    totalLength = 0
    for size in sorted(sizes, reverse=True):
        totalLength += size
        if totalLength >= targetLength:
            return size
    else:
        raise Exception("Target length never reached!\
                         \nN=%d, target=%d, total=%d" % (N,
                                                         targetLength,
                                                         totalLength))

######
# A set of methods for plotting the quality of reseq hits to a set of contigs


def plotHitStats(axes, sequenceFile, hitsFile,
                 referenceLengths=None,
                 sequenceFormat='fasta',
                 bins=20, hlog=False, lengthRange=None,
                 barcolor='b', baredgecolor='k', hcolor='r', params=None,
                 **kwargs):
    """
    Given two or three matplotlib.axes.AxesSubplot objects create plot in
    each binned by sequence length:

     * overlay a histogram of sequence lengths on the fraction of sequences
       in each bin that have a hit
     * same bins as above, but use total sequence bases on top of fraction
       of bases covered by hits
     * if fasta or lengths of reference hits given, plot (using same bins)
       fraction of reference bases used in hits

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
      * referenceLengths=None: if give, create 3rd plot using given
        dictionary of hits. It can also just be the fasta of the reference
        sequences and the code will look up the lengths.
     * plotting:
      * bins=20 Number of length bins to divide sequence data into
      * barcolor='b' Color of data bars
      * baredgecolor='k' Color of data bar edges
      * hcolor='r' Color of histogram line and axis labels
      * lengthRange=None Can be used to force the x axis to span a
        specific range
      * hlog=False If set to True, histogram data plotted in log scale
    """

    # get sequence lengths
    lengths = getSequenceLengths(sequenceFile, format=sequenceFormat)

    # parse hit file
    if params is None:
        params = edl.blastm8.FilterParams(**kwargs)
    hits = getSequenceHits(hitsFile, params)

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
    sequenceHits = {}
    hitCount = 0
    m8stream = edl.blastm8.M8Stream(hitsFile)
    for seqid, hits in edl.blastm8.filterM8Stream(m8stream,
                                                  params,
                                                  returnLines=False):
        hits = list(hits)
        if len(hits) == 0:
            continue
        hitCount += len(hits)
        sequenceHits[seqid] = hits
    logging.debug("Parsed %d hits for %d sequences fromm %d lines" %
                  (hitCount, len(sequenceHits), m8stream.lines))
    return sequenceHits


def getSequenceLengths(sequenceFile, format='fasta'):
    """
    Get the sequence sizes from a fasta or other file
    """
    if format == 'CAF':
        return getContigLengthsFromCAF(sequenceFile)
    sequenceLengths = {}
    for record in SeqIO.parse(sequenceFile, format=format):
        sequenceLengths[record.id] = len(record)
    logging.debug("Parsed lengths for %d sequences" % (len(sequenceLengths)))
    return sequenceLengths


def plotTranscriptHitRateByLengthBins(ax, lengths, hits,
                                      bins=20,
                                      lengthRange=None,
                                      barcolor='b',
                                      baredgecolor='k',
                                      hlog=False,
                                      hcolor='r'):
    """
    Given a dictionary of transcript lengths and a dictionary of hits,
    Produce a plot of hit rate by length bin.
    """

    # Don't try to plot empty data
    if len(lengths) == 0:
        raise Exception("Lengths cannot be empty!")

    # Draw counts as steppted histogram
    ax2 = ax.twinx()
    transcriptCounts, boundaries = ax2.hist(lengths.values(),
                                            bins=bins,
                                            range=lengthRange,
                                            histtype='step',
                                            log=hlog,
                                            color=hcolor)[:2]
    ax2.set_ylabel('counts', color=hcolor)
    for tl in ax2.get_yticklabels():
        tl.set_color(hcolor)

    # count hits by bin
    hitCounts = numpy.zeros(transcriptCounts.shape)
    for transcript in hits:
        try:
            index = getBin(lengths[transcript], boundaries)
        except ValueError:
            # length was outside range
            continue
        hitCounts[index] += 1

    # normalize hit counts by transcript counts
    hitRate = hitCounts / transcriptCounts
    # remove infinities
    hitRate[transcriptCounts == 0] = 0

    # Draw histogram bars
    lefts = boundaries[:-1]
    widths = [boundaries[i + 1] - boundaries[i]
              for i in range(len(boundaries) - 1)]
    ax.bar(lefts, hitRate, width=widths,
           color=barcolor, edgecolor=baredgecolor)
    ax.set_ylim([0, 1])
    ax.set_ylabel('hit rate')
    ax.set_xlabel('transcript length')


def plotTranscriptCoverageByLengthBins(ax, lengths, hits,
                                       bins=20,
                                       lengthRange=None,
                                       barcolor='b',
                                       baredgecolor='k',
                                       hlog=False,
                                       hcolor='r',
                                       includeMissed=True):
    """
    Given a dictionary of transcript lengths and a dictionary of hits,
    Produce a plot of coverate rate by length bin. IE: What fracton of
     total transcript bases were matched.
    """

    # Don't try to plot empty data
    if len(lengths) == 0:
        raise Exception("Lengths cannot be empty!")

    transcriptCounts, boundaries = numpy.histogram(
        lengths.values(), bins=bins, range=lengthRange)

    # count bases by bin
    hitBaseCounts = numpy.zeros(transcriptCounts.shape)
    totalBaseCounts = numpy.zeros(transcriptCounts.shape)
    for transcript, hitList in hits.iteritems():
        try:
            index = getBin(lengths[transcript], boundaries)
        except ValueError:
            # length was outside range
            continue
        totalBaseCounts[index] += lengths[transcript]
        hitBaseCounts[index] += longestHit(hitList)

    if includeMissed:
        for transcript, length in lengths.iteritems():
            if transcript not in hits:
                totalBaseCounts[index] += length

    # Simulate stepped histogram of total bases
    ax2 = ax.twinx()
    x, y = getSteppedBars(totalBaseCounts, boundaries)
    if hlog:
        ax2.set_yscale("log", nonposy='clip')
    ax2.plot(x, y, color=hcolor)
    ax2.set_ylabel('total bases', color=hcolor)
    for tl in ax2.get_yticklabels():
        tl.set_color(hcolor)

    # normalize hit counts by transcript counts
    hitRate = hitBaseCounts / totalBaseCounts
    # remove infinities
    hitRate[totalBaseCounts == 0] = 0

    # Draw histogram bars
    lefts = boundaries[:-1]
    widths = [boundaries[i + 1] - boundaries[i]
              for i in range(len(boundaries) - 1)]
    ax.bar(lefts, hitRate, width=widths,
           color=barcolor, edgecolor=baredgecolor)
    ax.set_ylim([0, 1])
    ax.set_ylabel('bases matched')
    ax.set_xlabel('transcript length')


def plotHitCoverageByLengthBins(ax, lengths, hits, referenceLengths,
                                bins=20,
                                lengthRange=None,
                                barcolor='b',
                                baredgecolor='k',
                                hlog=False,
                                hcolor='r',
                                includeMissed=False):
    """
    Given a dictionary of transcript lengths, a dictionary
     of hits, and a dict of reference sequence lengths...
    Produce a plot of reference coverate rate by length bin.
     IE: What fracton of total residues in the reference sequences were
         matched.

    The param referenceLengths can be a dictionary from hit names to
    lengths or a fasta file of sequences. The names in both should
    match the hit names in the "hits" dictionary.
    """
    # Don't try to plot empty data
    if len(lengths) == 0:
        raise Exception("Lengths cannot be empty!")

    transcriptCounts, boundaries = numpy.histogram(
        lengths.values(), bins=bins, range=lengthRange)

    get_hit_length = build_get_hit_length_function(referenceLengths)

    # count bases by bin
    hitBaseCounts = numpy.zeros(transcriptCounts.shape)
    referenceBaseCounts = numpy.zeros(transcriptCounts.shape)
    totalBaseCounts = numpy.zeros(transcriptCounts.shape)
    for transcript, hitList in hits.iteritems():
        try:
            index = getBin(lengths[transcript], boundaries)
        except ValueError:
            # length was outside range
            continue
        totalBaseCounts[index] += lengths[transcript]
        firstHit = hitList[0]
        hitLength = getHitLength(firstHit.hit)
        logger.debug("Hit of length %d goes from %d to %d" % (hitLength,
                                                              firstHit.hstart,
                                                              firstHit.hend))
        referenceBaseCounts[index] += hitLength
        hitBaseCounts[index] += numpy.abs(firstHit.hend - firstHit.hstart) + 1

    if includeMissed:
        for transcript, length in lengths.iteritems():
            if transcript not in hits:
                totalBaseCounts[index] += length

    # Simulate stepped histogram of total bases
    ax2 = ax.twinx()
    x, y = getSteppedBars(totalBaseCounts, boundaries)
    if hlog:
        ax2.set_yscale("log", nonposy='clip')
    ax2.plot(x, y, color=hcolor)
    ax2.set_ylabel('total bases', color=hcolor)
    for tl in ax2.get_yticklabels():
        tl.set_color(hcolor)

    # normalize hit counts by transcript counts
    hitRate = hitBaseCounts / referenceBaseCounts
    # remove infinities
    hitRate[totalBaseCounts == 0] = 0

    # Draw histogram bars
    lefts = boundaries[:-1]
    widths = [boundaries[i + 1] - boundaries[i]
              for i in range(len(boundaries) - 1)]
    ax.bar(lefts, hitRate, width=widths,
           color=barcolor, edgecolor=baredgecolor)
    ax.set_ylim([0, 1])
    ax.set_ylabel('% reference matched')
    ax.set_xlabel('transcript length')


def build_get_hit_length_function(referenceLengths):
    """
    Given the referenceLengths parameter return a lambda function that will
    map a reference sequence id to its sequence length

    The referenceLenths parameter may be either a python dict or a str name
    of a fasta file. In the latter case, the file is parsed to get lengths
    """
    if isinstance(referenceLengths, str):
        import screed
        # assume we have the path to a fasta file
        # has it been parsed by screed?
        if not os.path.exists("%s_screed" % (referenceLengths)):
            # TODO: just use Bio.SeqIO to get lengths if
            #   screed module or screed index is missing.
            #   screed is overkill here.
            screed.read_fasta_sequences(referenceLengths)
        refScreed = screed.ScreedDB(referenceLengths)

        return lambda h: len(refScreed[h]['sequence'])
    else:
        return lambda h: referenceLengths[h]


def plotSortedContigLengths(ax, lengths,
                            linecolor='b',
                            log=False,
                            minContigLength=500):
    """
    Given a list of contig lengths, create a stepped plot of ordered lengths.
    """

    # Don't try to plot empty data
    if len(lengths) == 0:
        raise Exception("Lengths cannot be empty!")

    x, y = getSteppedBars(sorted([l for l in lengths if l > minContigLength],
                                 reverse=True))

    if log:
        ax.set_yscale("log", nonposy='clip')
    objects = ax.plot(x, y, color=linecolor)
    ax.set_ylabel('contig length')
    ax.set_xlabel('contig number')
    return objects


def longestHit(hits):
    return max([numpy.abs(hit.qend - hit.qstart) + 1 for hit in hits])


def getSteppedBars(values, boundaries=None):
    """
    return lists of x and y coordinates for a stepped line/bar plot
    given N values (bar heights) and N+1 boundaries (bar edges)
    """
    x = []
    y = []

    def translate(b):
        return b if numpy.isfinite(b) else 0

    if boundaries is None:
        boundaries = range(len(values) + 1)

    x.append(boundaries[0])
    y.append(0)
    for i in range(len(values)):
        y.append(translate(values[i]))
        x.append(boundaries[i])
        y.append(translate(values[i]))
        x.append(boundaries[i + 1])
    x.append(boundaries[-1])
    y.append(0)

    return x, y


def getBin(value, binboundaries):
    """
    Given a value and a set of N+1 ordered values defining N segments or bins,
    return the index of the bin containing value.

    Note, this currently uses brute force and could be sped up with a
    simple dividing by two search
    """
    if binboundaries[0] > value:
        raise ValueError("Value too low")
    for i in range(len(binboundaries) - 1):
        if binboundaries[i] <= value and value <= binboundaries[i + 1]:
            return i
    raise ValueError("Value too high")

#

sequenceRE = re.compile(r'^Sequence\s+:\s+(\S+)')
dnaRE = re.compile(r'^DNA\s+:\s+(\S+)')
qualityRE = re.compile(r'^BaseQuality\s+:\s+(\S+)')


def getContigLengthsFromCAF(cafFile):
    """
    Return a dictionary of contig names and lengths from a CAF file
    """
    if (isinstance(cafFile, str)):
        cafHandle = open(cafFile)
    else:
        cafHandle = cafFile

    lengths = {}
    try:
        while True:
            line = cafHandle.next()
            m = sequenceRE.match(line)
            if m is None:
                continue

            # found a sequence header
            seqName = m.group()

            if cafHandle.next().strip() != 'Is_contig':
                continue

            # it is a contig
            contigName = seqName
            # jump to sequence
            while True:
                line = cafHandle.next()
                if dnaRE.match(line):
                    break

            # parse out sequence
            sequenceLength = 0
            while True:
                line = cafHandle.next()
                sequenceLength += len(line.strip())
                if qualityRE.match(line):
                    break

            lengths[seqName] = sequenceLength

    except StopIteration:
        pass

    if (isinstance(cafFile, str)):
        cafHandle.close()

    return lengths

if __name__ == '__main__':
    main()
