"""
A collection of methods and classes to simplify interactive parsing of
hit tables
"""
from edl.blastm8 import *
from edl.hits import *
import numpy as np
import logging
import pandas
logger = logging.getLogger(__name__)


def countHits(infile, **kwargs):
    """
    Count hits from a hit table.

    Calls edl.hits.parseM8FileIter with the following optional parameters:
        hitStringMap (None): dictionary (or file) mapping hit IDs to
        something else
        format (GENE): hit table format
        filterPct (0): only consider hits within this % of top score for
        each read
        parseStyle (ACCS): how to process hit data into an identifying string
        countMethod ('all'): how to resolve hits to multiple sequences
        taxonomy (None): An edl.taxon.Taxonomy object or directory
        conatining taxdmp
        rank (None): Maximum rank to resolve hits
    """

    # if taxonomy or hitStringMap are file names, parse them
    taxonomy = kwargs.pop('taxonomy', None)
    if isinstance(taxonomy, str):
        taxonomy = readTaxonomy(
            taxonomy, namesMap=kwargs.pop(
                'namesMap', False))
    hitStringMap = kwargs.pop('hitStringMap', None)
    if isinstance(hitStringMap, str):
        if taxonomy is not None:
            # the mapped hit ids will need to be ints
            valueType = kwargs.pop('valueType', int)
        else:
            valueType = kwargs.pop('valueType', None)
        hitStringMap = parseMapFile(hitStringMap, valueType=valueType)

    # if infile is name (and not handle), open as a handle
    if isinstance(infile, str):
        inhandle = open(infile)
    else:
        inhandle = infile

    # get iterator over reads that will parse hits
    hitIter = parseM8FileIter(inhandle,
                              hitStringMap,
                              kwargs.pop('format', GENE),
                              kwargs.pop('filterPct', 0),
                              kwargs.pop('parseStyle', ACCS),
                              kwargs.pop('countMethod', 'all'),
                              taxonomy=taxonomy,
                              rank=kwargs.pop('rank', None))

    # count the hits
    (total, counts) = countIterHits(hitIter,
                                    allMethod=kwargs.pop('allMethod', ALLEQ),
                                    returnMap=False)

    logger.info("Total hits: %s" % total)
    if isinstance(infile, str):
        inhandle.close()

    return counts


def countRefSeqHits(filename, **kwargs):
    """
    Counts hits from a hit table with defaults that make sense for a
    RefSeq search.


    Calls countHits with these defaults (that can be overriden):
        hitStringMap:
                /common/FASTA/NCBI/RefSeq/LATEST/acc.to.taxid.protein.plus
        taxonomy: read from /common/FASTA/NCBI/RefSeq/LATEST/taxdump

    You can also override the defaults with the parameter 'parseFilePrefix'
    which will replace '/common/FASTA/NCBI/RefSeq/LATEST/' in the defaults
    above.

    Plus the defaults set in countHits:
        format (GENE): hit table format
        filterPct (0): only consider hits within this % of top score for each
        read
        parseStyle (ACCS): how to process hit data into an identifying string
        countMethod ('all'): how to resolve hits to multiple sequences
        rank (None): Maximum rank to resolve hits
    """

    parseFilePrefix = kwargs.pop(
        'parseFilePrefix',
        '/common/FASTA/NCBI/RefSeq/LATEST/')
    kwargs.setdefault(
        'hitStringMap',
        parseFilePrefix +
        'acc.to.taxid.protein.plus')
    kwargs.setdefault('taxonomy', parseFilePrefix + 'taxdump')

    return countHits(filename, **kwargs)


def getCountDataFrame(*args, **kwargs):
    """
    For each file name (in *args), count hits (with options in **kwargs).
    Merge all counts into a dataframe and return.
    Pass useRefSeq=True to use countRefSeqHits defaults.
    """
    return getCountDataFramePanda(*args, **kwargs)


def getCountDataFrameNp(*args, **kwargs):
    """
    For each file name (in *args), count hits (with options in **kwargs).
    Merge all counts into a dataframe and return.
    Pass useRefSeq=True to use countRefSeqHits defaults.
    """
    if kwargs.pop('useRefSeq', False):
        countHitsFn = countRefSeqHits
    else:
        countHitsFn = countHits

    countsByFile = {}
    keys = []
    for filename in args:
        counts = countHitsFn(filename, **kwargs)
        keys.extend(counts.keys())
        countsByFile[filename] = counts

    index = pandas.Index(set(keys))
    arrays = []
    for filename in args:
        data = []
        for key in index:
            data.append(countsByFile[filename].get(key, 0))
        arrays.append(np.array(data))
    return pandas.DataFrame(arrays, index=args, columns=index)


def getCountDataFramePanda(*args, **kwargs):
    """
    For each file name (in *args), count hits (with options in **kwargs).
    Merge all counts into a dataframe and return.
    Pass useRefSeq=True to use countRefSeqHits defaults.
    """
    if kwargs.pop('useRefSeq', False):
        countHitsFn = countRefSeqHits
    else:
        countHitsFn = countHits

    countsByFile = {}
    keys = []
    for filename in args:
        counts = countHitsFn(filename, **kwargs)
        keys.extend(counts.keys())
        countsByFile[filename] = counts

    return pandas.DataFrame(countsByFile)
