from edl.hits import parseHits
from edl.util import tupleIteratorToMap
import logging
import pandas
logger = logging.getLogger(__name__)


def collapseDataFrame(frame, indices, dest='Other', axis=1, how='sum'):
    """
    Collapse all the matching rows or columns into a single one
    """
    if axis == 0:
        frame = frame.transpose()

    indices = list(indices)
    if dest in frame.columns and dest not in indices:
        indices.append(dest)

    newColumn = pandas.DataFrame(frame[indices].sum(axis=1), columns=[dest])
    newFrame = frame.drop(indices, axis=1).join(newColumn)

    if axis == 0:
        return newFrame.transpose()
    else:
        return newFrame


def crossTabulateHits(hitFile1, hitFile2,
                      skipHeader1=True, skipHeader2=True,
                      readCol1=0, readCol2=0,
                      hitCol1=-1, hitCol2=-1,
                      hitSep1=None, hitSep2=None):
    """
    Given two hit tables, cross reference reads to build a 2D DataFrame
    of hit counts
    """
    # get a map from reads to hits in file 1
    hitmap = tupleIteratorToMap(
        parseHits(
            hitFile1,
            readCol1,
            hitCol1,
            skipHeader1,
            hitSep1))
    # get an iterator over hits in file 2
    crossHits = parseHits(hitFile2, readCol2, hitCol2, skipHeader2, hitSep2)

    counts = {}
    types2 = set()
    readsOnlyIn2 = []
    readMatchCount = 0
    for read, hits in crossHits:
        try:
            hits2 = hitmap.pop(read)
            readMatchCount += 1
        except KeyError:
            readsOnlyIn2.append(read)
            continue

        for h2 in hits2:
            for h1 in hits:
                h1counts = counts.setdefault(h1, {})
                h1counts[h2] = h1counts.setdefault(h2, 0) + 1
            types2.add(h2)

    # generate dataframe
    logger.warn(
        "%d of %d reads only in file 1" %
        (len(hitmap), readMatchCount))
    logger.warn(
        "%d of %d reads only in file 2" %
        (len(readsOnlyIn2), readMatchCount))
    return pandas.DataFrame(counts)
