"""
Pick a best hit for each read based on abundances from singular hits
"""
import logging
from urllib.parse import quote_plus
import numpy as np
from edl.hits import binAndMapHits, getHitTranslator, HITID
from edl import blastm8
logger = logging.getLogger(__name__)


def pickBestHitByAbundance(m8stream,
                           filterParams=None,
                           returnLines=True,
                           returnTranslations=False,
                           organismCounts=None,
                           winnerTakeAll=False,
                           sequenceWeights=None,
                           **kwargs):
    """
    Given a hit table with (potentially) multiple hits for each read.
    Select the best hit for each read. Hits are parsed from given hit
    table (m8stream) if given a FilterParams object, otherwise it is
    assumed that m8stream is an iterator over Hit objects. Remaining
    keyword arguments are used to translate hits to accessions, organisms,
    or anything else using a HitTranslator.

    Ambiguous hits (multiple 'best' hits to one read) are resolved as follows:
        given a set of reads that all hit the same list of translated hits:
            divvy up reads so that the abundance ratios change minimally

    Abundance is recorded for whatever the Hittranslator returns. If a hit
    map and taxonomy are given, this will be organisms, if only the parseStyle
    is given and it's set to ACC, then accessions will be the currency. The
    default is HITID.

    Yields (read,hit) tuples, (read, [translated hits]) tuples, or hit
    table lines.
    """
    if returnLines and returnTranslations:
        returnLines = False
        logger.warn("returnTranslations overrides returnLines!")

    # filtered hits
    if filterParams is None:
        hitIter = m8stream
    else:
        hitIter = blastm8.filterM8Stream(
            m8stream, filterParams, returnLines=False)

    # custom function for pulling orgs from hits
    #  if no settings given, just use the hit ID as the 'organism'
    kwargs.setdefault("parseStyle", HITID)
    hitTranslator = getHitTranslator(**kwargs)

    # we need to keep track of lots of things
    orgCounts = {}
    totalReads = 0
    unambiguousReads = 0
    ambiguousReads = 0
    sameOrgCount = 0
    ambiguousHits = {}

    # Check to see if organism counts were given
    if organismCounts is not None:
        if isinstance(organismCounts, str):
            organismCounts = getOrganismCountsFromFile(organismCounts)

    # loop over hits and yield unambiguous ones
    # Save ambiguous hits and org abundances
    logger.debug(str(hitIter))
    for (read, hits) in hitIter:
        logger.debug("Read: %s" % (read))
        totalReads += 1
        hitByOrg = {}
        orgs = []
        count = 0
        for hit in hits:
            count += 1
            hitOrgs = hitTranslator.translateHit(hit)
            logger.debug("Hit: %s (%s), %s" % (hit.hit, hitOrgs, hit.score))
            orgs.extend(hitOrgs)
            for org in hitOrgs:
                if org in hitByOrg:
                    # This should be REALLY rare.
                    sameOrgCount += 1
                    sameOrgExample = (read, hit.hit, org)
                    logger.warn(
                        "Read (%s) has two best hits to same org (%s)!" %
                        (read, org))
                    # always keep the first alphabetically, for reproducibility
                    if hit.hit < hitByOrg[org].hit:
                        hitByOrg[org] = hit
                else:
                    hitByOrg[org] = hit
        orgs = tuple(sorted(set(orgs)))
        if count == 0:
            # This *should* never happen
            logger.error("No hits for %s!!!!!" % (read))
            raise Exception(
                "Read (%s) has not hits. This shouldn't happen." %
                (read))
        elif count == 1 or len(hitByOrg) == 1:
            logger.debug("Read is UNambiguous")
            unambiguousReads += 1
            for org in orgs:
                if sequenceWeights is not None:
                    increment = sequenceWeights.get(read, 1)
                else:
                    increment = 1
                orgCounts[org] = orgCounts.get(org, 0) + increment
            if returnLines:
                yield hit.line
            elif returnTranslations:
                yield (read, orgs)
            else:
                yield (read, hit)
        else:
            logger.debug("Read IS ambiguous")
            ambiguousReads += 1
            if organismCounts is None:
                # If we don't have count data to start, save these til the end
                ambiguousHits.setdefault(orgs, []).append(hitByOrg)
            else:
                # Use given counts to resolve
                for (
                    hit, org) in assignHits(
                    orgs, [
                        hitByOrg, ], organismCounts, winnerTakeAll):
                    yield formatReturn(hit,
                                       org,
                                       returnLines,
                                       returnTranslations)

    logger.info("Processed %d reads:" % (totalReads))
    logger.info(
        "Collected unambiguous counts for %d orgs from %d reads" %
        (len(orgCounts), unambiguousReads))

    # if we used given organism counts, then we are done
    if organismCounts is not None:
        return

    # otherwise, we have ambiguous reads to resolve
    logger.info(
        "Need to resolve %d ambiguous reads hitting %d orgs" %
        (ambiguousReads, len(ambiguousHits)))

    if sameOrgCount > 0:
        elements = list(sameOrgExample)
        elements.insert(0, sameOrgCount)
        logger.warn(
            "found %d cases where a read had an extra hit to the same "
            "organism. For Example: %s (%s,%s)" %
            tuple(elements))

    # loop over ambiguous hits (grouped by possible orgs) and pick one for
    # each read
    ambiguousReads = 0
    # for orgs, hits in ambiguousHits.items():
    for orgs in sorted(ambiguousHits.keys()):
        hits = ambiguousHits[orgs]
        for (hit, org) in assignHits(orgs, hits, orgCounts, winnerTakeAll):
            ambiguousReads += 1
            yield formatReturn(hit, org, returnLines, returnTranslations)

    logger.info(
        "Selected top hit for %d ambiguous reads for a total of %d "
        "returned hit assignments" %
        (ambiguousReads, ambiguousReads + unambiguousReads))


def formatReturn(hit, org, returnLines, returnTranslations):
    if returnLines:
        return hit.line
    elif returnTranslations:
        return (hit.read, [org, ])
    else:
        return (hit.read, hit)


def getOrganismCountsFromFile(orgCountFile):
    """
    Parse the given file into dictionary from organism name to count of hits
    """
    orgCounts = {}
    with open(orgCountFile) as f:
        for line in f:
            (organism, counts) = line.split(None, 1)
            orgCounts[organism] = int(counts)
    return orgCounts


def assignHits(orgs, hits, orgCounts, winnerTakeAll):
    """
    Returns a genererator over hits.

    Given a subset of organisms, a dictionary of organism counts,
     and a list of hit sets
      where each element is a dictionary from organsism to blastm8.Hit object

    Return for each set, a single hit.

    If winnerTakeAll is true, return the hit corresponding to the
    most abundant organism
    If not, return hits such that the returned proportion of orgnisms
    is close to the proportion in the global hit count dictionary.
    """
    if winnerTakeAll:
        return assignHitsToMostCommon(orgs, hits, orgCounts)
    else:
        return assignHitsByProportion(orgs, hits, orgCounts)


def assignHitsToMostCommon(orgs, hits, orgCounts):
    """
    For each ambiguous hit, return the hit to the most abundant organism.
    Tie-break alphabetically.
    """
    # Find the most abundant organism
    mostAbundantOrg = None
    highestCount = 0
    for org in orgs:
        ocount = orgCounts.get(org, .5)
        if ocount > highestCount:
            highestCount = ocount
            mostAbundantOrg = org

    # for each set of hits, return the hit to this organism
    for hitsByOrg in hits:
        yield (hitsByOrg[mostAbundantOrg], mostAbundantOrg)


def assignHitsByProportion(orgs, hits, orgCounts):
    """
    Deterministic version
    Given:
        tuple of organisms
        list of hit dicts:
            each is a map from organisms to hit objects
            map from organisms to unambiguous hit counts
    Yield for each hit dictionary:
        tuple (hit,org) so that count ratios of returned orgs match ratio
        of unambiguous counts
    """
    # calculate fraction that should go to each organism
    indices = range(len(orgs))
    ocounts = [orgCounts.get(org, .5) for org in orgs]
    ocounts = np.array(ocounts) / float(sum(ocounts))

    # expected # hits per org
    hits = list(hits)
    hitCounts = ocounts * len(hits)

    # for each org, pick n=floor(expected # hits) reads to assign to hit for
    # that org
    for i in sorted(indices, key=lambda i: hitCounts[i], reverse=True):
        org = orgs[i]
        while hitCounts[i] >= 1 and len(hits) > 0:
            hitCounts[i] -= 1
            hit = hits.pop()
            yield (hit[org], org)

    # sort orgs by remainder, and assign remaining hits in that order
    orgsByRemainder = [orgs[i]
                       for i in sorted(indices, key=lambda i: hitCounts[i])]
    while len(hits) > 0:
        hit = hits.pop()
        org = orgsByRemainder.pop()
        yield (hit[org], org)


###############
# Given LCA-like orgnaism assignments, push count back up to tree tips
#   (or max rank)
#  Total org/taxa counts should be good, but individual read assignments
#    may be suspect
###############
def redistributeHits(hitMap, rank):
    """
    given a map of read ids to assigned Taxa, redistribute any hits that
    have been
    dropped to LowestCommonAncestor to higher descendants based on
    proportion of
    unambiguous hits
    """
    logger.info("redistributing LCA assignments")
    if isinstance(hitMap, dict):
        hits = binHits(hitMap)
    else:
        (hits, hitMap) = binAndMapHits(hitMap)

    root = hits.keys()[0].getRootNode()
    redistributeHitsForNode(root, hits, rank)

    # generate new hitmap
    for (node, nhits) in hits.items():
        if isinstance(nhits, list):
            for read in nhits:
                hitMap[read] = [node, ]
        else:
            logger.warn("hits (%s) not a list for %s!" % (nhits, node.name))

    return hitMap


def redistributeHitsForNode(node, hits, rank):
    """
    recursive call used by redistributeHits
    if rank is not rank and any children have hits, redistribute
        hits to children
    """

    if rank is not None and rank == node.rank:
        logger.debug("Node %s has rank %s, skipping" % (node.name, rank))
        return

    nodeHits = hits.get(node, [])
    if not isinstance(nodeHits, list):
        nodeHits = [nodeHits, ]
        hits[node] = nodeHits
        logger.warn("Hit for %s was a list: %s" % (node.name, nodeHits[0]))
    nodeHitCount = len(nodeHits)

    # check children for hits
    childCounts = {}
    total = 0
    for child in node.children:
        if child is node:
            # root is sometimes a child of itself!
            continue
        # add up all hits to child and its children and on up the tree
        kidHits = getTotalHits(child, hits)
        if kidHits > 0:
            total += kidHits
            childCounts[child] = kidHits
    logger.debug("Child hits: %s" % (childCounts))

    if nodeHitCount != 0:
        if total > 0:
            # redistribute
            logger.debug(
                "Redistributing %d hits from %s to %s" %
                (nodeHitCount, node.name, [
                    n.name for n in childCounts.keys()]))
            logger.debug(str(childCounts))
            remainders = {}
            for child in sorted(
                    childCounts.keys(),
                    key=lambda kid: childCounts[kid],
                    reverse=True):
                logger.debug("---\n%s\n----" % (nodeHits))
                # calculate number of hits for this child (as a float)
                newKidHits = nodeHitCount * childCounts[child] / float(total)
                logger.debug(
                    "Add %f = %d * %d / %d new hits to %s" %
                    (newKidHits, nodeHitCount, childCounts[child],
                        total, child.name))

                if child not in hits:
                    # special case where child has children with hits, but no
                    # hits itself
                    hits[child] = []

                # move hits one at a time to child until remainder <1
                while newKidHits >= 1 and len(nodeHits) > 0:
                    nextHit = nodeHits.pop(0)
                    logger.debug(
                        "nkh: %f  child: %s hit: %s" %
                        (newKidHits, child.name, nextHit))
                    hits[child].append(nextHit)
                    newKidHits -= 1
                remainders[child] = newKidHits

            # sort children by remainder size and assign remaining hits in that
            # order
            logger.debug(
                "%d hits left. Remainders: %s" %
                (len(nodeHits), remainders))
            mostDeserving = sorted(
                remainders.keys(),
                key=lambda kid: remainders[kid],
                reverse=True)
            while len(nodeHits) > 0:
                hits.get(mostDeserving.pop(0), []).append(nodeHits.pop(0))

    # Now call this method on all tthe children recursively
    for child in childCounts.keys():
        redistributeHitsForNode(child, hits, rank)


def getTotalHits(node, hits):
    nodeHits = len(hits.get(node, ()))
    for child in node.children:
        nodeHits += getTotalHits(child, hits)
    return nodeHits


def multipleFileWrapper(m8files):
    """
    return a generator over all lines in the given list of files. The
    generator will have an object variable (.fileDict) that maps read
    names back to the source file

    Read names are altered to have the url-enoded file name as a prefix.
    To get file name or tag from altered reads, use:
        from urllib.parse import unquote_plus
        encoded_tag, read_name = read_name.split('/',1)
        tag = unquote_plus(encoded_tag)
    """
    return _multipleFileGeneratorPrefixed(m8files)


def _multipleFileGeneratorPrefixed(m8files):
    """
    Iterate over all lines in a given list of streams as a single stream
    of hits, but keep track of which read came from which file in the given
    dictionary by renaming reads (adding file name to start)
    """
    for m8file in m8files:
        # allow for list of (file,tag) tuples
        if isinstance(m8file, tuple):
            m8file, fileTag = m8file
        else:
            fileTag = None

        # get file stream
        m8stream = M8Stream(m8file, file_tag=fileTag)

        for line in m8stream:
            yield line


class M8Stream(blastm8.M8Stream):
    """
    Variant of M8Stream that encodes a file tag into the read names
    """

    def __init__(self, fileName, *args, file_tag=None):
        blastm8.M8Stream.__init__(self, fileName, *args)
        if file_tag:
            self.file_tag = quote_plus(file_tag)
        else:
            self.file_tag = quote_plus(self.fileName)

    def __next__(self):
        line = self.file_tag + "/" + blastm8.M8Stream.__next__(self)
        return line

# old method


def _multipleFileGenerator(
        m8files,
        filterParams,
        readFileDict,
        returnLines=True):
    """
    Iterate over all lines in a given list of streams as a single stream
    of hits, but keep track of which read came from which file in the given
    dictionary

    Note to self: this causes us to parse Hits twice if the resulting lines
    are parsed. It may speed things up to be able to return the Hit objects
    """
    for m8file in m8files:
        # allow for list of (file,tag) tuples
        if isinstance(m8file, tuple):
            m8file, fileTag = m8file
        else:
            fileTag = None

        # get file stream
        m8stream = blastm8.M8Stream(m8file)
        if fileTag is None:
            # let M8Stream class work out the file name for tagging
            fileTag = m8stream.fileName

        for hit in blastm8.getHitStream(m8stream, filterParams):
            # build map from reads to files/tags
            readFileDict[hit.read] = fileTag
            if returnLines:
                yield hit.line
            else:
                yield hit
