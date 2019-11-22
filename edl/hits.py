from edl.util import parseMapFile
from edl.taxon import getNodeFromHit, \
                      getAncestorClosestToRank, \
                      readTaxonomy, \
                      add_taxonomy_dir_argument
from edl.blastm8 import filterM8Stream, \
                        FilterParams,  \
                        formatsWithNoDescription, \
                        add_hit_table_arguments
from edl.expressions import accessionRE, nrOrgRE, koRE, giRE, pfamRE
import logging
logger = logging.getLogger(__name__)

#############
# Constants #
#############
ACCS = 'accs'
ORGS = 'orgs'
KEGG = 'kegg'
PFAM = 'pfam'
GIS = 'gis'
HITID = 'hitid'
HITDESC = 'hitdesc'
parsingREs = {
    ORGS: nrOrgRE,
    ACCS: accessionRE,
    KEGG: koRE,
    GIS: giRE,
    PFAM: pfamRE}
ALLEQ = 'all'
FIRST = 'first'
PORTION = 'portion'


def translateHits(hitMap, hitTranslation):
    for (read, hit) in hitMap.items():
        if isinstance(hit, type([])):
            newHits = []
            for h in hit:
                t = hitTranslation.get(h, None)
                if t is not None:
                    if isinstance(t, type([])):
                        newHits.extend(t)
                    else:
                        newHits.append(t)
                else:
                    newHits.append(h)
            hitMap[read] = list(set(newHits))
        else:
            t = hitTranslation.get(hit, None)
            if t is not None:
                hitMap[read] = t


def translateCounts(counts, translation):
    for key in counts.keys():
        newKey = translation.get(key, None)
        if newKey is not None:
            count = counts.pop(key)
            counts[newKey] = counts.setdefault(newKey, 0) + count


def binHits(hitMap):
    """
    return map of assignments to list of reads
    """
    hits = {}
    for (read, hit) in hitMap.items():
        if isinstance(hit, list):
            for h in hit:
                hits.setdefault(h, []).append(read)
        else:
            hits.setdefault(hit, []).append(read)

    return hits


def binAndMapHits(hitIter):
    """
    return map of assignments to list of reads
    """
    hits = {}
    hitMap = {}
    for (read, hit) in hitIter:
        hitMap[read] = hit
        if isinstance(hit, list):
            for h in hit:
                hits.setdefault(h, []).append(read)
        else:
            hits.setdefault(hit, []).append(read)

    return (hits, hitMap)


def loadSequenceWeights(weightFiles):
    """
    Load and merge list of sequence weight maps.
    """
    if len(weightFiles) > 0:
        sequenceWeights = {}
        for weightFile in weightFiles:
            sequenceWeights.update(parseMapFile(weightFiles, valueType=int))
    else:
        sequenceWeights = None
    return sequenceWeights


def add_weight_arguments(parser, multiple=False):
    action = 'store'
    default = None
    helpText = "File listing counting weights by sequence id. This is \
used for clustered or assembled data where each read (or contig) could \
represent any number of raw reads. The file should be a simple two-column \
tab-separated table with sequence-ids in the first column and integer \
weights in the second. "
    if multiple:
        action = 'append'
        default = []
        helpText += "For multiple files, supply the flag (-w or \
--sequenceWeights) for each file name. Concatenating all tables into \
one file will have the same net result."

    parser.add_argument("-w", "--sequenceWeights", dest='weights',
                        action=action, default=default, help=helpText)


def add_count_arguments(parser, defaults={}):
    default = defaults.get('cutoff', 0.01)
    parser.add_argument(
        "-c",
        "--cutoff",
        dest="cutoff",
        type=float,
        default=default,
        help="Cutoff for showing taxa. If a fractional count for a taxa "
             "is below this value, it will be folded up into its parent "
             "domain. Defaults to: %s" % default,
        metavar="CUTOFF")
    default = defaults.get('allMethod', ALLEQ)
    parser.add_argument(
        "-a",
        "--allMethod",
        dest="allMethod",
        default=default,
        choices=(
            FIRST,
            ALLEQ,
            PORTION),
        help="%r means +1 for every hit found for each read. %r means"
             " +1 to the first hit for each read. %r means +1/(nhits) for all"
             " hits of each read. Defaults to %r" % (ALLEQ,
                                                     FIRST,
                                                     PORTION,
                                                     default))


def getAllMethod(allMethod):
    return allMethods[allMethod]


def applyFractionalCutoff(counts, threshold=None, cutoff=None, label='Other'):
    """
    For any value in the dict below cutoff, remove and add to 'other' value
    """
    if threshold is None:
        if cutoff is None:
            logger.warn("Nothing to do for applyFractionalCutoff")
            return
        threshold = float(cutoff) * sum(counts.values())

    osum = 0
    for key in list(counts.keys()):
        if key == label:
            continue

        count = counts[key]
        if count < threshold:
            osum += count
            del counts[key]

    counts[label] = osum + counts.get(label, 0)
    return counts


def countIterHits(hitIter, allMethod=ALLEQ, weights=None, returnMap=True):
    """
    bin counts by hit and find total
    return map from assignments to number of reads
    and dict of original mappings
    """
    countHitsForRead = getAllMethod(allMethod)
    total = 0
    counts = {}
    if returnMap:
        hitMap = {}
    multiplier = 1
    for (read, hit) in hitIter:
        total += 1
        if returnMap:
            hitMap[read] = hit
        if weights is not None:
            multiplier = weights.get(read, 1)

        if isinstance(hit, type([])):
            countHitsForRead(hit, counts, multiplier=multiplier)
        else:
            counts[hit] = multiplier + counts.get(hit, 0)
    if returnMap:
        return (total, counts, hitMap)
    return (total, counts)


def _oneCountPerHit(hits, counts, multiplier=1):
    for hit in hits:
        counts[hit] = multiplier + counts.get(hit, 0)


def _portionHitCount(hits, counts, multiplier=1):
    multiplier = multiplier / float(len(hits))
    _oneCountPerHit(hits, counts, multiplier=multiplier)


def _countFirstHit(hits, counts, multiplier=1):
    counts[hits[0]] = multiplier + counts.get(hits[0], 0)


def countHits(hitMap):
    """
    bin counts by hit and find total
    return map from assignments to number of reads
    """
    total = 0
    counts = {}
    if isinstance(hitMap, dict):
        hitIter = hitMap.items()
    else:
        hitIter = hitMap
    for (read, hit) in hitIter:
        total += 1
        if isinstance(hit, type([])):
            for h in hit:
                counts[h] = 1 + counts.get(h, 0)
        else:
            counts[hit] = 1 + counts.get(hit, 0)
    return (total, counts)


def parseAndFilterM8Stream(inhandle, options):
    """
    runs the input stream through m8 filtering
    and then through parseM8Hits to get map from each read to all hits
    """
    inhandle = filterM8Stream(inhandle, options, return_lines=False)
    logger.info("Parsing hits")

    # since filter already parses hits, use that info
    infoInDescription = options.parseStyle in [KEGG, ORGS, PFAM]
    return parseM8Hits(inhandle, infoInDescription)


def parseM8File(inhandle,
                hitStringMap,
                options,
                parsingStyle,
                countMethod,
                taxonomy=None,
                rank=None,
                ignoreEmptyHits=True,
               ):
    """
    Wrapper method that combines filterM8, parseHits, and process hits to:
        filter hits using format and scorePct
        map reads to hits using parseHits
        translate hits using processHits

    If taxonomy is not None, hits will be TaxNode objects
    contMethod can only be LCA if taxonomy given

    Return a dict from read to hits
    """

    hitIter = parseM8FileIter(inhandle,
                              hitStringMap,
                              options,
                              parsingStyle,
                              countMethod,
                              taxonomy=taxonomy,
                              rank=rank,
                              ignoreEmptyHits=ignoreEmptyHits,
                             )

    hitMap = {}
    for (read, hits) in hitIter:
        hitMap[read] = hits

    logger.info("Done counting %d hits" % (len(hitMap)))

    return hitMap


def parseM8FileIter(inhandle,
                    hitStringMap,
                    options,
                    parsingStyle,
                    countMethod,
                    taxonomy=None,
                    rank=None,
                    ignoreEmptyHits=True,
                   ):
    """
    Wrapper method that combines filterM8, parseHits, and process hits to:
        filter hits using format and scorePct
        map reads to hits using parseHits
        translate hits using processHits

    If taxonomy is not None, hits will be TaxNode objects
    contMethod can only be LCA if taxonomy given

    Return an iterator over (read,hits) tuples.
    """

    # get map from reads to lists of hit strings
    logger.info("Parsing hits")
    # filters and parses
    # options.parseStyle = parsingStyle
    hitIter = filterM8Stream(inhandle, options, return_lines=False)

    # apply org or acc translation
    # apply map of hit names if given'
    # look up taxon node
    hitIter = processHits(
        hitIter,
        hitStringMap=hitStringMap,
        parseStyle=parsingStyle,
        taxonomy=taxonomy,
        rank=rank)

    # apply count method
    hitIter = applyCountMethod(hitIter, countMethod, ignoreEmptyHits)

    return hitIter


def parseHitsIter(
        hitIter,
        hitStringMap,
        parsingStyle,
        countMethod,
        taxonomy=None,
        rank=None,
        ignoreEmptyHits=None):
    """
    Same as parseM8FileIter, but takes in an iterator over Hit objects

    Simply runs processHits and applyCountMethod
    """
    # apply org or acc translation
    # apply map of hit names if given'
    # look up taxon node
    hitIter = processHits(
        hitIter,
        hitStringMap=hitStringMap,
        parseStyle=parsingStyle,
        taxonomy=taxonomy,
        rank=rank)

    # debugKey="F4UZ9WW02HMBZJ"
    # logger.debug("Hits for %s: %r" % (debugKey,hitMap[debugKey]))

    # apply count method
    hitIter = applyCountMethod(hitIter, countMethod, ignoreEmptyHits)

    return hitIter


def sortedHitIterator(hitMap):
    """
    Given a dictionary of reads to hits, return in order
    """
    for read in sorted(hitMap.keys()):
        yield (read, hitMap[read])


def applyCountMethod(hitIter, method, ignoreEmpty=True):
    # chose function that applies method
    if method == 'LCA' or method == 'rLCA':
        getBestHit = _findLeastCommonAncestor
    elif method == 'first':
        getBestHit = _takeFirstHit
    elif method == 'all':
        getBestHit = _returnAllHits
    elif method == 'consensus':
        getBestHit = _returnConsensus
    elif method == 'most':
        getBestHit = _returnMostCommon
    if ignoreEmpty:
        removeEmptyFunc = _removeEmpty
    else:
        removeEmptyFunc = _return_value

    # apply method to hit map
    hitsIn = 0
    hitsOut = 0
    reads = 0
    for (read, hits) in hitIter:
        reads += 1
        hitsIn += len(hits)
        hits = getBestHit(hits)
        hits = removeEmptyFunc(hits)
        if hits is not None:
            hitsOut += len(hits)
            yield (read, hits)
        logger.debug("%s=>%r" % (read, hits))

    logger.info(
        "Collected %d hits into %d hits for %d reads" %
        (hitsIn, hitsOut, reads))


def _findLeastCommonAncestor(hits):
    """
    Given a list of hits as TaxNode objects, find the least common ancestor.
    Hits that are not TaxNodes are ignored.
    """

    # check for hits not translated to TaxNode objects
    i = 0
    while i < len(hits):
        if hits[i] is None:
            hits.pop(i)
        elif isinstance(hits[i], type("")):
            logger.info(
                "Skipping hit: %s (cannot translate to taxon)" %
                (hits.pop(i)))
        else:
            i += 1

    # make sure there are some hits to process
    if len(hits) == 0:
        # sys.exit("No hits given!")
        return None

    # get LCA for these hits
    hit = hits[0]
    for i in range(1, len(hits)):
        hit = hit.getLCA(hits[i])

    return [hit, ]


def _returnMostCommon(hits):
    counts = {}
    for hit in hits:
        count = counts.get(hit, 0)
        count += 1
        counts[hit] = count

    logger.debug(repr(counts))

    bestCount = 0
    bestHit = None
    for (hit, count) in counts.items():
        if count > bestCount:
            bestHit = [hit, ]
            bestCount = count
        elif count == bestCount:
            bestHit.append(hit)

    return bestHit


def _takeFirstHit(hits):
    if len(hits) > 0:
        return hits[0:1]
    else:
        logger.debug("No hits!")
        return None


def _returnAllHits(hits):
    return list(set(hits))


def _returnConsensus(hits):
    hits = _returnAllHits(hits)
    if len(hits) == 1:
        return hits
    else:
        return None


def _return_value(value):
    return value


def _removeEmpty(hits):
    if hits is None:
        return hits

    while True:
        try:
            hits.remove(None)
        except ValueError:
            break

    while True:
        try:
            hits.remove('')
        except ValueError:
            break

    if len(hits) > 0:
        return hits
    else:
        return []


def parseHits(inhandle, readCol, hitCol, skipFirst, hitSep):
    """
    read over lines and pull out (read,[hits]) pairs given:
        inhandle: iterable set of strings (ie lines in a file)
        readCol: index of column with read name
        hitCol: index of column with hit name (-1 => every non-read column)
        skipFirst: skip first line if True
        hitSep: if not None, split data in hit column with this separator
    """
    logger.debug("BEGIN parseHits(in, %r, %r, %r, %r)" %
                 (readCol, hitCol, skipFirst, hitSep))

    # get line parsing function
    if hitSep == 'eval':
        extractReadHits = _getReadHitsEval
    else:
        hitCol = int(hitCol)
        if hitCol < 0:
            extractReadHits = _getReadHitsAll
        elif hitSep is not None:
            extractReadHits = _getReadHitsSep
        else:
            extractReadHits = _getReadHitsSimple

    if skipFirst:
        next(inhandle)

    hitCount = 0
    lineCount = 0
    lastRead = None
    for line in inhandle:
        lineCount += 1
        cells = line.rstrip('\n\r').split('\t')
        (read, hits) = extractReadHits(cells, readCol, hitCol, hitSep)
        if read != lastRead:
            if lastRead is not None:
                yield (lastRead, readHits)
            readHits = list(hits)
            lastRead = read
        else:
            readHits.extend(hits)
        hitCount += len(hits)
    if lastRead is not None:
        yield (lastRead, readHits)

    logger.info("Read %d hits from %d lines" % (hitCount, lineCount))


def parseM8Hits(hitIter, returnHitDescriptions):
    logger.debug("BEGIN parseM8Hits()")

    lastRead = None
    hitCount = 0
    readCount = 0
    for read, hits in hitIter:
        readCount += 1
        fields = []
        for hit in hits:
            hitCount += 1
            if returnHitDescriptions:
                fields.append(hit.hitDesc)
            else:
                fields.append(hit.hit)
        yield (read, fields)

    logger.info("Read %d hits from %d reads" % (hitCount, readCount))

# -- helpers for parseHits -- #
# the following functions take a line from a table and return a read name
# and an iterable collection of hits


def _getReadHitsEval(cells, readCol, hitCol, hitSep):
    """
    use eval to evaluate contents of hit cell. If resulting object is
    not iterable, put it into a tuple
    """
    read = cells[readCol]
    hit = cells[hitCol]

    # try to evaluate expression
    try:
        hit = eval(hit)
    except Exception:
        logger.warn("exception from 'eval(%r)'" % (hit))

    # make it iterable if it's not
    try:
        getattr(hit, '__iter__')
    except AttributeError:
        hit = (hit,)

    return (read, hit)


def _getReadHitsAll(cells, readCol, hitCol, hitSep):
    """
    every entry in cells (other than read) is a hit
    """
    read = cells.pop(readCol)
    return(read, cells)


def _getReadHitsSep(cells, readCol, hitCol, hitSep):
    """
    use hitSep to divide hit cell in to multipl hits
    """
    read = cells[readCol]
    hitCell = cells[hitCol]
    hits = hitCell.strip().split(hitSep)
    return (read, hits)


def _getReadHitsSimple(cells, readCol, hitCol, hitSep):
    read = cells[readCol]
    hit = cells[hitCol]
    return (read, (hit,))

# -- end helpers for parseHits -- #


class HitTranslator:
    """
    Given a list of (function,data,returnType) tuples ("mappings")
    Return an object with translateHit method that will apply the
    mappings to a hit
    """

    def __init__(self, mappings, useDesc=False, hitsAreObjects=True):
        self.mappings = mappings
        if mappings is None or len(mappings) == 0:
            self.applyMappings = self.returnSame
        if hitsAreObjects:
            if useDesc:
                self.getId = self.getDescription
        else:
            self.getId = self.returnSame

    def getId(self, hit):
        return hit.hit

    def getDescription(self, hit):
        return hit.hitDesc

    def returnSame(self, hit):
        return hit

    def translateHit(self, hit):
        return self.applyMappings([self.getId(hit), ])

    def applyMappings(self, hits):
        for (mapFunc, mapping, retType) in self.mappings:
            newHits = []
            for hit in hits:
                mapped = mapFunc(hit, mapping)
                if retType is list:
                    newHits.extend(mapped)
                elif retType is str:
                    newHits.append(mapped)
                else:
                    if isinstance(mapped, list) or isinstance(mapped, tuple):
                        newHits.extend(mapped)
                    else:
                        newHits.append(mapped)
            hits = newHits
        return hits


def getHitTranslator(
        hitStringMap=None,
        parseStyle=ORGS,
        taxonomy=None,
        rank=None,
        defaultToNone=True,
        hitsAreObjects=True):
    """
    Return a function that will return a list of organsims from a single hit.
        hitStringMap (None): dictionary mapping hit IDs to something else
        parseStyle (ORGS): how to process hit data into an identifying string
        taxonomy (None): An edl.taxon.Taxonomy object or directory
                          conatining taxdmp
        rank (None): Maximum rank to resolve hits
        hitsAreObjects: True if hits are edl.blastm8.Hit objects, else strings
    """
    parseRE = parsingREs.get(parseStyle, None)
    if logger.getEffectiveLevel() <= logging.INFO:
        if hitStringMap is None:
            mapstr = 'None'
        else:
            mapstr = '%d keys' % (len(hitStringMap))

        if parseRE is None:
            exprstr = 'None'
        else:
            exprstr = parseRE.pattern

        if taxonomy is None:
            taxstr = 'None'
        else:
            taxstr = '%d ids' % (len(taxonomy.idMap))

        logger.info(
            "Creating hit translator:\n default to None: %r\n map: %s\n "
            "parsing %s: %s\n taxa: %s\n rank: %s" %
            (defaultToNone, mapstr, parseStyle, exprstr, taxstr, rank))

    # set up variables
    infoInDescription = parseStyle in [KEGG, ORGS, PFAM]
    mappings = []
    if defaultToNone:
        mapFunction = _simpleMapNoneFunction
    else:
        mapFunction = _simpleMapFunction

    # initial parsing of hit id or description via regular expression
    if parseRE is not None:
        mappings.append((_findAllREfunctionSimpler, parseRE, list))

    # optional look up table
    if hitStringMap is not None:
        mappings.append((mapFunction, hitStringMap, None))

    # optional conversion to Taxon objects
    if taxonomy is not None:
        if parseStyle == ORGS:
            if defaultToNone:
                mappings.append((getNodeFromHit, taxonomy.nameMap, str))
            else:
                mappings.append((_getNodeHitFunction, taxonomy.nameMap, str))
        else:
            mappings.append((mapFunction, taxonomy.idMap, str))
        if rank is not None:
            mappings.append((getAncestorClosestToRank, rank, str))

    return HitTranslator(
        mappings,
        useDesc=infoInDescription,
        hitsAreObjects=hitsAreObjects)

# turn hit lines into organisms or KOs or anything else


def processHits(hitIter, **kwargs):
    """
    Take an in iterator over read,hits tuples and apply mappings using
    a HitTranslator
    """
    translator = getHitTranslator(**kwargs)

    # translate hits
    for (key, hits) in hitIter:
        logger.debug("%s => %s" % (key, hits))
        newHits = []
        for h in hits:
            newHits.extend(translator.translateHit(h))
        logger.debug(str(newHits))
        yield (key, newHits)


def processHitsOld(
        hitIter,
        mapping=None,
        expr=None,
        taxIdMap=None,
        taxNameMap=None,
        defaultToNone=True,
        rank=None):
    """
    Take a map of reads (or other keys) to lists of hits and translate hits.

    Can use the following steps in this order with any steps omitted:
        simpile dictionary translation using 'mapping'
        regular expression (where every captured group is returned as a hit)
        a translation to taxNode objects by one of:
            simple dictionary translation using taxIdMap
            name based look up using edl.taxon.getNodeFromHit() and taxNameMap

    if defaultToNone is changed to False, anything not found in one of
        the mappings
      (mapping, taxIdMap, or taxNameMap)
    """

    if logger.getEffectiveLevel() <= logging.DEBUG:
        if mapping is None:
            mapstr = 'None'
        else:
            mapstr = '%d keys' % (len(mapping))

        if expr is None:
            exprstr = 'None'
        else:
            exprstr = expr.pattern

        if taxIdMap is None:
            if taxNameMap is None:
                taxstr = 'None'
            else:
                taxstr = '%d names' % (len(taxNameMap))
        else:
            taxstr = '%d ids' % (len(taxIdMap))

        logger.debug(
            "Starting processHits:\n default to None: %r\n map: %s\n "
            "exp: %s\n taxa: %s\n rank: %s" %
            (defaultToNone, mapstr, exprstr, taxstr, rank))

    # set the functions to use:
    if mapping is None:
        mapFunction = _passFunction
    elif defaultToNone:
        mapFunction = _simpleMapNoneFunction
    else:
        mapFunction = _simpleMapFunction

    exprFunction = _findAllREfunction

    if taxIdMap is not None:
        taxMap = taxIdMap
        if defaultToNone:
            taxFunction = _simpleMapNoneFunction
        else:
            taxFunction = _simpleMapFunction
    elif taxNameMap is not None:
        taxMap = taxNameMap
        if defaultToNone:
            taxFunction = getNodeFromHit
        else:
            taxFunction = _getNodeHitFunction
    else:
        taxMap = None
        taxFunction = _passFunction

    if taxMap is None or rank is None:
        rankFunction = _passFunction
    else:
        rankFunction = getAncestorClosestToRank

    # translate hits
    for (key, hits) in hitIter:
        logger.debug("%s => %s" % (key, hits))
        newHits = []
        for h in hits:
            # find all matches to expr, may be more than one
            hs = exprFunction(h, expr)
            logger.debug("%s => %s" % (h, hs))
            for hit in hs:
                hts = mapFunction(hit, mapping)
                if not (isinstance(hts, list) or isinstance(hts, tuple)):
                    hts = [hts]
                for hit in hts:
                    hit = taxFunction(hit, taxMap)
                    hit = rankFunction(hit, rank)
                    newHits.append(hit)
        logger.debug(str(newHits))
        yield (key, newHits)

# helper functions for processHits
# each function takes a hit and something else, and then reutrns a
# translated hit


def _passFunction(hit, mapping):
    return hit


def _simpleMapFunction(hit, mapping):
    newHit = mapping.get(hit, hit)
    logger.debug("%s --> %r" % (hit, newHit))
    return newHit


def _simpleMapNoneFunction(hit, mapping):
    newHit = mapping.get(hit, None)
    logger.debug("%s --> %r" % (hit, newHit))
    return newHit


def _getNodeHitFunction(hit, taxMap):
    newHit = getNodeFromHit(hit, taxMap)
    if newHit is None:
        return hit
    else:
        return newHit


def _findAllREfunctionSimpler(hit, expr):
    hits = expr.findall(hit)
    if len(hits) == 0:
        return [hit, ]
    else:
        return hits


def _findAllREfunction(hit, expr):
    if expr is None:
        return (hit,)
    hits = expr.findall(hit)
    if len(hits) == 0:
        return [hit, ]
    else:
        return hits

# end helper functions for processHits


def add_taxon_arguments(parser, defaults={}, choices={}):
    # get format and filter_top_pct options from blastm8
    add_hit_table_arguments(parser, defaults,
                            flags=['format', 'filter_top_pct'])

    # specific to taxon parsing:
    parser.add_argument(
        "-m",
        "--mapFile",
        dest="mapFile",
        default=defaults.get(
            "mapFile",
            None),
        metavar="MAPFILE",
        help="Location of file containing table of with db hit name "
             "as first column and taxa or taxonids in second column. "
             "Defaults to '%s'" % (defaults.get("mapFile", None)))
    parser.add_argument(
        "-p",
        "--parseStyle",
        default=defaults.get(
            "parseStyle",
            ACCS),
        choices=[
            ACCS,
            GIS,
            ORGS,
            HITID,
            HITDESC],
        help="What should be parsed from the hit table: accessions('accs'), "
             "'gis', organsim names in brackets ('orgs'), the full hit "
             "name('hitid'), or the full hit description('hitdesc'). "
             "(defaults to '%s')" % (defaults.get("parseStyles", ACCS)))
    parser.add_argument(
        "-C",
        "--countMethod",
        dest="countMethod",
        default=defaults.get(
            "countMethod",
            "first"),
        choices=choices.get(
            'countMethod',
            ('first',
             'most',
             'all',
             'LCA',
             'consensus')),
        help="How to deal with counts from multiple hits. (first, most: "
             "can return multiple hits in case of a tie, LCA: MEGAN-like, "
             "all: return every hit, consensus: return None unless all "
             "the same). Default is %s" % (defaults.get("countMethod",
                                                        "first")),
        metavar="COUNTMETHOD")
    add_taxonomy_dir_argument(parser, defaults)


def readMaps(options, namesMap=False):
    """
    Load the taxonomy and id->to->taxid maps requested by user
    """
    return (readTaxonomyFiles(options, namesMap=namesMap), readIDMap(options))


def readTaxonomyFiles(options, namesMap=False):
    """
    load the taxonomy specififed by the user. Create a name lookup map if
    parseStyle is 'orgs'
    """
    # read taxonomy
    if options.taxdir is not None:
        getTaxNames = namesMap or options.parseStyle == ORGS
        taxonomy = readTaxonomy(options.taxdir, namesMap=getTaxNames)
        logging.info("Read %d nodes from tax dump" % (len(taxonomy.idMap)))
    else:
        taxonomy = None
        if options.countMethod == 'LCA' or options.countMethod == 'rLCA':
            raise Exception('Cannot use LCA without providng a taxonomy (-n)')
        logging.info("No taxonomy needed")

    return taxonomy


def readIDMap(options):
    """
    Load the specififed lookup table for hit IDs. If the parseStyle
    requested is 'gis', convert keys to integers. The values are always
    convereted to integeres since they are assumed to be taxids
    """
    # map reads to hits
    if options.parseStyle == GIS:
        keyType = int
    else:
        keyType = None
        if options.taxdir is not None:
            valueType = int
        else:
            valueType = None
    return parseMapFile(options.mapFile, valueType=valueType, keyType=keyType)


allMethods = {ALLEQ: _oneCountPerHit,
              FIRST: _countFirstHit,
              PORTION: _portionHitCount}

############
# Tests
############


def test():
    import sys
    global myAssertEq, myAssertIs
    from test import myAssertEq, myAssertIs

    if len(sys.argv) > 2:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.warn
    logging.basicConfig(stream=sys.stderr, level=loglevel)
    logger.setLevel(loglevel)

    hits = testParseHits(sys.argv[1])
    testTranslateAndCountHits(hits)


def testParseHits(testFile):
    # test line parsing methods
    cells = [1, 2, 3, 4, "(4,5)", "6,7"]
    (read, hitIter) = _getReadHitsSimple(cells, 0, 2, None)
    hits = []
    for h in hitIter:
        hits.append(h)

    myAssertEq(read, 1)
    myAssertEq(len(hits), 1)
    myAssertEq(hits[0], 3)

    (read, hitIter) = _getReadHitsSep(cells, 1, 5, ',')
    hits = []
    for h in hitIter:
        hits.append(h)
    myAssertEq(read, 2)
    myAssertEq(hits, ['6', '7'])

    (read, hitIter) = _getReadHitsAll(list(cells), 3, -1, None)
    hits = []
    for h in hitIter:
        hits.append(h)
    myAssertEq(read, 4)
    myAssertEq(len(hits), 5)
    myAssertEq(hits, [1, 2, 3, "(4,5)", "6,7"])

    # give it a test file
    hitIter = parseHits(open(testFile), 0, -1, True, None)
    hits = {}
    for r, h in hitIter:
        hits[r] = h
    logging.debug(repr(hits))
    myAssertEq(len(hits), 29)
    myAssertEq(hits['000023_2435_2174'], ['Prochlorococcus'])
    myAssertEq(hits['000178_2410_1152'], ['Bacteria <prokaryote>'])
    myAssertEq(hits['000093_2435_2228'], ['Candidatus Pelagibacter'])

    return hits


def testTranslateAndCountHits(hits):
    (total, counts) = countHits(hits)
    myAssertEq(total, 29)
    myAssertEq(counts["Prochlorococcus"], 10)
    myAssertEq(counts['root'], 7)

    translateHits(hits,
                  {'Bacteria <prokaryote>': 'other',
                   'root': 'other',
                   'Candidatus Pelagibacter': 'Pelagibacter'})
    myAssertEq(hits['000178_2410_1152'], ['other'])
    myAssertEq(hits['000093_2435_2228'], ['Pelagibacter'])


if __name__ == '__main__':
    test()
