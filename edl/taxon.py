import re
import logging
import numpy as np

logger = logging.getLogger(__name__)

##############
# Classes    #
##############


class Taxonomy:
    """
    A container for taxonomy data: contains two maps: id to Node and name
    to Node
    """

    def __init__(self, idMap, nameMap, realNameMap, path=None, rootNode=None):
        self.idMap = idMap
        self.nameMap = nameMap
        self.realNameMap = realNameMap
        self.path = path
        if rootNode is None:
            rootNode = next(iter(idMap.values())).getRootNode()
        self.root = rootNode

    def __str__(self):
        return "Taxonomy with %d nodes" % (len(self.idMap))

    def __repr__(self):
        return "Taxonomy(%s)" % (self.path)


class TaxNode:
    """
    A node in a pylogenetic tree. This has a parent and many children
    """
    domains = ["Bacteria", "Archaea", "Eukaryota", "Viruses", "Viroids"]
    namedNodes = {}

    def __init__(self, taxid, parentid, rank):
        """
        """
        self.id = taxid
        self.taxid = taxid
        self.parentid = parentid
        self.rank = rank
        self.parent = None
        self.children = []
        self.name = ""
        self.translation = None
        self.lineage = None
        self.lineage_strings = {}

    def __hash__(self):
        return hash(self.getLineageString(';'))

    def __repr__(self):
        return "TaxNode(%s,%s,%s)" % (
            repr(self.id), repr(self.parentid), repr(self.rank))

    def __key__(self):
        return self.getLineageString(';')

    def __lt__(self, other):
        return self.__key__() < other.__key__()

    def __eq__(self, other):
        return self.__key__() == other.__key__() if isinstance(
            other, self.__class__) else False

    def __str__(self):
        if self.name == "":
            return str(self.taxid)
        else:
            if (self.isNameGeneric()):
                if self is not self.parent:
                    return "%s(%s)" % (self.name, str(self.parent))
        return self.name

    def setParent(self, parent):
        self.parent = parent
        if parent is not None:
            parent.children.append(self)

    def isAncestorOf(self, node):
        """
        Return true if node contained within this node
        """
        return node is not None and self in node.getLineage()

    def getLCA(self, node):
        """
        Given another node in the same tree, find the most recent (lowest)
        common ancestor.

        Node will cache their lineages for faster retrieval
        """
        if self in node.getLineage():
            if logger.getEffectiveLevel() <= logging.DEBUG:
                lineageString = ""
                for n in node.getLineage():
                    lineageString += str(n) + ','
                logger.debug("%s found in [%s]" % (str(self), lineageString))
            return self
        else:
            if self.parent is None or self.parent is self:
                return self
            return self.parent.getLCA(node)

    def transmogrify(self, rank, taxList):
        """
        Given a rank string, and a list of taxa, return a single
        taxon name by the following rules:

        1) if taxon or an ancestor matches a name in taxList, return
            matching name.
        2) if taxon or an ancestor has the indicated rank, return
            matching taxon name.
        3) if taxon or a parent matches list of domains, return the domain
        4) print a warning and return None
        """

        if self.translation is None:
            if self.name in taxList:
                self.translation = self.name
                logger.debug("%s in tax list" % (self.name))
                return self.translation

            if self.parent is None or self.parent is self:
                self.translation = self.name
                logger.debug("map to self %s" % (self.name))
                return self.translation

            parentTranslation = self.parent.transmogrify(rank, taxList)

            if self.rank == rank and parentTranslation not in taxList:
                self.translation = self.name
                logger.debug("%s is rank %s" % (self.name, rank))
                return self.translation

            self.translation = parentTranslation
            logger.debug(
                "%s is using parent's translation: %s" %
                (self.name, self.translation))
        else:
            logger.debug(
                "%s already translated to %s" %
                (self.name, self.translation))
        return self.translation

    def getAncestorClosestToRank(self, rank, **kwargs):
        return getAncestorClosestToRank(self, rank, **kwargs)

    def getAncestorAtRank(self, rank):
        """
        return the TaxNode ancestor of this node which has the given rank.
        Return None if no suitable ancestor found.
        """
        if self.rank == rank:
            return self
        else:
            if self.parent is None:
                return None
            if self.parent is self:
                return None
            return self.parent.getAncestorAtRank(rank)

    def getRootNode(self):
        """
        go up through parents til we get to the root
        """
        logger.info("getting root node through %s" % (self.name))
        if self.parent is not None and self.parent is not self:
            return self.parent.getRootNode()
        else:
            return self

    def isNameGeneric(self):
        name = spaceRE.sub(',', self.name)
        if name[0:10] == 'uncultured':
            return True
        if name[0:13] == 'environmental':
            return True
        if metagenomeRE.search(name):
            return True
        if name[0:12] == 'endosymbiont':
            return True

        return False

    def getCollapsedCounts(self, counts, cutoff, hitTranslations):
        """
        walk the tree and fill hitTranslation hash with map from nodes
        under cutoff to node where counts should be aggregated
        """
        count = 0
        countedNodes = [self]
        if self in counts:
            # use value in counts if it's there
            count = counts[self]
            logger.debug("%s has %d hits" % (self, count))

        # add sum of children (those under cutoff)
        for child in self.children:
            if child is self:
                logger.warn(
                    "node %s(%s) is child of itself!" %
                    (repr(self), str(self)))
                continue
            (kidCount, kidsCountedNodes) = child.getCollapsedCounts(
                counts, cutoff, hitTranslations)
            if kidCount is not None:
                logging.debug(
                    "Adding %d to %s from %s" %
                    (kidCount, child, self))
                count += kidCount
                countedNodes.extend(kidsCountedNodes)

        # if this node has a generic name, defer to parent
        if self.isNameGeneric():
            logger.debug(
                "%s is too generic, giving %d hits to parent: %s" %
                (self.name, count, self.parent.name))
            return (count, countedNodes)

        # if this node is over cutoff, add to counts
        if count >= cutoff:
            logger.info(
                "keeping %d hits in %s (from %d nodes)" %
                (count, str(self), len(countedNodes)))
            name = self.name
            for node in countedNodes:
                hitTranslations[node] = self
            return (None, [])
        else:
            # otherwise return count for parent to use
            logger.debug("Passing %d hits to parent from %s" % (count, self))
            return (count, countedNodes)

    def getLineageString(self, sep):
        if sep not in self.lineage_strings:
            if self.parent is None or self.parent is self:
                self.lineage_strings[sep] = self.name
            else:
                self.lineage_strings[sep] = \
                    sep.join((self.parent.getLineageString(sep),
                              self.name))
        return self.lineage_strings[sep]

    def getLineage(self):
        """
        return a list of all the nodes in this node's ancestry
        """
        if self.lineage is None:
            if self.parent is None or self.parent is self:
                self.lineage = tuple([self, ])
            else:
                lineage = list(self.parent.getLineage())
                lineage.append(self)
                self.lineage = tuple(lineage)
        return self.lineage

    def compareRanks(self, comparisons):
        for kid in self.children:
            if kid is self or kid is None:
                continue
            if self.rank is not None and self.rank.strip() != "no rank":
                kid.compareRank(self.rank, comparisons)
            kid.compareRanks(comparisons)

    def compareRank(self, ancestorRank, comparisons):
        if self.rank is not None and self.rank.strip() != "no rank":
            compKey = (ancestorRank, self.rank)
            comparisons[compKey] = comparisons.get(compKey, 0) + 1

        for kid in self.children:
            if kid is self or kid is None:
                continue
            kid.compareRank(ancestorRank, comparisons)

    def generateMemberTaxids(node):
        for child in node.children:
            for taxid in generateMemberTaxids(child):
                yield taxid
        yield node.id

    @staticmethod
    def getNamedNode(name):
        if name not in TaxNode.namedNodes:
            node = TaxNode(name, None, None)
            TaxNode.namedNodes[name] = node
        return TaxNode.namedNodes[name]

    @staticmethod
    def addToTreeFromString(taxString, tree={}):
        if 'root' not in tree:
            if len(tree) > 0:
                raise Error('tree must have root node!')
            root = TaxNode('root', None, None)
            root.name = root.id
            tree['root'] = root

        lineage = scRE.split(taxString)
        logger.debug("parsing %s: %s" % (taxString, str(lineage)))
        lastNode = tree['root']
        for taxon in lineage:
            taxon = taxon.strip()
            taxon = removeSpaces(taxon)
            if (taxon in tree) and (tree[taxon].parent is lastNode):
                lastNode = tree[taxon]
            else:
                newNode = TaxNode(taxon, lastNode.id, None)
                newNode.name = newNode.id
                newNode.setParent(lastNode)
                tree[taxon] = newNode
                lastNode = newNode
        return lastNode

################
# compiled REs #
################
cladeRE = re.compile(r'clade')
parensRE = re.compile(r'\([^\(\)]+\)')
lastSemicolonRE = re.compile(r'^.*;([^;]+)$')
spaceRE = re.compile("\s")
dotRE = re.compile("\.")
scRE = re.compile(r';+')
metagenomeRE = re.compile(r'metagenome')

#############
# Functions #
#############

# this is a list (in order) of the ranks in the ncbi tax dump
ranks = [
    'forma',
    'varietas',
    'subspecies',
    'species',
    'species subgroup',
    'species group',
    'subgenus',
    'genus',
    'subtribe',
    'tribe',
    'subfamily',
    'family',
    'superfamily',
    'parvorder',
    'infraorder',
    'suborder',
    'order',
    'superorder',
    'infraclass',
    'subclass',
    'class',
    'superclass',
    'subphylum',
    'phylum',
    'superphylum',
    'subkingdom',
    'kingdom',
    'superkingdom']

# The next few things were an attempt to automatically determine the order
# of the ranks from the taxonomic tree.
"""
sortKey={}
def getSortKey(rank):
    return sortKey.get(rank,len(rank))
"""
comparisons = {}


def compareRanks(r1, r2):
    r1anc = comparisons.get((r1, r2), 0)
    r2anc = comparisons.get((r2, r1), 0)

    if r1anc > 0 and r2anc > 0:
        logger.warn(
            "ambiguos relationshp between %s and %s: (%d,%d)" %
            (r1, r2, r1anc, r2anc))

    if r1anc == 0 and r2anc == 0:
        logger.warn(
            "no information for %s and %s: (%d, %d)" %
            (r1, r2, r1anc, r2anc))

    return cmp(r1anc, r2anc)


def deduceRankOrder(taxMap):
    import numpy as np

    logger.info("figuring out ranks!")
    root = next(iter(taxMap.values())).getRootNode()

    # this generates a map of tuples to counts
    # comparisons[(ranka,rankb)] == 4  means ranka was an ancestor to rankb 4
    # times
    global comparisons
    comparisons = {}
    root.compareRanks(comparisons)

    logger.info(
        "There are %d entries in the comparison map!" %
        (len(comparisons)))

    # get list of all ranks
    ranks = []
    for key in comparisons:
        ranks.extend(key)
    ranks = set(ranks)
    logger.info("%d ranks: %s" % (len(ranks), ranks))

    """
    # some testing
    global sortKey
    sortKey={}
    aCounts={}
    dCounts={}
    for key in comparisons:
        v = comparisons[key]
        aCounts[key[0]]=aCounts.get(key[0],0) + v
        dCounts[key[1]]=dCounts.get(key[1],0) + v

        # just checking:
        if (key[1],key[0]) in comparisons:
            logger.warn("%s is an ancestor of %s %d times and "
                        "vice versa %d times!" %
                 (key[0],key[1],v,comparisons[(key[1],key[0])]))

    for rank in ranks:
        if rank not in aCounts:
            if rank not in dCounts:
                logger.warn("rank %s is not counted!" % (rank))
                continue
            sortKey[rank]=600.0
        elif rank not in dCounts:
            sortKey[rank]=-600.0

        else:
            sortKey[rank]=np.arctan(float(aCounts[rank])/dCounts[rank])

        logger.info("key for %s is %s (%d/%d)" % (rank,str(sortKey[rank]),
                    aCounts.get(rank,1),dCounts.get(rank,0)))
    """

    # sort ranks
    ranks = sorted(ranks, cmp=compareRanks)
    logger.info("sorted ranks: %s" % str(ranks))

    return ranks

_taxonomies = {}


def readTaxonomy(taxDir, namesMap=False):
    """
    read the names.dmp and nodes.dmp files in this directory and build a tree
    return a map from taxon id to each node in the tree
    """

    if taxDir in _taxonomies:
        # if this taxonomy has already been parse, just re-use it
        if not namesMap or len(_taxonomies[taxDir].nameMap) > 0:
            return _taxonomies[taxDir]

    logger.info("read nodes")

    taxMap = {}
    nameMap = {}
    realNameMap = {}
    # build tree from nodes.dmp
    for line in open("%s/nodes.dmp" % (taxDir)):
        cells = line.split(r'|')
        taxid = int(cells[0].strip())
        parentid = int(cells[1].strip())
        rank = cells[2].strip()
        node = TaxNode(taxid, parentid, rank)
        taxMap[taxid] = node

    logger.info("link nodes")

    # link parents and children
    for node in taxMap.values():
        parent = None
        try:
            parent = taxMap[node.parentid]
        except:
            logger.warn(
                "Can't find parent (%s) for %s" %
                (node.parentid, node.id))
        else:
            node.setParent(parent)

    logger.info("name nodes")

    # name nodes from names.dmp
    for line in open("%s/names.dmp" % (taxDir)):
        cells = line.split(r'|')
        taxid = int(cells[0].strip())
        name = cells[2].strip()
        name2 = cells[1].strip()
        if name == "":
            name = name2
            name2 = None

        quality = cells[3].strip()
        if quality == "scientific name":
            node = taxMap[taxid]
            node.name = name
            if namesMap:
                realNameMap[simplifyString(name)] = node
        elif namesMap:
            node = taxMap[taxid]

        if namesMap:
            if name2 is None or name2 == name:
                names = [name, ]
            else:
                names = [name, name2]

            for name in names:
                name = simplifyString(name)
                mapnode = nameMap.get(name, 0)
                if mapnode == 0:
                    # not in map, add it
                    nameMap[name] = node
                elif mapnode is not None:
                    # already in map
                    if mapnode is not node:
                        # already in map with different taxon
                        lca = node.getLCA(mapnode)
                        nameMap[name] = lca

    taxonomy = Taxonomy(taxMap, nameMap, realNameMap)
    _taxonomies[taxDir] = taxonomy
    return taxonomy


def simplifyString(string):
    return dotRE.sub("", removeSpaces(string.lower()))


def removeSpaces(string):
    return spaceRE.sub("", string)

nameTranslations = {'asaia lannensis': 'asaia lannaensis',
                    'uncultured haptophyte': 'haptophyta',
                    'acidisoma sibiricum': 'acidisoma sibirica'}
methodCount = {
    'sp': 0,
    'par': 0,
    'map': 0,
    'raw': 0,
    'none': 0,
    'pre': 0,
    'sub': 0}


def getNodeFromHit(hit, nameMap, exhaustive=True):
    """
    Use a number of tricks to map the organism name given by 'hit' to
    a taxon object in 'nameMap'
    """
    if hit is None:
        return hit

    # remove extra formatting
    hit = simplifyString(hit)

    #  remove an initial '|'
    if hit[0:1] == '|':
        hit = hit[1:]

    """
    First, try a simple look up:
        nameMap contains all the synonyms in the NCBI tax dump
        the synonyms and hit are lowercase have had all spaces remove
    """
    try:
        taxNode = nameMap[hit]
        if logger.getEffectiveLevel() <= logging.DEBUG:
            methodCount['raw'] += 1
        return taxNode
    except KeyError:
        pass

    # try to remove parens
    m = parensRE.search(hit)
    if m:
        hit = parensRE.sub('', hit)
        try:
            taxNode = nameMap[hit]
            if logger.getEffectiveLevel() <= logging.DEBUG:
                methodCount['par'] += 1
            return taxNode
        except KeyError:
            pass

    # hard coded translations
    if hit in nameTranslations:
        hit = nameTranslations[hit]
        try:
            taxNode = nameMap[hit]
            if logger.getEffectiveLevel() <= logging.DEBUG:
                methodCount['map'] += 1
            return taxNode
        except KeyError:
            pass

    # replace clade with cluster
    (newHit, count) = cladeRE.subn('cluster', hit)
    if count > 0:
        try:
            taxNode = nameMap[newHit]
            return taxNode
        except KeyError:
            pass

    if exhaustive:
        startingName = None
        # look for name that starts with this complete hit
        # or name that is found is start of hit
        hitLen = len(hit)
        for name in nameMap:
            nameLen = len(name)
            if hitLen < nameLen:
                # except for cases like 'alteromonas sp.', ...
                if hit[-2:] != 'sp':
                    # check to see if hit is substring of name
                    if name[0:hitLen] == hit:
                        logger.debug("%s changed to %s" % (hit, name))
                        if logger.getEffectiveLevel() <= logging.DEBUG:
                            methodCount['pre'] += 1
                        return nameMap[name]
            else:
                if hit[0:nameLen] == name:
                    if startingName is None:
                        startingName = name
                    else:
                        if nameLen > len(startingName):
                            startingName = name
        else:
            # there was no 'pre' match, so take longest 'sub' match
            if startingName is not None:
                logger.debug('%s changed to %s' % (hit, startingName))
                if logger.getEffectiveLevel() <= logging.DEBUG:
                    methodCount['sub'] += 1
                return nameMap[startingName]

    logger.warn("Can't translate name: %s" % (hit))
    if logger.getEffectiveLevel() <= logging.DEBUG:
        methodCount['none'] += 1
    return None


def getAncestorClosestToRank(node, rank, **kwargs):
    """
    This is an attempt to get something close to the requested rank even when
    the organism has no acnestral taxon with that exact rank

    The named ranks on either side of our target are found
    (eg prehaps kingdom and class if phylum is missing)
    Then based on how many unranked items are in the lineage between these AND
     the number of ranks skipped, interpolate which ancestral taxon is
     closest to
     the target rank
    """
    # Set the fall back (aka default) to starting node unless set by caller
    default = kwargs.pop('default', node)

    # Set behavior in special case:
    #  if the first ranked ancestor is beyond the target rank
    #  then use the child of that ancestor (ie the previous on in the lineage)
    #  if this is set to False, just return the 'default'
    useChildOfFirstRankedAncestor = kwargs.pop(
        'useChildOfFirstRankedAncestor', True)

    # This only works on TaxNode objects
    if not isinstance(node, TaxNode):
        logger.debug("Not a TaxNode (%r)" % (node))
        return default

    # Walk through the lineage and find ranked taxa as reference points
    # Get all the ancestors/parents of this Taxon
    lineage = list(node.getLineage())
    lineage.reverse()
    if rank == 'domain':
        rank = 'superkingdom'
    targetIndex = ranks.index(rank)
    logger.debug("looking for rank: %s (%d)" % (rank, targetIndex))
    lastIndex = -1
    lastIndexedAnc = None
    lastAnc = node
    # For each ancestor/parent in this TaxNode's lineage
    for anc in lineage:
        try:
            # If it has a rank, where is that rank in the heirarchy?
            ancRankIndex = ranks.index(anc.rank)
            logger.debug(
                "rank of %s is %s(%d)" %
                (anc, anc.rank, ancRankIndex))
        except ValueError:
            # hard coded special cases
            # family:
            # SAR11,SAR116,SAR324,SAR86,SAR406,SAR202,SUP05,SAR92,OMZ60,
            if anc.id in [
                    54526,
                    62654,
                    131190,
                    62672,
                    62680,
                    648176,
                    655184,
                    745004,
                    744996]:
                ancRankIndex = ranks.index('family')
            else:
                ancRankIndex = -1

        if ancRankIndex >= 0:
            # An exact match is easy, just return it
            if ancRankIndex is targetIndex:
                logger.debug("MATCH!")
                return anc
            elif ancRankIndex > targetIndex:

                # if we've hit the next rank without hitting any other ranks:
                if lastIndexedAnc is None:
                    # We've hit the a lower rank without hitting any other
                    # if this is an ancestor of the node itself and
                    # and it caller requests it, return previous ancestor
                    if anc is not node and useChildOfFirstRankedAncestor:
                        logger.debug("Take previous!")
                        return lastAnc

                    # Otherwise return the default
                    if useChildOfFirstRankedAncestor:
                        logger.debug("Node is already ranked. Won't do it.")
                    return default

                # We've hit the next rank and there was a previous rank:
                # try to interpolate
                ancIndex = lineage.index(anc)
                lastAncIndex = lineage.index(lastIndexedAnc)
                logger.debug(
                    "Trying to interpolate between %d and %d based "
                    "on %d between %d and %d" % (ancIndex, lastAncIndex,
                                                 targetIndex, ancRankIndex,
                                                 lastIndex))
                ranksSkipped = ancRankIndex - lastIndex
                ancsSkipped = ancIndex - lastAncIndex
                rankAdjustment = ancRankIndex - targetIndex
                ancAdjustment = (
                    float(rankAdjustment) / ranksSkipped) * ancsSkipped
                logger.debug("rolling back by %s" % (str(ancAdjustment)))
                return lineage[int(np.floor(ancIndex - ancAdjustment))]

        if ancRankIndex > -1:
            lastIndex = ancRankIndex
            lastIndexedAnc = anc
        lastAnc = anc

    logger.debug("Nothing found close to %s" % rank)
    return default

############
# Tests
############


def test():
    import sys
    global myAssertEq, myAssertIs
    from test import myAssertEq, myAssertIs

    ndir = sys.argv[1]
    if len(sys.argv) > 2:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.WARN
    logging.basicConfig(stream=sys.stderr, level=loglevel)

    test_root_node()
    test_get_lineage()
    test_collapse_counts()
    ncbiTree = test_read_ncbi(ndir)
    test_get_ancestor(ncbiTree)
    test_transmogrify(ncbiTree)


def test_transmogrify(tree):
    orgLists = []
    orgLists.append(['Bacteria <prokaryote>', 'Archaea',
                     'Prochlorales', 'Rickettsiales', 'Eukaryota'])
    orgLists.append(['Gammaproteobacteria',
                     'Alphaproteobacteria',
                     'Deltaproteobacteria'])
    ranks = ['phylum', 'superkingdom', 'genus']
    ids = [439493, 939841, 655186, 1046240, 333146]
    answers = {439493:
               {'phylum':
                {str(orgLists[0]): 'Rickettsiales',
                 str(orgLists[1]): 'Alphaproteobacteria'},
                'superkingdom':
                {str(orgLists[0]): 'Rickettsiales',
                 str(orgLists[1]): 'Alphaproteobacteria'},
                'genus':
                {str(orgLists[0]): 'Rickettsiales',
                 str(orgLists[1]): 'Alphaproteobacteria'}},
               939841:
               {'phylum':
                {str(orgLists[0]): 'Prochlorales',
                 str(orgLists[1]): 'Cyanobacteria'},
                'superkingdom':
                {str(orgLists[0]): 'Prochlorales',
                 str(orgLists[1]): 'Bacteria <prokaryote>'},
                'genus':
                {str(orgLists[0]): 'Prochlorales',
                 str(orgLists[1]): 'Prochlorococcus'}},
               655186:
               {'phylum':
                {str(orgLists[0]): 'Bacteria <prokaryote>',
                 str(orgLists[1]): 'Gammaproteobacteria'},
                'superkingdom':
                {str(orgLists[0]): 'Bacteria <prokaryote>',
                 str(orgLists[1]): 'Gammaproteobacteria'},
                'genus':
                {str(orgLists[0]): 'Bacteria <prokaryote>',
                 str(orgLists[1]): 'Gammaproteobacteria'}},
               1046240:
               {'phylum':
                {str(orgLists[0]): 'Bacteria <prokaryote>',
                 str(orgLists[1]): 'Gammaproteobacteria'},
                'superkingdom':
                {str(orgLists[0]): 'Bacteria <prokaryote>',
                 str(orgLists[1]): 'Gammaproteobacteria'},
                'genus':
                {str(orgLists[0]): 'Bacteria <prokaryote>',
                 str(orgLists[1]): 'Gammaproteobacteria'}},
               333146:
               {'phylum':
                {str(orgLists[0]): 'Archaea',
                 str(orgLists[1]): 'Euryarchaeota'},
                'superkingdom':
                {str(orgLists[0]): 'Archaea',
                 str(orgLists[1]): 'Archaea'},
                'genus':
                {str(orgLists[0]): 'Archaea',
                 str(orgLists[1]): 'Ferroplasma'}}
               }

    for orgs in orgLists:
        for rank in ranks:
            for node in tree.values():
                node.translation = None

            for taxid in ids:
                node = tree[taxid]
                try:
                    assert node.transmogrify(rank, orgs) == answers[
                        taxid][rank][str(orgs)]
                except AssertionError:
                    logging.warn(
                        "%d:%s at rank %s%s goes to: %s, not %s" %
                        (taxid, node.name, rank, str(orgs), node.transmogrify(
                            rank, orgs), answers[taxid][rank][
                            str(orgs)]))


def test_get_ancestor(tree):
    answers = ((112233, 'order', 30483), (654321, 'order', 4892),
               (1032926, 'phylum', 35493))
    try:
        for data in answers:
            assert tree[data[0]].getAncestorAtRank(data[1]).id == data[2]
    except AssertionError:
        for data in answers:
            node = tree[data[0]]
            ance = node.getAncestorAtRank(data[1])
            "Ancestor of %d:%s at rank %s is %d:%s, not %d:%s(%s)" % (
                data[0], node.name, data[1], ance.id, ance.name,
                data[2], tree[data[2]].name)


def test_collapse_counts():
    from hits import countHits, translateHits
    tree = {}
    node1 = TaxNode.addToTreeFromString(
        'Bacteria;Cyanobacteria;Prochlorococcus', tree)
    node2 = TaxNode.addToTreeFromString(
        'Bacteria;Gammaproteobacteria;SUP05', tree)
    node3 = TaxNode.addToTreeFromString(
        'Archaea;Thermoplasmata;Ferroplasma', tree)

    hits = {
        '1': node1,
        '2': node1,
        '3': node2,
        '4': node3,
        '5': node3,
        '6': node3,
        '7': node3,
        '8': node3,
        '9': node3}
    (total, counts) = countHits(hits)
    collapsedHits = {}
    node1.getRootNode().getCollapsedCounts(counts, .3 * 9, collapsedHits)
    translateHits(hits, collapsedHits)
    (ctotal, collapsedCounts) = countHits(hits)
    expectedCounts = {node3: 6, node1.parent.parent: 3}
    try:
        assert collapsedCounts == expectedCounts
    except AssertionError:
        print("total: %d, ctotal: %d" % (total, ctotal))
        print(str(counts))
        print(str(collapsedHits))
        print(str(collapsedCounts))
        print(str(expectedCounts))
        raise AssertionError


def test_root_node():
    tree = {}
    node1 = TaxNode.addToTreeFromString(
        'Bacteria;Cyanobacteria;Prochlorococcus', tree)
    node2 = TaxNode.addToTreeFromString(
        'Bacteria;Gammaproteobacteria;SUP05', tree)
    node3 = TaxNode.addToTreeFromString(
        'Archaea;Thermoplasmata;Ferroplasma', tree)

    root = node1.getRootNode()
    assert root.name == 'root'
    assert root.parent is None or root.parent is root


def test_get_lineage():
    tree = {}
    lineage = 'Bacteria;Cyanobacteria;Prochlorococcus'
    node1 = TaxNode.addToTreeFromString(lineage, tree)
    try:
        assert node1.getLineageString(';') == 'root;' + lineage
    except AssertionError:
        logging.warn(
            "%s is not %s" %
            (node1.getLineageString(';'), 'root;' + lineage))
        raise AssertionError
    lineage = 'Bacteria;Gammaproteobacteria;SUP05'
    node2 = TaxNode.addToTreeFromString(lineage, tree)
    try:
        assert node2.getLineageString(';') == 'root;' + lineage
    except AssertionError:
        logging.warn(
            "%s is not %s" %
            (node2.getLineageString(';'), 'root;' + lineage))
        raise AssertionError


def test_read_ncbi(ndir):
    taxNames = True
    taxonomy = readTaxonomy(ndir, taxNames)
    taxIds = taxonomy.idMap
    taxNames = taxonomy.nameMap
    myAssertEq(len(taxIds), 783145)
    myAssertEq(len(taxNames), 1101991)

    # pick some random things to check
    myAssertEq(taxIds[123456].name, 'Psammomoya choretroides')
    myAssertIs(
        taxNames[
            simplifyString('Psammomoya choretroides')],
        taxIds[123456])
    myAssertIs(taxNames[simplifyString('Psammomoya choretroides '
                                       '(F.Muell.) Diels & Loes.')],
               taxIds[123456])
    myAssertEq(taxIds[123499].parent.id, 50537)

    return taxIds


def add_taxonomy_dir_argument(parser, defaults={}):
    parser.add_argument(
        "-n",
        "--ncbiTaxDir",
        dest="taxdir",
        metavar="PATH",
        default=defaults.get(
            "taxdir",
            None),
        help="Directory with unpacked ncbi tax dump (specifically names.dmp "
             "and nodes.dmp) and use to translate org names in desc, "
             "otherwise try to find lineage info in desc. Default is: %s" %
             (defaults.get("taxdir", None)))


if __name__ == '__main__':
    test()
