#! /usr/bin/python
"""
Library of methods and regular expressions for parsing files from the
KEGG gene ontology
"""

import sys
import re
import logging
import os
logger = logging.getLogger(__name__)

##############
# Classes    #
##############
################
# compiled REs #
################
koRE = re.compile(r'\b(K\d{5})\b')
kokoRE = re.compile(r'^ENTRY\s+(K\d{5})\b')
endSectionRE = re.compile(r'^\S')
definitionRE = re.compile(r'^DEFINITION\s+(\S.*\S)\s*$')
classRE = re.compile(r'^CLASS\s+(\S.*)$')
ecRE = re.compile(r'^DEFINITION.*\[(EC:[-0-9\.]+)\]')
genesRE = re.compile(r'^GENES\s+(\S.*\S)\s*$')
trailingBracketRE = re.compile(r'\s*\[[^\[\]]+\]\s*$')
kegkoRE = re.compile(r'^D\s+.+\b(K\d{5})(?:<\/a>)?\s*(.*)$')
britekoRE = re.compile(r'^[A-Z]\s+(K\d{5})\s*(.*)$')
geneListRE = re.compile(r'(?<=\s)([a-zA-Z0-9_.-]+)\b')
orgRE = re.compile('^([A-Z]{3,4}):')
cogGroupRE = re.compile(r'(.+)\[(.+)\]')
cogMapRE = re.compile(r'^\[(\S+)\]\s+(\S+)\s+(\S.+)$')

#############
# Functions #
#############


def readSEEDTree(treeFile):
    """
    Return nested dictionary where first dict is map from levels (1,2,3)
    and next dict is map from role to name.

    This is a simply formatted file with 4 columns:
    "role\tsubsystem\tlevel 2\t level 1"
    """
    seedTree = {'1': {}, '2': {}, '3': {}}
    with open(treeFile) as f:
        for line in f:
            (role, l3, l2, l1) = line.rstrip().split('\t')
            seedTree['1'][role] = l1
            seedTree['3'][role] = l3
            seedTree['2'][role] = l2
    return seedTree


def readCogTree(mapFile):
    """
    return maps from CDD id to COGID, COG description, and COG category
    """
    cogMap = {'gene': {}, 'description': {}, 'group': {}}
    with open(mapFile) as f:
        for line in f:
            cdd, cog, gene, description, count = line.rstrip().split('\t')
            description, group = cogGroupRE.match(description).groups()
            cogMap['gene'][cdd] = cog
            cogMap['description'][cdd] = description
            groups = [re.sub(' +', ' ', g.strip()) for g in group.split("/")]
            cogMap['group'][cdd] = groups

    # hack to make things work with the previous methods (ie SEED)
    cogMap['3'] = cogMap['group']
    return cogMap


def readCogTreeFromWhog(mapFile):
    """
    Return a map from COG  id to category
    """
    cogMap = {'gene': {}, 'group': {}}
    with open(mapFile) as f:
        for line in f:
            line.rstrip()
            m = cogMapRE.maatch(line)
            if m:
                category = m.group(1)
                cog = m.group(2)
                description = m.group(3)
                cogMap['gene'][cog] = description
                cogMap['group'][cog] = category
    return cogMap


def parseSeedMap(mapFile):
    """
    Return a dictionary mapping from accession to subsystem.

    The SEED map (refseq2md52role.gz) starts with two summary lines followed
    by three columns (accession\thash/sum\tSubsystem):

    Mapped roles:9004(subsys.txt 2 subsystems2peg)
    Unmapped roles:521(s                ubsys.txt)
    YP_9218461b41910965945b806d5defc49ad1a224CO dehydrogenases      \
            maturation factor, CoxF family
    YP_0012863261e472ed51c0df8feb03ee296a0  e55de4CO dehydrogenases \
            maturation factor, CoxF family
    """
    accMap = {}
    with open(mapFile) as f:
        f.next()
        f.next()
        for line in f:
            # logger.debug(line)
            (acc, code, subsys) = line.rstrip('\r\n').split('\t', 2)
            # logger.debug("Mapped %s to %s (sum: %s)" % (acc,subsys,code))
            accMap.setdefault(acc.strip(), []).append(subsys.strip())

    return accMap


def _stripKeggKeyPrefix(key):
    return key.split(":", 1)[1]


def parseLinkFile(mapFile, stripKeys=False, stripVals=True):
    """
    Parse the gene_ko.list file from KEGG

    hsa:10001       ko:K15128
    hsa:10002       ko:K08546
    hsa:10003       ko:K01301

    with the possibility of duplicate records
    """
    if mapFile is None:
        return None

    logger.info("parsing map file: %s" % (mapFile))
    translation = {}
    rows = 0
    badRows = 0
    lastKey = None
    for line in open(mapFile):
        cells = line.split('\t')
        if len(cells) > 1:
            rows += 1
            key = cells[0].strip()
            if stripKeys:
                key = _stripKeggKeyPrefix(key)
            value = cells[1].strip()
            if stripVals:
                # strip 'ko:' from start of each value
                value = _stripKeggKeyPrefix(value)
            if key == lastKey:
                translation[key].append(value)
            else:
                translation[key] = [value, ]
        else:
            badRows += 1
    if badRows > 0:
        logger.warn("%d rows in map file had too few columns!" % (badRows))

    logger.info(
        "Read %d records from %d lines of %s" %
        (len(translation), rows, mapFile))
    return translation


def parseModuleMap(mapFile):
    """
    Parse module file to dict
    """
    return parseLinkFile(mapFile, stripKeys=True, stripVals=False)


def parseGeneKOMap(koFile):
    """
    scan ko file an build map from gene names to kos
    """
    koMap = {}
    koCount = 0
    inGenes = False
    logger.info("Reading kos and gene names from %s" % koFile)
    for line in open(koFile):
        # find KO first
        # looking for a KO line
        match = kokoRE.match(line)
        if match:
            ko = match.group(1)
            logger.debug("Start of %s" % (ko))
            koCount += 1
            if logger.getEffectiveLevel() >= logging.DEBUG:
                if koCount % 1000 == 0:
                    logging.debug("Parsed %d KOs" % koCount)
            continue

        # look for information on this KO
        if not inGenes:
            # looking for a GENE line
            match = genesRE.match(line)
            if match:
                inGenes = True
                geneString = match.group(1)
                logger.debug("found genes: %s" % (geneString))
                _mapGenes(koMap, ko, geneString)
            continue

        # everyline is a gene line until further notice
        elif not endSectionRE.match(line):
            # not the end of gene section: reading more genes
            geneString = line.strip()
            logger.debug("found genes: %s" % (geneString))
            _mapGenes(koMap, ko, geneString)

        else:
            # found all gene strings
            logger.debug("End of genes")
            inGenes = False
            ko = None

    logger.info("Mapped %d genes to %d kos" % (len(koMap), koCount))
    return koMap


def _mapGenes(koMap, ko, geneString):
    """
    geneString looks like "ORG: geneid(genename) geneid(genename) geneid geneid
    while ids in KeggGenes look like: "org:geneid"
    """
    org = orgRE.match(geneString).group(1).lower()
    genes = geneListRE.findall(geneString)
    for gene in genes:
        kGene = "%s:%s" % (org, gene)
        koMap.setdefault(kGene, []).append(ko)


def readKEGGFile(kFile, keggLevel):
    if len(kFile) > 4 and kFile[-4:] == ".keg":
        return readKeggFile(kFile, keggLevel)
    else:
        return readKOFile(kFile, keggLevel)


def readKOFileLevels(kofile, keggLevel):
    raise Exception(
        "Cannot parse KEGG levels form ko file. Sorry. Use the brite "
        "hierarchy file: brite/ko/ko00001.keg")


def readKOFile(kofile, keggLevel):
    """
    Scan ko file and build map from ko to names at the given level.
    The level can be one of:
        NAME, PATHWAY, EC, DEFINITION,
          or a level in the CLASS heirachy: 1, 2, or 3

    The returned dictionary maps KO strings to lists of names.
    """

    # based on kegg level, define what a 'name' line looks like
    if keggLevel in ['1', '2', '3', 1, 2, 3]:
        # return readKOFileLevels(kofile, keggLevel)
        logger.warn(
            "KEGG leveles can not be parsed from ko files starting "
            "around 2013. This may fail unless you are using an older file")
        reString = 'CLASS'
    else:
        reString = keggLevel

    koMap = {}
    reString = keggLevel

    nameRE = re.compile("^%s\s+(\S.*\S)\s*$" % (reString))
    nameSep = None
    nameIndex = None
    if keggLevel == 'NAME':
        nameSep = ','
    elif reString == 'CLASS':
        nameSep = ';'
        nameIndex = int(keggLevel) - 1
    elif keggLevel == 'EC':
        nameRE = ecRE

    logger.info(
        "Looking for %s with %s and (%s,%s)" %
        (repr(keggLevel),
         repr(
            nameRE.pattern),
            repr(nameSep),
            repr(nameIndex)))

    ko = None
    inName = False
    names = []
    for line in open(kofile):
        # find KO first
        # looking for a KO line
        match = kokoRE.match(line)
        if match:
            ko = match.group(1)
            logger.debug("Start of %s" % (ko))

        # look for information on this KO
        elif not inName:
            # looking for a name line
            match = nameRE.match(line)
            if match:
                inName = True
                nameString = match.group(1)
                # contents of name line processed differently based on sep and
                # index
                names.extend(_parseName(nameString, nameSep, nameIndex))
                logger.debug("found names: %s" % (nameString))
            continue

        elif not endSectionRE.match(line):
            # not the end of name section: reading more names
            nameString = line.strip()
            names.extend(_parseName(nameString, nameSep, nameIndex))
            logger.debug("found names: %s" % (nameString))

        else:
            # found all name strings
            logger.debug("End of names: %s" % (str(names)))
            koMap[ko.strip()] = names
            names = []
            inName = False
            ko = None

    if len(koMap) == 0:
        if keggLevel in ['1', '2', '3', 1, 2, 3]:
            raise Exception(
                "Cannot parse KEGG levels from newer ko files. Sorry. "
                "Use the brite heirarchy file: brite/ko/ko00001.keg")
        else:
            raise Exception(
                "No kegg levels were parsed! The newer ko files have a "
                "modified format. You may need to try the brite heirachy "
                "file: brite/ko/ko00001.keg")
    return koMap


def _parseName(string, sep, index):
    """
    parses string from name line based on separator char and index:
        if no sep given, return whole string
        if index given, return that element from split (using sep)
        otherwise, return all elements from split
    """
    if sep is None:
        return [_removeTrailingBrackets(string).strip(), ]
    else:
        bits = string.split(sep)
        if index is None:
            bits[-1] = _removeTrailingBrackets(bits[-1]).strip()
            return bits
        else:
            return [_removeTrailingBrackets(bits[index]).strip()]


def _removeTrailingBrackets(string):
    """
    strips bracketed strings from end
    """
    return trailingBracketRE.sub('', string)

##
# readKegFile(keggFile, keggLevel)
# map KOs to pathway given .keg file from KEGG
##


def readKeggFile(keggFile, keggLevel):
    logger.info("Parsing level %s from %s" % (keggLevel, keggFile))
    kmap = {}
    if keggLevel == 'PATHWAY':
        keggLevel = '3'
    kfile = open(keggFile)
    keggLevelMap = {
        '1': re.compile(r'^A\s*<[bB]>(.*)</[bB]>\s*$'),
        '2': re.compile(r'^B\s*<[bB]>(.*)</[bB]>\s*$'),
        '3': re.compile(r'^C\s*(\S.*)$')}
    pathRE = keggLevelMap.get(str(keggLevel), None)
    # We'll need to get ko lists or descriptsion from brite files
    briteRE = re.compile(r'^C\s*.+\[BR:(ko\d+)\]\s*$')
    detailsRE = None
    if keggLevel == 'NAME':
        detailsRE = re.compile(r'^(.*);')
    elif keggLevel == 'DEFINITION':
        detailsRE = re.compile(r';\s*(\S.*)$')
    elif keggLevel == 'EC':
        sys.exit(
            "Sorry, cannot currently parse EC from a .keg file, you'll "
            "need to use the 'ko' file")
        # detailsRE=re.compile(r'\[(EC:[0-9\.]+)\]\s*$')

    logger.debug("Looking for ko with: %s" % (kegkoRE.pattern))
    if pathRE is not None:
        logger.debug("Looking for path with: %s" % (pathRE.pattern))
    if detailsRE is not None:
        logger.debug("Looking for details with: %s" % (detailsRE.pattern))

    desc = ''
    if pathRE is not None:
        # if we're looking for a path in hierarchy...
        for line in kfile:
            # Is this a header in the level we cara about (1,2, or 3)
            m = pathRE.search(line)
            if m:
                desc = _removeTrailingBrackets(m.group(1))
                # Is this a BRITE place holder? if so, scan for KOs
                m = briteRE.search(line)
                if m:
                    processBriteFile(m.group(1), desc, keggFile, kmap)
                continue

            # Is this a BRITE place holder? if so, scan for KOs
            m = briteRE.search(line)
            if m:
                processBriteFile(m.group(1), desc, keggFile, kmap)
                continue

            # Is this a KO line, add to most recently seen level header
            m = kegkoRE.search(line)
            if m:
                ko = m.group(1)
                # make sure we have something to map to
                if desc == '':
                    die("Pathway info not found before ko: %s" % ko)
                try:
                    if desc not in kmap[ko]:
                        kmap[ko].append(desc)
                except KeyError:
                    kmap[ko] = [desc]

    else:
        # if we just want descriptions
        for line in kfile:
            # we still need to catch the BR:brite files for KO details
            m = briteRE.search(line)
            if m:
                getDescriptionsFromBriteFile(
                    m.group(1), keggFile, detailsRE, kmap)
                continue

            m = kegkoRE.search(line)
            if m:
                ko = m.group(1)
                details = m.group(2)
                if detailsRE is None:
                    kmap[ko] = [details]
                else:
                    m = detailsRE.search(details)
                    if m:
                        kmap[ko] = kmap.setdefault(ko, []).append(m.group(1))

    logger.info("Read %d KOs" % (len(kmap)))
    return kmap


def getDescriptionsFromBriteFile(pathway, file1, detailsRE, kmap):
    """
    Parse the given brite pathway for KO details
    """
    logging.debug("getting descriptions from BRITE:%s" % pathway)
    briteDir = os.path.split(file1)[0]
    briteFile = briteDir + os.path.sep + pathway + ".keg"
    try:
        with open(briteFile) as f:
            for line in f:
                m = britekoRE.search(line)
                if m:
                    ko = m.group(1)
                    details = m.group(2)
                    if detailsRE is None:
                        kmap[ko] = [details]
                    else:
                        m = detailsRE.search(details)
                        if m:
                            kmap[ko] = kmap.setdefault(
                                ko, []).append(m.group(1))
    except IOError:
        logging.warn(
            "Cannot parse brite pathway: %s from %s" %
            (pathway, briteFile))


def processBriteFile(pathway, desc, file1, kmap):
    logging.debug("getting KOs from BRITE:%s" % pathway)
    briteDir = os.path.split(file1)[0]
    briteFile = briteDir + os.path.sep + pathway + ".keg"
    koCount = 0
    try:
        with open(briteFile) as f:
            for line in f:
                m = britekoRE.search(line)
                if m:
                    ko = m.group(1)
                    kmap.setdefault(ko, []).append(desc)
                    koCount += 1
        logging.debug("Assigned %d KOs to pathway: %s" % (koCount, desc))
    except IOError:
        logging.warn(
            "Cannot parse brite pathway: %s from %s" %
            (pathway, briteFile))


def add_path_arguments(parser, defaults={}, choices={}, helps={}):
    # get format and filterPct arguments from blastm8
    from edl.hits import HITID, ACCS, GIS, KEGG, HITDESC, PFAM
    from edl.blastm8 import add_hit_table_arguments
    add_hit_table_arguments(parser, defaults, flags=['format', 
                                                     'filterPct',
                                                     'sort'
                                                    ])

    # specific to pathway parsing:
    pgroup = parser.add_argument_group(
        "Pathway Arguments",
        "These arguments control the mapping of hits to gene "
        "function heirarchies like KEGG or SEED""")
    pgroup.add_argument(
        "-m",
        "--mapFile",
        dest="mapFile",
        default=defaults.get(
            "mapFile",
            None),
        metavar="MAPFILE",
        help="Location of file containing table of with db hit name as "
             "first column and geneIDs (Knumber) in second column.")
    pgroup.add_argument(
        "-M",
        "--mapStyle",
        default='auto',
        choices=[
            'auto',
            'kegg',
            'tab',
            'seed'],
        help="What type of mapping file are you using: simple tab "
             "separated list of IDs and kos/subsystems/domains, the "
             "genes_ko.list file from KEGG (which adds ko: to the K "
             "numbers and can have multiple records for each gene id), "
             "or the 3 column file from SEED. By default, this script "
             "will inspect the file and guess, but you can force 'kegg', "
             "'seed' or 'tab' with this argument.")
    pgroup.add_argument(
        "-p",
        "--parseStyle",
        default=defaults.get(
            "parseStyle",
            HITID),
        choices=[
            ACCS,
            GIS,
            KEGG,
            HITID,
            HITDESC,
            PFAM],
        help="What should be parsed from the hit table: accessions('accs'), "
             "'gis', K numbers in description ('kegg'), the full hit "
             "name('hitid'), or the full hit description('hitdesc'). "
             "(defaults to '%s')" % (defaults.get("parseStyle",
                                                  HITID)))
    pgroup.add_argument(
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
             'consensus')),
        help=helps.get(
            "countMethod",
            "How to deal with counts from multiple hits. (first, most: "
            "can return multiple hits, all: return every hit, consensus: "
            "return None unless all the same). Do not use most or consensus "
            "with more than one level at a time. Default is %s" %
            (defaults.get(
                "countMethod",
                "first"))),
        metavar="COUNTMETHOD")
    if defaults.get("filterForPath", False):
        action = 'store_false'
        default = True
        helpstr = 'Consider all hits. By deafult, only hits with path \
assignments are used.'
    else:
        action = 'store_true'
        default = False
        helpstr = 'Ignore hits with no entry in pathway map (-m). By default \
all hits are used and if the best hit(s) is(are) to sequences with no path, \
then the read will not be assigned to a path'
    pgroup.add_argument(
        "-r",
        "--filterForPath",
        action=action,
        dest="mappedHitsOnly",
        default=default,
        help=helpstr)
    add_pathways_argument(pgroup, defaults)
    parser.add_argument_group(pgroup)


def add_pathways_argument(parser, defaults={}):
    parser.add_argument(
        "-T",
        "--heirarchyType",
        default=defaults.get(
            "heirarchyType",
            'kegg'),
        choices=[
            'kegg',
            'seed',
            'cazy',
            'cog',
            'kegg_module'],
        help="What kind of functional heirarchy to use. 'kegg', seed', "
             "or 'cazy'. Defaults to '%s'" % (defaults.get("heirarchyType",
                                                           'kegg')))
    parser.add_argument(
        "-H",
        "--heirarchyFile",
        metavar="HEIRARCHY_FILE",
        default=defaults.get(
            'heirarchyFile',
            None),
        help="File containing pathway/subsystem/genefaimly heirarchy "
             "(either ko or ko00001.keg for KEGG or susbsys.txt for SEED). "
             "Defaults to %s" % (defaults.get('heirarchyFile', None)))


############
# Tests
############
def test():
    import sys

    if len(sys.argv) > 3:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.WARN
    logging.basicConfig(stream=sys.stderr, level=loglevel)
    logger.setLevel(logging.INFO)

    kegg_nosetest(sys.argv[1], sys.argv[2])


def kegg_nosetest(ko_map, kegg_file):
    global myAssertEq, myAssertIs
    from edl.test import myAssertEq, myAssertIs
    testReadKeggFile(kegg_file)
    testParseGeneLink(ko_map)
    # testParseGeneKOMap(ko_map)
    # testReadKoFile(ko_map)


def testParseGeneLink(koFile):
    gkmap = parseLinkFile(koFile)
    myAssertEq(gkmap['ggo:101148121'], ['K16534'])
    myAssertEq(gkmap['olu:OSTLU_15108'], ['K11126'])
    myAssertEq(gkmap['ebt:EBL_c03070'], ['K02879'])
    myAssertEq(gkmap['pec:W5S_4205'], ['K00363'])
    myAssertEq(gkmap['buc:BU148'], ['K03101'])
    myAssertEq(gkmap['smaf:D781_0330'], ['K06925'])
    myAssertEq(gkmap['nkr:NKOR_05565'], ['K03524'])


def testReadKeggFile(keggFile):
    kDmap = readKeggFile(keggFile, 'DESCRIPTION')
    myAssertEq(kDmap['K09630'], ['PRSS36; protease, serine, 36 [EC:3.4.21.-]'])
    kPmap = readKeggFile(keggFile, 'PATHWAY')
    assert('K00397' in kPmap)
    myAssertEq(kPmap['K00399'], ['01200 Carbon metabolism',
                                 '00680 Methane metabolism', '01000 Enzymes'])
    k2map = readKeggFile(keggFile, 2)
    myAssertEq(k2map['K13810'][1].lower(), 'Carbohydrate Metabolism'.lower())
    myAssertEq(k2map['K13810'][0].lower(), 'Overview'.lower())
    myAssertEq(k2map['K00399'][1].lower(), 'Energy Metabolism'.lower())
    myAssertEq(k2map['K00399'][0].lower(), 'Overview'.lower())
    k3map = readKeggFile(keggFile, 3)
    myAssertEq(k3map['K13810'],
               ['01230 Biosynthesis of amino acids',
                '00010 Glycolysis / Gluconeogenesis',
                '00030 Pentose phosphate pathway',
                '00500 Starch and sucrose metabolism',
                '00520 Amino sugar and nucleotide sugar metabolism',
                '01000 Enzymes',
                '01000 Enzymes'])
    myAssertEq(k3map['K00399'], ['01200 Carbon metabolism',
                                 '00680 Methane metabolism', '01000 Enzymes'])
    myAssertEq(
        k3map['K03404'], [
            '00860 Porphyrin and chlorophyll metabolism', '01000 Enzymes'])
    myAssertEq(k3map['K01976'], ['01000 Enzymes'])
    myAssertEq(k3map['K07347'],
               ['02000 Transporters',
                '02044 Secretion system',
                '02035 Bacterial motility proteins',
                '05133 Pertussis'])
    myAssertEq(k3map['K09630'], ['01000 Enzymes', '01002 Peptidases'])
    k3mapQ = readKeggFile(keggFile, '3')
    for k in k3map.keys():
        try:
            myAssertEq(k3map[k], k3mapQ[k])
        except AsserionError:
            raise AssertionError(
                "level 3 classes for %s do not match:\n%s\n%s" %
                (k, k3map[k], k3mapQ[k]))


def testReadKoFile(koFile):
    kPmap = readKOFile(koFile, 'PATHWAY')
    assert('K00397' not in kPmap)
    myAssertEq(kPmap['K00399'],
               ['ko00680  Methane metabolism',
                'ko01200  Carbon metabolism'])

    kEmap = readKOFile(koFile, 'EC')
    myAssertEq(kEmap['K00397'], ['EC:1.8.99.-'])
    myAssertEq(kEmap['K00399'], ['EC:2.8.4.1'])

if __name__ == '__main__':
    test()
