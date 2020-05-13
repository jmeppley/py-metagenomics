#! /usr/bin/python
"""
Library of methods and regular expressions for parsing files from the
KEGG gene ontology
"""

import logging
import os
import re
import sys
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


def parse_KEGG_file(k_file, kegg_level):
    """ Checks filename and runs:
        parse_KO_file if base filename is ko
        parse_keg_file if file extension is .keg """
    if os.path.basename(k_file) == 'ko':
        return parse_ko_file(k_file, kegg_level)
    elif len(k_file) > 4 and k_file[-4:] == ".keg":
        return parse_keg_file(k_file, kegg_level)
    else:
        raise Exception("I don't know what to do with file: %s"
                        % (os.path.basename(k_file)))


def parse_ko_file(ko_file, level):
    results = {}

    # synonyms
    if level in ['PATHWAYS', 'PATH', 'PATHS']:
        level = 'PATHWAY'
    if level in ['DESCRIPTION', 'FUNCTION']:
        level = 'DEFINITION'
    if level == 'EC':
        level = 'ko01000:4'

    if re.match(r'(ko\d\d\d\d\d:)?(\d+)', str(level)):
        # these are in the BRITE heirachy
        #  eg: k00001:2 for level to of the main hierarchy
        brite_hier, brite_level = \
            re.match(r'(?:(ko\d\d\d\d\d):)?(\d+)', str(level)).groups()
        brite_level = int(brite_level)
        if brite_hier is None:
            # if its just a number, assume k00001
            brite_hier = 'ko00001'
        logger.debug(f"Looking for level {brite_level} in {brite_hier}")

    with open(ko_file) as ko_handle:

        # look for single line per entry
        if level in ['NAME', 'DEFINITION']:
            kw_expr = re.compile(r'^(ENTRY|{})(\s+)(\S.*)?'.format(level))
            try:
                for i, line in enumerate(ko_handle):
                    m = kw_expr.match(line)
                    if m:
                        keyword, spaces, value = m.groups()
                        if keyword == 'ENTRY':
                            ko = value.split()[0].strip()
                        elif keyword == level:
                            results[ko] = value.strip()
            except Exception as exc:
                print(f'Error on line {i}:\n{line}')
                raise exc

        # there can be multiple pathways after and including the PATHWAY line
        elif level == 'PATHWAY':
            kw_expr = re.compile(r'^(ENTRY|{})(\s+)(\S.*)?'.format(level))

            def skip(line, indent, pathways):
                return

            def add_pathway(line, indent, pathways):
                pathways.append(line[indent:-1])

            pathways, indent = None, 0
            for line in ko_handle:
                m = kw_expr.match(line)
                if m:
                    keyword, spaces, value = m.groups()
                    if keyword == 'ENTRY':
                        ko = value.split()[0].strip()
                        indent = 5 + len(spaces)
                        process_line = skip
                        continue
                    elif keyword == level:
                        process_line = add_pathway
                        pathways = results.setdefault(ko, [])
                    else:
                        process_line = skip
                        continue

                process_line(line, indent, pathways)

        else:
            # BRITE
            entry_rexp = re.compile(r'^ENTRY\s+(K\d+)')
            brite_rexp = \
                re.compile(r'^((?:BRITE)?\s+)(\S.+\S)\s*\[BR:(ko\d+)\]')
            end_brite_rexp = re.compile(r'^\S')
            level_rexp = re.compile(r'^(\s+)(\S.+)')

            lines = iter(enumerate(ko_handle))
            try:
                # outer loop looping over Entries
                while True:

                    # find next Entry line
                    for i, line in lines:
                        m = entry_rexp.match(line)
                        if m:
                            ko = m.group(1)
                            break
                    else:
                        # no more entries
                        break

                    # find start of BRITE
                    for i, line in lines:
                        m = brite_rexp.match(line)
                        if m:
                            spaces, name, hierarchy = m.groups()
                            if hierarchy == brite_hier:
                                brite_indent = len(spaces)
                                brite_levels = results.setdefault(ko, [])
                                break

                    # process BRITE lines
                    for i, line in lines:
                        if end_brite_rexp.match(line) or \
                                brite_rexp.match(line):
                            # start of next hierarchy or next keyword section
                            break

                        spaces, level_name = level_rexp.match(line).groups()
                        # level is number of spaces beyond original indent
                        if len(spaces) - brite_indent == brite_level:
                            brite_levels.append(level_name)

                    # end while outer loop
            except StopIteration:
                # I don't think we ever get here
                pass
            except Exception as exc:
                print(f"error on line {i}:\n{line}")
                print(f"found {len(results)} kos so far")
                raise exc

    return results


def parse_keg_file(keg_file, level):
    """ Parse KEGG metadata from brite .keg files
        level: one of

            * PATH, PATHWAY, or PATHWAYS
            * 1 - 6 or A - F
            * DEFINITION, DESCRIPTION, or FUNCITON

    """

    # synonyms
    if level in ['PATHWAYS', 'PATHWAY', 'PATHS']:
        level = 'PATH'
    if level in ['DESCRIPTION', 'FUNCTION']:
        level = 'DEFINITION'
    if str(level) in {'1', '2', '3', '4', '5', '6'}:
        level = 'ABCDEF'[int(level) - 1]

    ko_def_rexp = re.compile(r'^[B-F]\s+(K\d\d\d\d\d)\s+(\S.+\S)\s*$')
    level_rexp = re.compile(r'^([A-F])\s*(\S.+)')
    path_rexp = re.compile(r'\s*\[PATH:\s*ko\d+\s*\]')
    html_rexp = re.compile(r'</?[a-z]+/?>')

    results = {}
    with open(keg_file) as keg_handle:

        # two types of parsing
        if level == 'DEFINITION':
            # just looking for the line with the K# and description
            #  (ignore hierarchy)
            for line in keg_handle:
                m = ko_def_rexp.match(line)
                if m:
                    ko, desc = m.groups()
                    results[ko] = desc

        elif level in ['A', 'B', 'C', 'D', 'E', 'F', 'PATH']:
            # looking for level, and all KOs after that
            print(f"looking for {level}")
            level_name = None
            for line in keg_handle:

                # check for ko first, because it also looks like a level
                m = ko_def_rexp.match(line)
                if m:
                    if level_name is not None:
                        # if we've seen a level we like
                        # save ko level
                        ko = m.group(1)
                        results.setdefault(ko, []).append(level_name)
                        continue

                m = level_rexp.match(line)
                if m:
                    letter, name = m.groups()
                    if letter == level:
                        # found a header at the target level, remember name
                        level_name = html_rexp.sub('', name)
                    elif (level == 'PATH' and path_rexp.search(name)):
                        # found a header at the target level, remember name
                        level_name = path_rexp.sub('', name)
                    elif letter < level:
                        # we've gone back up a level, don't label anything
                        level_name = None
                    continue

            # remove duplicates and ensure order
            results = {ko: sorted(set(paths))
                       for ko, paths in results.items()}

        else:
            if level == 'NAME':
                raise Exception("I can't parse name from a .keg file. "
                                "Please use the file: ko/ko")
            if level == 'EC':
                raise Exception("For EC use the ko/ko file or the EC brite: "
                                "brite/ko/ko01000.keg and choose a level")
            raise Exception(f"I don't know what level {level} is!")

    return results


def add_path_arguments(parser, defaults={}, choices={}, helps={}):
    # get format and filter_top_pct arguments from blastm8
    from edl.hits import HITID, ACCS, GIS, KEGG, HITDESC, PFAM
    from edl.blastm8 import add_hit_table_arguments
    add_hit_table_arguments(parser, defaults, flags=['format',
                                                     'filter_top_pct',
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
    default = defaults.get('tab_map_delim', None)
    pgroup.add_argument("--tab_map_delim",
                        default=default,
                        help=("Delimiter to parse multiple assignments in "
                              "map from ids to ko/path/fam. Only used for "
                              "tabular mapping tables. Defaults to {}"
                              .format(str(default))))
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
    if defaults.get("filter_for_path", False):
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
        "--filter_for_path",
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
    kDmap = parse_keg_file(keggFile, 'DESCRIPTION')
    myAssertEq(kDmap['K01623'],
               'ALDO; fructose-bisphosphate aldolase, class I [EC:4.1.2.13]')
    kPmap = parse_keg_file(keggFile, 'PATHWAY')
    assert('K04519' in kPmap)
    assert('K15634' in kPmap)
    myAssertEq(kPmap['K03011'],
               ['00230 Purine metabolism',
                '00240 Pyrimidine metabolism',
                '03020 RNA polymerase',
                "05016 Huntington's disease",
                '05169 Epstein-Barr virus infection'])
    k2map = parse_keg_file(keggFile, 2)
    myAssertEq(k2map['K13810'][0].lower(), 'Carbohydrate Metabolism'.lower())
    myAssertEq(k2map['K13810'][1].lower(), 'Overview'.lower())
    myAssertEq(k2map['K00399'][0].lower(), 'Energy Metabolism'.lower())
    myAssertEq(k2map['K00399'][1].lower(), 'Overview'.lower())
    k3map = parse_keg_file(keggFile, 3)
    myAssertEq(k3map['K13810'], [
        '00010 Glycolysis / Gluconeogenesis [PATH:ko00010]',
        '00030 Pentose phosphate pathway [PATH:ko00030]',
        '00500 Starch and sucrose metabolism [PATH:ko00500]',
        '00520 Amino sugar and nucleotide sugar metabolism [PATH:ko00520]',
        '01230 Biosynthesis of amino acids [PATH:ko01230]',
    ])
    myAssertEq(k3map['K00399'], [
        '00680 Methane metabolism [PATH:ko00680]',
        '01200 Carbon metabolism [PATH:ko01200]',
    ])
    myAssertEq(k3map['K03404'], [
        '00860 Porphyrin and chlorophyll metabolism [PATH:ko00860]',
    ])
    k3mapQ = parse_keg_file(keggFile, '3')
    for k in k3map.keys():
        try:
            myAssertEq(k3map[k], k3mapQ[k])
        except AssertionError:
            raise AssertionError(
                "level 3 classes for %s do not match:\n%s\n%s" %
                (k, k3map[k], k3mapQ[k]))


def testReadKoFile(koFile):
    kPmap = parse_ko_file(koFile, 'PATHWAY')
    assert('K00397' not in kPmap)
    myAssertEq(kPmap['K00399'],
               ['ko00680  Methane metabolism',
                'ko01200  Carbon metabolism'])

    kEmap = parse_ko_file(koFile, 'EC')
    myAssertEq(kEmap['K00397'], ['EC:1.8.99.-'])
    myAssertEq(kEmap['K00399'], ['EC:2.8.4.1'])


if __name__ == '__main__':
    test()
