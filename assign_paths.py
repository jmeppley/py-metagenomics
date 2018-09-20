#! /usr/bin/env python
"""
Takes a single hit table file and generates a table (or tables) of
pathway/gene family assignments for the query sequences (aka 'reads').
Assignments can be for gene families, gene classes, or pathways. Multiple
pathway or classification levels can be given. If they are, an assignment
will be made at each level.

This differs from assignPathsToReadsFromBlast.py in that: (1) it can
handle CAZy and SEED, (2) it will output multiple levels in one file,
(3) multiple assignments are always printed on multiple lines.

This script will work with KEGG, SEED, or CAZy. CAZy only has one level
of heirarchy, the others have 3. The CAZy heirarchy is apparent from the
hit name and needs no supporting files. KEGG and SEED require mapping files
to identify gene families and heirachy files to report levels other than
the gene family or ortholog level. Both SEED and KEGG have three levels
of classifications that can be indicated with a 1, 2, or 3. The words
"subsystem" and "pathway" are synonyms for level 3.

If a count method is selected that can produce multiple assignments per
read, each assignment will be printed on a new line.

NOTE: in KEGG (and SEED) a single ortholog (role) may belong to multiple
pathways (subsystems). A hit to such an ortholog will result in extra
assignment values for that query sequence (1 for each pathway it belongs to).
"""

import argparse
import logging
import re
from edl import hits, kegg
from edl.util import add_IO_arguments, add_universal_arguments, \
        inputIterator, parseMapFile, setup_logging


def main():
    description = __doc__
    parser = argparse.ArgumentParser(description)
    add_IO_arguments(parser)
    parser.add_argument("-l", "--level", dest="levels", default=None,
                        metavar="LEVEL", action="append",
                        help=""" Level(s) to collect counts on. Use flag
                      multiple times to specify multiple levels. If multiple
                      values given, one table produced for each with rank
                      name appended to file name. Levels can be an integer
                      (1-3) for KEGG or SEED levels, any one of 'gene',
                      'role', 'family',
                      'ko', or 'ortholog' (which are all synonyms), or
                      anything not synonymous with 'gene' to
                      get CAZy groups. Defaults to ortholog/role and
                      levels 1, 2, and 3 for KEGG and SEED
                      and gene and group for CAZy and COG.""")
    parser.add_argument(
        '-S',
        '--squash',
        dest='splitForLevels',
        default=True,
        action='store_false',
        help="Don't split assignment rows if gene maps to multiple pathways, "
             "just squash them into one row using python list syntax")

    # format, ortholog heirarchy, and more
    kegg.add_path_arguments(parser)

    # log level and help
    add_universal_arguments(parser)
    arguments = parser.parse_args()
    setup_logging(arguments)

    # Set defaults and check for some conflicts
    if arguments.levels is None and arguments.heirarchyFile is None:
        # using hit names only
        arguments.levels = [None]
    else:
        if arguments.heirarchyFile is None \
                and arguments.heirarchyType != 'cazy':
            logging.warn("Type: %s" % (arguments.heirarchyType))
            parser.error("Cannot select levels without a heirarchy (ko) file")
        if arguments.levels is None:
            # set a default
            if arguments.heirarchyType is 'kegg':
                arguments.levels = ['ko', '1', '2', 'pathway']
            if arguments.heirarchyType is 'seed':
                arguments.levels = ['role', '1', '2', 'subsystem']
            else:
                arguments.levels = ['gene', 'group']

        try:
            # Make sure the level list makes sense
            arguments.levels = cleanLevels(arguments.levels)
        except Exception as e:
            parser.error(str(e))

    # map reads to hits
    if arguments.mapFile is not None:
        if arguments.mapStyle == 'auto':
            with open(arguments.mapFile) as f:
                firstLine = next(f)
                while len(firstLine) == 0 or firstLine[0] == '#':
                    firstLine = next(f)
            if koMapRE.search(firstLine):
                arguments.mapStyle = 'kegg'
            elif seedMapRE.search(firstLine):
                arguments.mapStyle = 'seed'
            elif tabMapRE.search(firstLine):
                arguments.mapStyle = 'tab'
            elif cogMapRE.search(firstLine):
                arguments.mapStyle = 'cog'
            else:
                raise Exception(
                    "Cannot figure out map type from first line:\n%s" %
                    (firstLine))

        logging.info("Map file seems to be: %s" % (arguments.mapStyle))
        if arguments.mapStyle == 'kegg':
            valueMap = kegg.parseLinkFile(arguments.mapFile)
        elif arguments.mapStyle == 'seed':
            valueMap = kegg.parseSeedMap(arguments.mapFile)
        elif arguments.mapStyle == 'cog':
            valueMap = kegg.parseCogMap(arguments.mapFile)
        else:
            if arguments.parseStyle == hits.GIS:
                keyType = int
            else:
                keyType = None
            valueMap = parseMapFile(
                arguments.mapFile,
                valueType=None,
                valueDelim=arguments.tab_map_delim,
                keyType=keyType)
        if len(valueMap) > 0:
            logging.info("Read %d items into map. EG: %s" %
                         (len(valueMap), next(iter(valueMap.items()))))
        else:
            logging.warn("Read 0 items into value map!")
    else:
        valueMap = None

    # set up level mapping
    levelMappers = [getLevelMapper(l, arguments) for l in arguments.levels]

    # parse input files
    for (inhandle, outhandle) in inputIterator(arguments):
        logging.debug(
            "Reading from %s and writing to %s" %
            (inhandle, outhandle))
        hitMapIter = hits.parseM8FileIter(
            inhandle,
            valueMap,
            hits.FilterParams.create_from_arguments(arguments),
            #arguments.hitTableFormat,
            #arguments.filter_top_pct,
            arguments.parseStyle,
            arguments.countMethod,
            ignoreEmptyHits=arguments.mappedHitsOnly)

        if arguments.levels == [None]:
            arguments.levels = ['Hit']
        outhandle.write("Read\t%s\n" % ('\t'.join(arguments.levels)))
        for read, hitIter in hitMapIter:
            assignments = []
            for hit in hitIter:
                logging.debug("Hit: %s" % (hit))
                assignment = []
                for levelMapper in levelMappers:
                    assignment.append(levelMapper(hit))
                assignments.append(assignment)
            logging.debug("Read %s has %d hits" % (read, len(assignments)))
            for assignment in assignments:
                for assignmentList in handleMultipleMappings(
                        assignment, arguments):
                    outhandle.write(
                        "%s\t%s\n" %
                        (read, "\t".join(assignmentList)))


def cleanLevels(levelList):
    # don't allow duplicates
    newList = list(set(levelList))
    newList.sort(key=lambda l: levelList.index(l))

    # return levels
    return newList


def getCazyGroup(gene):
    m = cazyRE.search(gene)
    if m:
        cazygroup = m.group(1)
        logging.debug("Mapping %s to %s" % (gene, cazygroup))
    else:
        logging.debug(
            "Could not parse group from %s with r'%s'" %
            (gene, cazyRE.pattern))
        cazygroup = gene
    return cazygroup


# levels equivalent to returning just the ko/gene
koSyns = [None, 'ko', 'gene', 'ortholog', 'family', 'role']
level3Syns = ['subsystem', 'pathway', 'group']
koMapRE = re.compile(r'\sko:K\d{5}')
seedMapRE = re.compile(r'^Mapped roles:')
cogMapRE = re.compile(r'^\d+\t[KC]OG\d+\t')
tabMapRE = re.compile(r'^[^\t]+\t[^\t+]')
cazyRE = re.compile(r'([a-zA-Z]+)\d+')


def getLevelMapper(level, arguments):
    if level in koSyns:
        return lambda h: h
    if arguments.heirarchyType == 'cazy':
        return getCazyGroup

    lookupLevel = level if level not in level3Syns else '3'

    if arguments.heirarchyType == 'kegg':
        # Ideally, we'd be able to parse the heirachy once, but the current
        # KEGG code just retuns simple mappings
        logging.info(
            "Reading KEGG level %s assignments from %s" %
            (level, arguments.heirarchyFile))
        geneTranslation = kegg.readKEGGFile(
            arguments.heirarchyFile, lookupLevel)
    else:
        # SEED or COG/KOG
        if arguments.heirarchyType == 'seed':
            logging.info(
                "Reading SEED subsystem assignments from %s" %
                (arguments.heirarchyFile))
            seedTree = kegg.readSEEDTree(arguments.heirarchyFile)
        elif arguments.heirarchyType == 'cog':
            logging.info(
                "Reading COG subsystem assignments from %s" %
                (arguments.heirarchyFile))
            seedTree = kegg.readCogTree(arguments.heirarchyFile)

        geneTranslation = seedTree[lookupLevel]
    return lambda gene: geneTranslation.get(gene, gene)


def handleMultipleMappings(assignmentList, arguments):
    """
    given list of assignments where each element coresponds to a hierarchy
    level and can contain a sigle assignment or a list of assignments
    return a list of data "rows"
    """
    if arguments.splitForLevels:
        # return a line for each possibility
        newAssignmentListArray = [assignmentList]
        for i in range(len(assignmentList)):
            item = assignmentList[i]
            if isinstance(item, list) or isinstance(item, tuple):
                logging.debug("Splitting %r" % (item))
                itemList = set(item)
                tempList = []
                for al in newAssignmentListArray:
                    for item in sorted(itemList):
                        al[i] = item
                        tempList.append(sorted(list(al)))
                newAssignmentListArray = tempList
        for al in newAssignmentListArray:
            for i in range(len(al)):
                al[i] = str(al[i])
        return newAssignmentListArray
    else:
        assignment = [
            str(a) if isinstance(
                a, str) or a is None else str(
                sorted(
                    list(
                        set(a)))) for a in assignmentList]
        return [assignment, ]


if __name__ == '__main__':
    main()
