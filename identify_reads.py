#!/usr/bin/env python
"""
identify_reads.py

Given two lists of taxids and one or more hit tables,
 print out (for each table) the reads that:
     have their best hits in taxid list 1
     have all other hits in either list
"""
import argparse
import sys
import re
import os
import logging
from edl.taxon import *
from edl.blastm8 import filterM8Stream
from edl.hits import *
from edl.util import add_universal_arguments, setup_logging, \
    parseMapFile, add_IO_arguments, inputIterator, parse_list_to_set
from edl.expressions import nrOrgRE


def main():
    description = """
    Given two lists of taxids and one or more hit tables, identify reads that:
     (1) have their best hits in taxid list 1
     (2) have all other hits in either list

    Finally, print out either the hits (that match the target group) for
    these reads or just read names (-r). The -F filter limits which hits
    are used in part (2) as well as which are printed.

    The countMethod (-C) option is not used.
    """
    parser = argparse.ArgumentParser(description=description)
    add_IO_arguments(parser)
    add_taxon_arguments(
        parser,
        defaults={
            'mapFile': None,
            'parseStyle': ACCS,
            'filter_top_pct': -1,
            'countMethod': 'all',
            'taxdir': None})
    parser.add_argument(
        "-g",
        "--targetTaxonGroup",
        dest="group1",
        default=None,
        metavar="TAXON",
        action='append',
        help="Taxon to identify reads in. Top hits (as defined by "
             "--topHitPct) must be in this group. It can be a taxid, "
             "a name, or a file listing taxids. Use multiple times to "
             "specify a list of organisms. Use -a to specify whether "
             "all or at least one of the top hits must match.")
    parser.add_argument(
        "-a",
        "--any",
        default=False,
        action="store_true",
        help="If specified, accept reads where any top hit is to an organism "
             "in the target taxon/taxa. By default, all top hits must be "
             "in the target group.")
    parser.add_argument(
        '-t',
        '--topHitPct',
        default=0,
        type=float,
        help="How close(as a percentage to the best score a hit must be "
             "to qualify as a top hit. Default is 0, ie must have the best "
             "score. Use 100 to get all hits.")
    parser.add_argument(
        "-G",
        "--outerTaxonGroup",
        dest="group2",
        default=None,
        metavar="TAXON",
        action="append",
        help="Broader taxon to limit reads. All hits (use -F to limit "
             "these hits) must be in the target group or this group. Again, "
             "it can be a taxid, a name, or a file listing taxids. "
             "It can also be inkoved multiple times to choose multiple "
             "groups.")
    parser.add_argument(
        '-r',
        '--reads',
        default=False,
        action="store_true",
        help="Output just read names. By default, print the relevant hit "
             "lines for each read")

    # log level and help
    add_universal_arguments(parser)
    arguments = parser.parse_args()
    setup_logging(arguments)

    # check args
    if arguments.group1 is None:
        parser.error("Please use -g to specify a target taxonomic group")

    if arguments.taxdir is not None:
        taxonomy = readTaxonomy(arguments.taxdir, namesMap=True)
    else:
        taxonomy = None

    group_1_set = get_group_set(arguments.group1, taxonomy)
    group_2_set = get_group_set(arguments.group2, taxonomy)
    logging.debug(
        "Group 1 has %d entries and 439482 in group1 is %s" %
        (len(group_1_set), 439482 in group_1_set))
    if group_2_set is not None:
        logging.debug(
            "Group 2 has %d entries and 439482 in group2 is %s" %
            (len(group_2_set), 439482 in group_2_set))

    # map reads to hits
    if arguments.parseStyle == GIS:
        keyType = int
    else:
        keyType = None
    accToTaxMap = parseMapFile(
        arguments.mapFile,
        valueType=int,
        keyType=keyType)

    # set up some function pointers
    global hitRE
    hitRE = parsingREs.get(arguments.parseStyle, None)
    if arguments.parseStyle == ORGS:
        getTaxid = _getOrgTaxid
    elif arguments.parseStyle == HITID:
        getTaxid = _getHitidTaxid
    elif arguments.parseStyle == HITDESC:
        getTaxid = _getHitdescTaxid
    else:
        getTaxid = _getExprTaxid

    # for filtering:
    filterParams = FilterParams.create_from_arguments(arguments)
    logging.debug(repr(filterParams))

    # loop over hit tables
    for (inhandle, outhandle) in inputIterator(arguments):
        readCount = 0
        goodReadCount = 0
        printCount = 0

        # parse file
        for (
                read,
                hits) in filterM8Stream(
                inhandle,
                filterParams,
                returnLines=False):
            readCount += 1
            bestScore = 0
            hitTaxids = {}
            for hit in hits:
                score = hit.score
                taxids = []
                # does this hit have at least one associated taxid in group2?
                for taxid in getTaxid(hit, accToTaxMap, taxonomy):
                    if taxid is None:
                        break
                    if group_2_set is not None and taxid not in group_2_set:
                        break
                    taxids.append(taxid)
                if len(taxids) == 0:
                    # nothing matched in the wider group
                    break
                hitTaxids[hit] = taxids

                # find the top score
                if score > bestScore:
                    bestScore = score
            else:
                # if we get here, then every hit was in wider taxon list
                logging.debug(
                    "Checking best hits for %s (top score: %.1f)" %
                    (read, bestScore))
                all = True
                recognized = []
                for hit, taxids in _getBestHitTaxids(
                        hitTaxids, bestScore, arguments.topHitPct):
                    if _anyTaxidInGroup(taxids, group_1_set):
                        logging.debug("%s (%r)  is in group 1" % (hit, taxids))

                        recognized.append(hit)
                    else:
                        logging.debug(
                            "%s (%r) is not in group 1" %
                            (hit, taxids))
                        all = False
                if len(recognized) == 0:
                    # if none of the best are in our target list, next read
                    logging.debug(
                        "No best hits for %s are in group 1" %
                        (read))
                    continue
                if (not arguments.any) and (not all):
                    # next read unless user said any or all hits are in list
                    logging.debug(
                        "Not all best hits for %s are in group 1" %
                        (read))
                    continue

                # if we get here, then the read is a match
                goodReadCount += 1
                if arguments.reads:
                    logging.debug("Keeping %s" % (read))
                    outhandle.write(read)
                    outhandle.write('\n')
                else:
                    logging.debug(
                        "Keeping %d hits for %s" %
                        (len(recognized), read))
                    for hit in sorted(
                        recognized,
                        key=lambda h: (
                            h.score,
                            h.hit)):
                        outhandle.write(hit.getLine(filterParams))
                        printCount += 1

        if arguments.reads:
            logging.info("Printed %d of %d reads" % (goodReadCount, readCount))
        else:
            logging.info(
                "Printed %d lines for %d of %d reads" %
                (printCount, goodReadCount, readCount))


def _anyTaxidInGroup(taxids, group_set):
    for taxid in taxids:
        if taxid in group_set:
            return True
    return False


def _getBestHitTaxids(hits, bestScore, pct):
    cutoff = bestScore - pct * bestScore
    for (hit, taxids) in hits.items():
        if hit.score >= cutoff:
            yield hit, taxids


def _getOrgTaxid(hit, taxmap, taxonomy):
    for org in nrOrgRE.findall(hit.hitDesc):
        taxon = getNodeFromHit(org, taxonomy.nameMap)
        if taxon is not None:
            yield taxon.id


def _getExprTaxid(hit, taxmap, taxonomy):
    for hit.hit in hitRE.findall(hit.hit):
        taxon = taxmap.get(hit.hit, None)
        if taxon is not None:
            yield taxon


def _getHitidTaxid(hit, taxmap, taxonomy):
    taxon = taxmap.get(hit.hit, None)
    if taxon is not None:
        yield taxon


def _getHitdescTaxid(hit, taxmap, taxonomy):
    taxon = taxmap.get(hit.hitDesc, None)
    if taxon is not None:
        yield taxon


def get_group_set(groups, taxonomy):
    """
    Return set of all indicated taxids
    input may be a taxid, taxon name, or file listing taxids
    """
    if groups is None or len(groups) == 0:
        return None
    taxidset = set()
    for group in groups:
        if os.path.isfile(group):
            taxidset.update(parse_list_to_set(group, keyType=int))
        else:
            taxon = getTaxonFromArg(taxonomy, group)
            taxidset.update(create_taxid_set(taxon))

    return taxidset


def generateMemberTaxids(node):
    for child in node.children:
        for taxid in generateMemberTaxids(child):
            yield taxid
    yield node.id


def create_taxid_set(node):
    idset = set()
    for taxid in generateMemberTaxids(node):
        idset.add(taxid)
    return idset


def getTaxonFromArg(taxonomy, arg):
    """
    Given a string from the user, attempt to translate it into a taxon object
    """
    try:
        taxon = taxonomy.idMap[int(arg)]
    except Exception:
        taxon = getNodeFromHit(arg, taxonomy.nameMap)
        if taxon is None:
            raise "Cannot find taxon: %s"
    return taxon


if __name__ == '__main__':
    main()
