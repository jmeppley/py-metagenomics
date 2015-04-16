#! /usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
"""
identifyReads.py

Given two lists of taxids and one or more hit tables,
 print out (for each table) the reads that:
     have their best hits in taxid list 1
     have all other hits in either list
"""
from optparse import OptionParser
import sys, re, os, logging
from edl.taxon import *
from edl.blastm8 import filterM8Stream
from edl.hits import *
from edl.util import parseMapFile, addUniversalOptions, setupLogging, parseListToMap, inputIterator, addIOOptions
from edl.expressions import nrOrgRE

def main():
    usage = "usage: %prog -O ORTHOLOGY [OPTIONS] BLAST_M8_FILES"
    description = """
    Given two lists of taxids and one or more hit tables, identify reads that:
     (1) have their best hits in taxid list 1
     (2) have all other hits in either list

    Finally, print out either the hits (that match the target group) for these reads or just read names (-r). The -F filter limits which hits are used in part (2) as well as which are printed.

    The countMethod (-C) option is not used.
    """
    parser = OptionParser(usage, description=description)
    addIOOptions(parser)
    addTaxonOptions(parser,defaults={'mapFile':None,'parseStyle':ACCS,'filterPct':-1,'countMethod':'all','taxdir':None})
    parser.add_option("-g", "--targetTaxonGroup", dest="group1", default=None, metavar="TAXON", action='append',
                      help="Taxon to identify reads in. Top hits (as defined by --topHitPct) must be in this group. It can be a taxid, a name, or a file listing taxids. Use multiple times to specify a list of organisms. Use -a to specify whether all or at least one of the top hits must match.")
    parser.add_option("-a","--any", default=False, action="store_true", help="If specified, accept reads where any top hit is to an organism in the target taxon/taxa. By default, all top hits must be in the target group.")
    addUniversalOptions(parser)
    parser.add_option('-t','--topHitPct', default=0, type='float',
                      help='How close (as a %) to the best score a hit must be to qualify as a top hit. Default is 0, ie must have the best score. Use 100 to get all hits.')
    parser.add_option("-G", "--outerTaxonGroup", dest="group2", default=None, metavar="TAXON", action="append",
                      help="Broader taxon to limit reads. All hits (use -F to limit these hits) must be in the target group or this group. Again, it can be a taxid, a name, or a file listing taxids. It can also be inkoved multiple times to choose multiple groups.")
    parser.add_option('-r','--reads', default=False, action="store_true",
                      help="Output just read names. By default, print the relevant hit lines for each read")

    (options, args) = parser.parse_args()

    if options.about:
        print description
        exit(0)

    # check args
    setupLogging(options,description)
    if options.group1 is None:
        parser.error("Please use -g to specify a target taxonomic group")

    if options.taxdir is not None:
        taxonomy = readTaxonomy(options.taxdir, namesMap=True)
    else:
        taxonomy = None

    group1Map=getGroupMap(options.group1,taxonomy)
    group2Map=getGroupMap(options.group2,taxonomy)
    logging.debug("Group 1 has %d entries and 439482 in group1 is %s" % (len(group1Map),group1Map.get(439482,False)))
    if group2Map is not None:
        logging.debug("Group 2 has %d entries and 439482 in group2 is %s" % (len(group2Map),group2Map.get(439482,False)))

    # map reads to hits
    if options.parseStyle==GIS:
        keyType=int
    else:
        keyType=None
    accToTaxMap = parseMapFile(options.mapFile,valueType=int,keyType=keyType)

    # set up some function pointers
    global hitRE
    hitRE=parsingREs.get(options.parseStyle,None)
    if options.parseStyle == ORGS:
        getTaxid=_getOrgTaxid
    elif options.parseStyle == HITID:
        getTaxid=_getHitidTaxid
    elif options.parseStyle == HITDESC:
        getTaxid=_getHitdescTaxid
    else:
        getTaxid=_getExprTaxid

    # for filtering:
    filterParams = FilterParams.createFromOptions(options)
    logging.debug(repr(filterParams))

    # loop over hit tables
    for (inhandle,outhandle) in inputIterator(args,options):
        readCount=0
        goodReadCount=0
        printCount=0

        # parse file
        for (read,hits) in filterM8Stream(inhandle, filterParams, returnLines=False):
            readCount+=1
            bestScore=0
            hitTaxids={}
            for hit in hits:
                score=hit.score
                taxids=[]
                # does this hit have at least one associated taxid in group2?
                for taxid in getTaxid(hit,accToTaxMap,taxonomy):
                    if taxid is None:
                        break
                    if group2Map is not None and not group2Map.get(taxid,False):
                        break
                    taxids.append(taxid)
                if len(taxids)==0:
                    # nothing matched in the wider group
                    break
                hitTaxids[hit]=taxids

                # find the top score
                if score>bestScore:
                    bestScore=score
            else:
                # if we get here, then every hit was in wider taxon list
                logging.debug("Checking best hits for %s (top score: %.1f)" % (read,bestScore))
                all=True
                recognized=[]
                for hit,taxids in _getBestHitTaxids(hitTaxids,bestScore,options.topHitPct):
                    if _anyTaxidInGroup(taxids,group1Map):
                        logging.debug("%s (%r)  is in group 1" % (hit,taxids))

                        recognized.append(hit)
                    else:
                        logging.debug("%s (%r) is not in group 1" % (hit,taxids))
                        all=False
                if len(recognized)==0:
                    # if none of the best are in our target list, next read
                    logging.debug("No best hits for %s are in group 1" % (read))
                    continue
                if (not options.any) and (not all):
                    # next read unless user said any or all hits are in list
                    logging.debug("Not all best hits for %s are in group 1" % (read))
                    continue

                # if we get here, then the read is a match
                goodReadCount+=1
                if options.reads:
                    logging.debug("Keeping %s" % (read))
                    outhandle.write(read)
                    outhandle.write('\n')
                else:
                    logging.debug("Keeping %d hits for %s" % (len(recognized),read))
                    for hit in sorted(recognized,key=lambda h: (h.score,h.hit)):
                        outhandle.write(hit.getLine(filterParams))
                        printCount+=1

        if options.reads:
            logging.info("Printed %d of %d reads" % (goodReadCount,readCount))
        else:
            logging.info("Printed %d lines for %d of %d reads" % (printCount,goodReadCount, readCount))

def _anyTaxidInGroup(taxids, groupMap):
    for taxid in taxids:
        if groupMap.get(taxid,False):
            return True
    return False

def _getBestHitTaxids(hits, bestScore, pct):
    cutoff=bestScore - pct*bestScore
    for (hit,taxids) in hits.iteritems():
        if hit.score >= cutoff:
            yield hit,taxids

def _getOrgTaxid(hit,taxmap,taxonomy):
    for org in nrOrgRE.findall(hit.hitDesc):
        taxon = getNodeFromHit(org,taxonomy.nameMap)
        if taxon is not None:
            yield taxon.id

def _getExprTaxid(hit,taxmap,taxonomy):
    for hit.hit in hitRE.findall(hit.hit):
        taxon = taxmap.get(hit.hit,None)
        if taxon is not None:
            yield taxon

def _getHitidTaxid(hit,taxmap,taxonomy):
    taxon = taxmap.get(hit.hit,None)
    if taxon is not None:
        yield taxon

def _getHitdescTaxid(hit,taxmap,taxonomy):
    taxon = taxmap.get(hit.hitDesc,None)
    if taxon is not None:
        yield taxon

def getGroupMap(groups,taxonomy):
    """
    Return map from all indicated taxids to True
    input may be a taxid, taxon name, or file listing taxids
    """
    if groups is None or len(groups)==0:
        return None
    taxidmap={}
    for group in groups:
        if os.path.isfile(group):
            taxidmap.update(parseListToMap(group,keyType=int))
        else:
            taxon = getTaxonFromArg(taxonomy, group)
            taxidmap.update( createTaxMap(taxon) )

    return taxidmap

def generateMemberTaxids(node):
    for child in node.children:
        for taxid in generateMemberTaxids(child):
            yield taxid
    yield node.id

def createTaxMap(node):
    idmap={}
    for taxid in generateMemberTaxids(node):
        idmap[taxid]=True
    return idmap

def getTaxonFromArg(taxonomy, arg):
    """
    Given a string from the user, attempt to translate it into a taxon object
    """
    try:
        taxon=taxonomy.idMap[int(arg)]
    except:
        taxon=getNodeFromHit(arg,taxonomy.nameMap)
        if taxon is None:
            raise("Cannot find taxon: %s" % arg)
    return taxon

if __name__=='__main__':
    main()
