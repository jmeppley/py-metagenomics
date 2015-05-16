#!/usr/bin/env python
"""
Build ncbi-like taxonomy files for Silva or PR2.

The files nodes.dmp and names.dmp and lastdb.tax are written to teh specified output directory. The first two define the taxonomy (in NCBI format) and the last is a mapping from sequence IDs to taxids.

For silva, the taxonomy text file is used to get the ranks for each taxon. Then the fasta headers are parsed to find where to place each sequence on the tree. 

For PR2, the fasta file is all that's needed, but you can optionanly choose to simplify the hidids. If you specify an output fasta file, the first element (accession.start.end) is used, but if you don't, the mapping will use the original ID's which include the full taxonomy string. 
"""
from optparse import OptionParser
import sys, logging, os, re
import edl.util, edl.taxon, edl.silva

## Some ranks need to be renamed to work with existing scripts
# superkingdom->major_clade: NCBI uses superkingdom for domain, so this has to be changed to a reasonable alternative
rankMapping={"superkingdom":"major_clade"}

def main():
    usage = "usage: %prog [options] FASTA_FILE OUTPUT_DIR"
    description = __doc__

    parser = OptionParser(usage, description=description)
    parser.add_option("-t", "--taxfile", dest="taxfile",
                      metavar="FILE", help="Read Silva ranks from FILE")
    parser.add_option("-d","--dbType",default="silva",choices=['silva','pr2'],
            help="Which database are we importing: 'silva' (default) or 'pr2'")
    parser.add_option("-f","--fastaout",default=None, metavar='FILE',
            help="Write fasta with modified headers to FILE. Only used for pr2")
    parser.add_option("-i","--idStart", default=0, type=int,
            help="Start taxid counting with this number")

    # logging and help
    edl.util.addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    edl.util.setupLogging(options, description)
    logger=logging.getLogger(sys.argv[0])

    # command line arguments (this could be more robust...)
    logger.debug(args)
    if len(args)==3:
        #hack to support old method for just Silva
        options.taxfile = args.pop(0)
    (fastafile, dumpDir) = args

    # parse the input files
    if options.dbType=='pr2':
        rootNode, taxmap = buildPR2Tree(fastafile, 
                                        fastaout=options.fastaout,
                                        nextId=options.idStart,
                                        )
    else:
        if options.taxfile is None:
            parser.error("You must provide the tax file for Silva (-t)")
        rootNode, taxmap = buildSilvaTree(options.taxfile, fastafile, logger)

    nodesFile = os.path.sep.join((dumpDir, 'nodes.dmp'))
    namesFile = os.path.sep.join((dumpDir, 'names.dmp'))
    taxFile = os.path.sep.join((dumpDir, 'lastdb.tax'))

    logger.info("Writing taxonomy to %s and %s" % (nodesFile, namesFile))
    with open(nodesFile,'w') as nodeFile:
        with open(namesFile,'w') as nameFile:
            writeDumpFiles(rootNode, nodeFile, nameFile)

    logger.info("Writing taxid table to %s" % (taxFile))
    with open(taxFile,'w') as taxMapFile:
        for (hitid, taxid) in taxmap.iteritems():
            taxMapFile.write("%s\t%s\n" % (hitid, taxid))

def treeGenerator(node, kidsFirst=False, **kwargs):
    if not kidsFirst:
        yield node
    for kid in sorted(node.children,**kwargs):
        if kid is not node:
            for n in treeGenerator(kid, kidsFirst=kidsFirst, **kwargs):
                yield n
    if kidsFirst:
        yield node

def writeDumpFiles(rootNode, nodeStream, nameStream):
    for node in treeGenerator(rootNode):
        nid = node.id
        nname = node.name
        if node==rootNode:
           nparent = node.id
        else:
           nparent = node.parent.id
        if node.rank == 'domain':
            nrank = 'superkingdom'
        elif node.rank in edl.taxon.ranks:
            nrank = node.rank
        else:
            nrank = "no rank"
        nodeStream.write("%s\t|\t%s\t|\t%s\t\n" % (nid,nparent,nrank))
        nameStream.write("%s\t|\t%s\t|\t\t|\tscientific name\t\n" % (nid,nname))

def getOrgsFromSilvaFasta(fasta):
    with open(fasta) as f:
        for line in f:
            if len(line)>0 and line[0]=='>':
                (read,desc) = line[1:].strip().split(None,1)
                yield (read,desc)

ranks=['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
def buildPR2Tree(fastaFile, fastaout=None, nextId=0):
    """
    Given a fasta file from PR2, parse the headerlines into a taxonomy and a map from hit-ids to taxon ids
    example line

    >10-044.1.1773|Eukaryota|Stramenopiles|Stramenopiles_X|Oomycota|Oomycota_X|Oomyc
    ota_XX|Oomycota_XXX|Oomycota_XXX+sp.

    The fasta headers are too long. We'll just save the first element here.
    NOTE: This means the fasta file will have to be modified to match!

    """
    
    # create core of tree from taxonomy text file
    taxMap={}
    prTree={}
    taxaByNameAndRank={r:{} for r in ranks}
    rootNode=edl.taxon.TaxNode('root',None,None)
    prTree[rootNode.id]=rootNode
    if fastaout:
        outhandle = open(fastaout,'w')
    with open(fastaFile) as f:
        for line in f:
            # parse header line
            if len(line)>0 and line[0]=='>':
                names=line[1:].rstrip('\n\r').split("|")
                hitid=names.pop(0)

                # build tree from header
                lastNode=rootNode
                for rank, taxName in zip(ranks,names):
                    if taxName not in taxaByNameAndRank[rank]:
                        nextId+=1
                        newNode=edl.taxon.TaxNode(nextId, lastNode.id, None)
                        newNode.name=taxName
                        newNode.rank=rank
                        newNode.setParent(lastNode)
                        prTree[newNode.id]=newNode
                        taxaByNameAndRank[rank][taxName]=newNode
                        lastNode=newNode
                    else:
                        lastNode=taxaByNameAndRank[rank][taxName]

                # save hitid->taxid mapping
                taxMap[hitid]=lastNode.id

                # print modified header
                if fastaout:
                    outhandle.write('>%s %s\n' % (hitid,";".join(names)))
            else:
                if fastaout:
                    outhandle.write(line)

    # There should only be one child of root: Eukaryota
    if len(rootNode.children) != 1:
        logging.warn("There should only be one child of root! Not %d:\n%s" % \
                (len(rootNode.children),
                 ", ".join([c.name for c in rootNode.children])))
    rootNode=rootNode.children[0]
    logging.info("New root is: %s" % (rootNode.name))

    if fastaout:
        outhandle.close()

    return (rootNode, taxMap)


def buildSilvaTree(taxFile, fastaFile, logger):
    """
    Given a text taxonomy file (lineage <tab> id <tab> rank) and a fasta file with full lineages as the description:
    Return the root node from a taxonomy of edl.taxon.Node objects and a mapping from fasta record IDs to taxids.
    """
    rankMap=edl.util.parseMapFile(taxFile, keyCol=0, valueCol=2, skipFirst=0)
    silvaTaxidMap=edl.util.parseMapFile(taxFile, keyCol=0, valueCol=1, valueType=int, skipFirst=0)
    
    # create core of tree from taxonomy text file
    silvaTree={}
    maxTaxid=max(silvaTaxidMap.values())
    for (lineage, rank) in rankMap.iteritems():
        node=edl.silva.SilvaTaxNode.addToTreeFromString(lineage.strip("; "), silvaTree)
        node.rank = rankMapping.get(rank,rank)
        node.ncbi_tax_id = silvaTaxidMap[lineage]
        if not isinstance(node.ncbi_tax_id,int):
            logger.warn("NCBI taxid is not an int: %s (%s)" % (node.ncbi_tax_id, node.name))

    logger.info("Built tree of %d taxa with the largest ID of %d" % (len(silvaTree),maxTaxid))

    # Add leaves to tree from lineages in fasta file and build mapping
    taxmap={}
    for (hitid,lineage) in getOrgsFromSilvaFasta(fastaFile):
        node = edl.silva.SilvaTaxNode.addToTreeFromString(lineage, silvaTree)
        taxmap[hitid]=node.id

    logger.info("Added nodes from fasta file for a total of %d" % (len(silvaTree)))

    rootNode=silvaTree.itervalues().next().getRootNode()
    # make sure everything is OK
    for node in edl.util.treeGenerator(rootNode):
        if not isinstance(node.id,int):
            if "ncbi_tax_id" in dir(node):
                node.id = int(node.ncbi_tax_id)
            else:
                maxTaxid+=1
                node.id=maxTaxid

    return (rootNode, taxmap)

if __name__ == '__main__':
    main()

