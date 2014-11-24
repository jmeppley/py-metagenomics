import sys, logging, os, re
import edl.util, edl.taxon, edl.silva

logging.basicConfig(stream=sys.stdout)
logger=logging.getLogger(sys.argv[0])
logger.setLevel(logging.INFO)
#logger.setLevel(logging.DEBUG)

## Some ranks need to be renamed to work with existing scripts
# superkingdom->major_clade: NCBI uses superkingdom for domain, so this has to be changed to a reasonable alternative
rankMapping={"superkingdom":"major_clade"}

def main():
    # command line arguments (this could be more robust...)
    logger.debug(sys.argv)
    (taxfile, fastafile, dumpDir) = sys.argv[1:]

    # parse the input files
    rootNode, taxmap = buildSilvaTree(taxfile, fastafile)

    nodesFile = os.path.sep.join((dumpDir, 'nodes.dmp'))
    namesFile = os.path.sep.join((dumpDir, 'names.dmp'))
    taxFile = os.path.sep.join((dumpDir, 'lastdb.tax'))

    logger.info("Writing taxonomy to %s and %s" % (nodesFile, namesFile))
    with open(nodesFile,'w') as nodeFile:
        with open(namesFile,'w') as nameFile:
            writeDumpFiles(rootNode, nodeFile, nameFile)

    logger.info("Writing taxid table to %s" % (taxFile))
    with open(taxFile,'w') as taxMapFile:
        for (hitid, taxNode) in taxmap.iteritems():
            taxMapFile.write("%s\t%d\n" % (hitid, taxNode.id))

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

def buildSilvaTree(taxFile, fastaFile):
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
        taxmap[hitid]=node

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

