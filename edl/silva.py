import edl.taxon, re, logging, edl.util, os

ignoreList=['marinemetagenome','uncultured','incertaesedis']
translations={'RoseobactercladeNAC11-7lineage':'uncultured Roseobacter NAC11-7',
              'OCS116clade':'uncultured OCS116 cluster bacterium',
              'OrderIIIIncertaeSedis':'Bacteroidetes Order III. Incertae sedis',
              'mitochondria':'Eukaryota',
              'Mitochondria':'Eukaryota',
              'chloroplast':'Eukaryota',
              'Chloroplast':'Eukaryota',
              'Cosenzaeamyxofaciens':'Proteus myxofaciens',
              'Escherichiablattae':'Escherichia blattae',
              'Psychrobacillusinsolitus':'Bacillus insolitus',
              'Psychrobacilluspsychrotolerans':'Bacillus psychrotolerans',
              'Psychrobacilluspsychrodurans':'Bacillus psychrodurans',
              'Kandleriavitulina':'Lactobacillus vitulinus',
              'Eggerthiacatenaformis':'Lactobacillus catenaformis',
             }

def getNodeFromHit(hit, taxonomy):
    hit = edl.taxon.removeSpaces(hit)
    hit=translations.get(hit,hit)
    hit=edl.taxon.simplifyString(hit)

    # Some strings don't map well to NCBI. Return none and let code try parent.
    if hit in ignoreList:
        return None
    if len(hit)>10 and hit[:10]=='uncultured':
        return None

    node=taxonomy.realNameMap.get(hit,None)
    if node is not None:
        return node

    return edl.taxon.getNodeFromHit(hit, taxonomy.nameMap, exhaustive=False)

def filterStream(stream,expression):
    """
    Simple generator to filer an input stream
    """
    for string in stream:
        if expression.search(string):
            yield string

class SilvaTaxNode(edl.taxon.TaxNode):
    """
    A node in a pylogenetic tree. This has a parent and many children
    """

    @staticmethod
    def addToTreeFromString(taxString, tree={}, removeSpaces=True):
        if 'root' not in tree:
            if len(tree)>0:
                raise Error('tree must have root node!')
            root =  SilvaTaxNode('root',None,None)
            root.name = root.id
            tree['root']=root

        if removeSpaces:
            lineage = scRE.split(edl.taxon.removeSpaces(taxString))
        else:
            lineage = scRE.split(taxString)
        logging.debug("parsing %s: %s"  % (taxString,str(lineage)))
        lastNode = tree['root']
        for i in range(len(lineage)):
            taxon = ";".join(lineage[:i+1])
            if (taxon in tree):
                lastNode = tree[taxon]
            else:
                #logging.debug("Adding %s as child to %s" % (taxon, lastNode.name))
                newNode = SilvaTaxNode(taxon, lastNode.id, None)
                newNode.name = lineage[i]
                newNode.setParent(lastNode)
                tree[taxon] = newNode
                lastNode=newNode
        return lastNode


    # override versoin of this method in library
    def getCollapsedCounts2(self, counts, cutoff, collapsedCounts):
        """
        walk the tree and fill collapseCounts hash with counts for nodes
        over the cutoff value
        """
        count=0
        if self in counts:
            # use value in counts if it's there
            count=counts[self]
        else:
            # otherwise use sum of children (those under cutoff)
            for child in self.children:
                if child is self:
                    continue
                kidCount = child.getCollapsedCounts(counts, cutoff, collapsedCounts)
                if kidCount is not None:
                    count += kidCount

        # if this node has a generic name, defer to parent
        if self.isNameGeneric():
            logging.debug("%s is too generic, giving %d hits to parent: %s" % (self.name,count,self.parent.name))
            return count

        # if this node is over cutoff, add to counts
        if count>=cutoff:
            name = self.name
            #name = self.getLineageString('; ')
            if name not in collapsedCounts:
                collapsedCounts[name]=count
            else:
                logging.warn("Duplicate name in collapsed count: %s" % name)
                collapsedCounts[name]+=count
            return None
        else:
            # otherwise return count for parent to use
            return count


def translateSSUTaxToNCBIDump(ssuFasta,ssuRanks,outDir):
    """
    Given an SSU FASTA file and ranks file form Silva, generate a
    direcory with names.dmp and nodes.dmp in the NCBI style plus a
    mapping from SSU hit IDs to taxids in the new files.
    """

    # Create ranked nodes from rank file
    rankMap=edl.util.parseMapFile(ssuRanks, valueCol=2, skipFirst=1)
    silvaTree={}
    for (lineage, rank) in rankMap.iteritems():
        node=SilvaTaxNode.addToTreeFromString(lineage, silvaTree, removeSpaces=False)
        node.rank = rank

    logging.debug("Generated %d taxa from ranks file" % (len(silvaTree)))

    # Add nodes from fasta file (and keep hit associations)
    taxmap={}
    for (hitid,lineage) in getOrgsFromSSUFasta(ssuFasta):
        node = SilvaTaxNode.addToTreeFromString(lineage, silvaTree, removeSpaces=False)
        taxmap[hitid]=node

    logging.debug("Generated %d taxa from fasta file" % (len(silvaTree)))

    # Assign taxids to nodes (parsing expects integers)
    i=0
    rootNode = silvaTree.itervalues().next().getRootNode()
    for node in edl.util.treeGenerator(rootNode):
        i+=1
        node.id=i

    # Write dump files
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    with open(os.path.sep.join((outDir,'nodes.dmp')),'w') as nodeFile:
        with open(os.path.sep.join((outDir, 'names.dmp')),'w') as nameFile:
            writeDumpFiles(rootNode, nodeFile, nameFile)


    # Write hit->tax mapping file
    with open(os.path.sep.join((outDir,'ssu.hitid.to.taxid.silva')),'w') as taxMapFile:
        for (hitid, taxNode) in taxmap.iteritems():
            taxMapFile.write("%s\t%d\n" % (hitid, taxNode.id))

def writeDumpFiles(rootNode, nodeStream, nameStream):
    for node in edl.util.treeGenerator(rootNode):
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

def getOrgsFromSSUFasta(ssuFasta):
    with open(ssuFasta) as f:
        for line in f:
            if len(line)>0 and line[0]=='>':
                (read,desc) = line[1:].strip().split(None,1)
                yield (read,desc)




################
# compiled REs #
################
silvaFastaRE=re.compile(r'^>(\S+_[LS]SU)\s+(\S.*)$')
parensRE = re.compile(r'\([^\(\)]+\)')
lastSemicolonRE = re.compile(r'^.*;([^;]+)$')
semicolonRE=re.compile(r';')
spaceRE=re.compile("\s")
dotRE=re.compile("\.")
scRE=re.compile(r';+')
metagenomeRE=re.compile(r'metagenome')
silvaLineRE=re.compile(r'_(SSU|LSU)\b')
