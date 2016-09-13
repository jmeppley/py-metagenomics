"""
Tools for generating sunburst figures.
convertToNX: create a networkx representation of the figure
convertToJSON: create (from nx) a JSON rep for use with d3.js
plotSunburst: use matplotlib to create a figure
"""

import networkx as nx
import numpy as np
import re, sys, logging
from edl.util import treeGenerator
try:
    import simplejson as json
except:
    import json as json
from edl.taxon import ranks as rankOrder

logger=logging.getLogger(__name__)

# Constants and expressions
nameCleanRE=re.compile(r'\s*\<[^<>]+\>\s*$')
pi=np.pi
NAME='name'
VALUE='size'
COLOR='color'
KIDS='children'

def convertToNx(counts, leaves=True, ranks=['superkingdom','class','genus'], colorMap=None,**kwargs):
    """
    Take a dict mapping TaxNodes (from NCBI taxonomy) to counts.
    Return a netowrkx directed graph (tree) limited to selected ranks
    """
    if ranks is None:
        raise Exception("This does not work unless you specify a rank. Try going straing to JSON with convertTaxToJSON()")

    valueKey=kwargs.get(VALUE,VALUE)
    nameKey=kwargs.get(NAME,NAME)
    colorKey=kwargs.get(COLOR,COLOR)

    # merge counts for None in to counts for root
    root = next(iter(counts.keys())).getRootNode()
    counts[root]=counts.get(root,0)+counts.pop(None,0)

    # set up graph
    G=nx.DiGraph()
    _add_node(G,root,colorMap,nameKey,colorKey)
    G.node[root][valueKey]=counts[root]

    ranks = sorted(ranks, key=rankOrder.index, reverse=False)

    #for node, count in _nodeCountGenerator(root,counts):
    #for node in _treeGenerator(root):
    used = set([root,])
    for node, count in counts.items():
        if node in used:
            # We've already seen this node
            continue

        # Is this an outermost node?
        if leaves and isLeaf(node,counts):
            nodeIsLeaf = True
            logger.info("%s(%d) is leaf!" % (node,count))
        else:
            nodeIsLeaf = False
            logger.debug("%s (%d) is not leaf!" % (node,count))

        # Walk towards root from node until we hit a node we've already done
        parent=node
        countSum=0
        lastNode = None
        while True:
            logger.debug("Parent: %s(%s)" % (parent.name,parent.id))
            # is this node one of the ranks to display
            if parent.rank in ranks or nodeIsLeaf or parent is root:
                if parent not in used:
                    # Add node for this leaf or ranked node
                    _add_node(G, parent, colorMap, nameKey, colorKey)
                    G.node[parent][valueKey] = counts.get(parent,0)
                    logger.debug("Adding parent: %s (%d)" % (parent,G.node[parent][valueKey]))
                    breakAfter=False
                    used.add(parent)
                else:
                    breakAfter=True
                if parent is not node:
                    logger.debug("adding counts to %s (%d + %d)" % (parent, G.node[parent][valueKey],countSum))
                    G.node[parent][valueKey] = G.node[parent].get(valueKey,0) + countSum
                    if lastNode is not None:
                        logger.debug("linking %s to %s" % (lastNode, parent))
                        G.add_edge(parent, lastNode)
                if breakAfter:
                    break
                lastNode=parent
                nodeIsLeaf=False
                countSum=0
            else:
                # save counts for next level
                if not nodeIsLeaf:
                    if parent not in used:
                        countSum += counts.get(parent,0)
                        used.add(parent)
            parent = parent.parent
        used.add(node)

    #print (G.nodes())
    #print (G.edges())
    return (G,root)

def isLeaf(taxon,counts):
    """
    This is a leaf if none of it's descendants have any counts
    (and if it has counts, but we know that already)
    """
    return _sumChildren(taxon, counts)==0

def _sumChildren(node, counts):
    try:
        return node._childrenTotal
    except:
        total=0
        for kid in node.children:
            total+=counts.get(kid,0)
            total+=_sumChildren(kid, counts)
        node._childrenTotal = total
        return total

def _clearKidSums(root):
    for kid in node.children:
        _clearKidSums(kid)
    # I wonder if this will work...
    del root._childrenTotal

def _add_node(graph, taxNode, colorMap, nameKey, colorKey):
    graph.add_node(taxNode)
    name=taxNode.name
    graph.node[taxNode][nameKey]=nameCleanRE.sub('',name)

    if colorMap is not None:
        colorMapKey=None
        for key in colorMap.keys():
            if key.isAncestorOf(taxNode):
                if colorMapKey is None:
                    colorMapKey=key
                else:
                    if colorMapKey.isAncestorOf(key):
                        colorMapKey=key
        graph.node[taxNode][colorKey]=colorMap[colorMapKey]

def _nodeCountGenerator(root,counts):
    for node in treeGenerator(root, kidsFirst=True, key=lambda n: n.name):
        if node in counts:
            yield (node, counts[node])

def plotSunburst(tree, root=None, total=None, sort=None, **kwargs):
    """
    Take a tree (as NetworkX Directed Graph) and plot as sunburst

    keyword arguments:
        root: node to put at center (default: assume nodes are taxon.TaxNodes)
        total: sum of all values (default: walk the tree first to calculate)
        sort: key (in node dicts) to sort on within each node
        figsize: passed to matplotlib.figure()
        exterior_labels: True: (default) labels outside outrer ring, False: all labels over figure
        polar: True: (default) 'sunburst', False: 'Icicle'
        value: key in node dicts to get value from (default: %r)
        color: key in node dicts to get color from (default: %r)
        name: key in node dicts to get label from (default: %r)
        all others passed to matplolib.bars()
    """ % (VALUE, COLOR, NAME)
    logger.info("Generating sunburst")

    if root is None:
        root=next(iter(tree.nodes())).getRootNode()

    if total is None:
        total=sumCounts(tree, kwargs.get(VALUE,VALUE))
        logger.debug("Total value=%s" % total)

    from matplotlib.pyplot import figure
    fig = figure(figsize=kwargs.pop('figsize',(10,10)))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=kwargs.get('polar',True))
    ax.set_frame_on(False)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    exteriorLabels=kwargs.get('exterior_labels',True)
    _plotArc(tree, root, ax, total, sort, **kwargs)

    if exteriorLabels:
        ax.set_ylim(top=ax.get_ylim()[1]+1)

def _plotArc(tree, node, ax, total, sort, level=1, theta=0, usedNodes=None, **kwargs):
    nodeSum=0
    childTheta=theta
    if usedNodes is None:
        usedNodes=set()
    usedNodes.add(node)
    kids=0
    for child in getChildren(node, tree, sort):
        kids+=1
        if child in usedNodes:
            continue
        if level>30:
            raise Exception("Too deep!")
        childValue=_plotArc(tree, child, ax, total, sort, level=level+1, theta=childTheta, usedNodes=usedNodes, **kwargs)
        nodeSum+=childValue
        childTheta+=2.0*pi*childValue/total

    polar=kwargs.pop('polar',True)
    exteriorLabels=kwargs.pop('exterior_labels',True)
    valueKey=kwargs.pop(VALUE,VALUE)
    colorKey=kwargs.pop(COLOR,COLOR)
    nameKey=kwargs.pop(NAME,NAME)

    nodeSum+=tree.node[node].get(valueKey,0)
    width=2.0*pi*nodeSum/total
    color=tree.node[node].get(colorKey,'r')

    rects = ax.bar([theta,], [1,], width=width, bottom=level, color=color, **kwargs)

    # Labels
    x=theta+width/2.0
    if exteriorLabels and kids==0:
        y=level+1.5
    else:
        y=level+.5

    if polar:
        if x>pi/2.0 and x<3*pi/2.0:
            rotation=radToDeg(x+pi)
        else:
            rotation=radToDeg(x)
    else:
        rotation=90
    ax.text(x,y,tree.node[node].get(nameKey,""),ha='center',va='center',rotation=rotation)
    return nodeSum

def radToDeg(angle):
    return 180*angle/pi

def getChildren(node,tree,sortType):
    if sortType is None:
        return tree.edge[node].keys()
    else:
        return sorted(tree.edge[node].keys(),key=lambda n: tree.node[n].get(sortType,n))

def sumCounts(tree, valueKey=VALUE):
    """
    Sum all the values in a networkX tree
    """
    treesum=0
    for node, data in tree.nodes(data=True):
        treesum+=data.get(VALUE,0)
    return treesum

def convertTaxToJSON(node, counts, **kwargs):
    """
    Walk the entire taxonomic tree and build JSON tree of all lineages with any counts. Returns (total,tree) tuple.
    """

    valueKey=kwargs.get(VALUE,VALUE)
    nameKey=kwargs.get(NAME,NAME)
    kidsKey=kwargs.get(KIDS,KIDS)

    tree={}
    count=counts.get(node,0)
    tree[valueKey]=count
    tree[nameKey]=node.name
    children=[]
    for kid in node.children:
        (total,child) = convertTaxToJSON(kid, counts, **kwargs)
        if total>0:
            count+=total
            children.append(child)
    if len(children)>0:
        tree[kidsKey]=children
    return (count,tree)

def convertToJSON(nxTree,node):
    """
    Take a NetworkX tree and convert to dicts and lists for easy export to JSON (via repr())
    """
    tree=nxTree.node[node]
    logger.debug("NX->JSON: %s(%d)" % (tree[NAME],tree[VALUE]))
    children=[]
    for child in nxTree.edge[node].keys():
        children.append(convertToJSON(nxTree,child))
    if len(children)>0:
        tree[KIDS]=children
    #logger.debug(json.dumps(tree,indent=1))
    return tree

def plotSunburstJSON(tree, total=None, sort=None, **kwargs):
    """
    Take in nested dicts/lists where each node is defined by a dict and connection are defined by the list in the key: CHILDREN.
    optional arguments:
        keys in JSON tree:
            %r: Value to plot (default to %r, value default to 0 if not in node dict)
            %r: Color of node (default to %r, color defaults to 'r' if not in node dict)
            %r: Children of node (default to %r)
        total: total value for sizing segments (if missing, tree is crawled to calculate)
        sort: key in node dicts to sort on
        figsize: passed to matplotlib.figure to set figure size
        exterior_labels: True: (default) labels outside outrer ring, False: all labels over figure
        polar: True: (default) 'sunburst', False: 'Icicle'
        all others passed to matplolib.Axes.bars()
    """ % (VALUE,VALUE,COLOR,COLOR,KIDS,KIDS)
    logger.info("Genreating sunburst from JSON")

    # Calculate total if not given
    if total is None:
        total=sumCountsJSON(tree,**kwargs)
        logger.debug("Tree total=%s" % total)

    # create figure
    from matplotlib.pyplot import figure
    fig = figure(figsize=kwargs.pop('figsize',(10,10)))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=kwargs.get('polar',True))
    ax.set_frame_on(False)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    exteriorLabels=kwargs.get('exterior_labels',True)

    # recursive plotting function
    plotArcJSON(tree, ax, total, sort, **kwargs)

    # adjust figure range to accomodate labels outside bars
    if exteriorLabels:
        ax.set_ylim(top=ax.get_ylim()[1]+1)

def _getJSONChildren(tree,kidsKey,sort):
    children=tree.get(kidsKey,[])
    if sort is not None:
        if isinstance(sort,list) or isinstance(sort,tuple):
            if len(sort)==2 and isinstance(sort[1],bool):
                reverse=sort[1]
                if isinstance(sort[0],tuple) or isinstance(sort[0],list):
                    keyList=sort[0]
                else:
                    keyList=(sort[0],)
            else:
                keyList=sort
                reverse=False
        else:
            keyList=(sort,)
            reverse=False

        sortKey = lambda d: tuple(d.get(k,None) for k in keyList)
        children.sort(key=sortKey)
    return children

def plotArcJSON(tree, ax, total, sort, level=1, theta=0, usedNodes=None, **kwargs):
    nodeSum=0
    childTheta=theta
    kids=0

    for child in _getJSONChildren(tree,kwargs.get(KIDS,KIDS),sort):
        kids+=1
        if level>40:
            # Watch out for endless recursion
            raise Exception("Too deep!")
        childValue=plotArcJSON(child, ax, total, sort, level=level+1, theta=childTheta, usedNodes=usedNodes, **kwargs)
        nodeSum+=childValue
        childTheta+=2.0*pi*childValue/total

    polar=kwargs.pop('polar',True)
    exteriorLabels=kwargs.pop('exterior_labels',True)
    valueKey=kwargs.pop(VALUE,VALUE)
    colorKey=kwargs.pop(COLOR,COLOR)
    nameKey=kwargs.pop(NAME,NAME)
    kidsKey=kwargs.pop(KIDS,KIDS)

    nodeValue=tree.get(valueKey,0)
    nodeSum+=nodeValue
    width=2.0*pi*nodeSum/total
    color=tree.get(colorKey,'r')

    # plot the bar (arc) for this node
    rects = ax.bar([theta,], [1,], width=width, bottom=level, color=color, **kwargs)

    # Labels
    if nodeValue>0 or kids>1:
        x=theta+width/2.0
        if exteriorLabels and kids==0:
            y=level+1.5
        else:
            y=level+.5

        if polar:
            if x>pi/2.0 and x<3*pi/2.0:
                rotation=radToDeg(x+pi)
            else:
                rotation=radToDeg(x)
        else:
            rotation=90
        ax.text(x,y,tree.get(NAME,""),ha='center',va='center',rotation=rotation)
    return nodeSum

def sumCountsJSON(tree,**kwargs):
    treesum=0
    treesum+=tree.get(kwargs.get(VALUE,VALUE),0)
    for kid in tree.get(kwargs.get(KIDS,KIDS),[]):
        treesum+=sumCountsJSON(kid,**kwargs)
    #print ("total=%s" % treesum)
    return treesum

def applyCutoff(tree, cutoff, **kwargs):
    """
    Walk the tree and apply a cutoff. Moving values from underrepresented nodes into their parents.
    Return the sum of all nodes.
    """
    # cutoff assumed to be a minimum value unless it's under 1
    if cutoff<1:
        total = sumCountsJSON(tree,**kwargs)
        if cutoff<=0:
            # if cutoff is less than or equal to zero, it has no effect
            return total

        #if cutoff between 1 and 0, multiply by total to get absolute cutoff
        cutoff = cutoff * total

    return cutoffTree(tree, cutoff, **kwargs)

def cutoffTree(tree, cutoff, **kwargs):
    valueKey=kwargs.get(VALUE,VALUE)
    kidsKey=kwargs.get(KIDS,KIDS)
    treeSum=tree.get(valueKey,0)
    childrenToKeep=[]

    for child in tree.get(kidsKey,[]):
        # get the sum for each child
        childSum = cutoffTree(child, cutoff, **kwargs)
        treeSum+=childSum
        if childSum>=cutoff:
            # keep if over cutoff
            childrenToKeep.append(child)
        else:
            # add sum to value if under cutoff
            tree[valueKey]=tree.get(valueKey,0)+childSum

    # replace kids list with updated list
    tree[kidsKey]=childrenToKeep
    return treeSum

def test(projectDir='.'):
    logging.basicConfig(logLevel=logging.DEBUG)
    import os,sys
    taxDir = os.path.sep.join([projectDir,"test","data"])
    import edl.taxon
    global taxonomy
    taxonomy=edl.taxon.readTaxonomy(taxDir)
    global myAssertEq, myAssertIs
    from edl.test import myAssertEq, myAssertIs
    logging.info("Starting Tests")
    
    logging.warn("Testing NX conversion")
    testNXConversion()
    logging.warn("Testing JSON conversion")
    testNXJSON()

def testNXConversion():
    # check to see if sum stays the same
    # small test
    testCounts={}
    logger.debug("TEST: small nx")
    testCounts[taxonomy.idMap[1291540]]=4
    testCounts[taxonomy.idMap[1236689]]=10
    usum=0
    for count in testCounts.values():
        usum+=count
    logger.debug("..without species")
    (nxtreeo,rooto) = convertToNx(testCounts,ranks=['superkingdom','phylum','order'])
    myAssertEq(usum,nxSum(nxtreeo))
    logger.debug("..with species")
    (nxtreeos,rootos) = convertToNx(testCounts,ranks=['superkingdom','phylum','order','species'])
    myAssertEq(usum,nxSum(nxtreeos))
    # slightly larger test
    logger.debug("TEST: medium nx")
    testCounts={}
    testCounts[taxonomy.idMap[1291540]]=4
    testCounts[taxonomy.idMap[1236689]]=10
    testCounts[taxonomy.idMap[439481]]=26
    testCounts[taxonomy.idMap[115531]]=3
    testCounts[taxonomy.idMap[379547]]=10
    usum=0
    for count in testCounts.values():
        usum+=count
    (nxtreeo,rooto) = convertToNx(testCounts,ranks=['superkingdom','phylum','order'])
    (nxtreeos,rootos) = convertToNx(testCounts,ranks=['superkingdom','phylum','order','species'])
    myAssertEq(usum,nxSum(nxtreeos))
    myAssertEq(usum,nxSum(nxtreeo))
    # a little larger:
    logger.debug("TEST: large nx")
    testCounts = {taxonomy.idMap[115531]: 3,
                  taxonomy.idMap[379547]: 10,
                  taxonomy.idMap[439481]: 26,
         taxonomy.idMap[673860]: 36,
         taxonomy.idMap[913335]: 73,
         taxonomy.idMap[913336]: 89,
         taxonomy.idMap[1032475]: 520,
         taxonomy.idMap[1130317]: 2,
         taxonomy.idMap[1130318]: 4,
         taxonomy.idMap[1131268]: 290,
         taxonomy.idMap[1236689]: 10,
         taxonomy.idMap[1291540]: 4}
    usum=0
    for count in testCounts.values():
        usum+=count
    (nxtreeo,rooto) = convertToNx(testCounts,ranks=['superkingdom','phylum','order'])
    (nxtreeos,rootos) = convertToNx(testCounts,ranks=['superkingdom','phylum','order','species'])
    myAssertEq(usum,nxSum(nxtreeos))
    myAssertEq(usum,nxSum(nxtreeo))

def testNXJSON():
    # check to see if sum stays the same
    # small test
    logger.debug("TEST: small NXJS")
    testCounts={}
    testCounts[taxonomy.idMap[1291540]]=4
    testCounts[taxonomy.idMap[1236689]]=10
    testCounts[taxonomy.idMap[1032475]]=520
    usum=0
    for count in testCounts.values():
        usum+=count
    (nxtreeo,rooto) = convertToNx(testCounts,ranks=['superkingdom','phylum','order'])
    (nxtreeos,rootos) = convertToNx(testCounts,ranks=['superkingdom','phylum','order','species'])
    jtreeo = convertToJSON(nxtreeo,rooto)
    jtreeos = convertToJSON(nxtreeos,rootos)
    myAssertEq(usum,jSum(jtreeos))
    myAssertEq(usum,jSum(jtreeo))
    # slightly larger test
    testCounts={}
    testCounts[taxonomy.idMap[1291540]]=4
    testCounts[taxonomy.idMap[1236689]]=10
    testCounts[taxonomy.idMap[439481]]=26
    testCounts[taxonomy.idMap[115531]]=3
    testCounts[taxonomy.idMap[379547]]=10
    usum=0
    for count in testCounts.values():
        usum+=count
    (nxtreeo,rooto) = convertToNx(testCounts,ranks=['superkingdom','phylum','order'])
    (nxtreeos,rootos) = convertToNx(testCounts,ranks=['superkingdom','phylum','order','species'])
    jtreeo = convertToJSON(nxtreeo,rooto)
    jtreeos = convertToJSON(nxtreeos,rootos)
    myAssertEq(usum,jSum(jtreeos))
    myAssertEq(usum,jSum(jtreeo))
    # a little larger:e
    testCounts = {taxonomy.idMap[115531]: 3,
                  taxonomy.idMap[379547]: 10,
                  taxonomy.idMap[439481]: 26,
         taxonomy.idMap[673860]: 36,
         taxonomy.idMap[913335]: 73,
         taxonomy.idMap[913336]: 89,
         taxonomy.idMap[1032475]: 520,
         taxonomy.idMap[1130317]: 2,
         taxonomy.idMap[1130318]: 4,
         taxonomy.idMap[1131268]: 290,
         taxonomy.idMap[1236689]: 10,
         taxonomy.idMap[1291540]: 4}
    usum=0
    for count in testCounts.values():
        usum+=count
    (nxtreeo,rooto) = convertToNx(testCounts,ranks=['superkingdom','phylum','order'])
    (nxtreeos,rootos) = convertToNx(testCounts,ranks=['superkingdom','phylum','order','species'])
    jtreeo = convertToJSON(nxtreeo,rooto)
    jtreeos = convertToJSON(nxtreeos,rootos)
    myAssertEq(usum,jSum(jtreeos))
    myAssertEq(usum,jSum(jtreeo))

def nxSum(nxtree):
    tsum=0
    for n,data in nxtree.nodes(data=True):
        tsum+=data['size']
    return tsum

def jSum(jtree):
    tsum=jtree['size']
    if 'children' in jtree:
        for kid in jtree['children']:
            tsum+=jSum(kid)
    return tsum
 
if __name__=='__main__':
    test()

