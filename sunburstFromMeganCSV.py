#!/usr/bin/python
"""
Take in a csv of taxids and counts and generate a sunburst plot.
"""
from optparse import OptionParser
import sys, json
import matplotlib
from edl.util import addIOOptions, addUniversalOptions, setupLogging, inputIterator as inputIteratorNormal
from edl.taxon import readTaxonomy
from edl.sunburst import *
from sunburstFromJSON import setColors, findNode, inputIterator as inputIteratorFig, ID

def main():
    usage = "usage: %prog OPTIONS JSON_FILE(s)"
    description = __doc__
    parser = OptionParser(usage, description=description)
    addIOOptions(parser)
    addUniversalOptions(parser)

    parser.add_option('-r', "--root", default=None, help="Plot a subset of the tree by choosing a root node for the subtree")
    parser.add_option('-c', "--colors", default=None, help="Set colors by mapping taxon names to color strings. Value should be a comma-separated list of id=color pairs (Bacteria=g,Archaea=r). The subtree of each mapped node will get the given color unless overridden by another entry.")
    parser.add_option('-C', "--cutoff", default=.025, type='float',
                      help="Trim nodes below this value. Interpreted as an absolute threshold if >1 and as fractional if <1. Set to 0 (or less) to turn off.")

    parser.add_option("-R", "--ranks", default='superkingdom,phylum,family',
                      help="Ranks to inclued in sunburst, default: %default")
    parser.add_option("-n", "--ncbiTaxDir", dest="taxdir", metavar="PATH", default='/common/FASTA/NCBI/RefSeq/LATEST/taxdump',
                      help="Directory with unpacked ncbi tax dump (specifically names.dmp and nodes.dmp) and use to translate taxids into taxa. Default is %default")
    parser.add_option('-i', '--icicle', default=False, action='store_true',
                      help="Print stacked bars in rectangular coordinates, not polar.")
    parser.add_option('-e', '--exterior_labels', default=False, action='store_true', help="Print labels for outermost nodes outside image")
    parser.add_option('-S', '--figsize', default=None,
                      help="Comma separated pair of numbers (in inches) for figure size")

    parser.add_option("-f", "--format", dest="format", default='pdf', choices=['png','ps','pdf','svg'],
		  help="Format for output image", metavar="FORMAT")
    parser.add_option("-J", "--JSON", default=False, action="store_true",
            help="output JSON tree of counts instead of figure")

    (options, args) = parser.parse_args()

    # check arguments
    setupLogging(options, description)

    if not options.JSON:
        # setup matplotlib
        backend = options.format
        if backend=='png':
            backend='agg'
        matplotlib.use(backend)
        import matplotlib.pyplot as plt

    # load taxonomy
    taxonomy = readTaxonomy(options.taxdir)

    # build rank list
    ranks=options.ranks.split(',')

    if options.JSON:
        # STandard iterator that returns handles
        inputIterator = inputIteratorNormal
    else:
        # version that defaults to adding format as suffix and returns name
        inputIterator = inputIteratorFig

    # process input files
    for (inhandle, outfile) in inputIterator(args, options):
        # load counts
        counts={}
        for line in inhandle:
            (taxid,count)=line.rstrip('\n\r').split(',')
            if taxid=="None":
                tax=None
            else:
                tax=taxonomy.idMap.get(int(taxid),None)
            counts[tax]=counts.get(tax,0)+int(count)


        # convert to JSON
        (nxtree,root) = convertToNx(counts, leaves=True, ranks=ranks)
        tree=convertToJSON(nxtree,root)

        # proecss user selected options
        kwargs=processOptions(options)

        # process JSON
        if options.colors is not None:
            setColors(tree, options.colors, **kwargs)
        if options.root is not None:
            newRoot=findNode(tree, options.root, **kwargs)
            if newRoot is not None:
                tree=newRoot

        total=applyCutoff(tree, options.cutoff, **kwargs)
        if options.JSON:
            putNodeCountsInOther(tree)
            outfile.write(json.dumps(tree, indent=2))
        else:
            # some of the matplotlib functions don't like extra arguments
            kwargs.pop(ID)

            # create figure
            plotSunburstJSON(tree,**kwargs)

            # save to file
            plt.savefig(outfile, format=options.format)

def putNodeCountsInOther(tree, **kwargs):
    """
    D3 tree layouts ignore node counts for nodes with children (non-leaves). So we need to add a leaf called "Other" for each non-leaf node with non-zero counts
    """
    valueKey=kwargs.get(VALUE,VALUE)
    kidsKey=kwargs.get(KIDS,KIDS)
    nameKey=kwargs.get(NAME,NAME)
    colorKey=kwargs.get(COLOR,COLOR)

    kids = tree.get(kidsKey,[])
    if len(kids)>0:
        value = tree.get(valueKey,0)
        if value > 0:
            newNode={nameKey:"Others",kidsKey:[],valueKey:value}
            if colorKey in tree:
                newNode[colorKey]=tree[colorKey]
            kids.append(newNode)
        for kid in kids:
            putNodeCountsInOther(kid, **kwargs)

#############
# Functions #
#############
def processOptions(options):
    """
    Create kwargs object from options
    """
    kwargs={}
    kwargs[ID]=NAME
    kwargs['sort']="color,size"
    if options.figsize is not None:
        figsize = [ int(bit) for bit in options.figsize.split(",") ]
        kwargs['figsize']=tuple(figsize[:2])
    kwargs['polar']=not options.icicle
    kwargs['exterior_labels']=options.exterior_labels

    return kwargs

if __name__ == '__main__':
    main()

