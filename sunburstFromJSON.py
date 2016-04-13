#!/usr/bin/python
"""
Take in a JSON tree and create a suburst PNG
"""
from optparse import OptionParser
import sys, re, json, os
import matplotlib
from edl.util import addIOOptions, addUniversalOptions, setupLogging
from edl.sunburst import *

ID='id'

def main():
    usage = "usage: %prog OPTIONS JSON_FILE(s)"
    description = """
    Generates a sunburst plot for each input JSON tree.
    """
    parser = OptionParser(usage, description=description)
    addIOOptions(parser)
    addUniversalOptions(parser)

    parser.add_option('-r', "--root", default=None, help="Plot a subset of the tree by choosing a root node for the subtree")
    parser.add_option('-c', "--colors", default=None, help="Set colors by mapping node IDs to color strings. Value should be a comma-separated list of id=color pairs (Bacteria=g,Archaea=r). The subtree of each mapped node will get the given color unless overridden by another entry. If omitted, colors pulled from JSON (using colorkey) with red as the default. If present without --colorkey setting, colors in JSON will be ignored.")
    parser.add_option('-s','--sort', default=None,
                      help="List of keys to sort on for plotting, NOTE: sorting on the value key will give suprising results for lower level nodes as sum of nested values will not be included. To get desired behavior, add a total value key to your tree and sort on that.")

    parser.add_option('-I','--idkey', default='name', help="String to use as key for node IDs. Default: %default")
    parser.add_option('-L','--labelkey', default='name', help="String to use as key for node labels. Default: %default")
    parser.add_option('-C','--colorkey', default='color', help="String to use as key for node colors. Default: %default")
    parser.add_option('-V','--valuekey', default='size', help="String to use as key for node sizes. Default: %default")
    parser.add_option('-K','--kidskey', default='children', help="String to use as key for list of child nodes. Default: %default")

    parser.add_option('-i', '--icicle', default=False, action='store_true',
                      help="Print stacked bars in rectangular coordinates, not polar.")
    parser.add_option('-e', '--exterior_labels', default=False, action='store_true', help="Print labels for outermost nodes outside image")
    parser.add_option('-S', '--figsize', default=None,
                      help="Comma separated pair of numbers (in inches) for figure size")

    parser.add_option("-f", "--format", dest="format", default='pdf', choices=['png','ps','pdf','svg'],
		  help="Format for output image", metavar="FORMAT")

    (options, args) = parser.parse_args()

    # check arguments
    setupLogging(options, description)

    # setup matplotlib
    backend = options.format
    if backend=='png':
        backend='agg'
    matplotlib.use(backend)
    import matplotlib.pyplot as plt

    for (inhandle, outfile) in inputIterator(args, options):
        # import JSON
        tree=json.load(inhandle)

        # proecss user selected options
        kwargs=processOptions(options)

        # process JSON
        if options.colors is not None:
            setColors(tree, options.colors, **kwargs)
        if options.root is not None:
            newRoot=findNode(tree, options.root, **kwargs)
            if newRoot is not None:
                tree=newRoot

        # some of the matplotlib functions don't like extra arguments
        kwargs.pop(ID)

        # create figure
        plotSunburstJSON(tree,**kwargs)

        # save to file
        plt.savefig(outfile)

#############
# Functions #
#############
def processOptions(options):
    """
    Create kwargs object from options
    """
    kwargs={}
    if options.sort is not None:
        kwargs['sort']=options.sort
    if options.figsize is not None:
        figsize = [ int(bit) for bit in options.figsize.split(",") ]
        kwargs['figsize']=tuple(figsize[:2])
    kwargs[ID]=options.idkey
    kwargs[NAME]=options.labelkey
    kwargs[COLOR]=options.colorkey
    kwargs[VALUE]=options.valuekey
    kwargs[KIDS]=options.kidskey
    kwargs['polar']=not options.icicle
    kwargs['exterior_labels']=options.exterior_labels

    return kwargs

def setColors(tree, colorString, **kwargs):
    """
    Given the color mapping string, add color entries to nodes of tree
    """
    idKey=kwargs.get(ID,ID)
    colorkey=kwargs.get(COLOR,COLOR)
    kidskey=kwargs.get(KIDS,KIDS)

    logger.debug("Color String: %s" % colorString)
    colorMap={}
    for pair in colorString.split(','):
        (key,value) = pair.split('=')
        logger.debug("Coloring %s: %s" % (key,value))
        colorMap[key]=value

    _setColor(tree, colorMap, idKey, colorkey, kidskey, 'r')

def _setColor(tree, colorMap, idkey, colorkey, kidskey, default):
    color=colorMap.get(tree.get(idkey,None),default)
    tree[colorkey]=color
    for kid in tree.get(kidskey,[]):
        _setColor(kid, colorMap, idkey,colorkey, kidskey, color)

def findNode(tree, nodeId, **kwargs):
    if tree.get(kwargs.get(ID,ID),None) == nodeId:
        return tree
    for kid in tree.get(kwargs.get(KIDS,KIDS),[]):
        node=findNode(kid, nodeId, **kwargs)
        if findNode(kid, nodeId, **kwargs) is not None:
            return node
    return None

def inputIterator(infileNames, options):
    """
    take list of input files from infileName (which can be None) and infileNames (which is a list of any length)
    if no file names give, read from stdin
    if outfile is not given, add extension to filename
    if multiple infiles, treat outfile as suffix
    if cwd set to False (for multipleinpus) create output files in same dir as inputs. Otherwise, create files in current dir withinput names plus suffix
    """
    outfileName=options.outfile
    if 'infile' in dir(options) and options.infile is not None:
        infileNames.append(options.infile)
    if len(infileNames)==0:
        inhandle=sys.stdin
        if outfileName==None:
            outfileName="sunburst." + options.format
            logger.info("IO: STDIN -> %s" % outFileName)
            yield (inhandle,outfileName)
        else:
            logger.info("IO: STDIN -> %s" % outfileName)
            yield(inhandle,outfileName)
    elif len(infileNames)==1:
        infileName=infileNames[0]
        inhandle=open(infileName)
        if outfileName==None:
            # use format as suffix
            if options.cwd:
                # strip path info first
                (infilePath,infileName)=os.path.split(infileName)
                infileName="./"+infileName
            outfileName="%s.%s" % (infileName,options.format)
            logger.info("IO: %s -> %s" % (infileName,outfileName))
            yield (inhandle,outfileName)
        else:
            logger.info("IO: %s -> %s" % (infileName,outfileName))
            yield(inhandle,outfileName)
        inhandle.close()
    else:
        if outfileName==None:
            outfileName="."+options.format
        for infileName in infileNames:
            inhandle=open(infileName)
            # use outfileName as suffix
            if options.cwd:
                # strip path info first
                (infilePath,infileName)=os.path.split(infileName)
                infileName="./"+infileName
            thisoutfile="%s%s" % (infileName,outfileName)
            logger.info("IO: %s -> %s" % (infileName,thisoutfile))
            yield (inhandle,thisoutfile)
            inhandle.close()

if __name__ == '__main__':
    main()
