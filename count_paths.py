#! /usr/bin/python
"""
    Takes hit tables (eg from blast, last, or hmmer) and generates a table (or tables) of hit counts for gene families,
gene classes, or pathways. Multiple pathway or classification levels can be given. If they are, multiple output tables 
will be generated with the level as a suffix.

    This script will work with KEGG, SEED, or CAZy. CAZy only has one level of heirarchy, the others have 3. The CAZy 
heirarchy is apparent from the hit name and needs no supporting files. KEGG and SEED require mapping files to identify
gene families and heirachy files to report levels other than the gene family or ortholog level. Both SEED and KEGG
have three levels of classifications that can be indicated with a 1, 2, or 3. The words "subsystem" and "pathway" are
synonyms for level 3.

    The 'topHit' count method is applied so the any ambiguous hits are given to the hit that is most common (in all
samples). In the end each read is mapped to a single hit which may or may not be part of the gene heirarchy. 
    
    NOTE: in KEGG (and SEED) a single ortholog (role) may belong to multiple pathways (subsystems). A hit to such an
ortholog will result in a +1 for the count values of ALL the pathways it belongs to. Therefore, the pathway counts
will almost always add up to more than the total number of reads.
"""
import sys, re, logging, argparse
from urllib.parse import unquote_plus
from edl import redistribute, kegg
from edl.hits import add_weight_arguments, loadSequenceWeights, add_count_arguments, GIS, FilterParams, getHitTranslator
from edl.util import add_universal_arguments, setup_logging, parseMapFile
from edl.expressions import accessionRE, nrOrgRE

def main():
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("input_files", nargs="*", 
                        type=argparse.FileType('r'),
                        default=[sys.stdin,],
                        metavar="INFILE",
                        help="List of input files. Omit to read from STDIN")
    parser.add_argument("-o", "--outfile", dest="output_file",
                      metavar="OUTFILE", help="Write count table to OUTFILE")
    parser.add_argument("-l", "--level", dest="levels", default=None,
                      metavar="LEVEL", action="append",
                      help=""" Level(s) to collect counts on. Use flag 
                      multiple times to specify multiple levels. If multiple 
                      values given, one table produced for each with rank 
                      name appended to file name. Levels can be an integer 
                      (1-3) for KEGG or SEED levels, any one of 'gene', 'role', 'family', 
                      'ko', or 'ortholog' (which are all synonyms), or  
                      anything not synonymous with 'gene' to 
                      get CAZy groups. Defaults to ortholog/role and 
                      levels 1, 2, and 3 for KEGG and SEED
                      and gene and group for CAZy and COG.""")

    # option for deconvoluting clusters or assemblies
    add_weight_arguments(parser, multiple=True)

    # cutoff options
    add_count_arguments(parser)

    # format, ortholog heirarchy, and more
    kegg.add_path_arguments(parser,defaults={'countMethod':'tophit'},
                                   choices={'countMethod':('tophit',)},
                                   helps={'countMethod':"Currently the only supported method is tophit: \
                                                         This breaks any ties by choosing the most abundant hit \
                                                         based on other unambiguous assignments."})

    # log level and help
    add_universal_arguments(parser)
    arguments = parser.parse_args()
    setup_logging(arguments)

    if len(arguments.input_files)==0:
        parser.error("Must supply at least one m8 file to parse")

    # Set defaults and check for some conflicts
    if arguments.levels is None and arguments.heirarchyFile is None:
        # using hit names only
        arguments.levels=[None]
    else:
        if arguments.heirarchyFile is None and arguments.heirarchyType != 'cazy':
            logging.warn("Type: %s" % (arguments.heirarchyType))
            parser.error("Cannot select levels without a heirarchy (ko) file")
        if arguments.levels is None:
            # set a default
            if arguments.heirarchyType is 'kegg':
                arguments.levels=['ko','1','2','pathway']
            if arguments.heirarchyType is 'seed':
                arguments.levels=['role','1','2','subsystem']
            else:
                arguments.levels=['gene','group']

        try:
            # Make sure the rank lists make sense
            arguments.levels=cleanLevels(arguments.levels)
        except Exception as e:
            parser.error(str(e))

    # load weights file
    sequenceWeights = loadSequenceWeights(arguments.weights)

    # only print to stdout if there is a single level
    if len(arguments.levels)>1 and arguments.output_file is None:
        parser.error("STDOUT only works if a single level is chosen!")

    cutoff=arguments.cutoff

    # map reads to hits
    if arguments.mapFile is not None:
        if arguments.mapStyle == 'auto':
            with open(arguments.mapFile) as f:
                firstLine=next(f)
                while len(firstLine)==0 or firstLine[0]=='#':
                    firstLine=next(f)
            if koMapRE.search(firstLine):
                arguments.mapStyle='kegg'
            elif seedMapRE.search(firstLine):
                arguments.mapStyle='seed'
            elif tabMapRE.search(firstLine):
                arguments.mapStyle='tab'
            #elif cogMapRE.search(firstLine):
            #    arguments.mapStyle='cog'
            else:
                raise Exception("Cannot figure out map type from first line:\n%s" % (firstLine))

        logging.info("Map file seems to be: %s" % (arguments.mapStyle))
        if arguments.mapStyle=='kegg':
            valueMap=kegg.parseLinkFile(arguments.mapFile)
        elif arguments.mapStyle=='seed':
            valueMap=kegg.parseSeedMap(arguments.mapFile)
        #elif arguments.mapStyle=='cog':
        #    valueMap=kegg.parseCogMap(arguments.mapFile)
        else:
            if arguments.parseStyle == GIS:
                keyType=int
            else:
                keyType=None
            valueMap = parseMapFile(arguments.mapFile,valueType=None,keyType=keyType)
        if len(valueMap)>0:
            logging.info("Read %d items into map. EG: %s" % (len(valueMap),next(iter(valueMap.items()))))
        else:
            logging.warn("Read 0 items into value map!")
    else:
        valueMap=None

    # parse input files
    fileCounts={}
    totals={}
    fileLabels={}
    sortedLabels=[]

    # Allow for file names to be preceded with TAG=
    for file_handle in arguments.input_files:
        filename=file_handle.name
        bits=filename.split("=",1)
        if len(bits) > 1:
            (filetag,filename)=bits
        else:
            filetag=filename
        fileLabels[filename]=filetag
        # keep order so that column order matches arguments
        sortedLabels.append(filetag)
        fileCounts[filetag]={}
        totals[filetag]=0

    #TODO: incorporate weights into tophit algorithm!
    if arguments.countMethod == 'tophit':
        # Process all files at once and use overall abundance to pick best hits
        from edl import redistribute
        params = FilterParams.create_from_arguments(arguments)
        multifile = redistribute.multipleFileWrapper(fileLabels.keys())

        # don't give any hit translation, just use hit ids for redistribution
        readHits = redistribute.pickBestHitByAbundance(multifile,
                                              filterParams=params,
                                              returnLines=False,
                                              winnerTakeAll=True,
                                              parseStyle=arguments.parseStyle,
                                              sequenceWeights=sequenceWeights)
        # define method to turn Hits into Genes (kos, families)
        hitTranslator=getHitTranslator(parseStyle=arguments.parseStyle,
                                       hitStringMap=valueMap)
        translateHit = lambda hit: hitTranslator.translateHit(hit)[0]

        # use read->file mapping and hit translator to get file based counts
        #  from returned (read,Hit) pairs
        increment=1
        for (read_name, hit) in readHits:
            file_tag, read_name = read_name.split("/",1)
            file_tag = unquote_plus(file_tag)
            gene = translateHit(hit)
            if gene is None:
                gene = "None"
            logging.debug("READ: %s\t%s\t%s\t%s" % (file_tag,read_name,hit.hit,gene))
            genecount = fileCounts[file_tag].setdefault(gene,0)
            if sequenceWeights is not None:
                increment = sequenceWeights.get(read_name,1)
            fileCounts[file_tag][gene]=genecount+increment
            totals[file_tag]+=increment
        logging.debug(str(totals))

    else:
        # TODO: support for other methods
        parser.error("Only tophit counting is supported at the moment. Sorry.")

    logging.debug(repr(fileCounts))
    printCountTablesByLevel(fileCounts, totals, sortedLabels, arguments)

def cleanLevels(levelList):
    # don't allow duplicates
    levelList=list(set(levelList))

    # return levels
    return levelList

def printCountTablesByLevel(fileCounts, totals, fileNames, options):
    """
    Create a new file for each level with a tab separated table of counts
    """
    cutoff=options.cutoff

    if options.heirarchyType == 'seed':
        logging.info("Reading SEED subsystem assignments from %s" % (options.heirarchyFile))
        seedTree = kegg.readSEEDTree(options.heirarchyFile)
    elif options.heirarchyType == 'cog':
        logging.info("Reading COG subsystem assignments from %s" % (options.heirarchyFile))
        seedTree = kegg.readCogTree(options.heirarchyFile)

    # create an output table for each requested level
    for level in options.levels:
        logging.debug("Processing level %s" % (level))
        translateToPaths = level not in koSyns
        descString = None
        if translateToPaths:
            if options.heirarchyType == 'cazy':
                geneTranslator = getCazyGroup
            else:
                lookupLevel = level if level not in level3Syns else '3'

                if options.heirarchyType == 'kegg':
                    # Ideally, we'd be able to parse the heirachy once, but the current KEGG code just retuns simple mappings
                    logging.info("Reading KEGG level %s assignments from %s" % (level,options.heirarchyFile))
                    geneTranslation = kegg.readKEGGFile(options.heirarchyFile, lookupLevel)
                else:
                    # SEED or COG/KOG
                    geneTranslation = seedTree[lookupLevel]
                geneTranslator = lambda gene: geneTranslation.get(gene,gene)

        elif level is not None and options.heirarchyType == 'kegg':
            # return descriptions if level explicitly set to ko (or syn.)
            descString = "Description"
            logging.info("Reading KO descriptions from %s" % options.heirarchyFile)
            geneTranslation = kegg.readKEGGFile(options.heirarchyFile, "DESCRIPTION")
            geneTranslator = lambda gene: "%s\t%s" % (gene,
                                             geneTranslation.get(gene,gene))
        elif level is not None and options.heirarchyType == 'cog':
            # return descriptions if level explicitly set to ko (or syn.)
            descString = "Description\tCategories"
            geneTranslator = lambda gene: "%s\t%s\t%s" % (
                                     seedTree['gene'].get(gene, gene),
                                     seedTree['description'].get(gene,"None"),
                                     seedTree['group'].get(gene,"None"))
        else:
            # just return gene if no level set or not KEGG/COG/KOG
            geneTranslator = lambda gene: gene
        

        # For each level, try to force all counts to be at that level
        fileLevelTotals={}
        levelCounts = {}
        levelPaths = {}
        thresholds = {}
        for (filename, counts) in fileCounts.items():
            fileLevelTotals[filename]=0
            thresholds[filename] = totals[filename]*cutoff
            fileLevelCounts = levelCounts.setdefault(filename, {})

            fileTotal=0
            for gene in sorted(counts.keys()):
                # get the counts from this node
                geneCount = counts[gene]
                fileTotal+=geneCount

                # translate gene to pathway (or not depending on above code)
                pathway = geneTranslator(gene)

                # update counts
                # Some KOs will map to multiple pathways,
                #  so... allow for multiple translated values
                if not(isinstance(pathway, list) or isinstance(pathway, tuple)):
                    pathway = [pathway,]
                for indPathway in pathway:
                    fileLevelCounts[indPathway] = fileLevelCounts.get(indPathway,0) + geneCount
                    levelPaths[indPathway]=True

            logging.debug("File %s has %d hits (had %d)" % (filename, fileTotal, totals[filename]))

        #logging.debug(repr(levelPaths))
        #logging.debug(repr(levelCounts))
        if logging.getLogger().level <= logging.DEBUG:
            for (filename,counts) in levelCounts.items():
                logging.debug("File %s hs %d counts" % (filename,sum(counts.values())))

        # apply cutoff
        for pathway in list(levelPaths.keys()):
            # check to see if pathway is over cutoff in any file
            over=False
            for (filename,fileLevelCount) in levelCounts.items():
                flPathCount=fileLevelCount.get(pathway,0)
                fileLevelTotals[filename]+=flPathCount
                if flPathCount > thresholds[filename]:
                    over=True
            if not over:
                # this pathway is not over the cutoff for any file
                levelPaths.pop(pathway)
                other='Other'
                levelPaths[other]=True
                for (filename, fileLevelCount) in levelCounts.items():
                    fileLevelCount[other]=fileLevelCount.get(other,0) + fileLevelCount.pop(pathway,0)

        if logging.getLogger().level <= logging.DEBUG:
            for (filename,counts) in levelCounts.items():
                logging.debug("File %s has %d counts" % (filename,sum(counts.values())))
                missed=False
                for path in counts.keys():
                    if path not in levelPaths:
                        missed=True
                        logging.debug("Missing pathway %s has %d counts for %s" % (path, counts[path], filename))
                if not missed:
                    logging.debug("There are no missing pathways from %s" % (filename))

        logging.debug("Final file counts: %r" % (fileLevelTotals))

        # output file
        if options.output_file is None:
            outs = sys.stdout
        else:
            if len(options.levels)>1:
                outfile = "%s.%s" % (options.output_file, level)
            else:
                outfile = options.output_file
            outs = open(outfile,'w')

        # write to file(s?)
        # header
        if level in koSyns:
            # Header for when level is the gene
            if descString is not None:
                outs.write("Gene\t%s\t%s\n" % (descString,
                                               '\t'.join(fileNames)))
            else:
                outs.write("Gene\t%s\n" % ('\t'.join(fileNames)))
        else:
            # Header for when level is a pathway or group
            outs.write("Pathway\t%s\n" % ('\t'.join(fileNames)))

        for pathway in sorted(levelPaths.keys()):
            outs.write(str(pathway))
            for filename in fileNames:
                outs.write("\t")
                outs.write(str(levelCounts[filename].get(pathway,0)))
            outs.write("\n")

        # close out stream
        if options.output_file is not None:
            outs.close()

def getCazyGroup(gene):
    m=cazyRE.search(gene)
    if m:
        cazygroup = m.group(1)
        logging.debug("Mapping %s to %s" % (gene, cazygroup))
    else:
        logging.debug("Could not parse group from %s with r'%s'" % (gene, cazyRE.pattern))
        cazygroup=gene
    return cazygroup

# levels equivalent to returning just the ko/gene
koSyns = [None, 'ko', 'gene', 'ortholog', 'family', 'role']
level3Syns = ['subsystem','pathway','group']
koMapRE = re.compile(r'\sko:K\d{5}')
seedMapRE = re.compile(r'^Mapped roles:')
cogMapRE = re.compile(r'^\d+\t[KC]OG\d+\t')
tabMapRE = re.compile(r'^[^\t]+\t[^\t+]')
cazyRE = re.compile(r'([a-zA-Z]+)\d+')

if __name__ == '__main__':
    main()

