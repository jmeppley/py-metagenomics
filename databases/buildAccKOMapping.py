#!/usr/bin/python
"""
Given a list of fasta files, and the KEGG gene links directory, build a map from accession to KO

Step 1: read in map from KeggGene IDs to KOs

Step 2: read map from KeggGene IDs to GIs:
    keep only mappings from GIs to KOs

Step 3: read in Accession->GIlists from fasta headers:
    Build map from accessions to KOs and delete GIs from GI->KO map as we go
"""

from optparse import OptionParser
import sys, re
from edl.util import *

def main():
    usage = "usage: %prog [OPTIONS] BLAST_M8_FILE[S]"
    description = __doc__
    parser = OptionParser(usage, description=description)
    parser.add_option("-o", "--outfile", dest="outfile",
                      metavar="OUTFILE", help="Write count table to OUTFILE")
    parser.add_option("-d", "--gisInDesc", default=False, action="store_true", help="Look for GIs in description")
    parser.add_option("-l", "--keggGeneLinkDir", dest="keggGeneLinkDir", default="/common/FASTA/KEGG/LATEST/links",
                      metavar="DIR",
                      help="Directory containing lists of kegg gene ids mapped to other identifiers. Defaults to: /common/FASTA/KEGG/LATEST/links")

    # log level and help
    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options, description)

    # map kos to gis
    gi2kos = processKeggGeneLinks(options.keggGeneLinkDir)
    logging.info("Parsed ko map for %d gis" % len(gi2kos))
    if logging.getLogger().level <= logging.DEBUG:
        logging.debug(repr(["%s:\t%s" % (key,gi2kos[key]) for key in gi2kos.keys()[:10]]))

    # pull acc/gi associations from FASTA and write out the acc->ko map
    for (inhandle, outhandle) in inputIterator(args, options):
        logging.info("Pulling acc->gi map from fasta headers in %s" % (inhandle))
        logging.info("Writing acc->ko to %s" % (outhandle))
        printAcc2KOMap(inhandle, outhandle, gi2kos, options.gisInDesc)

## END MAIN ##

giRE=re.compile(r'gi\|(\d+)\|')
rsGIRE=re.compile(r'^>(gi\|(\d+)\S+\s.*)$')
accessionRE = re.compile(r'(?:edl|dbj|emb|gb|pdb|pir|prf|ref|sp|tpd|tpe|tpg)\|{1,2}([-0-9a-zA-Z_]+)\.\d{1,2}(?:\||$)')
ncbiGeneIdRE = re.compile(r'^(\S+)\s+ncbi-gi:(\d+)')

def processKeggGeneLinks(linkdir):
    gid2kos={}
    # get gene to ko map
    kg2koFile=linkdir + os.path.sep + "genes_ko.list"
    with open(kg2koFile) as f:
        for line in f:
            (gid,ko) = line.split()[:2]
            gid2kos.setdefault(gid,[]).append(ko[3:])

    gi2kos={}
    # merge with gene to gi map
    kg2giFile=linkdir + os.path.sep + "genes_ncbi-gi.list"
    with open(kg2giFile) as f:
        for line in f:
            #(gid,gi) = line.split()[:2]
            try:
                (gid,gi) = ncbiGeneIdRE.match(line).groups()
            except:
                logging.error(line)
                logging.error(ncbiGeneIdRE.pattern)
                raise
            try:
                kos=gid2kos.pop(gid)
            except KeyError:
                continue
            gi2kos.setdefault(int(gi),[]).extend(kos)

    gid2kos.clear()
    return gi2kos

def printAcc2KOMap(inhandle, outhandle, gi2kos, giInDesc):
    """
    Parse accession from  ID and all GIs from ID and description in input FASTA
    Print map from accessions to KOs (removing GIs form gi2ko as we go)
    """
    for line in inhandle:
        # only consider header lines
        if len(line)<1 or line[0]!='>':
            continue

        logging.debug(line)
        (header,description)=line.split(None,1)
        m=accessionRE.search(header)
        if m is None:
            acc=header
            logging.warn("No accession in %s" % (header))
        else:
            acc=m.group(1)
        logging.debug("ACC: %s" % (acc))
        kos=set()
        if giInDesc:
            gistring=description
        else:
            gistring=header
        for gi in giRE.findall(gistring):
            logging.debug("GI: %s" % (gi))
            try:
                for ko in gi2kos.pop(int(gi)):
                    kos.add(ko)
            except:
                continue
        if len(kos)>0:
            outhandle.write("%s\t%s\n" % (acc,",".join(kos)))

if __name__ == '__main__':
    main()
