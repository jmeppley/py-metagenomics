#!/usr/bin/python
"""
Read RefSeq catalog from stdin:
    get taxid, name, accession from first three columns
Read taxonomy from taxdump dir
For each RefSeq entry that is a protein and is not (in)vertebrate:
    if taxid not in taxonomy:
        use name to find updated taxid
    print accesion (without suffix) and taxid
"""

from optparse import OptionParser
import sys, logging, re
sys.path.append('/common/bin/scripts')
from edl.taxon import readTaxonomy, getNodeFromHit

def main():
    usage='%prog [OPTIONS] TAXDUMP_PATH'
    description='reduce full RefSea catalog (from STDIN) to acc->taxid map using TAXDUMP to verify taxids'

    parser = OptionParser(usage, description=description)
    parser.add_option("-g", "--genomic", default=False, action="store_true", help="output genoic accessions instead of proteins")
    parser.add_option("-v", "--verbose",
                      action="count", dest="verbose", default=1,
                      help="Print log messages. Use twice for debugging")
    parser.add_option("-q", '--quiet', dest='verbose',
                      action="store_const", const=0,
                      help="Suppress warnings. Only print fatal messages")
    parser.add_option("-A", "--about",
              action="store_true", dest="about", default=False,
              help="Print description")

    (options, args) = parser.parse_args()

    if options.about:
        print description
        exit(0)

    # check args
    if options.verbose==0:
        loglevel=logging.ERROR
    elif options.verbose==1:
        loglevel=logging.WARN
    elif options.verbose==2:
        loglevel=logging.INFO
    elif options.verbose>=3:
        loglevel=logging.DEBUG
    logging.basicConfig(stream=sys.stderr, level=loglevel)
    logging.info("Log level set to %r(%d)" % (loglevel,options.verbose))

    if len(args) != 1:
        parser.error("Please supply TAXDUMP path in command line")

    # parse catalog
    accToOrg={}
    accRE=re.compile(r'\b([A-Z]{2}_[A-Z]*\d+)(\.\b)?\b')
    protRE=re.compile(r'^([ANXWYZ]P_[A-Z]*\d+)$')
    logging.info("reading catalog from STDIN")
    for line in sys.stdin:
        cells=line.rstrip('\r\n').split('\t')
        (taxid,name,acc)=cells[0:3]
        try:
            taxid=int(taxid)
        except:
            pass
        logging.debug(acc)
        acc=accRE.match(acc).group(1)
        logging.debug("'%s'" % acc)
        m=protRE.match(acc)
        if (m==None)==options.genomic:
            # will match if acc matches exp and genomic is false OR
            #            if acc doesn't match and genomic is true
            accToOrg[acc]=(taxid,name)
            logging.debug("USing acc: %s" % (acc))
        else:
            logging.debug("Skipping acc: %s" % (acc))

    # load taxonomy
    logging.info("loading taxonomy from %s" % (args[0]))
    taxonomy = readTaxonomy(args[0],namesMap=True)

    # print table
    changes=0
    for (acc,(taxid,name)) in accToOrg.iteritems():
        if taxid not in taxonomy.idMap:
            node = getNodeFromHit(name, taxonomy.nameMap)
            logging.debug("Changing %s to %s" % (taxid,node.id))
            taxid = node.id
            changes+=1
        print '\t'.join([acc,str(taxid)])
    logging.info("Changed %d taxon ids" % changes)

if __name__ =='__main__':
    main()
