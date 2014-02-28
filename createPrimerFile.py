#!/usr/bin/python
import re
from edl.util import *
from optparse import OptionParser

def main():
    usage = "usage: %prog [PRIMER_TEMPLATE BARCODE [BARCODE]]"
    description="""
Generates a primer file for use iwth trimmomatic (and illuminPrep.py)
Given a template file and barcode(s), replaces placeholders with these barcodes

If no arguments are given, it returns (via STDOUT) a dummy file with all N's
    """

    parser = OptionParser(usage, description=description)
    addUniversalOptions(parser)
    (options, args) = parser.parse_args()
    setupLogging(options, description)

    # check/adjust options
    if len(args)==0:
        dummy=True
        primers=getDummyPrimers()
    else:
        dummy=False
        if len(args)==2 or len(args)==3:
            primers=getPrimers(*args)
        else:
            parser.error("Wrong number of arguments, should be 0, 2, or 3")

    print primers

#############
# Functions
def getDummyPrimers():
    return """>Prefix/1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>Prefix/2
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"""

def getPrimers(*args):
    """
    Create a primer file from a template with one or two barcodes
    replace "{BARCODE}" in file
    """
    if len(args)==2:
        (template,barcode)=args
        barcodes=[barcode,barcode]
    else:
        (template,barcode1,barcode2)=args
        barcodes=(barcode1,barcode2)

    logging.info("Processing %s with %s" % (template, barcodes))
    subst=0
    primers=""
    with open(template) as f:
        for line in f:
            if subst>0 and not line.startswith(">"):
                logging.debug(line)
                line=barCodeRcRE.sub(reverseComplement(barcodes[subst-1]),line)
                line=barCodeRE.sub(barcodes[subst-1],line)
                logging.debug(line)
            elif line.startswith(">Prefix"):
                try:
                    subst=int(line.split(None,1)[0][-1:])
                except ValueError:
                    logging.error("Prefix Id must end with /1 or /2. Offending line:\n %s" % (line))
                    sys.exit(-1)
            primers+=line
    return primers

def reverseComplement(sequence):
    newSeq=""
    for i in xrange(len(sequence)-1,-1,-1):
        c=sequence[i]
        if c=='A':
            c='T'
        elif c=='a':
            c='t'
        elif c=='t':
            c='a'
        elif c=='T':
            c='A'
        elif c=='G':
            c='C'
        elif c=='g':
            c='c'
        elif c=='C':
            c='G'
        elif c=='c':
            c='g'
        newSeq+=c
    return newSeq


#############
# Expressions
barCodeRcRE=re.compile(r'\{EDOCRAB\}')
barCodeRE=re.compile(r'\{BARCODE\}')

if __name__ == '__main__':
    main()
