#!/usr/bin/env python
import re, argparse
from edl.util import *

def main():
    description="""
Generates a primer file for use iwth trimmomatic (and illuminPrep.py)
Given a template file and barcode(s), replaces placeholders with these barcodes

    """

    parser = argparse.ArgumentParser(description=description)
    add_universal_arguments(parser)
    parser.add_argument('template',type=argparse.FileType('r'),
                    default=sys.stdin,
                    help="primer template file with {BARCODE} placeholders")
    parser.add_argument('barcodes',nargs='+', 
                    help="1 or 2 barcode sequences")
    arguments = parser.parse_args()
    setup_logging(arguments)

    # check/adjust options
    if len(arguments.barcodes)>2:
        parser.error("Wrong number of arguments, there should be either 1 or 2 barcode sequences")

    primers=get_primers(**vars(arguments))

    print (primers)

#############
# Functions
def getDummyPrimers():
    return """>Prefix/1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>Prefix/2
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"""

def get_primers(template="dummy", barcodes=["NNNNNN"], **kwargs):
    """
    Create a primer file from a template with one or two barcodes
    replace "{BARCODE}" in file
    """
    if template=='dummy':
        return getDummyPrimers()

    if len(barcodes)==1:
        barcodes.append(barcodes[0])

    logging.info("Processing %s with %s" % (template, barcodes))
    subst=0
    primers=""
    for line in template:
        if subst>0 and not line.startswith(">"):
            logging.debug(line)
            line=barCodeRcRE.sub(reverse_complement(barcodes[subst-1]),line)
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

def reverse_complement(sequence):
    newSeq=""
    for i in range(len(sequence)-1,-1,-1):
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
