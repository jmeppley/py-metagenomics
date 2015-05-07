#!/usr/bin/env python
#$ -S /usr/bin/python
#$ -V
#$ -cwd
#$ -l greedy=8
#$ -l rps=1

import os, sys, tempfile, time
from edl.illumina import runPandaseq, fakePairs, pyTrim
from edl.tmatic import TMOptionsPE, TMOptionsSE
from edl.util import *
from Bio import SeqIO
from optparse import OptionParser

def main():
    usage = "usage: %prog [OPTIONS] -f FORWARD_FASTQ -r REVERSE_FASTQ -o OUTPUT"
    description="""
Runs preprocessing, quailty filtering, and read-pairing on illumina data.

Given:
    a pair of raw fastq or fastq.gz files from MiSeq
    a primers fasta file (optional)
    a minimum quality score for ends (optional)

Runs:
    TrimmomaticPE to remove primers
    Pandaseq to combine
    Simple qual based trimming

Produce:
    a single fastq (or fasta) file
    """

    parser = OptionParser(usage, description=description)
    parser.add_option("-f", "--forward",
                      metavar="FASTQ", help="Forward reads in fastq file")
    parser.add_option("-r", "--reverse",
                      metavar="FASTQ", help="reverse reads in fastq file")
    parser.add_option("-o", "--outfile", dest="outfile",
                      metavar="OUTFILE", help="Write joined pairs to OUTFILE.")
    parser.add_option("-e", "--errFile", default=None,
                       help="append logs to this file")
    parser.add_option("-6", "--phred64", default=False, action='store_true',
                      help="Input files will have their quality encoded as PHRED + 64 instead of the PHRED + 33. PHRED + 64 was used originally in the CASAVA pipeline from version 1.3 through 1.7. In CASAVA 1.8, the score is encoded as PHRED + 33, the default.")
    parser.add_option("-F","--fastq", default=True, action='store_true',
                      help="output in fastq format (default)")
    parser.add_option("-a","--fasta", dest='fastq', action='store_false',
                      help="output in fasta format")
    parser.add_option('-g', "--gapJoin", type='int', default=20,
                      help="Number of N's to insert between unjoined pairs. Set to -1 to leave unjoined.")
    parser.add_option('-p', "--primersFile", dest='primers', metavar='PRIMERS_FILE',
                      help="File containing primer sequences that Trimomatic should remove")
    parser.add_option('-C', "--trimPrimerCustom", default=None, metavar='TMATIC_SETTINGS',
                      help="Look for missed primer sequence with Trimmomatic and custom settings. Value should be a string like 2:40:15 where the numbers are <seed mismatches>:<palindrome clip threshold>:<simple clip threshold> as defined in http://www.usadellab.org/cms/index.php?page=trimmomatic. ")
    parser.add_option('-Q', "--endQuality", default=5, type='int', metavar='QUAL',
                      help='Minimum quality for bases at either end of final sequences, anything under this will be dropped. Default is 5.')
    parser.add_option('-P','--polyFrac', default=.9, type='float', metavar='FRAC',
                      help="maximum fraction of A's or T's allowed")
    parser.add_option('-m', "--minLength", default=45, type='int', metavar='LEN',
                      help='Minimum length of final reads. Defaults to 45. ')
    parser.add_option('-T','--noTrim',default=False, action='store_true',
                      help='skip all post-join trimming')
    parser.add_option('-t','--pandaseqThreshold',default=.32,type='float',metavar='FRAC',
                      help='Threshold for pandaseq quailty filter (between 0 and 1). Defaults to .32')
    parser.add_option("-L","--loglevel", type='int', dest="verbose",
                     help="Shortcut to set verbosity directly")
    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options, description)

    # input files
    forward=options.forward
    reverse=options.reverse

    pandaSeqLoc='pandaseq'

    if options.phred64:
        inputFormat='fastq-illumina'
    else:
        inputFormat='fastq'

    # working dir
    tempdir=tempfile.mkdtemp(dir='.')

    if options.errFile is None:
        options.errFile="%s.log" % (options.outfile)
    logstream=open(options.errFile,'a')

    ##
    # remove missed primers with tmatic
    if options.primers is not None:
        tmOptions=TMOptionsPE(forward,reverse,primers=options.primers,pref=tempdir,primerSettings=options.trimPrimerCustom)
        tmOptions.runTmatic()
        logging.debug("Trimmomatic messages:\n%s" % (tmOptions.stdout))
        logstream.write("TMATIC_PE:STDOUT:\n%s\n" % (tmOptions.stdout))
        logstream.write("TMATIC_PE:STDERR:\n%s\n" % (tmOptions.stderr))

        if tmOptions.exitcode!=0:
            logging.warn("TMATIC:STDERR:\n%s" % (tmOptions.stderr))
            sys.exit(tmOptions.exitcode)

        forward=tmOptions.outfiles['1p']
        reverse=tmOptions.outfiles['2p']
        singles=tmOptions.outfiles['1u']

    ##
    # Matepair joining
    # setup intermediate files if doing final triming
    if not options.noTrim:
        fastq=True
        pandaTmp="%s/pandaJoined" % tempdir
        pandaOut=open(pandaTmp,'w')
        #(pandaOut,pandaTmp)=tempfile.mkstemp(suffix='fastq')
        # save list of unpaired reads to trim before joining
        unpairedTmpList="%s/pandaUnpaired" % tempdir
        unpairedDest=open(unpairedTmpList,'w')
        #(unpairedDest,unpairedTmpList)=tempfile.mkstemp()
    else:
        fastq=options.fastq
        pandaOut=open(options.outfile,'w')
        # append unpaired reads to joined file
        unpairedDest=None

    # join ends
    runPandaseq(forward,reverse,pandaOut,unpairedList=unpairedDest,catchUnpaired=True,gap=options.gapJoin,errFile=logstream,fastq=fastq,inputFormat=inputFormat,threshold=options.pandaseqThreshold,exe=pandaSeqLoc)

    # add singes from tmatic
    if options.primers is not None:
        addSingles(singles,pandaOut,fastq)

    # close output stream
    pandaOut.close()
    if unpairedDest is not None:
        unpairedDest.close()

    ##
    # final trimming if requested
    if not options.noTrim:
        # run trimming
        logstream.write("Starting assembled pair Trimming:\n")
        pyTrim(pandaTmp,options.outfile,options,errFile=logstream,inputFormat=inputFormat,fastq=options.fastq)

        # close and delete temp file

        if logging.getLogger().level > logging.DEBUG:
            os.remove(pandaTmp)

        # trim unpaired reads and append to outfile as faked joins
        trimUnpairedReads(unpairedTmpList,forward,reverse, options, options.outfile, inputFormat, options.gapJoin, fastq=options.fastq, errFile=logstream)
        if logging.getLogger().level > logging.DEBUG:
            os.remove(unpairedTmpList)

    if logging.getLogger().level <= logging.DEBUG:
        logging.debug("temp files saved to: %s" % tempdir)
    else:
        for f in os.listdir(tempdir):
            os.remove(os.sep.join([tempdir,f]))
        try:
            time.sleep(2)
            os.rmdir(tempdir)
        except:
            logging.warn("Could not remove temporary directory: %s " % (tempdir))

def trimUnpairedReads(urList, fwd, rev, trimLimits, out, inputFormat, gap, fastq=True,errFile=None):
    """
    for read in list:
        pull out of fwd and rev fastq files
        run through tmatic (or use biopython?)
         trim is a dict with two keys: 'qual' and 'len' which define the lowest acceptable starting/ending quality and the minimum read length.
        fake join
        append to out file
    """
    reads=parseListToMap(urList)
    #if fake:
    logging.info("Submitting %d reads to be joined with %d N's" % (len(reads),gap))
    outstream = open(out,'a')
    fakePairs(reads,fwd,rev,outstream,trim=trimLimits,inputFormat=inputFormat,fastq=fastq,errFile=errFile,gap=gap)
    outstream.close()
    #else:
    #    out_fwd=out+".unpaired.fwd"
    #    out_rev=out+'.unpaired.rev'
    #    logging.debug("Trimming %d unpaired forward reads to %s" % (len(reads),out_fwd))
    #    pyTrim(fwd,out_fwd,trimLimits,errFile=errFile,inputFormat=inputFormat,fastq=fastq,reads=reads)
    #    logging.debug("Trimming %d unpaired reverse reads to %s" % (len(reads),out_rev))
    #    pyTrim(rev,out_rev,trimLimits,errFile=errFile,inputFormat=inputFormat,fastq=fastq,reads=reads)

def addSingles(singlesFile,outStream,fastq):
    """
    append data from given singles file to output stream
    if fastq is false, convert to fasta
    """

    if fastq:
        with open(singlesFile) as f:
            count=0
            for line in f:
                outStream.write(line)
                count+=1
            logging.info("Appended %d lines of fastq to %s" % (count,outStream))
    else:
        records=SeqIO.parse(singlesFile,format='fastq')
        count=SeqIO.write(records,outStream,'fasta')
        logging.info("Appended %d fasta records to %s" % (count,outStream))

if __name__ == '__main__':
    main()
