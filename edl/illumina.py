import re, logging, sys, gzip
import subprocess
from Bio import SeqIO, Seq, SeqRecord

logger=logging.getLogger(__name__)

overlapRE=re.compile(r'^\S+\s+INFO\s+BESTOLP\s+(\S+?)(?:\:(?:\d{1,2})|(?:\:[NACTG]+))?\s+(\S+)')
errorRE=re.compile(r'^\S+\s+ERR\s+(\S+)\s+(\S+)(?:\:\d)?\s+')
dbgRE=re.compile(r'^\S+\s+DBG')
def scanPandaseqLog(stream,errFile=None):
    """
    Scan the output of pandaseq and find the unpaired reads. Can also look for errors.
    Each requested type is returned as a dictionary object with read names as keys.
    The error dictionary has the error code as the value.
    """

    logger.info("Scaning pandaseq logging for unpaired reads")

    if errFile is not None:
        if isinstance(errFile,str):
            logger.debug("Opening errFile: %s" % errFile)
            errStream=open(errFile,'w')
        else:
            logger.debug("Copying err output to %s" % errFile)
            errStream=errFile

    keepNone = logger.getEffectiveLevel() > logging.WARN
    keepDBG = logger.getEffectiveLevel() <= logging.DEBUG

    unpaired={}
    #failed={}
    lineCount=0
    olpCount=0
    #unpCount=0
    errCount={}
    for line in stream:
        lineCount+=1

        # logging
        if errFile is not None and not keepNone:
            # only print line if...
            if keepDBG:
                # debugging is on
                errStream.write(line)
            else:
                # or this is NOT a DBG line.
                if dbgRE.match(line) is None:
                    errStream.write(line)
                else:
                    # may as well move to next line since we know this is DBG
                    continue

        # is this an BESTOLP line?
        m=overlapRE.match(line)
        if m:
            olpCount+=1
            read=m.group(1)
            score=m.group(2)
            #logger.debug("Score=%r" % score)
            if score == '-1':
                #unpCount+=1
                unpaired[read]=True
            continue

        # is this an ERR line?
        m=errorRE.match(line)
        if m:
            read=m.group(2)
            code=m.group(1)
            errCount[code]=1+errCount.get(code,0)
            logger.debug("%s Failed with %s" % (read,code))
            #failed[read]=code
            continue

    logger.debug("%d overlap lines in %d total lines" % (olpCount, lineCount))
    unpCount=len(unpaired)
    msg="""
#===============================
# PandaSeq Complete:
#  Processed: %d
#  Paired:    %d
#  Unpaired:  %d
#  Errors (LOWQ => under quality threshold. These are normal):
%s#===============================\n""" % (olpCount,olpCount-unpCount-sum(errCount.itervalues()),unpCount,errCountString(errCount))
    logger.info(msg)
    if errFile is not None:
        errFile.write(msg)
        if isinstance(errFile,str):
            errStream.close()


    errCount['reads']=olpCount
    errCount['paired']=olpCount-unpCount
    errCount['unpaired']=unpCount

    logger.debug("%d errors, %d unpaired" % (sum(errCount.values()),len(unpaired)))
    return (unpaired,errCount)

def errCountString(errCount):
    retStr=""
    for (code,count) in errCount.iteritems():
        retStr="%s#   %s:   %d\n" %(retStr, code, count)
    return retStr

def runPandaseq(forward, reverse, output, unpairedList=None, catchUnpaired=True, exe='/common/bin/pandaseq', fastq=False, gap=20, errFile=None, inputFormat='fastq', threshold=0.6):
    """
    Runs pandaseq on the two given (forward,reverse) fastq files, printing the paired
    reads to the file or stream named in the third argument. By default, unparied reads are
    joined with 20 N's and added to the output file.

    if unpairedList is set to a filename, unpaired read names will be printed there
    instead of joining them and appending to output.

    if catchUnpaired is set to False, the output will not be scanned.

    If errFile is set, the pandaseq logs will be saved to this file, plus some summary counts. Pandaseq DBG lines will be dropped unless logging level is DEBUG.

    Note: input files are assumed to be fastq, even if output is fasta.
    """

    counts={}

    # start by building up a command
    #  debug: S gives summary stats, B gives us overlaps
    command=[exe, '-f', forward, '-r', reverse, '-d', 'rBfkmS', '-t', str(threshold)]
    if fastq:
        command.append('-F')
    if inputFormat=='fastq-solexa' or inputFormat=='fastq-illumina':
        command.append('-6')

    # set up streams and pipes
    if isinstance(output,str):
        pairsOut=open('output','w')
    else:
        pairsOut=output
    if catchUnpaired:
        logOut=subprocess.PIPE
    else:
        if errFile is not None:
            if isinstance(errfile,str):
                logOut=open(errFile)
            else:
                logOut=errFile
        else:
            logOut=None

    # start program
    logger.info("Running: %s" % (command))
    p=subprocess.Popen(command,stdout=pairsOut,stderr=logOut)

    if not catchUnpaired:
        # just let it run
        exitcode=p.wait()
        if exitcode != 0:
            logger.warn("pandaseq failed with code: %d" % exitcode)
            sys.exit(exitcode)
    else:
        # catch list of unpaired reads
        # orig approach didn't iterate over lines
        #(stdout,stderr)=p.communicate()
        #unpaired=scanPandaseqLog(stderr,singles=True,errFile=errFile)
        # this works,but apparently may deadlock
        (unpaired,counts)=scanPandaseqLog(p.stderr,errFile=errFile)

        logger.info("PAndaseq left %d read pairs unassembled" % (len(unpaired)))
        exitcode=p.wait()
        if exitcode != 0:
            logger.warn("pandaseq failed with code: %d" % exitcode)
            sys.exit(exitcode)

        if unpairedList is None:
            # append
            logger.info("appending unpaired reads with gap=%d" % gap)
            fakePairs(unpaired,forward,reverse,pairsOut,gap=gap,fastq=fastq,inputFormat=inputFormat,errFile=errFile)
        else:
            # write names to file
            if isinstance(unpairedList,str):
                singleFile=open(unpairedList,'w')
            else:
                singleFile=unpairedList
            for single in unpaired:
                singleFile.write("%s\n" % (single))
            if isinstance(unpairedList,str):
                singleFile.close()

    if isinstance(output,str):
        pairsOut.close()

def fakePairs(singles,forward,reverse,paired,gap=20,fastq=True,trim=None,inputFormat='fastq',errFile=None,batchSize=10000):
    """
    Given:
        a list of read names,
        forward and revers fastq files,
        file to append to
        gap size
    Pull out forward and reverse reads and join with N's as specified by gap
    WRite combined reads to 'paired' (can be filename or handle)
    """
    trimmingCounts={}

    # if singles is not a dict, make it so
    if not isinstance(singles,dict):
        tempDict={}
        for s in singles:
            tempDict[s]=True
        singles=tempDict

    # Open input files. Use gzip if needed
    if len(reverse)>3 and reverse[-3:]=='.gz':
        rrecords=SeqIO.parse(gzip.open(reverse), format=inputFormat)
    else:
        rrecords=SeqIO.parse(reverse, format=inputFormat)
    if len(forward)>3 and forward[-3:]=='.gz':
        frecords=SeqIO.parse(gzip.open(forward), format=inputFormat)
    else:
        frecords=SeqIO.parse(forward, format=inputFormat)

    joinCount=0
    fwdCount=0
    revCount=0
    nulCount=0
    totalJoined=0
    totalWritten=0
    if gap>=0:
        gapStr='N'*gap
        gapQ=[0]*gap
    fakedJoins=[]
    if fastq:
        outputFormat=inputFormat
    else:
        outputFormat='fasta'
    for frec in frecords:
        # get matching reverse record and make sure they are still synced
        try:
            rrec=rrecords.next()
        except StopIteration:
            logger.warn("Too few records in reverse fastq: %s" % reverse)
            sys.exit(1)
            if rrec.id[:-1] != frec.id[:-1]:
                logger.warn("Mismatched IDs in forward and reverse fastq (%s :: %s)\n%s\n%s" % (frec.id,rrec.id,forward,reverse))

        # ignore all but singles
        if not singles.pop(rrec.id,False):
            continue

        # optional trim
        if trim is not None:
            (frec,fwhy)=trimEnds(frec,trim)
            trimmingCounts[fwhy]=trimmingCounts.get(fwhy,0)+1
            (rrec,rwhy)=trimEnds(rrec,trim)
            trimmingCounts[rwhy]=trimmingCounts.get(rwhy,0)+1

        # join seqs
        newRecs=[]
        if frec is None:
            if rrec is None:
                logger.debug("Both ends trimmed to oblivion: (%s,%s)" % (fwhy,rwhy))
                nulCount+=1
                continue
            logger.debug("Forward seq trimmed to oblivion (%s)" % (fwhy))
            newRecs.append(revComp(rrec,qual=fastq,suffix=".rev"))
            revCount+=1
        elif rrec is None:
            logger.debug("Reverse seq trimmed to oblivion (%s)" % (rwhy))
            newRecs.append(frec)
            fwdCount+=1
        elif gap>=0:
            newSeq=frec.seq + Seq.Seq(gapStr,frec.seq.alphabet) + rrec.seq.reverse_complement()
            newRec=SeqRecord.SeqRecord(newSeq,id=frec.id,name=frec.name,description="Faked join")
            joinCount+=1
            # join quality
            if fastq:
                newRec.letter_annotations['phred_quality'] = \
                    frec.letter_annotations['phred_quality'] + \
                    gapQ + \
                    list(reversed(rrec.letter_annotations['phred_quality']))
            newRecs.append(newRec)
        else:
            # gap < 0 means don't join...add separately
            newRecs.append(frec)
            newRecs.append(revComp(rrec,qual=fastq,suffix=".rev"))

        fakedJoins.extend(newRecs)
        if len(fakedJoins)>=batchSize:
            numWritten =  SeqIO.write(fakedJoins, paired, format=outputFormat)
            if numWritten != len(fakedJoins):
                logger.warn("Only %d of %d faked joins written!" % (numWritten, len(fakedJoins)))
            totalJoined+=len(fakedJoins)
            totalWritten+=numWritten
            del fakedJoins[:]

    # Confirm that the reverse record iterator is also finished
    try:
        rrec=rrecords.next()
        # should not get here, iterator should be done
        logger.warn("Extra records in reverse fastq (%s):\n%s" % (rrec.id,reverse))
    except StopIteration:
        # this is what we expect
        pass

    numWritten =  SeqIO.write(fakedJoins, paired, format=outputFormat)
    if numWritten != len(fakedJoins):
        logger.warn("Only %d of %d faked joins written!" % (numWritten, len(fakedJoins)))
    totalWritten+=numWritten
    totalJoined+=len(fakedJoins)

    # Report some counts
    msg="#======================\n# Faked joins\n#  Total: %s\n" % (totalJoined)
    if trim is not None:
        msg+="# Joined: %d\n# FwdOnly: %d\n# RevOnly: %d\n# Dropped: %d\n" % (joinCount, fwdCount, revCount, nulCount)
        msg+="#======================\n# End trimming:\n"
        for status,count in trimmingCounts.iteritems():
            msg+="# %s: %d\n" % (status,count)
    msg+="#======================\n"

    if errFile is not None:
        if isinstance(errFile,str):
            with open(errFile,'a') as errstream:
                errstream.write(msg)
        else:
            errFile.write(msg)

    if trim is not None:
        logger.info(msg)
        logger.debug("# Total written: %d" % (totalWritten))
    else:
        logger.info("# Faked joins: %d" % (totalJoined))

    if len(singles)>0:
        raise Exception("Some (%d) unpaired read names were not found! EG: %s" % (len(singles),singles.keys()[0]))

    return numWritten

def revComp(record, qual=False, suffix=''):
    newRec=SeqRecord.SeqRecord(record.seq.reverse_complement(),id=record.id+suffix,name=record.name,description=record.description)
    if qual:
        newRec.letter_annotations['phred_quality']=list(reversed(record.letter_annotations['phred_quality']))
    return newRec

def trimEnds(record,trimOptions):
    """
    Given a sequence record, remove low quality bases from ends
    """
    minlen=trimOptions.minLength
    cutoff=trimOptions.endQuality
    polyfrac=trimOptions.polyFrac
    rlen=len(record)

    if cutoff>0:
        scores=record.letter_annotations['phred_quality']
        max=rlen-minlen+1
        for i in xrange(max):
            if scores[i]>=cutoff:
                start=i
                break
        else:
            # no scores over cutoff.
            return (None, "Dropped for qulaity/length")

        for i in xrange(max-start):
            if scores[-1-i]>=cutoff:
                if i==0:
                    end=None
                else:
                    end=-i
                break
        else:
            # no scores over cutoff before we hit the min length
            return (None, "Dropped for qulaity/length")
    else:
        if rlen<minlen:
            return (None, "Dropped for length")
        start=0
        end=None

    if polyfrac<1.0:
        #logger.debug("Checking for poly T or A")
        if checkForPoly1(record.seq,rlen,polyfrac):
            return (None, "Dropped for polyA")

    return (record[start:end],"Trimmed")

def checkForPoly1(seq,rlen,frac):
    otherLim=(1.0-frac)*rlen
    #logger.debug("max=%s" % otherLim)
    tcount=0
    acount=0
    for c in seq:
        if c!='T' and c!='t':
            tcount+=1
        elif c!='A' and c!='a':
            acount+=1
        if acount>otherLim and tcount>otherLim:
            #logger.debug("not t: %d   not a: %d" % (tcount,acount))
            return False
    #logger.debug("not t: %d   not a: %d" % (tcount,acount))
    return True

tRE=re.compile(r'[tT]')
aRE=re.compile(r'[aA]')
def checkForPoly(seq,rlen,frac):
    otherLim=(1.0-frac)*rlen
    tcount=0
    acount=0
    for c in seq:
        if tRE.match(c) is None:
            tcount+=1
        elif aRE.match(c) is None:
            acount+=1
        if acount>otherLim and tcount>otherLim:
            return False
    return True

def checkForPoly2(seq,rlen,frac):
    counts={}
    for c in seq:
        counts[c]=counts.get(c,0)+1
    lim=frac*rlen
    for c in ('a','t'):
        if counts.get(c,0)+counts.get(c.upper(),0)>lim:
            return False
    return True

def pyTrim(infile, outfile, trimOptions, errFile=None, inputFormat='fastq', fastq=True,reads=None):
    """
    Apply simple trimming using biopython
    trimOptions is an object with the vaules: endQuality, polyFrac, minLength
    if reads is set to a dict of read names, only process named reads
    if fastq is set to False, ouput fasta
    """
    records=SeqIO.parse(infile,format=inputFormat)
    if reads is not None:
        records = screenRecords(records,reads)

    if isinstance(errFile,str):
        errstream=open(errFile,'w')
    else:
        errstream=errFile
    trimmedRecords=recordTrimmer(records,trimOptions,errstream)

    if fastq:
        outputFormat=inputFormat
    else:
        outputFormat='fasta'
    count=SeqIO.write(trimmedRecords, outfile, format=outputFormat)
    logger.debug("Wrote %d records to %s" % (count,outfile))

def screenRecords(records, reads):
    for record in records:
        if record.id in reads:
            yield record

def recordTrimmer(records, trimOptions, errstream):
    recordCount=0
    trimmedCounts={}
    for record in records:
        recordCount+=1
        (trimmed,status) = trimEnds(record, trimOptions)
        trimmedCounts[status]=trimmedCounts.get(status,0)+1
        if trimmed is not None:
            yield trimmed
        else:
            if errstream is not None and logger.getEffectiveLevel()<=logging.WARN:
                errstream.write("Sequence %s trimmed to 0 length\n" % (record.id))
            logger.debug("Sequence %s trimmed to 0 length" % (record.id))
    logger.debug("Trimmer processed %d records: %s" % (recordCount, trimmedCounts))

    msg="#==================\n# Assembled pair Trimming:\n"
    for status in trimmedCounts:
        msg+="# %s: %d\n" % (status,trimmedCounts[status])
    msg+="#==================\n"

    logger.info(msg)
    if errstream is not None:
        errstream.write(msg)

