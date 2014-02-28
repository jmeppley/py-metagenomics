#! /usr/bin/python
import sys, logging, re
from optparse import OptionGroup
from edl.util import parseExp, LineCounter, openInputFile
logger=logging.getLogger(__name__)

#############
# Constants #
#############
GENE='gene'
LIZ='liz'
YANMEI='yanmei'
BLASTPLUS='blast'
FRHIT='frhit'
LAST0='last'
HMMSCANDOM='hmmscandom'
HMMSEARCHDOM='hmmsearchdom'
HMMSCAN='hmmscan'
HMMSEARCH='hmmsearch'
SAM='sam'
formatsWithNoDescription=[LAST0, FRHIT, BLASTPLUS, SAM]
cigarRE=re.compile(r'\d+[^\d]')

#############
# Classes   #
#############
class M8Stream(LineCounter):
    """ Inherits line counting from edl.util.LineCounter plus ability to open anything from edl.util.openIntpuFile. The interface is a simple next() method and a close() method that only does anything if the first argument is a file name (and not a file-like object)
    """
    def __init__(self, fileName, *args):
        LineCounter.__init__(self, openInputFile(fileName, *args))
        self.fileName=fileName

    def close(self):
        """
        Close any file objects that were opened during __init__
        """
        if isinstance(self.fileName,str):
            self.rawStream.close()

class FilterParams:
    @staticmethod
    def createFromOptions(options, ignore=[], translate={}):
        """
        Translate an options object created using the addHitTableOptions below into a FilterParams object. The attributes format, sort, and sortReads are expected to be in the options object with the 'hitTable' prefix. The attributes bits, evalue, pctid, aln, length, hitsPerRead, hspsPerHit, and nonoverlapping with the prefix 'filter'. For example, options.hitTableFormat will be copied to params.format and options.filterTopPct to params.topPct.

        translate= should be set to a dictionary mapping attributes of the options object to the desired attribuets of the FilterParams object.

        ignore= should be a list of the standard params to be skipped (eg: filterTopPct).
        """
        params=FilterParams()

        # check the nonstandard one(s)
        for (oname,pname) in translate.iteritems():
            if hasattr(options,oname):
                setattr(params,pname,getattr(options,oname))

        # get the standard ones
        for param in ['bits','evalue','pctid','length', 'aln', 'hitsPerRead','hspsPerHit','nonoverlapping','topPct']:
            oparam = 'filter' + param[0].upper() + param[1:]
            if oparam not in ignore:
                if hasattr(options,oparam):
                    setattr(params,param,getattr(options,oparam))
        for param in ['format','sort','sortReads']:
            oparam = 'hitTable' + param[0].upper() + param[1:]
            if oparam not in ignore:
                if hasattr(options,oparam):
                    setattr(params,param,getattr(options,oparam))

        logging.debug("%r" % (options))
        logging.debug("%r" % (params))
        return params

    def __init__(self, format=GENE, topPct=-1, bits=0.0, evalue=None, pctid=0.0, length=0, aln=None,hitsPerRead=0,hspsPerHit=0,nonoverlapping=False, sort=None, sortReads=None):
        self.format=format
        self.topPct=topPct
        self.bits=bits
        self.evalue=evalue
        self.pctid=pctid
        self.length=length
        self.aln=aln
        self.hitsPerRead=hitsPerRead
        self.hspsPerHit=hspsPerHit
        self.nonoverlapping=nonoverlapping
        self.sortReads=sortReads
        self.sort=sort

    def __repr__(self):
        return "FilterParams(format=%r, topPct=%r, bits=%r, evalue=%r, pctid=%r, length=%r, aln=%r, hitsPerRead=%r, hspsPerHit=%r, nonoverlapping=%r sort=%r, sortReads=%r)" % (self.format, self.topPct, self.bits, self.evalue, self.pctid, self.length, self.aln, self.hitsPerRead, self.hspsPerHit, self.nonoverlapping, self.sort, self.sortReads)

class EmptyHitException(Exception):
    pass

class Hit:
    @staticmethod
    def getHit(line, options):
        try:
            return Hit(line, options)
        except EmptyHitException:
            return None
        except:
            logger.warn("Error parsing line:\n%s" % (line))
            raise

    """
    Object representing a single hit from a search program like blast or fr-hit
    """
    def __init__(self, line, options):
        """
        Creates Hit object from a line in a hit table. Options can be any class with a 'format' property or a string. The string value of options or options.format must be one of the recognized formats: 'gene','liz','yanmei','last','frhit'
        """
        self.line = line
        self.setFormat(options)
        self.parseLine(line)

    def __repr__(self):
        return "Hit(%r,%r)" % (self.line.rstrip('\n\r'), self.format)

    def setFormat(self, options):
        if type(options) is type(GENE):
            # just to make this work if instantiated from repr
            self.format=options
        else:
            self.format=options.format

        if self.format == GENE:
            self.parseLine=self.parseGeneLine
        elif self.format == LIZ:
            self.parseLine = self.parseLizLine
        elif self.format == YANMEI:
            self.parseLine = self.parseYanmeiLine
        elif self.format == BLASTPLUS:
            self.parseLine = self.parseBlastPlusLine
        elif self.format == LAST0:
            self.parseLine = self.parseLastalLine
        elif self.format == FRHIT:
            self.parseLine = self.parseFrHitLine
        elif self.format == HMMSEARCHDOM:
            self.parseLine = self.parseHmmSearchDomLine
        elif self.format == HMMSCANDOM:
            self.parseLine = self.parseHmmScanDomLine
        elif self.format == HMMSEARCH:
            self.parseLine = self.parseHmmSearchLine
        elif self.format == HMMSCAN:
            self.parseLine = self.parseHmmScanLine
        elif self.format == SAM:
            self.parseLine = self.parseSamLine
        else:
            sys.exit("Unknown format: %s" % (self.format))

    def parseGeneLine(self, line):
        cells=line.rstrip('\n\r').split('\t')
        self.read=cells[0]
        self.readDesc=cells[1]
        self.hit=cells[2]
        self.hitDesc=cells[3]
        self.pctid=float(cells[4])
        self.mlen=int(cells[5])
        self.qstart=int(cells[6])
        self.qend=int(cells[7])
        self.hstart=int(cells[8])
        self.hend=int(cells[9])
        self.score=float(cells[10])
        self.evalue=parseExp(cells[11])
        self.aln=float(cells[12])

    def parseLizLine(self, line):
        cells=line.rstrip('\n\r').split('\t')
        self.read=cells[0]
        self.hit=cells[1]
        self.hitDesc=cells[2]
        try:
            self.pctid=float(cells[3])
        except:
            self.pctid=0
        self.mlen=int(cells[4])
        self.qstart=int(cells[5])
        self.qend=int(cells[6])
        self.hstart=int(cells[7])
        self.hend=int(cells[8])
        self.score=float(cells[9])
        self.evalue=parseExp(cells[10])
        self.aln=float(cells[11])

    def parseSamLine(self, line):
        if line[0]=='@':
            raise EmptyHitException("reference sequence line")
        cells=line.rstrip('\n\r').split('\t')
        self.read=cells[0]
        self.hit=cells[2]
        if self.hit=='*':
            raise EmptyHitException("No match")
        (alen,alenh,alenq,qstart,qend,pctid)=parseCigarString(cells[5])
        self.mlen=alen
        self.qstart=qstart
        self.qend=qend
        self.astart=int(cells[3])
        self.aend=self.astart+alenh-1
        if pctid is not None:
            self.pctid=pctid
        else:
            self.pctid=0
        for tagstr in cells[11:]:
            if len(tagstr)>2 and tagstr[:2]=='AS':
                self.score = float(tagstr.split(":")[2])
                # for now, we only carea bout the score tag
                break
        else:
            raise Exception("No score (AS tag) found in line:\n%s" % (line))

    #target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
    def parseHmmScanDomLine(self, line):
        if line[0]=='#':
            raise EmptyHitException("Comment line")
        cells=line.rstrip('\n\r').split()
        self.read=cells[3]
        qlen=int(cells[5])
        self.hit=cells[0]
        hitLen=int(cells[2])
        self.evalue=float(cells[6])
        self.score=float(cells[7])
        self.hstart=int(cells[15])
        self.hend=int(cells[16])
        self.qstart=int(cells[17])
        self.qend=int(cells[18])
        self.hitDesc=cells[22]
        self.pctid=0
        self.mlen=1+self.qend-self.qstart
        self.aln=self.mlen/float(qlen)

    #target name        accession  query name           accession    E-value  score  bias ue  score  bias   exp reg clu  ov env dom rep inc description of target
    def parseHmmScanLine(self, line):
        if line[0]=='#':
            raise EmptyHitException("Comment line")
        cells=line.rstrip('\n\r').split()
        self.read=cells[2]
        self.hit=cells[0]
        self.evalue=float(cells[4])
        self.score=float(cells[5])
        self.hitDesc=cells[18]
        self.pctid=0
        self.mlen=0

    #target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
    def parseHmmSearchDomLine(self, line):
        if line[0]=='#':
            raise EmptyHitException("Comment line")
        cells=line.rstrip('\n\r').split()
        self.read=cells[0]
        qlen=int(cells[2])
        self.hit=cells[3]
        hitLen=int(cells[5])
        self.evalue=float(cells[6])
        self.score=float(cells[7])
        self.hstart=int(cells[15])
        self.hend=int(cells[16])
        self.qstart=int(cells[17])
        self.qend=int(cells[18])
        self.readDesc=cells[22]
        self.pctid=0
        self.mlen=1+self.qend-self.qstart
        self.aln=self.mlen/float(qlen)

    #target name        accession  query name           accession    E-value  score  bias ue  score  bias   exp reg clu  ov env dom rep inc description of target
    def parseHmmSearchLine(self, line):
        if line[0]=='#':
            raise EmptyHitException("Comment line")
        cells=line.rstrip('\n\r').split()
        self.read=cells[0]
        self.hit=cells[2]
        self.evalue=float(cells[4])
        self.score=float(cells[5])
        self.readDesc=cells[18]
        self.pctid=0
        self.mlen=0

    # score name1   start1  alnSize1        strand1 seqSize1        name2   start2  alnSize2   strand2 seqSize2        blocks
    def parseLastalLine(self, line):
        logger.debug(line)
        if line[0]=='#':
            raise EmptyHitException("Comment line")
        cells=line.rstrip('\n\r').split('\t')
        self.score=float(cells[0])
        self.hit=cells[1]

        hmlen=int(cells[3])
        hlen=int(cells[5])
        if cells[4] == '+':
            self.hstart=int(cells[2]) + 1
            self.hend=self.hstart + hmlen - 1
        else:
            # start was in reverse strand
            self.hstart = hlen - int(cells[2])
            self.hend = self.hstart - hmlen + 1

        self.read=cells[6]
        qmlen=int(cells[8])
        qlen=int(cells[10])
        if cells[9] == '+':
            self.qstart=int(cells[7]) + 1
            self.qend=self.qstart + qmlen - 1
        else:
            # start was in reverse strand
            self.qstart = qlen - int(cells[7])
            self.qend = self.qstart - qmlen + 1
        self.aln=qmlen/float(qlen)

        self.evalue=None
        self.pctid=None
        self.mlen=computeLastHitValues(cells[11])
        self.hitDesc=None
        logger.debug("Span: %d-%d" % (self.qstart,self.qend))

    def parseFrHitLine(self, line):
        cells=line.rstrip('\n\r').split('\t')
        self.read=cells[0]
        self.evalue=parseExp(cells[2])
        self.score=None
        self.mlen=cells[3]
        self.qstart=int(cells[4])
        self.qend=int(cells[5])
        self.pctid=float(cells[7])
        self.hit=cells[8]
        self.hstart=int(cells[9])
        self.hend=int(cells[10])

    def parseBlastPlusLine(self, line):
        cells=line.rstrip('\n\r').split('\t')
        self.read=cells[0]
        self.hit=cells[1]
        self.pctid=float(cells[2])
        self.mlen=int(cells[3])
        self.mismatch=int(cells[4])
        self.gaps=int(cells[5])
        self.qstart=int(cells[6])
        self.qend=int(cells[7])
        self.hstart=int(cells[8])
        self.hend=int(cells[9])
        self.evalue=parseExp(cells[10])
        self.score=float(cells[11])

    def parseYanmeiLine(self, line):
        cells=line.rstrip('\n\r').split('\t')
        self.read=cells[0]
        self.hit=cells[1]
        self.pctid=float(cells[2])
        self.mlen=int(cells[3])
        self.mismatch=int(cells[4])
        self.gaps=int(cells[5])
        self.qstart=int(cells[6])
        self.qend=int(cells[7])
        self.hstart=int(cells[8])
        self.hend=int(cells[9])
        self.evalue=parseExp(cells[10])
        self.score=float(cells[11])
        self.hitDesc=cells[12]

    def getAln(self):
        if 'aln' in dir(self):
            return self.aln
        else:
            sys.exit("Cannot calculate alignment percentage from data in this m8 format")

    def checkForOverlap(self,regions):
        """
        check to see if this hit overlaps any already hit region. The regions array should be a list of (Start,end) pairs indicating resiongs already hit.
        """
        logger.debug("Looking for [%d,%d] in %s" % (self.qstart,self.qend,regions))

        # make sure start is smaller than end
        start = min(self.qstart,self.qend)
        end = max(self.qstart,self.qend)

        # compare with all used ranges
        for i in range(len(regions)):
            occupiedRange = regions[i]
            # hit cannot intersect an used range
            if (start>=occupiedRange[1] or end<=occupiedRange[0]):
                # does not overlap this range (try next range)
                continue
            else:
                # overlaps. We are done here
                return None

        # we get here if there was no overlap
        regions.append((start,end))
        return regions

    def getLine(self,options):
        return self.line

#############
# Functions #
#############

def computeLastHitValues(blocks):
    """
    given the query length and 'blacks' sting from a last hit
    return the:
        match length

    the blocks string looks something like this:
        "73,0:1,15,0:1,13,0:1,9"
        where integer elements indicate lenghts of matches and
        colon separated elements indicate lengths of gaps
    """
    matchLen=0
    for segment in blocks.split(','):
        try:
            matchLen+=int(segment)
        except ValueError:
            (hml,qml) = segment.split(":")
            mml=max(int(hml),int(qml))
            matchLen+=mml

    return matchLen

def getReadCol(format):
    logger.info("Getting read col for format: %s" % format)
    if format == LAST0:
        return 6
    else:
        return 0
def getHitCol(format, useDesc=False):
    logger.info("Getting hit col for format: %s" % format)
    if format == GENE:
        if useDesc:
            hitCol=3
        else:
            hitCol=2
    elif format == LAST0:
        if useDesc:
            raise Exception("lastal does not report the hit description, sorry")
        else:
            hitCol=1
    elif format == LIZ:
        if useDesc:
            hitCol=2
        else:
            hitCol=1
    elif format == YANMEI:
        if useDesc:
            hitCol=12
        else:
            hitCol=1
    elif format == BLASTPLUS:
        if useDesc:
            raise Exception("The default blast table does not keep the hit description")
        else:
            hitCol=1
    elif format == SAM:
        if useDesc:
            raise Exception("The SAM format does not keep the hit description")
        else:
            hitCol=2
    else:
        sys.exit("I'm sorry, I don't understand the format: %s" % (format))

    return hitCol

def filterM8(instream, outstream, options):
    """
    Filter instream and write output to outstream
    """
    #if options.sortReads:
    #    instream=sortLines(instream)
    for line in filterM8Stream(instream, options):
        outstream.write(line)

def sortLines(instream):
    logger.info("Sorting input lines")
    lines=[]
    for line in instream:
        lines.append(line)
    lines.sort()
    logger.debug("Done sorting")
    return lines.__iter__()

def getHitStream(instream, options):
    if options.sortReads:
        return getSortedHits(instream, options)
    else:
        return getUnsortedHitStream(instream, options)

def getSortedHits(instream, options):
    """
    Read in every line to a list of hits, sort by read and return individually
    """
    hits=[]
    for line in instream:
        if len(line)==0:
            continue
        hit = Hit.getHit(line,options)
        if hit is not None:
            hits.append(hit)

    sortKey = lambda h: (h.read)

    hits.sort(key=sortKey)
    return hits

def getUnsortedHitStream(instream, options):
    """
    Simply parse lines into Hits one at a time and return them
    """
    for line in instream:
        if len(line)==0:
            continue
        hit = Hit.getHit(line,options)
        if hit is not None:
            yield hit

def filterM8Stream(instream, options, returnLines=True):
    """
    return an iterator over the lines in given input stream that pass filter
    """

    currentRead=None
    hits=[]
    needsFilter=doWeNeedToFilter(options)
    for hit in getHitStream(instream, options):
        if hit.read != currentRead:
            if currentRead is not None:
                logging.debug("processing %d hits for %s" % (len(hits), currentRead))
                if options.sort is not None:
                    sortHits(hits,options.sort)
                if needsFilter:
                    hits=filterHits(hits, options, returnLines=returnLines)
                if returnLines:
                    for line in hits:
                        yield line
                else:
                    yield (currentRead, hits)
                hits=[]
            currentRead=hit.read
        #logging.debug("%s -> %s" % (hit.read, hit.hit))
        hits.append(hit)

    if currentRead is not None:
        if options.sort is not None:
            sortHits(hits,options.sort)
        if needsFilter:
            hits=filterHits(hits, options,returnLines=returnLines)
        if returnLines:
            for line in hits:
                yield line
        else:
            yield (currentRead, hits)

def sortHits(hits, sortType):
    """
    take all the hits for a read as a list and sort them by score and hit name
    """

    logger.debug("Sorting hits (byScore = %s)" % (sortType))

    # set up sort key
    if sortType=='evalue':
        sortKey = lambda h: (h.evalue, h.hit)
    else:
        sortKey = lambda h: (h.score*-1, h.hit)

    # sort in place
    hits.sort(key=sortKey)

def doWeNeedToFilter(options):
    if options.topPct>=0:
        return True
    if options.bits > 0:
        return True
    if options.evalue is not None:
        return True
    if options.pctid>0:
        return True
    if options.length>0:
        return True
    if options.hitsPerRead>0:
        return True
    if options.hspsPerHit>0:
        return True
    return False

def filterHits(hits, options, returnLines=True):
    # A topPct cutoff, requires finding the top score first
    if options.topPct >= 0:
        # get the best score in this set of hits
        bestScore=0
        if options.sort=='score':
            if len(hits)>0:
                bestScore=hits[0].score
        else:
            for hit in hits:
                if hit.score>bestScore:
                    bestScore=hit.score

        tpScore = bestScore - (options.topPct * bestScore / 100.0)
        minScore = max((tpScore,options.bits))
        logger.debug("Cutoff (%s) is max of bits(%s) and %d%% less than max(%s)" % (minScore, options.bits,options.topPct,bestScore))
    else:
        minScore = options.bits

    # apply filters
    hitCount=0
    hspCounts={}
    hitRegions=[]
    for hit in hits:
        hspCount=hspCounts.get(hit.hit,0)
        logger.debug("hit: %s::%s - score:%s" % (hit.read,hit.hit,hit.score))

        # Simple comparison tests
        if options.format != FRHIT and hit.score<minScore:
            logger.debug("score too low: %r" % hit.score)
            continue
        if options.format != LAST0 and options.evalue is not None and hit.evalue>options.evalue:
            logger.debug("evalue too high: %r" % hit.evalue)
            continue
        if hit.pctid > 0 and hit.pctid<options.pctid:
            logger.debug("pct ID too low: %r < %r" % (hit.pctid,options.pctid))
            continue
        if hit.mlen<options.length:
            logger.debug("hit too short: %r" % hit.mlen)
            continue
        if options.aln is not None and hit.getAln()<options.aln:
            logger.debug("aln fraction too low: %r" % hit.getAln())
            continue
        if options.hitsPerRead>0 and hitCount>=options.hitsPerRead:
            logger.debug("Too many hits")
            continue
        if options.hspsPerHit>0 and hspCount>=options.hspsPerHit:
            logger.debug("Too many HSPs")
            continue

        if options.nonoverlapping:
            # check for previous overlapping hits
            newHitRegions = hit.checkForOverlap(hitRegions)
            if newHitRegions is not None:
                hitRegions = newHitRegions
            else:
                continue

        # increment hit counts
        hspCounts[hit.hit]=hspCount+1
        hitCount+=1

        # print hit
        if returnLines:
            yield hit.getLine(options)
        else:
            yield hit

def addHitTableOptions(parser, defaults={}, flags=['format','filterPct','sortReads']):
    """
    Set up command line arguments for parsing and filtering an m8 file. By default, only the --format and --filterPct options are added, but any of the following can be enabled using the flags= keyword.

    general hit table handling:
        format (GENE, LIZ, BLASTPLUS, LAST, ...)
        sort (sort hits by 'score' or 'evalue')
        sortReads (sort lines by read name)

    filtering options:
        filterPct (aka topPct: minimum pct of best score for other scores. 0, for best score only, 100 for all hits)
        bits (minimum bit score)
        evalue
        pctid
        length (of the alignment)
        aln (fraction of query sequence aligned)
        hitsPerRead
        hspsPerHit
        nonoverlapping

For example, flags=['format','bits','filterPct'] would enable filtering on bit score by cutoff or by percent of the top hit. flags='all', will turn everything on.

    Default values can be changed by passing a dict mapping flags to new default values. Anything in this dict will be added to the flags list.
"""
    # merge explicit defaulst in to flags list
    if flags != 'all':
        flags=set(flags)
        for key in defaults.iterkeys():
            flags.add(key)

    ogroup = OptionGroup(parser, "Hit Table Options", """These options control the parsing and filtering of hit tables (from blast or lastal)""")
    if flags=='all' or 'format' in flags:
        ogroup.add_option('-f','--format', dest='hitTableFormat',
                          default=defaults.get("format",GENE),
                          choices=[GENE,
                                   LIZ,
                                   YANMEI,
                                   LAST0,
                                   BLASTPLUS,
                                   SAM,
                                   HMMSCANDOM,
                                   HMMSCAN,
                                   HMMSEARCHDOM,
                                   HMMSEARCH],
                          help="Format of input table: blast, last, hmmer, gene, yanmei, or liz. Default is %s" % (defaults.get("format",GENE)))
    if flags=='all' or 'filterPct' in flags:
        ogroup.add_option('-F','--filterPct',dest='filterTopPct',
                          default=defaults.get("filterPct",-1),
                          type='int',
                          help="If a positive number is given, only allow hits within this percent of the best hit. Use -1 for no filtering. Use 0 to just take the hit(s) with the best score. If hit table is from last, it will be automatically sorted by read. Sorting by read is slow and should be avoided if you already have a sorted and filtered hit table (e.g. skip this option). Default is %s" % (defaults.get("filterPct",-1)))
    if flags=='all' or 'bits' in flags:
        ogroup.add_option('-B','--bitScore', dest='filterBits',
                          type='int', default=defaults.get("bits",0),
                          help="Minimum bit score to allow. Default: %default")
    if flags=='all' or 'evalue' in flags:
        ogroup.add_option('-E','--evalue', dest='filterEvalue',
                          type='float', default=defaults.get("evalue",None),
                          help="Maximum evalue to allow. Default: %default")
    if flags=='all' or 'pctid' in flags:
        defVal=defaults.get("pctid",0)
        ogroup.add_option('-I','--pctid', dest='filterPctid',
                          type='float', default=defVal,
                          help="Minimum percent identity to allow. Default: %s" % (defVal))
    if flags=='all' or 'length' in flags:
        default=defaults.get("length",0)
        ogroup.add_option('-L','--length', dest='filterLength',
                          type='int', default=default,
                          help="Minimum alignment length to allow. Default: %s" % (default))
    if flags=='all' or 'aln' in flags:
        default=defaults.get("aln",None)
        ogroup.add_option('-N','--aln', dest='filterAln',
                          type='float', default=default,
                          help="Minimum aligned fraction to allow. Default: %s" % (default))
    if flags=='all' or 'hitsPerRead' in flags:
        default=defaults.get("hitsPerRead",0)
        ogroup.add_option('-H','--hitsPerRead', dest='filterHitsPerRead',
                          type='int', default=default,
                          help="Maximum number of hits to allow per read. 0 for all hits. Default: %s" % (default))
    if flags=='all' or 'hspsPerHit' in flags:
        default=defaults.get("hspsPerHit",1)
        ogroup.add_option('-P','--hspsPerHit', dest='filterHspsPerHit',
                          type='int', default=default,
                          help="Maximum number of HSPs to keep per hit. 0 for all HSPs. Default: %s" % (default))
    if flags=='all' or 'nonoverlapping' in flags:
        default=defaults.get("nonoverlapping",False)
        ogroup.add_option('-O','--nonoverlapping', action='store_true',
                          default=default, dest='filterNonoverlapping',
                          help="Ignore hits which overlap higher scoring hits.Default: %s" % (default))
    if flags=='all' or 'sort' in flags:
        default=defaults.get("sort",None)
        ogroup.add_option('-s', '--sort', dest='hitTableSort',
                          default=default, choices=['evalue','score'],
                          help="sort hits for each read by 'evalue' or 'score' before filtering. Secondarily sorted by hit id to make output more deterministic")
    if flags=='all' or 'sortReads' in flags:
        default=defaults.get("sortReads", False)
        ogroup.add_option('-S', '--sortReads', dest='hitTableSortReads',
                          default=default, action='store_true',
                          help="Sort input lines by read before parsing. Only needed if merging multiple results or when processing raw LAST output.  For REALLY large files, this may fill up available memory and be slow. If that happens, run through 'sort' first. (for last: sort -t $'\\t' -k 7,7 -k1rn,1 FILE)")

    parser.add_option_group(ogroup)

def parseCigarString(cigar):
    """
    M alignment match (can be a sequence match or mismatch)
    I insertion to the reference
    D deletion from the reference
    N skipped region from the reference
    S soft clipping (clipped sequences present in SEQ)
    H hard clipping (clipped sequences NOT present in SEQ)
    P padding (silent deletion from padded reference)
    = sequence match
    X sequence mismatch

    So, 5S6M1I4M means the first 5 bases were soft masked, the next 6 match, then an insertion in the query, then 4 matches.
    """
    alen=0
    alenh=0
    alenq=0
    qstart=1
    matches=0
    mismatches=0
    pctid=0

    for cigarBit in cigarRE.findall(cigar):
        bitlen=int(cigarBit[:-1])
        bittype=cigarBit[-1]

        if bittype=='S' or bittype=='H':
            if alenq==0:
                qstart+=bitlen
            else:
                #done reading matches, don't care about unmatched end
                break
        elif bittype=='M':
            alen+=bitlen
            alenh+=bitlen
            alenq+=bitlen
        elif bittype=='I':
            alenq+=bitlen
        elif bittype=='D' or bittype=='N':
            alenh+=bitlen
        elif bittype=='=':
            alen+=bitlen
            alenh+=bitlen
            alenq+=bitlen
            matches+=bitlen
        elif bittype=='X':
            alen+=bitlen
            alenh+=bitlen
            alenq+=bitlen
            mismatches+=bitlen
        else:
            raise Exception("Unrecognized CIGAR tag (%s) in string: %s" % (bittype, cigar))

    if matches>0:
        # (if no matches found, leave pctid at "0" so it is ignored)
        pctid = matches/float(matches+mismatches)
    qend = qstart + alenq - 1
    return (alen, alenh, alenq, qstart, qend, pctid)

def test():
    import sys
    global myAssertEq, myAssertIs
    from edl.test import myAssertEq, myAssertIs

    if len(sys.argv)>1:
        loglevel=logging.DEBUG
    else:
        loglevel=logging.WARN
    logging.basicConfig(stream=sys.stderr, level=loglevel)

    m8data=["001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3	gi|91763278|ref|ZP_01265242.1|	Peptidase family M48 [Candidatus Pelagibacter ubique HTCC1002]	61.7021276595745	94	282	1	57	150	134	4e-30	0.992957746478873\n",
        "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3	gi|71083682|ref|YP_266402.1|	M48 family peptidase [Candidatus Pelagibacter ubique HTCC1062]	61.7021276595745	94	282	1	40	133	134	4e-30	0.992957746478873\n",
        "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3	gi|262277211|ref|ZP_06055004.1|	peptidase family M48 family [alpha proteobacterium HIMB114]	65.9090909090909	88	264	1	63	150	132	9e-30	0.929577464788732\n",
        "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3	gi|254456035|ref|ZP_05069464.1|	peptidase family M48 [Candidatus Pelagibacter sp. HTCC7211]	66.2790697674419	86	258	1	65	150	132	2e-29	0.908450704225352\n",
        "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3	gi|118581678|ref|YP_902928.1|	peptidase M48, Ste24p [Pelobacter propionicus DSM 2379]	51.6129032258064	93	282	4	53	144	108	2e-22	0.982394366197183\n",
        "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3	gi|255534285|ref|YP_003094656.1|	zn-dependent protease with chaperone function [Flavobacteriaceae bacterium 3519-10]	47.3118279569892	93	282	4	55	146	102	1e-20	0.982394366197183\n",
        "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3	gi|317502588|ref|ZP_07960709.1|	M48B family peptidase [Prevotella salivae DSM 15606]	51.6129032258064	93	279	1	85	176	100	7e-20	0.982394366197183\n",
        "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3	gi|325104752|ref|YP_004274406.1|	peptidase M48 Ste24p [Pedobacter saltans DSM 12145]	48.3870967741936	93	279	1	59	150	100	7e-20	0.982394366197183\n",
        "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3	gi|256425464|ref|YP_003126117.1|	peptidase M48 Ste24p [Chitinophaga pinensis DSM 2588]	48.8888888888889	90	273	4	58	146	99.8	9e-20	0.950704225352113\n",
        "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3	gi|299142895|ref|ZP_07036022.1|	peptidase, M48 family [Prevotella oris C735]	50	94	282	1	58	150	99.4	1e-19	0.992957746478873\n"]

    #test1 passthrough
    logging.info("Starting passthrough test")
    m8stream=m8data.__iter__()
    params=FilterParams()
    outs = filterM8Stream(m8stream,params)
    myAssertEq(outs.next(),m8data[0])
    myAssertEq(outs.next(),m8data[1])
    myAssertEq(outs.next(),m8data[2])
    myAssertEq(outs.next(),m8data[3])
    myAssertEq(outs.next(),m8data[4])

    logging.info("Starting best test")
    m8stream=m8data.__iter__()
    params=FilterParams(topPct=0.)
    outs = filterM8Stream(m8stream,params)
    myAssertEq(outs.next(),m8data[0])
    myAssertEq(outs.next(),m8data[1])
    try:
        outs.next()
        sys.exit("There should only be 2 elements!")
    except StopIteration:
        pass

    logging.info("Starting n1 test")
    m8stream=m8data.__iter__()
    params = FilterParams(hitsPerRead=1)
    outs = filterM8Stream(m8stream,params)
    myAssertEq(outs.next(),m8data[0])
    try:
        outs.next()
        sys.exit("There should only be 1 elements!")
    except StopIteration:
        pass

if __name__ == '__main__':
    test()

