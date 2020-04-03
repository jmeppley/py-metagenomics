#! /usr/bin/python
"""
COllection of functions and classes to parse and filter  hit tables from:
    blast
    hmmer
    last
    diamond
    infernal
    SAM
    GFF

"""
import logging
import re
import sys
from edl.util import InputFile, parseExp
logger = logging.getLogger(__name__)

#############
# Constants #
#############
GENE = 'gene'
LIZ = 'liz'
YANMEI = 'yanmei'
BLASTPLUS = 'blast'
FRHIT = 'frhit'
LAST0 = 'last'
HMMSCANDOM = 'hmmscandom'
HMMSEARCHDOM = 'hmmsearchdom'
HMMSCAN = 'hmmscan'
HMMSEARCH = 'hmmsearch'
CMSEARCH = 'cmsearch'
CMSCAN = 'cmscan'
SAM = 'sam'
GFF = 'gff'
formatsWithNoDescription = [LAST0, FRHIT, BLASTPLUS, SAM]
cigarRE = re.compile(r'\d+[^\d]')

#############
# Classes   #
#############


class FilterParams:

    @staticmethod
    def create_from_arguments(arguments, ignore=[], translate={}):
        """
        Translate a NameSpace object created using the
        add_hit_table_arguments below into a FilterParams object.
        The attributes format, sort are expected
        to be in the arguments object with the 'hitTable' prefix.
        The attributes bits, evalue, pctid, aln, length, hits_per_read,
        hsps_per_hit, and nonoverlapping with the prefix 'filter'.
        For example, arguments.hitTableFormat will be copied to
        params.format and arguments.filter_top_pct to params.top_pct.

        translate= should be set to a dictionary mapping attributes
        of the arguments object to the desired attribuets of the
        FilterParams object.

        ignore= should be a list of the standard params to be skipped
        (eg: filter_top_pct).
        """
        params = FilterParams()

        # check the nonstandard one(s)
        for (oname, pname) in translate.items():
            if hasattr(arguments, oname):
                setattr(params, pname, getattr(arguments, oname))

        # get the standard ones
        for param in [
            'bits',
            'evalue',
            'pctid',
            'length',
            'aln',
            'hits_per_read',
            'hsps_per_hit',
            'nonoverlapping',
                'top_pct']:
            oparam = 'filter_' + param
            if oparam not in ignore:
                if hasattr(arguments, oparam):
                    setattr(params, param, getattr(arguments, oparam))
        for param in ['format', 'sort']:
            oparam = 'hitTable' + param[0].upper() + param[1:]
            if oparam not in ignore:
                if hasattr(arguments, oparam):
                    setattr(params, param, getattr(arguments, oparam))

        logging.debug("%r" % (arguments))
        logging.debug("%r" % (params))
        return params

    def __init__(
            self,
            format=GENE,
            top_pct=-1,
            bits=0.0,
            evalue=None,
            pctid=0.0,
            length=0,
            aln=None,
            hits_per_read=0,
            hsps_per_hit=0,
            nonoverlapping=-1,
            sort=None,
            bad_refs=None):
        self.format = format
        self.top_pct = top_pct
        self.bits = bits
        self.evalue = evalue
        self.pctid = pctid
        self.length = length
        self.aln = aln
        self.hits_per_read = hits_per_read
        self.hsps_per_hit = hsps_per_hit
        self.nonoverlapping = nonoverlapping
        self.sort = sort
        self.bad_refs = None

    def __repr__(self):
        return ("FilterParams(format=%r, top_pct=%r, bits=%r, evalue=%r, "
                "pctid=%r, length=%r, aln=%r, hits_per_read=%r, "
                "hsps_per_hit=%r, nonoverlapping=%r sort=%r)" % (
                        self.format, self.top_pct, self.bits, self.evalue,
                        self.pctid, self.length, self.aln, self.hits_per_read,
                        self.hsps_per_hit, self.nonoverlapping, self.sort,
                ))

    """
    def __setattr__(self, name, value):
        if name == 'pctid' and value > 1:
            logger.warning("requested PCTID filter is grater than 1, "
                           "so we're scaling it down by 100x")
            value = value / 100
        super().__setattr__(name, value)
    """


class EmptyHitException(Exception):
    pass


class Hit:
    """
    Object representing a single hit from a search program like blast
    or fr-hit
    """

    @staticmethod
    def getHit(line, options):
        try:
            return Hit(line, options)
        except EmptyHitException:
            return None
        except Exception:
            logger.warn("Error parsing line:\n%s" % (line))
            raise

    def __init__(self, line, options):
        """
        Creates Hit object from a line in a hit table. Options can
        be any class with a 'format' property or a string. The string
        value of options or options.format must be one of the recognized
        formats: 'gene','liz','yanmei','last','frhit'
        """
        self.line = line
        self.setFormat(options)
        self.parseLine(line)

    def __repr__(self):
        return "Hit(%r,%r)" % (self.line.rstrip('\n\r'), self.format)

    def setFormat(self, options):
        if isinstance(options, type(GENE)):
            # just to make this work if instantiated from repr
            self.format = options
        else:
            self.format = options.format

        if self.format == GENE:
            self.parseLine = self.parseGeneLine
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
        elif self.format == CMSEARCH:
            self.parseLine = self.parseCmSearchLine
        elif self.format == CMSCAN:
            self.parseLine = self.parseCmScanLine
        elif self.format == SAM:
            self.parseLine = self.parseSamLine
        elif self.format == GFF:
            self.parseLine = self.parseGFFLine
            self.to_gff = lambda self: self.line
        else:
            sys.exit("Unknown format: %s" % (self.format))

    def parseGeneLine(self, line):
        cells = line.rstrip('\n\r').split('\t')
        self.read = cells[0]
        self.readDesc = cells[1]
        self.hit = cells[2]
        self.hitDesc = cells[3]
        self.pctid = float(cells[4])
        self.mlen = int(cells[5])
        self.qstart = int(cells[6])
        self.qend = int(cells[7])
        self.hstart = int(cells[8])
        self.hend = int(cells[9])
        self.score = float(cells[10])
        self.evalue = parseExp(cells[11])
        self.aln = float(cells[12])

    def parseLizLine(self, line):
        cells = line.rstrip('\n\r').split('\t')
        self.read = cells[0]
        self.hit = cells[1]
        self.hitDesc = cells[2]
        try:
            self.pctid = float(cells[3])
        except ValueError:
            # leave unset if it's not a float
            pass
        self.mlen = int(cells[4])
        self.qstart = int(cells[5])
        self.qend = int(cells[6])
        self.hstart = int(cells[7])
        self.hend = int(cells[8])
        self.score = float(cells[9])
        self.evalue = parseExp(cells[10])
        self.aln = float(cells[11])

    def parseSamLine(self, line):
        if line[0] == '@':
            raise EmptyHitException("reference sequence line")
        cells = line.rstrip('\n\r').split('\t')
        self.read = cells[0]
        self.hit = cells[2]
        if self.hit == '*':
            raise EmptyHitException("No match")
        self.cigar = cells[5]
        (alen, alenh, alenq, qstart, qend, pctid) = \
            parseCigarString(self.cigar)
        self.mlen = alen
        self.qstart = qstart
        self.qend = qend
        self.hstart = int(cells[3])
        self.hend = self.hstart + alenh - 1
        self.qlen = len(cells[9])
        self.aln = float(self.mlen) / float(self.qlen)
        self.score = None
        for tagstr in cells[11:]:
            if tagstr.startswith('AS'):
                self.score = float(tagstr.split(":")[2])
            if tagstr.startswith('MD:Z:'):
                if pctid == 0:
                    pctid = get_alignment_percent_identity(tagstr[4:])
        if self.score is None:
            raise Exception("No score (AS tag) found in line:\n%s" % (line))
        if pctid != 0:
            self.pctid = pctid

    # target name        accession   tlen query name           accession
    # qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias
    # from    to  from    to  from    to  acc description of target
    def parseHmmScanDomLine(self, line):
        if line[0] == '#':
            raise EmptyHitException("Comment line")
        cells = line.rstrip('\n\r').split()
        self.read = cells[3]
        qlen = int(cells[5])
        self.hit = cells[0]
        # self.hlen = int(cells[2])
        self.evalue = float(cells[6])
        self.score = float(cells[7])
        self.hstart = int(cells[15])
        self.hend = int(cells[16])
        self.qstart = int(cells[17])
        self.qend = int(cells[18])
        self.hitDesc = cells[22]
        self.mlen = 1 + self.qend - self.qstart
        self.aln = self.mlen / float(qlen)

    # target name        accession  query name           accession    E-value
    # score  bias ue  score  bias   exp reg clu  ov env dom rep inc
    # description of target
    def parseHmmScanLine(self, line):
        if line[0] == '#':
            raise EmptyHitException("Comment line")
        cells = line.rstrip('\n\r').split()
        self.read = cells[2]
        self.hit = cells[0]
        self.evalue = float(cells[4])
        self.score = float(cells[5])
        self.hitDesc = cells[18]

    # target name        accession   tlen query name           accession
    # qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias
    # from    to  from    to  from    to  acc description of target
    def parseHmmSearchDomLine(self, line):
        if line[0] == '#':
            raise EmptyHitException("Comment line")
        cells = line.rstrip('\n\r').split()
        self.read = cells[0]
        qlen = int(cells[2])
        self.hit = cells[3]
        self.hitDesc = cells[4]
        # self.hlen = int(cells[5])
        self.evalue = float(cells[6])
        self.score = float(cells[7])
        self.hstart = int(cells[15])
        self.hend = int(cells[16])
        self.qstart = int(cells[17])
        self.qend = int(cells[18])
        self.readDesc = cells[22]
        self.mlen = 1 + self.qend - self.qstart
        self.aln = self.mlen / float(qlen)

    # target name        accession  query name           accession    E-value
    # score  bias ue  score  bias   exp reg clu  ov env dom rep inc
    # description of target
    def parseHmmSearchLine(self, line):
        if line[0] == '#':
            raise EmptyHitException("Comment line")
        cells = line.rstrip('\n\r').split()
        self.read = cells[0]
        self.hit = cells[2]
        self.evalue = float(cells[4])
        self.score = float(cells[5])
        self.readDesc = cells[18]

    # target name           accession query name           accession mdl mdl
    # from   mdl to seq from   seq to strand trunc pass   gc  bias  score
    # E-value inc description of target
    def parseCmSearchLine(self, line):
        if line[0] == '#':
            raise EmptyHitException("Comment line")
        cells = line.rstrip('\n\r').split()
        self.read = cells[0]
        self.hit = cells[3]
        self.hitDesc = cells[2]
        self.hstart = int(cells[5])
        self.hend = int(cells[6])
        self.qstart = int(cells[7])
        self.qend = int(cells[8])
        self.strand = cells[9]
        self.evalue = float(cells[15])
        self.score = float(cells[14])
        self.readDesc = cells[17]
        self.mlen = self.hend - self.hstart + 1

    # target name         accession query name           accession mdl mdl
    # from   mdl to seq from   seq to strand trunc pass   gc  bias  score
    # E-value inc description of target
    def parseCmScanLine(self, line):
        if line[0] == '#':
            raise EmptyHitException("Comment line")
        cells = line.rstrip('\n\r').split()
        self.read = cells[2]
        self.hit = cells[1]
        self.hitDesc = cells[0]
        self.hstart = int(cells[5])
        self.hend = int(cells[6])
        self.qstart = int(cells[7])
        self.qend = int(cells[8])
        self.strand = cells[9]
        self.evalue = float(cells[15])
        self.score = float(cells[14])
        self.readDesc = cells[17]
        self.mlen = self.hend - self.hstart + 1

    # score name1   start1  alnSize1        strand1 seqSize1        name2
    # start2  alnSize2   strand2 seqSize2        blocks
    def parseLastalLine(self, line):
        # logger.debug(line)
        if line[0] == '#':
            raise EmptyHitException("Comment line")
        cells = line.rstrip('\n\r').split('\t')
        self.score = float(cells[0])
        self.hit = cells[1]

        hmlen = int(cells[3])
        hlen = int(cells[5])
        if cells[4] == '+':
            self.hstart = int(cells[2]) + 1
            self.hend = self.hstart + hmlen - 1
        else:
            # start was in reverse strand
            self.hstart = hlen - int(cells[2])
            self.hend = self.hstart - hmlen + 1

        self.read = cells[6]
        qmlen = int(cells[8])
        qlen = int(cells[10])
        if cells[9] == '+':
            self.qstart = int(cells[7]) + 1
            self.qend = self.qstart + qmlen - 1
        else:
            # start was in reverse strand
            self.qstart = qlen - int(cells[7])
            self.qend = self.qstart - qmlen + 1
        self.aln = qmlen / float(qlen)
        if cells[9] == cells[4]:
            self.strand = "+"
        else:
            self.strand = "-"

        # some versions have evalues in the last few spots (eg: E=2.1e-09)
        self.evalue = float(cells[13][2:].strip()) if len(cells) > 13 else None
        self.mlen = computeLastHitValues(cells[11])
        self.hitDesc = None
        logger.debug("Span: %d-%d" % (self.qstart, self.qend))

    def parseFrHitLine(self, line):
        cells = line.rstrip('\n\r').split('\t')
        self.read = cells[0]
        self.evalue = parseExp(cells[2])
        self.score = None
        self.mlen = cells[3]
        self.qstart = int(cells[4])
        self.qend = int(cells[5])
        self.pctid = float(cells[7])
        self.hit = cells[8]
        self.hstart = int(cells[9])
        self.hend = int(cells[10])

    def parseBlastPlusLine(self, line):
        if line[0] == '#':
            raise EmptyHitException("Comment line")
        cells = line.rstrip('\n\r').split('\t')
        self.read = cells[0]
        self.hit = cells[1]
        self.pctid = float(cells[2])
        self.mlen = int(cells[3])
        self.mismatch = int(cells[4])
        self.gaps = int(cells[5])
        self.qstart = int(cells[6])
        self.qend = int(cells[7])
        self.hstart = int(cells[8])
        self.hend = int(cells[9])
        self.evalue = parseExp(cells[10])
        self.score = float(cells[11])

    def parseYanmeiLine(self, line):
        cells = line.rstrip('\n\r').split('\t')
        self.read = cells[0]
        self.hit = cells[1]
        self.pctid = float(cells[2])
        self.mlen = int(cells[3])
        self.mismatch = int(cells[4])
        self.gaps = int(cells[5])
        self.qstart = int(cells[6])
        self.qend = int(cells[7])
        self.hstart = int(cells[8])
        self.hend = int(cells[9])
        self.evalue = parseExp(cells[10])
        self.score = float(cells[11])
        self.hitDesc = cells[12]

    def parseGFFLine(self, line):
        if line[0] == '#':
            raise EmptyHitException("Comment")
        cells = line.rstrip('\n\r').split('\t')
        self.read = cells[0]
        self.source = cells[1]
        self.hit_type = cells[2]
        self.qstart = int(cells[3])
        self.qend = int(cells[4])
        self.score = float(cells[5])
        self.strand = cells[6]
        hit_data = dict([kv.split('=')
                         for kv in cells[8].strip(';').split(';')])
        if "ID" in hit_data:
            self.hit = hit_data['ID']
        elif "Name" in hit_data:
            self.hit = hit_data['Name']
        elif "Target" in hit_data:
            self.hit = hit_data['Target'].split()[0]
        self.hitDesc = hit_data.get('product', cells[8])
        self.evalue = self.score
        self.mlen = self.qend + 1 - self.qstart

    def getAln(self):
        try:
            return self.aln
        except AttributeError:
            sys.exit("Cannot calculate alignment percentage from"
                     "data in this m8 format")

    def checkForOverlap(self, regions, buffer):
        # make sure start is smaller than end
        start = min(self.qstart, self.qend)
        end = max(self.qstart, self.qend)

        # adjust start/end to account for buffer
        buf_st = start + buffer
        buf_en = end - buffer

        for i in range(len(regions)):
            occupiedRange = regions[i]
            # hit cannot intersect an used range
            if (buf_st >= occupiedRange[1] or buf_en <= occupiedRange[0]):
                # does not overlap this range (try next range)
                continue
            else:
                # overlaps. We are done here
                return ((start, end), occupiedRange)
        return ((start, end), None)

    def checkForOverlapAndAdd(self, regions, buffer):
        """
        check to see if this hit overlaps any already hit region.

        Regions must overlap by at least "buffer" bases to count.

        The regions array should be a list of (Start,end) pairs
        indicating regions already hit.
        """
        logger.debug(
            "Looking for [%d,%d] in %s" %
            (self.qstart, self.qend, regions))

        # compare with all used ranges
        hit_span, overlap_region = self.checkForOverlap(regions, buffer)
        if overlap_region is not None:
            return None

        # we get here if there was no overlap
        regions.append(hit_span)
        return regions

    def getLine(self, options):
        return self.line

    def to_gff(self):
        """
        An attempt to export any hit to GFF. Only partially tested
        """

        # The core data in a  GFF file:
        sequence_id = self.read
        source = self.format
        feature = self.hitDesc
        start = self.qstart
        end = self.qend
        score = self.score
        try:
            strand = self.strand
        except AttributeError:
            strand = "+" if start <= end else "-"
        # Ignore phase for now
        phase = '.'

        # create attribute field for attributes we have
        attributes = {}
        attributes['Target'] = " ".join(
            str(v) for v in [
                self.hit,
                self.hstart,
                self.hend])
        try:
            attributes['Cigar'] = self.cigar
        except AttributeError:
            pass
        attributes = ';'.join(["{0}={1}".format(k, attributes[k])
                               for k in attributes])
        # create line
        gff_line = '\t'.join(str(v) for v in
                             [sequence_id,
                              source,
                              feature,
                              start,
                              end,
                              score,
                              strand,
                              phase,
                              attributes]) \
            + '\n'
        return gff_line

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
    matchLen = 0
    for segment in blocks.split(','):
        try:
            matchLen += int(segment)
        except ValueError:
            (hml, qml) = segment.split(":")
            mml = max(int(hml), int(qml))
            matchLen += mml

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
            hitCol = 3
        else:
            hitCol = 2
    elif format == LAST0:
        if useDesc:
            raise Exception(
                "lastal does not report the hit description, sorry")
        else:
            hitCol = 1
    elif format == LIZ:
        if useDesc:
            hitCol = 2
        else:
            hitCol = 1
    elif format == YANMEI:
        if useDesc:
            hitCol = 12
        else:
            hitCol = 1
    elif format == BLASTPLUS:
        if useDesc:
            raise Exception(
                "The default blast table does not keep the hit description")
        else:
            hitCol = 1
    elif format == SAM:
        if useDesc:
            raise Exception("The SAM format does not keep the hit description")
        else:
            hitCol = 2
    else:
        sys.exit("I'm sorry, I don't understand the format: %s" % (format))

    return hitCol


def filterM8(instream, outstream, params, to_gff=False):
    """
    Filter instream and write output to outstream
    """
    logger.debug("blastm8.filterM8: reading from {}".format(instream.name))
    line_count = 0
    if to_gff and params.format != GFF:
        for read, hits in filterM8Stream(instream, params, return_lines=False):
            for hit in hits:
                line_count += 1
                outstream.write(hit.to_gff())
    else:
        for line in filterM8Stream(instream, params, return_lines=True):
            line_count += 1
            outstream.write(line)
    logger.debug("blastm8.filterM8: wrote {} lines".format(line_count))


def sortLines(instream):
    logger.info("Sorting input lines")
    lines = []
    for line in instream:
        lines.append(line)
    lines.sort()
    logger.debug("Done sorting")
    return lines.__iter__()


def getHitStream(instream, options):
    """
    Simply parse lines into Hits one at a time and return them
    """
    for line in instream:
        if len(line) == 0:
            continue
        try:
            hit = Hit.getHit(line, options)
        except Exception:
            print("ERROR parsing line")
            print(line)
            raise
        if hit is not None:
            yield hit


def generate_hits(hit_table, format=BLASTPLUS, **filter_args):
    """
    yeilds read,hit_iterator tuples

    default format is BLAST, change with format=format
    See class FilterParams for other arguments
    """
    with InputFile(hit_table) as m8stream:
        params = FilterParams(format=format, **filter_args)
        for read, hits in filterM8Stream(m8stream,
                                         params,
                                         return_lines=False,
                                         ):
            yield read, hits
        logger.debug("Read %d lines from %s table", m8stream.lines, format)


def filterM8Stream(instream, options, return_lines=True):
    """
    return an iterator over the lines in given input stream that pass filter
    """
    logger.debug("Processing stream: %s" % instream)

    current_read = None
    hits = []
    logger.debug(repr(options))

    def build_items_to_generate(current_read, hits,
                                sort=options.sort,
                                needs_filter=doWeNeedToFilter(options),
                                return_lines=return_lines):
        if current_read is not None:
            logger.debug(
                "processing %d hits for %s" %
                (len(hits), current_read))
            if options.sort is not None:
                sortHits(hits, options.sort)
            if needs_filter:
                hits = filterHits(hits, options)
            if return_lines:
                for hit in hits:
                    yield hit.line
            else:
                if needs_filter:
                    hits = list(hits)
                if len(hits) > 0:
                    yield (current_read, hits)

    for line_hit in getHitStream(instream, options):
        if line_hit.read != current_read:
            for item in build_items_to_generate(current_read, hits):
                yield item
            hits = []
            current_read = line_hit.read

        hits.append(line_hit)

    for item in build_items_to_generate(current_read, hits):
        yield item


def sortHits(hits, sortType):
    """
    take all the hits for a read as a list and sort them by score and hit name
    """

    logger.debug("Sorting hits (byScore = %s)" % (sortType))

    # set up sort key
    if sortType == 'evalue':
        sort_key = get_sort_by_evalue_key
    elif sortType == 'pctid':
        sort_key = get_sort_by_pctid_key
    else:
        sort_key = get_sort_by_score_key

    # sort in place
    hits.sort(key=sort_key)


def get_sort_by_evalue_key(hit):
    return (hit.evalue, hit.hit)


def get_sort_by_pctid_key(hit):
    return (hit.pctid * -1, hit.hit)


def get_sort_by_score_key(hit):
    return (hit.score * -1, hit.hit)


def doWeNeedToFilter(options):
    """ we can skip the filter step if we're going to let everything through
    """
    if options.top_pct >= 0:
        return True
    if options.bits > 0:
        return True
    if options.evalue is not None:
        return True
    if options.pctid > 0:
        return True
    if options.length > 0:
        return True
    if options.hits_per_read > 0:
        return True
    if options.hsps_per_hit > 0:
        return True
    if options.nonoverlapping >= 0:
        return True
    if options.bad_refs:
        return True
    return False


def filterHits(hits, options):
    if options.bad_refs:
        logger.debug("REmoving bad_refs from hits")
        hits = [h for h in hits if h.hit not in options.bad_refs]

    # A top_pct cutoff, requires finding the top score first
    if options.top_pct >= 0:
        # get the best score in this set of hits
        bestScore = 0
        if options.sort == 'score':
            if len(hits) > 0:
                bestScore = hits[0].score
        else:
            for hit in hits:
                if hit.score > bestScore:
                    bestScore = hit.score

        tpScore = bestScore - (options.top_pct * bestScore / 100.0)
        minScore = max((tpScore, options.bits))
        logger.debug(
            "Cutoff (%s) is max of bits(%s) and %d%% less than max(%s)" %
            (minScore, options.bits, options.top_pct, bestScore))
    else:
        minScore = options.bits

    # apply filters
    counted_hits = set()
    hit_count = 0
    hsp_counts = {}
    hit_regions = []
    for hit in hits:
        hsp_count = hsp_counts.get(hit.hit, 0)
        logger.debug("hit: %s::%s - score:%s" % (hit.read, hit.hit, hit.score))

        # Simple comparison tests
        if options.format != FRHIT and hit.score < minScore:
            logger.debug("score too low: %r" % hit.score)
            continue
        if options.format != LAST0\
                and options.evalue is not None\
                and hit.evalue > options.evalue:
            logger.debug("evalue too high: %r" % hit.evalue)
            continue

        # PCTID
        try:
            if hit.pctid < options.pctid:
                logger.debug("pct ID too low: %r < %r" %
                             (hit.pctid, options.pctid))
                continue
        except AttributeError:
            if options.pctid > 0:
                raise Exception(
                    "This hit type (%s) does not have a PCTID defined."
                    "You cannot filter by PCTID" %
                    (hit.format))

        if abs(hit.mlen) < options.length:
            logger.debug("hit too short: %r" % hit.mlen)
            continue
        if options.aln is not None and hit.getAln() < options.aln:
            logger.debug("aln fraction too low: %r" % hit.getAln())
            continue
        if options.hits_per_read > 0 \
           and hit_count >= options.hits_per_read \
           and hit.hit not in counted_hits:
            logger.debug("Too many hits")
            continue
        if options.hsps_per_hit > 0 and hsp_count >= options.hsps_per_hit:
            logger.debug("Too many HSPs")
            continue

        if options.nonoverlapping >= 0:
            # check for previous overlapping hits
            new_hit_regions = hit.checkForOverlapAndAdd(hit_regions,
                                                        options.nonoverlapping)
            if new_hit_regions is not None:
                hit_regions = new_hit_regions
            else:
                continue

        # increment hit counts
        hsp_counts[hit.hit] = hsp_count + 1
        counted_hits.add(hit.hit)
        hit_count = len(counted_hits)

        # print hit
        yield hit


def add_hit_table_arguments(parser,
                            defaults={},
                            flags=['format', 'filter_top_pct']):
    """
    Set up command line arguments for parsing and filtering an m8 file.
    By default, only the --format and --filter_top_pct options are added, but
    any of the following can be enabled using the flags= keyword.

    general hit table handling:
        format (GENE, LIZ, BLASTPLUS, LAST, ...)
        sort (sort hits by 'score', 'pctid', or 'evalue')

    filtering options:
        filter_top_pct (aka top_pct: minimum pct of best score for other
                     scores. 0, for best score only, 100 for all hits)
        bits (minimum bit score)
        evalue
        pctid
        length (of the alignment)
        aln (fraction of query sequence aligned)
        hits_per_read
        hsps_per_hit
        nonoverlapping

    For example, flags=['format','bits','filter_top_pct'] would enable
    filtering on bit score by cutoff or by percent of the top hit.
    flags='all', will turn everything on.

    Default values can be changed by passing a dict mapping flags to
    new default values. Anything in this dict will be added to the flags list.
    """
    # merge explicit defaulst in to flags list
    if flags != 'all':
        flags = set(flags)
        for key in defaults.keys():
            flags.add(key)

    agroup = parser.add_argument_group(
        "Hit Table Options",
        """These options control the parsing and filtering of hit
           tables (from blast or lastal)""")
    if flags == 'all' or 'format' in flags:
        agroup.add_argument(
            '-f',
            '--format',
            dest='hitTableFormat',
            default=defaults.get(
                "format",
                GENE),
            choices=[
                GENE,
                LIZ,
                YANMEI,
                LAST0,
                BLASTPLUS,
                SAM,
                GFF,
                CMSEARCH,
                CMSCAN,
                HMMSCANDOM,
                HMMSCAN,
                HMMSEARCHDOM,
                HMMSEARCH],
            help="Format of input table: blast, last, hmmer, gene, "
                 "yanmei, or liz. Default is %s" % (defaults.get("format",
                                                                 GENE)))
    if flags == 'all' or 'filter_top_pct' in flags:
        agroup.add_argument(
            '-F',
            '--filter_top_pct',
            dest='filter_top_pct',
            default=defaults.get(
                "filter_top_pct",
                -1),
            type=int,
            help="If a positive number is given, only allow hits within "
                 "this percent of the best hit. Use -1 for no filtering. "
                 "Use 0 to just take the hit(s) with the best score. "
                 "Default is {}"
                 .format(defaults.get("filter_top_pct", -1)))
    if flags == 'all' or 'bits' in flags:
        agroup.add_argument('-B', '--bitScore', dest='filter_bits',
                            type=int, default=defaults.get("bits", 0),
                            help="Minimum bit score to allow. Default: \
                          {}".format(defaults.get('bits', 0)))
    if flags == 'all' or 'evalue' in flags:
        agroup.add_argument('-E', '--evalue', dest='filter_evalue',
                            type=float, default=defaults.get("evalue", None),
                            help="Maximum evalue to allow. Default: \
                          {}".format(defaults.get('evalue', None)))
    if flags == 'all' or 'pctid' in flags:
        defVal = defaults.get("pctid", 0)
        agroup.add_argument(
            '-I',
            '--pctid',
            dest='filter_pctid',
            type=float,
            default=defVal,
            help=("Minimum percent identity to allow in range 0 to 100. "
                  "Default: %s" % (defVal))
        )
    if flags == 'all' or 'length' in flags:
        default = defaults.get("length", 0)
        agroup.add_argument(
            '-L',
            '--length',
            dest='filter_length',
            type=int,
            default=default,
            help="Minimum alignment length to allow. Default: %s" %
            (default))
    if flags == 'all' or 'aln' in flags:
        default = defaults.get("aln", None)
        agroup.add_argument(
            '-N',
            '--aln',
            dest='filter_aln',
            type=float,
            default=default,
            help="Minimum aligned fraction to allow. Default: %s" %
            (default))
    if flags == 'all' or 'hits_per_read' in flags:
        default = defaults.get("hits_per_read", 0)
        agroup.add_argument(
            '-H',
            '--hits_per_read',
            dest='filter_hits_per_read',
            type=int,
            default=default,
            help="Maximum number of hits to allow per read. 0 for all hits. "
                 "Default: %s" % (default))
    if flags == 'all' or 'hsps_per_hit' in flags:
        default = defaults.get("hsps_per_hit", 0)
        agroup.add_argument(
            '-P',
            '--hsps_per_hit',
            dest='filter_hsps_per_hit',
            type=int,
            default=default,
            help="Maximum number of HSPs to keep per hit. 0 for all HSPs. "
                 "Default: %s" % (default))
    if flags == 'all' or 'nonoverlapping' in flags:
        default = defaults.get("nonoverlapping", -1)
        agroup.add_argument(
            '-U',
            '--nonoverlapping',
            nargs="?",
            type=int,
            const=0,
            default=default,
            action='store',
            dest='filter_nonoverlapping',
            help="Ignore hits which overlap higher scoring hits. Optional "
                 "integer value may be specified to ignore small overlaps. "
                 "A negative value allows all overlaps. "
                 "Default: %s" % (default))
    if flags == 'all' or 'sort' in flags:
        default = defaults.get("sort", None)
        agroup.add_argument(
            '-s',
            '--sort',
            dest='hitTableSort',
            default=default,
            choices=[
                'evalue',
                'pctid',
                'score'],
            help="sort hits for each read by 'evalue', 'pctid' or 'score' "
                 "before "
                 "filtering. Secondarily sorted by hit id to make output "
                 "more deterministic")


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

    So, 5S6M1I4M means the first 5 bases were soft masked, the next 6 match,
    then an insertion in the query, then 4 matches.
    """
    alen = 0
    alenh = 0
    alenq = 0
    qstart = 1
    matches = 0
    mismatches = 0
    pctid = 0

    for cigarBit in cigarRE.findall(cigar):
        bitlen = int(cigarBit[:-1])
        bittype = cigarBit[-1]

        if bittype == 'S' or bittype == 'H':
            if alenq == 0:
                qstart += bitlen
            else:
                # done reading matches, don't care about unmatched end
                break
        elif bittype == 'M':
            alen += bitlen
            alenh += bitlen
            alenq += bitlen
        elif bittype == 'I':
            alenq += bitlen
        elif bittype == 'D' or bittype == 'N':
            alenh += bitlen
        elif bittype == '=':
            alen += bitlen
            alenh += bitlen
            alenq += bitlen
            matches += bitlen
        elif bittype == 'X':
            alen += bitlen
            alenh += bitlen
            alenq += bitlen
            mismatches += bitlen
        else:
            raise Exception(
                "Unrecognized CIGAR tag (%s) in string: %s" %
                (bittype, cigar))

    if matches > 0:
        # (if no matches found, leave pctid at "0" so it is ignored)
        pctid = matches / float(matches + mismatches)
    qend = qstart + alenq - 1
    return (alen, alenh, alenq, qstart, qend, pctid)


capturing_digits_re = re.compile(r'(\d+)')


def get_alignment_percent_identity(mdx_string):
    """
    Use the MD:Z string returned by BWA to get percent ID
    """
    matches = 0
    mismatches = 0
    for chunk in capturing_digits_re.split(mdx_string):
        if len(chunk) == 0:
            # first and last elements are often empty strings. Ignore them
            continue
        if chunk.startswith('^'):
            # this is an deletion, irrelevant for pctid
            continue
        try:
            # is it an integer? it's the nuimber of matches
            matches += int(chunk)
        except ValueError:
            # otherwise, it's the string of the reference that was mismatched
            mismatches += len(chunk)

    return 100 * float(matches) / float(matches + mismatches)


def setup_tests():
    import sys
    global myAssertEq, myAssertIs
    from edl.test import myAssertEq, myAssertIs

    if len(sys.argv) > 1:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.WARN
    logging.basicConfig(stream=sys.stderr, level=loglevel)


def test():
    setup_tests()
    m8data = ["001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3"
              "	gi|91763278|ref|ZP_01265242.1|	Peptidase family"
              " M48 [Candidatus Pelagibacter ubique HTCC1002]"
              "	61.7021276595745	94	282	1	57	150	134	4e-30"
              "	0.992957746478873\n",
              "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3"
              "	gi|71083682|ref|YP_266402.1|"
              "	M48 family peptidase [Candidatus Pelagibacter ubique"
              " HTCC1062]	61.7021276595745	94	282	1	40	133	134"
              "	4e-30	0.992957746478873\n",
              "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3"
              "	gi|262277211|ref|ZP_06055004.1|	peptidase family M48 family"
              " [alpha proteobacterium HIMB114]	65.9090909090909	88	264"
              "	1	63	150	132	9e-30	0.929577464788732\n",
              "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3"
              "	gi|254456035|ref|ZP_05069464.1|	peptidase family M48"
              " [Candidatus Pelagibacter sp. HTCC7211]	66.2790697674419"
              "	86	258	1	65	150	132	2e-29	0.908450704225352\n",
              "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3"
              "	gi|118581678|ref|YP_902928.1|	peptidase M48, Ste24p"
              " [Pelobacter propionicus DSM 2379]	51.6129032258064	93"
              "	282	4	53	144	108	2e-22	0.982394366197183\n",
              "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3"
              "	gi|255534285|ref|YP_003094656.1|	zn-dependent protease"
              " with chaperone function [Flavobacteriaceae bacterium"
              " 3519-10]	47.3118279569892	93	282	4	55	146	102"
              "	1e-20	0.982394366197183\n",
              "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3"
              "	gi|317502588|ref|ZP_07960709.1|	M48B family peptidase"
              " [Prevotella salivae DSM 15606]	51.6129032258064	93"
              "	279	1	85	176	100	7e-20	0.982394366197183\n",
              "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3"
              "	gi|325104752|ref|YP_004274406.1|	peptidase M48 Ste24p"
              " [Pedobacter saltans DSM 12145]	48.3870967741936	93"
              "	279	1	59	150	100	7e-20	0.982394366197183\n",
              "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3"
              "	gi|256425464|ref|YP_003126117.1|	peptidase M48 Ste24p"
              " [Chitinophaga pinensis DSM 2588]	48.8888888888889"
              "	90	273	4	58	146	99.8	9e-20	0.950704225352113\n",
              "001598_1419_3101	H186x25M  length=284 uaccno=E3N7QM101DQXE3"
              "	gi|299142895|ref|ZP_07036022.1|	peptidase, M48 family"
              " [Prevotella oris C735]	50	94	282	1	58	150	99.4"
              "	1e-19	0.992957746478873\n"]

    # test1 passthrough
    logging.info("Starting passthrough test")
    m8stream = m8data.__iter__()
    params = FilterParams()
    outs = filterM8Stream(m8stream, params, return_lines=True)
    myAssertEq(next(outs), m8data[0])
    myAssertEq(next(outs), m8data[1])
    myAssertEq(next(outs), m8data[2])
    myAssertEq(next(outs), m8data[3])
    myAssertEq(next(outs), m8data[4])

    logging.info("Starting best test")
    m8stream = m8data.__iter__()
    params = FilterParams(top_pct=0.)
    outs = filterM8Stream(m8stream, params)
    myAssertEq(next(outs), m8data[0])
    myAssertEq(next(outs), m8data[1])
    try:
        next(outs)
        sys.exit("There should only be 2 elements!")
    except StopIteration:
        pass

    logging.info("Starting n1 test")
    m8stream = m8data.__iter__()
    params = FilterParams(hits_per_read=1)
    outs = filterM8Stream(m8stream, params)
    myAssertEq(next(outs), m8data[0])
    try:
        next(outs)
        sys.exit("There should only be 1 elements!")
    except StopIteration:
        pass


def test_gff():
    setup_tests()

    line = 'KM282-20-02b-5_c283151\tcsearch\ttRNA\t303\t233\t' + \
           '40.6\t-\t.\tTarget=RF00005 2 70\n'
    hit = Hit(line, 'gff')
    myAssertEq(hit.read, 'KM282-20-02b-5_c283151')
    myAssertEq(hit.score, 40.6)
    myAssertEq(hit.hit, 'RF00005')

    line = 'KM282-20-02a-100_c12273\tbarrnap:0.7\trRNA\t9\t772\t' + \
           '6.8e-41\t+\t.\tName=12S_rRNA;product=12S ribosomal RNA\n'
    hit = Hit(line, 'gff')
    myAssertEq(hit.read, 'KM282-20-02a-100_c12273')
    myAssertEq(hit.evalue, 6.8e-41)
    myAssertEq(hit.hit, '12S_rRNA')

    line = 'KM282-20-02a-100_c1\tProdigal_v2.6.2\tCDS\t309\t686\t' + \
           '53.3\t-\t0\tID=1_2;partial=00;start_type=ATG;rbs_motif=' + \
           'AGGA;rbs_spacer=5-10bp;gc_cont=0.381;conf=100.00;score=' + \
           '54.54;cscore=45.68;sscore=8.86;rscore=5.50;uscore=-0.63;' + \
           'tscore=2.76;\n'
    hit = Hit(line, 'gff')
    myAssertEq(hit.read, 'KM282-20-02a-100_c1')
    myAssertEq(hit.score, 53.3)
    myAssertEq(hit.hit, '1_2')


if __name__ == '__main__':
    test()
