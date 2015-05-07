import re, logging, sys, os
from edl.expressions import accessionRE
logger=logging.getLogger(__name__)

class LineCounter():
    """
    File-like wrapper for a file-like object that counts the number of lines returned.
    """
    def __init__(self, stream):
        self.rawStream=stream
        self.lines=0

    def __iter__(self):
        return self

    def next(self):
        line = self.rawStream.next()
        self.lines+=1
        return line

def checkNoneOption(value):
    """
    Check if a given option has been set to the string 'None' and return None if so. Otherwise return the value unchanged
    """
    if value is None:
        return value
    if value=='None':
        return None
    if isinstance(value,list):
        if len(value)==1 and value[0]=='None':
            return None
    return value

def countBasesInFasta(fastaFile):
    """
    Given a fasta file, return a dict where the number of records and the total number of bases are given by 'records' and 'bases' respectively.
    """
    recordRE=re.compile(r'^>')
    whiteSpaceRE=re.compile(r'\s+')
    totalBases=0
    totalSeqs=0
    with open(fastaFile) as f:
        for line in f:
            if recordRE.match(line):
                totalSeqs+=1
                continue
            totalBases+=len(whiteSpaceRE.sub('',line))

    return {'records':totalSeqs,'bases':totalBases}

urlRE=re.compile(r'[a-z]+\:\/\/')
def openInputFile(infile, *args):
    """
    return input stream. Allow for text, gzipped text, or standard input if None given.
    """
    if infile is None:
        logging.info("Reading input from STDIN")
        return sys.stdin

    if isinstance(infile, str):
        if urlRE.match(infile):
            import urllib2
            return urllib2.urlopen(infile)
        if len(infile)>3 and infile[-3:]=='.gz':
            import gzip
            return gzip.GzipFile(infile,'rb')
        elif len(infile)>4 and infile[-4:]=='.bz2':
            import bz2
            return bz2.BZ2File(infile,'rb')
        else:
            return open(infile,'rU')
    else:
        return infile

def parseExp(string):
    try:
        return float(string)
    except ValueError:
        return float('1'+string)

class Stack():
    """
    Created a linked list for fun
    """
    class Item():
        def __init__(self, item, prev, next):
            self.item = item
            self.prev = prev
            self.next = next

        def getRepr(self):
            """ return a repr style string of this item and everything after it """
            thisRepr = repr(self.item)
            if self.next is not None:
                return thisRepr + ", " + self.next.getRepr()
            else:
                return thisRepr

        def __repr__(self):
            return repr(self.item)

        def __str__(self):
            return str(self.item)

    def __init__(self,data=None):
        self._first = None
        self._last = None
        if data is not None:
            for val in data:
                self.append(val)

    def pop(self):
        ret = self._last
        if ret is not None:
            self._last = ret.prev
            self._last.next = None
            if self._last is None:
                self._first = None
        return ret.item

    def append(self, item):
        new = Stack.Item(item,self._last,None)
        if self._last is None:
            self._first = new
        else:
            self._last.next = new
        self._last = new

    def shift(self, item):
        new = Stack.Item(item,None,self._first)
        if self._first is None:
            self._last=new
        else:
            self._first.prev = new
        self._first = new

    def unshift(self):
        ret = self._first
        if ret is not None:
           self._first = ret.next
           self._first.prev = None
           if self._first is None:
               self._last = None
        return ret.item

    def __len__(self):
       count = 0
       item = self._first
       while item is not None:
           count+=1
           item=item.next
       return count

    def isEmpty(self):
       return self._first is None

    def __getitem__(self, index):
       if index>=0:
           count = 0
           item = self._first
           if item is None:
               raise IndexError("list index out of range")
           while count<index:
               try:
                   count+=1
                   item = item.next
               except AttributeError:
                   raise IndexError("list index out of range")
           return item.item
       if index<0:
           count = -1
           item = self._last
           if item is None:
               raise IndexError("list index out of range")
           while count>index:
               try:
                   count-=1
                   item = item.prev
               except AttributeError:
                   raise IndexError("list index out of range")
           return item.item

    def __repr__(self):
       contents = ""
       if self._first is not None:
           contents = self._first.getRepr()
       return "[%s]" % (contents)

    def __getslice__(self,i,j):
       #TODO
       pass

    def __iter__(self):
        return Stack.Iter(self._first)

    class Iter():
       def __init__(self, item):
           self._next = item

       def __iter__(self):
           return self

       def next(self):
           ret = self._next
           if ret is None:
               raise StopIteration
           self._next = ret.next
           return ret

def testStack():
    s1=Stack()
    s1.shift(1)
    s1.shift(2)
    s1.shift(3)
    assert(len(s1)==3)
    s2=Stack((3,2,1))
    assert(len(s2)==3)
    assert(s2[2]==s1[-1])
    assert(1==s1.pop())
    s1.append(1)
    assert(3==s2.unshift())
    s2.shift(3)
    assert(s1 is not s2)
    assert(len(s2)==3)

def passThrough(x):
    return x

def parseMapFile(mapFile, delim="\t", keyType=None, valueType=None, keyCol=0, valueCol=1, skipFirst=0):
    """
    Parse a tabular file into a dictionary. The parameters are:

    delim:
        String to split lines on (using str.split()). None => whitespace. Defaults to tab
    keyType:
    valueType:
        Functions to apply to parsed strings (can be int(), float(), or anything else). Defaults to None for no action.
    keyCol:
		column index (starting with 0) of the key values (default: 0)
    valueCol:
		column index (starting with 0) of the values (default: 1)
    skipFirst:
        How many lines to skip at start of file. Defaults to 0.
    """
    if mapFile is None:
        return None

    if keyType is None:
        keyType=passThrough
    if valueType is None:
        valueType=passThrough

    logger.info("parsing map file: %s" % (mapFile))
    translation={}
    badRows=0
    minWidth=max(keyCol,valueCol)
    for line in open(mapFile):
        if skipFirst>0:
            skipFirst-=1
            continue
        cells = line.split(delim)
        if len(cells)>minWidth:
            key=keyType(cells[keyCol].strip())
            value=valueType(cells[valueCol].strip())
            translation[key]=value
        else:
            badRows += 1
    if badRows>0:
        logger.warn("%d rows in map file had too few columns!" % (badRows))

    logger.info("Read %d records from %s" % (len(translation), mapFile))
    return translation

def parseListToMap(listFile,delim=None,keyType=None,col=None):
    """
    Given a file, return a map where each line is a key and the value is True
    """

    if listFile is None:
        return None

    logger.info("parsing list: %s" % listFile)

    if keyType is None:
        keyType=passThrough

    if col is None and delim is None:
        processLine=lambda l: keyType(l.strip())
    else:
        if delim is not None and col is None:
            col=0
        logger.debug("List col: %d" % col)
        processLine=lambda l: keyType(l.split(delim)[col].strip())
        if delim is not None:
            logger.debug("List delim: %r" % delim)

    listHash={}
    with open(listFile) as f:
        for line in f:
            listHash[processLine(line)]=True

    logger.debug("Read %d items from list: %s" % (len(listHash),listFile))

    return listHash

##
# readClusterFile(file)
# return dict from read name to cluster size
##
readRE=re.compile(r'>(\S+)\.\.\.')
def readClusterFile(file):
    cfile = open(file)
    clusters = {}
    clusterread=""
    for line in cfile:
        if len(line)==0:
            continue
        if line[0]=='>':
            clusterread=''
            continue

        # get read (use regex?)
        thisread=""
        m = readRE.search(line)
        if m:
            thisread=m.group(1)
        else:
            sys.stderr.write("Cannot parse line: %s" % line)
            continue

        if clusterread=="":
            clusterread = thisread
            clusters[clusterread]=1
        else:
            clusters[clusterread]+=1

    if logger.getEffectiveLevel() <= logging.DEBUG:
        text=""
        for read, count in clusters.iteritems():
            text+="%s::%s:%d" % (text,read,count)
        logger.debug( "Parsed cluster file:\n %s" % text)

    return clusters

def tupleIteratorToMap(iterator):
    retMap={}
    for (key,value) in iterator:
        retMap[key]=value
    return retMap

def addIOOptions(parser, defaults={}):
    parser.add_option("-o", "--outfile", dest="outfile",
                      metavar="OUTFILE", help="Write count table to OUTFILE. Interpreted as a suffix if multiple input files given. Defaults to STDOUT.")
    parser.add_option("--cwd", default=False, action='store_true', help="If creating multiple output files from multiple inputs, create output files in current directory, not in directory with input files. (By default, a suffix is appended to the full path of the input file)")

#def inputIterator(infileNames, outfileName, infileName=None, cwd=True):
def inputIterator(infileNames, options):
    """
    take list of input files from infileName (which can be None) and infileNames (which is a list of any length)
    if no file names give, read from stdin
    if outfile is not given, pritn to stdout
    if multiple infiles, treat outfile as suffix
    if cwd set to False (for multipleinpus) create output files in same dir as inputs. Otherwise, create files in current dir withinput names plus suffix
    """
    outfileName=options.outfile
    if 'infile' in dir(options) and options.infile is not None:
        infileNames.append(options.infile)
    if len(infileNames)==0:
        inhandle=sys.stdin
        if outfileName==None:
            logger.info("IO: STDIN -> STDOUT")
            yield (inhandle,sys.stdout)
        else:
            outhandle=open(outfileName,'w')
            logger.info("IO: STDIN -> %s" % outfileName)
            yield(inhandle,outhandle)
            outhandle.close()
    elif len(infileNames)==1:
        infileName=infileNames[0]
        inhandle=open(infileName)
        if outfileName==None:
            logger.info("IO: %s -> STDOUT" % infileName)
            yield (inhandle,sys.stdout)
        else:
            outhandle=open(outfileName,'w')
            logger.info("IO: %s -> %s" % (infileName,outfileName))
            yield(inhandle,outhandle)
            outhandle.close()
        inhandle.close()
    else:
        for infileName in infileNames:
            inhandle=open(infileName)
            if outfileName==None:
                logger.info("IO: %s -> STDOUT" % (infileName))
                yield (inhandle,sys.stdout)
            else:
                # use outfileName as suffix
                if options.cwd:
                    # strip path info first
                    (infilePath,infileName)=os.path.split(infileName)
                    infileName="./"+infileName
                outhandle=open("%s%s" % (infileName,outfileName),'w')
                logger.info("IO: %s -> %s%s" % (infileName,infileName,outfileName))
                yield (inhandle,outhandle)
                outhandle.close()
            inhandle.close()

def addUniversalOptions(parser,addQuiet=True):
    parser.add_option("-v", "--verbose",
                      action="count", dest="verbose", default=1,
                      help="Print log messages. Use twice for debugging")
    if addQuiet:
        parser.add_option("-q", '--quiet', dest='verbose',
                          action="store_const", const=0,
                          help="Suppress warnings. Only print fatal messages")
    parser.add_option("-A", "--about",
              action="store_true", dest="about", default=False,
              help="Print description")

DEFAULT_LOGGER_FORMAT=':%(asctime)s::%(levelname)s:%(name)s:%(funcName)s:\n%(message)s'
def setupLogging(options, description, stream=sys.stderr, format=DEFAULT_LOGGER_FORMAT):
    """
    Do some basic setup common to all scripts.

    Given:
        an options object with:
            an integer logLevel value between 0(silent) and 3(debug)
            a boolean called 'about'
        a description string
    Print the description string and exit if about is True
    Set up a logger otherwise. Accepts stream and format key word arguments
    """
    if options.about:
        print description
        exit(0)

    verbose=options.verbose
    if verbose==0:
        loglevel=logging.ERROR
    elif verbose==1:
        loglevel=logging.WARN
    elif verbose==2:
        loglevel=logging.INFO
    elif verbose>=3:
        loglevel=logging.DEBUG
    logging.basicConfig(stream=stream, level=loglevel, format=format)
    logging.info("Log level set to %r(%d)" % (loglevel,verbose))

###############
# Methods for parsing tables to list (e.g. for screening)
##############
def addScreenOptions(parser, defaults={}, accs=False):
    parser.add_option("-l", "--listFile", dest="listFile",
                      metavar="LISTFILE", help="List of names to filter with"),
    parser.add_option("-k", "--keep",
                      action="store_true", dest="keep", default=False,
                      help="Keep listed reads instead of removing")
    parser.add_option("-D", "--listDelim", dest="listDelim", default=None,
                      help="list delimiter. If listColumn set, default is any whitespace, otherwise, the whole line (stripped of whietspace at the ends) is used. '\\t' will split on tab characters.", metavar="DELIM")
    parser.add_option("-C", "--listColumn", dest="listColumn", default=None, type='int',
                      help="Column in listFile to get names from. Defaults to 0 if a delimiter is set.")
    parser.add_option("-G", "--galaxy", default=False, action="store_true", help="Column indices should start with 1")
    if accs:
        parser.add_option("-a", "--accs",
                          action="store_true", dest="accs", default=False,
                          help="parse accession from read name in fasta")

dotRE=re.compile(r'(\.\d+)$')
def getScreenList(options, accs=False):
    """
    Use the listFile, listDelim, and listColumn fields in options to
    parse a text file into a dict of names all mapped to True

    If accs==True, remove version (s/\.\d+$//) from accessions.

    Uses parseListToMap to do the actual parsing
    """
    if accs:
        # hack the keyType to do our regex subs
        logging.debug("Stripping versions")
        translator=lambda r: dotRE.sub('',r)
    else:
        translator=None

    if options.listDelim is not None or options.listColumn is not None:
        if options.listColumn is None:
            options.listColumn=0
        if options.galaxy:
            options.listColumn-=1

    return parseListToMap(options.listFile, col=options.listColumn, delim=options.listDelim, keyType=translator)

def parseAcc(read):
    #logging.info("looking for acc in %s" % (repr(read)))
    m=accessionRE.search(read)
    if m:
        return m.group(1)
    return read

def treeGenerator(node, kidsFirst=False, **kwargs):
    """
    Yields the given node and children recursively.
    Starts with this node unless kidsFirst==True
    """
    if not kidsFirst:
        yield node
    for kid in sorted(node.children,**kwargs):
        if kid is not node:
            for n in treeGenerator(kid, kidsFirst=kidsFirst, **kwargs):
                yield n
    if kidsFirst:
        yield node

returnSelf=lambda x: x
def pairwise(items, sortKey=returnSelf):
    """
    iterate over unique pairs of items
    """
    itemList = sorted(items, key=sortKey)
    for i in range(len(itemList)):
        for j in range(i+1,len(itemList)):
            yield itemList[i],itemList[j]
from numpy import ceil
from numpy import log as nplog, exp as npexp
def asciiHistogram(histogram, log=False, width=60, label='length', maxLabelWidth=10):
    (values,edges)=histogram[:2]
    
    maxValue=max(values)
    
    centers=[int(float(sum(edges[i:i+2]))/2.) for i in range(len(values))]
    largestLabel = max(max([len(str(c)) for c in centers]),len(label))
    if largestLabel<6:
        largestLabel=6
    elif largestLabel>maxLabelWidth:
        largestLabel=maxLabelWidth
    
    plotWidth=width-largestLabel+1
    
    midPoint = npexp((nplog(maxValue)-nplog(.5))/2) if log else maxValue/2
    output="%s|count%s%s|%s%s|\n" % (rightPad(label,largestLabel),
                                     "".join([" " for i in range(plotWidth/2 - len(str(int(midPoint))) - len("count"))]),
                                     str(int(midPoint)),
                                     "".join([" " for i in range(int(ceil(plotWidth/2.)) - 1 - len(str(int(maxValue))))]),
                                     str(int(maxValue)),
                                     )
    #output+="%s|%s\n" % ("".join(["_" for i in range(largestLabel)]),
    #                     "".join(["_" for i in range(plotWidth)]),
    #                     )
    for i, v in enumerate(values):
        output+="%s|%s\n" % (rightPad(str(centers[i]),largestLabel),getBarString(v, maxValue, plotWidth, log))
    return output

logChars=['-','~','=','#']
def getBarString(value, maxValue, maxWidth, log):
    """
    return string of various signs (-,~,=,#) based on value and scale
    """
    if log:
        value=nplog(value)-nplog(.5) if value>0 else 0
        maxValue=nplog(maxValue)-nplog(.5)
    width=maxWidth*(value/float(maxValue))
    if width<1:
        return ''
    char=logChars[0]
    s=char
    while len(s)<width:
        if log:
            #print "s: %s, mw: %s, w: %s" % (s, maxWidth, width)
            char=logChars[int(ceil(len(logChars)*len(s)/float(maxWidth))-1)]
        s+=char
    return s

def rightPad(name, width):
    if width<6:
        width=6
        logger.warn("Can't force names to be fewer than 6 characters")
    if len(name)>width:
        # remove middle and insert elipsis to get to -width- characters
        return name[:width-4]+'***'+name[-1:]
    while len(name)<width:
        # pad with trailing space to get to 13 characters
        name+=' '
    return name

