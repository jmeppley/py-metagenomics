import re, logging, sys, os, numpy, argparse, random
from edl import __version__ 
from edl.expressions import accessionRE
logger=logging.getLogger(__name__)
VERSION='py-metagenomics-{}'.format(__version__)

class LineCounter():
    """
    File-like wrapper for a file-like object that counts the number of lines returned.
    """
    def __init__(self, stream):
        self.rawStream=stream
        self.lines=0

    def __iter__(self):
        return self

    def __next__(self):
        line = next(self.rawStream)
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

def parseMapFile(mapFile, delim="\t", keyType=None, valueType=None, keyCol=0, valueCol=1, valueDelim=None, skipFirst=0):
    """
    Parse a tabular file into a dictionary. The parameters are:

    delim:
        String to split lines on (using str.split()). None => whitespace.
        Defaults to tab
    keyType:
    valueType:
        Functions to apply to parsed strings (can be int(), float(), or
        anything else). Defaults to None for no action.
    keyCol:
        column index (starting with 0) of the key values (default: 0)
    valueCol:
        column index (starting with 0) of the values (default: 1)
    valueDelim:
        String to split value cells into multiple values (default: None => don't split)
    skipFirst:
        How many lines to skip at start of file. Defaults to 0.
    """
    if mapFile is None:
        return None

    if keyType is None:
        keyType=passThrough
    if valueType is None:
        valueType=passThrough
    if valueDelim is not None:
        baseValueType=valueType
        valueType = lambda value_cell: [baseValueType(v) for v in value_cell.split(valueDelim)]

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

def parse_list_to_set(listFile,delim=None,keyType=None,col=None):
    """
    Given a file, return a set where each line or cell (if delim and col are provided) is an entry in the set
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

    items=set()
    with open(listFile) as f:
        for line in f:
            items.add(processLine(line))

    logger.debug("Read %d items from list: %s" % (len(items),listFile))

    return items

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
        for read in clusters:
            count=clusters[read]
            text+="%s::%s:%d" % (text,read,count)
        logger.debug( "Parsed cluster file:\n %s" % text)

    return clusters

def tupleIteratorToMap(iterator):
    retMap={}
    for (key,value) in iterator:
        retMap[key]=value
    return retMap

def add_IO_arguments(parser, defaults={}):
    parser.add_argument("-o", "--outfile", dest="output_file",
                        metavar="OUTFILE", help="Write output to OUTFILE. Interpreted as a suffix if multiple input files given. Defaults to STDOUT.")
    parser.add_argument("--cwd", default=False, action='store_true', help="If creating multiple output files from multiple inputs, create output files in current directory, not in directory with input files. (By default, a suffix is appended to the full path of the input file)")
    parser.add_argument("input_files", nargs="*", 
                        type=argparse.FileType('r'),
                        default=[sys.stdin,],
                        metavar="INFILE",
                        help="List of input files. Omit to read from STDIN")

def inputIterator(arguments):
    """
    processes arguments from add_IO_arguments:
     input_files: list of file(...,'r') objects
     output_file: file name or prefix
    if no input file names give, read from stdin
    if outfile is not given, pritn to stdout
    if multiple infiles, treat outfile as suffix
    if cwd set to False (for multipleinpus) create output files in same dir as inputs. Otherwise, create files in current dir withinput names plus suffix

    yields pairs of input and output streams as 2-element tuples
    """
    outfile_name=arguments.output_file
    infiles=arguments.input_files
    if len(infiles)==1:
        inhandle=infiles[0]
        infile_name=inhandle.name
        if outfile_name==None:
            logger.info("IO: %s -> STDOUT" % infile_name)
            yield (inhandle,sys.stdout)
        else:
            outhandle=open(outfile_name,'w')
            logger.info("IO: %s -> %s" % (infile_name,outfile_name))
            yield(inhandle,outhandle)
            outhandle.close()
        inhandle.close()
    else:
        for inhandle in infiles:
            infile_name = inhandle.name
            if outfile_name==None:
                logger.info("IO: %s -> STDOUT" % (infile_name))
                yield (inhandle,sys.stdout)
            else:
                # use outfile_name as suffix
                if options.cwd:
                    # strip path info first
                    (infile_path,infile_name)=os.path.split(infile_name)
                    infile_name="./"+infile_name
                outhandle=open("%s%s" % (infile_name,outfile_name),'w')
                logger.info("IO: %s -> %s%s" % (infile_name,infile_name,outfile_name))
                yield (inhandle,outhandle)
                outhandle.close()
            inhandle.close()

def add_universal_arguments(parser,addQuiet=True):
    parser.add_argument('-V','--version', 
            action="version", version=VERSION)
    parser.add_argument("-v", "--verbose",
                      action="count", dest="verbose", default=1,
                      help="Print log messages. Use twice for debugging")
    if addQuiet:
        parser.add_argument("-q", '--quiet', dest='verbose',
                          action="store_const", const=0,
                          help="Suppress warnings. Only print fatal messages")

DEFAULT_LOGGER_FORMAT=\
        ':%(asctime)s::%(levelname)s:%(name)s:%(funcName)s:\n%(message)s'
def setup_logging(parsed_args, stream=sys.stderr, \
        format=DEFAULT_LOGGER_FORMAT):
    """
    Do some basic setup common to all scripts.

    Given:
        an parsed_arguments object with:
            an integer "verbose" value between 0(silent) and 3(debug)
    Set up a logger. Accepts stream and format key word arguments
    """
    verbose=parsed_args.verbose
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
def add_screen_arguments(parser, defaults={}, accs=False):
    parser.add_argument("-l", "--listFile", dest="listFile",
                      metavar="LISTFILE", help="List of names to filter with"),
    parser.add_argument("-k", "--keep",
                      action="store_true", dest="keep", default=False,
                      help="Keep listed reads instead of removing")
    parser.add_argument("-D", "--listDelim", dest="listDelim", default=None,
                      help="list delimiter. If listColumn set, default is any whitespace, otherwise, the whole line (stripped of whietspace at the ends) is used. '\\t' will split on tab characters.", metavar="DELIM")
    parser.add_argument("-C", "--listColumn", dest="listColumn", default=None, type='int',
                      help="Column in listFile to get names from. Defaults to 0 if a delimiter is set.")
    parser.add_argument("-G", "--galaxy", default=False, action="store_true", help="Column indices should start with 1")
    if accs:
        parser.add_argument("-a", "--accs",
                          action="store_true", dest="accs", default=False,
                          help="parse accession from read name in fasta")

dotRE=re.compile(r'(\.\d+)$')
def get_screen_list(arguments, accs=False):
    """
    Use the listFile, listDelim, and listColumn fields in arguments to
    parse a text file into a dict of names all mapped to True

    If accs==True, remove version (s/\.\d+$//) from accessions.

    Uses parse_list_to_set() to do the actual parsing
    """
    if accs:
        # hack the keyType to do our regex subs
        logging.debug("Stripping versions")
        translator=lambda r: dotRE.sub('',r)
    else:
        translator=None

    if arguments.listDelim is not None or arguments.listColumn is not None:
        if arguments.listColumn is None:
            arguments.listColumn=0
        if arguments.galaxy:
            arguments.listColumn-=1

    return parse_list_to_set(arguments.listFile, col=arguments.listColumn, delim=arguments.listDelim, keyType=translator)

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
                                     "".join([" " for i in range(int(plotWidth/2) - len(str(int(midPoint))) - len("count"))]),
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

def indexed_sample_generator(records, N, P=None):
    """
    Generate a random sample of N records from interator "records".

    If the total numer of records in terator (P) is not specified, 
    then reservoir sampling (storing entire sample in RAM) will be used!
    """
    if P>0:
        # use random.sample to get list of indices ahead of time
        indexes_to_return = set(random.sample(range(P),N))

        # return just records with tose indices
        for i, record in enumerate(records):
            if i in indexes_to_return:
                yield record
    else:
        # We don't know how many records to expect, so use reservoir
        for record in reservoir_sample(records, N=N, return_count=False):
            yield record

def reservoir_sample(iterator, N=100, return_count=False):
    """
    Randomly sample N items from an iterator of unknown size. Returns a list of the
    selected items. 
    
    If the iterator returns fewer than N items, all will be in the sample, but the sample will not have N elements.
    
    IF return_count set to true, return value is tuple: (sample, totalItemsFromIterator)
    """
    
    sample = []

    numerator = float(N)
    for i,item in enumerate(iterator):
        n=i+1
        val = numpy.random.rand()
        if n <= N:
            # Fill up sample with first N elements
            sample.append(item)
        elif val < numerator/n:
            # Replace random item in sample 
            # with next element with probabliity N/i
            sample[numpy.random.random_integers(0,N-1)]=item
            
    if return_count:
        return (sample, n)
    else:
        return sample

class ReservoirSamplingList(list):
    """
    Behaves mostly as a generic python container with one primary difference: it has a maximum size and
    only keeps a random subsample of added items once that size is reached.
    
    Other properties:
    
     * Insert() at an index is not supported. 
     * Deleting items is not supported.
     * Only append() is allowed to change the contents
     * The total_added veriable tracks how many iterms were added
     * The order of items is not guaranteed to match the order in which they were added. 
     * The order is nevertheless stable and elements can be retrieved by index or slice
    """
    
    def __init__(self, sample_size=100, iterable=[], preserve_order=False):
        """
        Create a new reservoirSampleingList with the given sample size (defaults to 100).
        
        If an iterable object is given as a second object,
         the reservoirSamplingList will be initialized with the items returned by iterating ofer the iterable.
        """
        self.N=int(sample_size)
        self.total_added=0
        if preserve_order:
            self.append=self._append_preserve_order
        else:
            self.append=self._append
        for item in iterable:
            self.append(item)
        
    def __delitem__(self, key):
        raise Exception("This container can only changed by appending items. Deleting is not supported")

    def _append(self, item):
        """
        Add an item to the reservoirSamplingList. 
        If the current size of the list is less than the sample_size, this item will simply be added.
        As the total added surpasses the sample_size, each new item will have a decreasing probability
         of replacing an existing item.
        This will simulate a subsample of the total set of added items.
        """
        self.total_added+=1
        if self.total_added <= self.N:
            # Fill up container with first N elements
            super(ReservoirSamplingList,self).append(item)
            return

        val = numpy.random.rand()
        if val < self.N/float(self.total_added):
            # Replace random item in sample
            # with next element with probabliity N/i
            self[numpy.random.random_integers(0,self.N-1)]=item

    def _append_preserve_order(self, item):
        """
        Add an item to the reservoirSamplingList. 
        If the current size of the list is less than the sample_size, this item will simply be added.
        As the total added surpasses the sample_size, each new item will have a decreasing probability
         of replacing an existing item.
        This will simulate a subsample of the total set of added items.
        """
        self.total_added+=1
        if self.total_added <= self.N:
            # Fill up container with first N elements
            super(ReservoirSamplingList,self).append(item)
            return

        val = numpy.random.rand()
        # Replace random item in sample
        # with next element with probabliity N/i
        if val < self.N/float(self.total_added):
            # pick random element to replace
            drop_index=numpy.random.random_integers(0,self.N-1)
            # shift following items over
            for i in xrange(drop_index,len(self)-1):
                self[i]=self[i+1]
            # add new one at the end
            self[i+1]=item

    def extend(self, iterable):
        for item in iterable:
            self.add(item)


