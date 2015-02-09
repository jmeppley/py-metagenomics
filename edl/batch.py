import tempfile, re, os, logging, sys
import numpy as np
from edl.util import openInputFile

def checkTmpDir(tmpDir,jobName):
    """
    Make sure tmp dir is empty.
    Create it if necessary.
    If name is None, create dir name based on job name
    """
    if tmpDir is None or tmpDir=="":
        tmpDir=tempfile.mkdtemp(suffix=jobName,dir=".")
    else:
        if os.path.exists(tmpDir):
            raise Exception("Temporary directory already exists! (%s)" % (tmpDir))
        else:
            os.makedirs(tmpDir)

    logging.debug("Created temporary directory: %s" % tmpDir)
    return tmpDir

def fragmentInput(infile, options, tmpdir, fragmentBase, suffix='.in'):
    """
    Wraps the following methods into one:
        getFileType(options, infile)
        getSizePerChunk(infile, options.splits, fileType, splitOnSize=options.splitOnSize)
        fragmentInputBySize(infile, tmpdir, options.chunk, fileType, fragmentBase,   splitOnSize=options.splitOnSize, suffix=)
    """
    fileType=getFileType(options, infile)
    if options.chunk is None:
        if options.splits is None:
            sys.exit("Please tell me how many chunks (-N) or how big (-C)!")
        options.chunk=getSizePerChunk(infile, options.splits, fileType, splitOnSize=options.splitOnSize)
    return fragmentInputBySize(infile, tmpdir, options.chunk, fileType, fragmentBase, splitOnSize=options.splitOnSize, suffix=suffix)

def getFileType(options, infile):
    """
    figure out how to split up the input file and return a FileType object with the necessary info.

    The options paramter is assumed to be an object with the following values:
        .infileType (None or a string indicating the input file type)
        .pattern (None or a regex pattern to find a record boundary)
        .numLines (None or an integer number of lines per record)
    if all of the above are None, the extension of the name in infile is used to guess
    """
    if options.pattern is not None:
        if options.numLines is not None:
            logging.error("Cannot use BOTH pattern and number of lines to split records!")
            raise Exception("Conflicting options: pattern and numLines")


    # if there is a file type, use that (And override if other options present)
    if options.infileType is not None:
        fileType = fileTypeMap[options.infileType]
        if fileType.sepRE is not None:
            logging.info('Using user chosen "%s" pattern to split files: /%s/' % (options.infileType, fileType.sepRE.pattern))
        else:
            logging.info('Using user chosen "%s" number of lines (%d) to split files' % (options.infileType, fileType.numLines))
    else:
        fileType=None

    if options.pattern is not None:
        sepRE = re.compile(options.pattern)
        logging.info('Using user supplied pattern to split files: /%s/' % (options.pattern))
        if fileType is None:
            fileType=FragmentableFileType(name=options.pattern,sepRE=sepRE)
        else:
            fileType.sepRE=sepRE
            fileType.numLines=None
    elif options.numLines is not None:
        logging.info('Using user supplied number of lines to split files: /%s/' % (options.numLines))
        if fileType is None:
            fileType=FragmentableFileType(str(options.numLines), numLines=options.numLines)
        else:
            logging.info("Overriding record line count (%d)!" % (options.numLines))
            fileType.sepRE=None
            fileType.numLines=options.numLines
    elif fileType is None:
        # try to guess from the extension
        if infile is None:
            sys.exit("We cannot infer the record separator from the filename when using standar input. Please specify the input file type with -T or a pattern to match the first line with -p")
        fileType = getTypeFromFileName(infile)
        if fileType.sepRE is not None:
            logging.info('File looks like %s, using pattern to split files: /%s/' % (fileType.name, fileType.sepRE.pattern))
        else:
            logging.info('File looks like %s, using number of lines (%d) to split files' % (fileType.name, fileType.numLines))
    return fileType

def getTypeFromFileName(filename):
    """
    find the file's extension in the fileExtensionMap
    return the separator RE for that file type

    If file not recognized, inform user how to set the pattern manually
    """
    ext = os.path.splitext(filename)[1]
    fileType = fileExtensionMap.get(ext,None)
    if fileType is not None:
        return fileType
    else:
        logging.warn("""Cannot guess type of %s!

The input file does not have a recognized extension:
%s

You must manually set the file type with '-T TYPE' where TYPE is in:
%s
OR set the record separator pattern with '-p PATTERN' where PATTERN is a regular expression. E.G. '^>" for fasta, or '^LOCUS' for gbk.
""" % (filename, fileExtensionMap.keys(), fileExtensionMap.keys()))
        sys.exit(65)

def getSizePerChunk(infile, splits, fileType, splitOnSize=False):
    """
    Get total size of all records and return target size for each chunk to end up with number of chunks specified by 'splits'
    """
    if infile is None:
        raise Exception("We cannot determine chunk size from STDIN!")

    if splitOnSize:
        # get a custom function that returns the size of this type of record
        recordSizer=fileType.sizer
    else:
        # just return 1 for each record
        recordSizer=recordCounter

    # loop through records
    inhandle = openInputFile(infile)
    totalSize = 0
    for record in fileType.recordStreamer(inhandle):
        totalSize+=recordSizer(record)
    inhandle.close()

    return calculateChunkSize(totalSize,splits)

def calculateChunkSize(size,splits):
    """
    how big should the fragments be?
    """
    chunk = int(np.ceil(size/float(splits)))
    logging.info("Setting chunk to: %d=ceil(%d/%d)" % (chunk,size,splits))
    return chunk

# Simply count records
recordCounter=lambda x: 1

def defaultRecordSizer(recordLines):
    """
    return the total number of characters in the record text
    """
    size=0
    for line in recordLines:
        size+=len(line)
    return size

def fastaRecordSizer(recordLines):
    """
    Returns the number of charcters in every line excluding:
        the first (header) line
        whitespace at the start and end of lines
    """
    size=0
    for i in xrange(1,len(recordLines)):
        size+=len(recordLines[i].strip())
    return size

def fastqRecordSizer(recordLines):
    """
    Returns the number of charaters in the lines between the sequence header (@) and the quality header (@) excluding whitespace at the start and end of lines
    """
    size=0
    for i in xrange(1,len(recordLines)):
        line=recordLines[i]
        if len(line)>0 and line[0]=='+':
            return size
        size+=len(line.strip())

    logging.warn("Did not find quality line in record: %s" % recordLines)
    return size

def gbRecordSizer(recordLines):
    """
    uses biopython to parse record
    """
    from Bio import SeqIO
    record = SeqIO.read(recordLines,'genbank')
    return len(record)

def fragmentInputBySize(infile, tmpdir, chunk, fileType, fragmentBase, splitOnSize=True, suffix='.in'):
    """
    Break up input into files of size chunk in tmpdir. Return number of fragments.
    """
    logging.debug("Fragmenting input: %r" % ({'infile':infile,'tmpDir':tmpdir,'chunk':chunk,'base':fragmentBase}))
    inhandle = openInputFile(infile)
    num = fragmentInputStreamBySize(inhandle, tmpdir, chunk, fileType, fragmentBase, splitOnSize=splitOnSize, suffix=suffix)
    if infile is not None:
        inhandle.close()
    return num

def fragmentInputStreamBySize(inhandle, tmpdir, chunk, fileType, fragmentBase, splitOnSize=True, suffix='.in'):
    if splitOnSize:
        # get a custom function that returns the size of this type of record
        recordSizer=fileType.sizer
    else:
        # just return 1 for each record
        recordSizer=lambda x: 1

    count=0
    num=1
    tmpFileName=getFragmentPath(tmpdir,fragmentBase,num,suffix)
    #logging.debug('Writing fragment (%d,%d,%d): %s' % (chunk,count,num,tmpFileName))
    tmpFile = open(tmpFileName, 'w')
    for record in fileType.recordStreamer(inhandle):
        recordSize=recordSizer(record)
        count+=recordSize

        if count>chunk:
            # close previous chunk and open new one
            tmpFile.close
            num+=1
            tmpFileName=getFragmentPath(tmpdir,fragmentBase,num,suffix)
            #logging.debug('Writing fragment (%d,%d,%d): %s' % (chunk,count,num,tmpFileName))
            tmpFile = open(tmpFileName, 'w')
            count=recordSize

        # write record
        tmpFile.writelines(record)

    tmpFile.close()

    return num

def getFragmentPath(directory, base, index, suffix='.in'):
    return "%s%s%s" % (directory, os.sep, getFragmentName(base,index,suffix=suffix))

def getFragmentName(base, index, suffix='.in'):
    return "%s%05d%s" % (base, index, suffix)

def getFragmentPrefix(base,index):
    return getFragmentName(base,index,suffix='.pre')

def formatCommand(command):
    """
    given a list of command elements, print string approximating what you'd type at a shell prompt
    """
    cmdstr=""
    logging.debug(repr(command))
    for arg in command:
        if " " in arg:
            cmdstr=cmdstr+" \""+arg+"\""
        else:
            cmdstr=cmdstr+" "+arg
    return cmdstr

def regexRecordGenerator(fileType, stream):
    """
    Using the sepRE setting in fileType, break the input stream into records
    """
    lastRecord=[]
    for line in stream:
        if fileType.sepRE.match(line):
            if len(lastRecord)>0:
                yield lastRecord
                del lastRecord[:]
        lastRecord.append(line)

    if len(lastRecord)>0:
        yield lastRecord

def linedRecordGenerator(fileType, stream):
    """
    Using the numLines setting in fileType, break the input stream into records
    """
    lastRecord=[]
    for index,line in enumerate(stream):
        if index%fileType.numLines==0:
            if len(lastRecord)>0:
                yield lastRecord
                del lastRecord[:]
        lastRecord.append(line)

    if len(lastRecord)>0:
        yield lastRecord

def addFragmentingOptions(parser,defaults={"splits":400}):
    parser.add_option("-L", "--recordLines", metavar="NUMLINES", dest='numLines', default=None, type="int",
                       help="Number of lines per record")
    parser.add_option("-P", "--pattern", metavar="PATTERN", dest='pattern', default=None,
                       help="Regular expression to split records")
    parser.add_option("-T","--infileType", dest='infileType', default=None,
                      choices=fileTypeMap.keys(),
                      help='Type of input file. Otherwise, choose by extension. Known types are: %choices')
    parser.add_option("-C", "--chunkSize", type="int", dest='chunk',  metavar="FRAG_SIZE",
                      help="The number of records per fragment. Overrides NUM_FRAGS")
    default=defaults.get("splits",None)
    parser.add_option("-N", "--numChunks", dest='splits', type='int', metavar="NUM_FRAGS", default=default,
                      help="The number of fragments to create (defaults to %default)")
    parser.add_option("-s", "--splitOnSize",  default=False, action='store_true',
                      help="create chunks based on record size, not number of records. For known sequence types (fasta, fastq, gb), sequence length is used, otherwize the full size of the record text is counted")

#################
# Classes
#################
class FragmentableFileType:
    def __init__(self, name, sepRE=None, numLines=None, sizer=None):
        self.name=name
        self.sepRE=sepRE
        if sepRE is not None:
            self._recordStreamer=regexRecordGenerator
        self.numLines=numLines
        if numLines is not None:
            self._recordStreamer=linedRecordGenerator
        if sizer==None:
            self.sizer=defaultRecordSizer
        else:
            self.sizer=sizer

    def recordStreamer(self, stream):
        return self._recordStreamer(self, stream)

##################
# Constants
##################
FASTA=FragmentableFileType('fasta',sizer=fastaRecordSizer,sepRE=re.compile(r'^>(\S+)'))
FASTQ=FragmentableFileType('fastq',sizer=fastqRecordSizer,numLines=4)
GENBANK=FragmentableFileType('gb',sizer=gbRecordSizer,sepRE=re.compile(r'^LOCUS'))
TABLE=FragmentableFileType('table',sepRE=re.compile(r'^'))
fileTypeMap={FASTA.name:FASTA,
             FASTQ.name:FASTQ,
             GENBANK.name:GENBANK,
             TABLE.name:TABLE,
            }
fileExtensionMap={'.fa':FASTA,
                  '.fna':FASTA,
                  '.faa':FASTA,
                  '.ffn':FASTA,
                  '.fasta':FASTA,
                  '.fastq':FASTQ,
                  '.gb':GENBANK,
                  '.gbk':GENBANK,
                  '.gbff':GENBANK,
                  '.gpff':GENBANK,
                  '.tab':TABLE,
                  '.tsv':TABLE,
                  '.csv':TABLE,
                  '.m8':TABLE,
                 }
