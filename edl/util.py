import argparse
import glob
import io
import logging
import os
import re
import sys
from edl import __version__
from edl.expressions import accessionRE
from numpy import ceil, random
from numpy import exp as npexp, log as nplog
logger = logging.getLogger(__name__)
VERSION = 'py-metagenomics-{}'.format(__version__)


def find_matching_files(search_string, wildcard_constraints=None):
    """
    Looks for files matching a pattern and captures the variable bits and
    yeilds
    (file_path, wildcard_values) tuples.

    For example:

    dict(find_matching_files("/path/{subdir}/{key}.static.{extension}"))

    returns

    {'/path/dir1/sample1.static.txt': {'subdir':'dir1',
                                       'key':'sample1',
                                       'extension':'txt'},
    {'/path/dir2/sample1.static.txt': {'subdir':'dir2',
                                       ...}
    ...
    }

    Like snakemake, you can specify constraints to limit what matches each
    wildcard.

    Unlike snakemake, a wildcard cannot cross the path separator. Each matched
    item must be a folder or file name or a fragment of either.
    """
    wildcard_re = re.compile(r'\{(\w+)\}')

    # replace wildcards with *'s for glob search
    glob_string = wildcard_re.sub('*', search_string)

    # get wildcard names
    wildcard_keys = wildcard_re.findall(search_string)

    # replace wildcards with (.+) or their constraints for rexp
    if wildcard_constraints is None:
        wildcard_constraints = {}

    def wc_repl(match):
        return "(" + wildcard_constraints.get(match.group(1), r'.+') + ")"
    search_string = re.sub(r'\.', r'\.', search_string)
    re_pattern = wildcard_re.sub(wc_repl, search_string)
    glob_re = re.compile(re_pattern)

    # use glob to filter files
    for full_path in glob.glob(glob_string):
        match = glob_re.match(full_path)
        if match:
            # It might not match if there are constraints
            wildcards = {}
            for wildcard, value in zip(wildcard_keys, match.groups()):
                if wildcard in wildcards:
                    # if a wc name is repeated, both values must match
                    if value != wildcards[wildcard]:
                        break
                else:
                    wildcards[wildcard] = value
            else:
                # only return if we didn't abort
                yield full_path, wildcards


class InputFile():
    """
    Simple line iterator from file or file handle that will count lines

    USes open_input_file() to handle:
        file names (incl gzipped or bzipped)
        already open handles
        None (implying stdin)
    """

    def __init__(self, input_data, *args):
        # use open_input_file to handle different input types
        self.raw_data = open_input_file(input_data, *args)
        # setup counting
        self.lines = 0
        self.enum = enumerate(self.raw_data, start=1)

        # save file name and whether we need to close anything
        self.close_me = False
        try:
            # if the file was already open, get name from property
            self.name = input_data.name
        except AttributeError:
            # we had to open it
            self.name = input_data
            if input_data is not None:
                # None would have given us stdin which we don't want to close
                self.close_me = True

    # define iter/next for line iteration
    def __iter__(self):
        return self

    def __next__(self):
        self.lines, line = next(self.enum)
        return line

    def next(self):
        return self.__next__()

    # define enter/exit to behave as context manager
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def close(self):
        if self.close_me:
            self.raw_data.close()

    def readline(self, *args):
        try:
            return next(iter(self))
        except StopIteration:
            return ""


def iterable_to_stream(iterable,
                       str_to_bytes=True,
                       buffer_size=io.DEFAULT_BUFFER_SIZE):
    """
    Lets you use an iterable (e.g. a generator)
    that yields bytestrings as a read-only
    input stream.

    The stream implements Python 3's newer I/O API
    (available in Python 2's io module).
    For efficiency, the stream is buffered.
    """

    if str_to_bytes:
        # encode strings as bytes
        iterable = (s.encode('utf-8') for s in iterable)

    class IterStream(io.RawIOBase):
        def __init__(self):
            self.leftover = None

        def readable(self):
            return True

        def readinto(self, b):
            try:
                length = len(b)  # We're supposed to return at most this much
                chunk = self.leftover or next(iterable)
                output, self.leftover = chunk[:length], chunk[length:]
                b[:len(output)] = output
                return len(output)
            except StopIteration:
                return 0    # indicate EOF
    return io.BufferedReader(IterStream(), buffer_size=buffer_size)


def checkNoneOption(value):
    """
    Check if a given option has been set to the string 'None' and return
    None if so. Otherwise return the value unchanged
    """
    if value is None:
        return value
    if value == 'None':
        return None
    if isinstance(value, list):
        if len(value) == 1 and value[0] == 'None':
            return None
    return value


def countBasesInFasta(fastaFile):
    """
    Given a fasta file, return a dict where the number of records and
    the total number of bases are given by 'records' and 'bases' respectively.
    """
    recordRE = re.compile(r'^>')
    whiteSpaceRE = re.compile(r'\s+')
    total_bases = 0
    total_seqs = 0
    with open(fastaFile) as f:
        for line in f:
            if recordRE.match(line):
                total_seqs += 1
                continue
            total_bases += len(whiteSpaceRE.sub('', line))

    return {'records': total_seqs, 'bases': total_bases}


urlRE = re.compile(r'[a-z]+\:\/\/')


def open_input_file(infile, *args):
    """
    return input stream. Allow for text, gzipped text, or standard input
    if None given.
    """
    if infile is None:
        logging.info("Reading input from STDIN")
        return sys.stdin

    if isinstance(infile, str):
        if urlRE.match(infile):
            import urllib2
            return urllib2.urlopen(infile)
        if len(infile) > 3 and infile[-3:] == '.gz':
            import gzip
            return gzip.GzipFile(infile, 'rb')
        elif len(infile) > 4 and infile[-4:] == '.bz2':
            import bz2
            return bz2.BZ2File(infile, 'rb')
        else:
            return open(infile, 'rt')
    else:
        return infile


def parseExp(string):
    try:
        return float(string)
    except ValueError:
        return float('1' + string)


def dict_lookup_default_to_query(dictionary):
    return lambda x: dictionary.get(x, x)


def passThrough(x):
    return x


def get_value_type_function(valueType, valueDelim):
    logger.debug("get_value_type_function(%r, '%s')", valueType, valueDelim)
    if valueType is None:
        if valueDelim is not None:
            valueType = str
        else:
            return passThrough
    if valueDelim is not None:
        logger.debug("mapped value delim is set to '%s'", valueDelim)
        baseValueType = valueType
        return lambda value_cell: [
            baseValueType(v) for v in value_cell.split(valueDelim)]
    return valueType


def get_process_line_function(col, delim, keyType):
    if col is None and delim is None:
        return lambda l: keyType(l.strip())
    else:
        if delim is not None and col is None:
            col = 0
        logger.debug("List col: %d" % col)
        if delim is not None:
            logger.debug("List delim: %r" % delim)
        return lambda l: keyType(l.split(delim)[col].strip())


def parseMapFile(
        mapFile,
        delim="\t",
        keyType=None,
        valueType=None,
        keyCol=0,
        valueCol=1,
        valueDelim=None,
        skipFirst=0):
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
        String to split value cells into multiple values
        (default: None => don't split)
    skipFirst:
        How many lines to skip at start of file. Defaults to 0.
    """
    if mapFile is None:
        return None

    if keyType is None:
        keyType = passThrough
    valueType = get_value_type_function(valueType, valueDelim)

    logger.info("parsing map file: %s" % (mapFile))
    translation = {}
    badRows = 0
    minWidth = max(keyCol, valueCol)
    with InputFile(mapFile) as lines:
        for line in lines:
            if skipFirst > 0:
                skipFirst -= 1
                continue
            cells = line.split(delim)
            if len(cells) > minWidth:
                key = keyType(cells[keyCol].strip())
                value = valueType(cells[valueCol].strip())
                translation[key] = value
                logger.debug("mapped %s to %s", key, value)
            else:
                badRows += 1
        num_lines = lines.lines
    if badRows > 0:
        logger.warn("%d of %d rows in map file had too few columns!" %
                    (badRows,
                     num_lines))

    logger.info("Read %d records from %d lines of %s" % (len(translation),
                                                         num_lines,
                                                         mapFile))
    return translation


def parse_list_to_set(listFile, delim=None, keyType=None, col=None):
    """
    Given a file, return a set where each line or cell (if delim
    and col are provided) is an entry in the set
    """

    if listFile is None:
        return None

    logger.info("parsing list: %s" % listFile)

    if keyType is None:
        keyType = passThrough

    processLine = get_process_line_function(col, delim, keyType)

    items = set()
    with open(listFile) as f:
        for line in f:
            items.add(processLine(line))

    logger.debug("Read %d items from list: %s" % (len(items), listFile))

    return items


##
# readClusterFile(file)
# return dict from read name to cluster size
##
readRE = re.compile(r'>(\S+)\.\.\.')


def readClusterFile(file):
    cfile = open(file)
    clusters = {}
    clusterread = ""
    for line in cfile:
        if len(line) == 0:
            continue
        if line[0] == '>':
            clusterread = ''
            continue

        # get read (use regex?)
        thisread = ""
        m = readRE.search(line)
        if m:
            thisread = m.group(1)
        else:
            sys.stderr.write("Cannot parse line: %s" % line)
            continue

        if clusterread == "":
            clusterread = thisread
            clusters[clusterread] = 1
        else:
            clusters[clusterread] += 1

    if logger.getEffectiveLevel() <= logging.DEBUG:
        text = ""
        for read in clusters:
            count = clusters[read]
            text += "%s::%s:%d" % (text, read, count)
        logger.debug("Parsed cluster file:\n %s" % text)

    return clusters


def tupleIteratorToMap(iterator):
    retMap = {}
    for (key, value) in iterator:
        retMap[key] = value
    return retMap


def add_IO_arguments(parser, defaults={}):
    parser.add_argument(
        "-o",
        "--outfile",
        dest="output_file",
        metavar="OUTFILE",
        help="Write output to OUTFILE. Interpreted as a suffix if "
             "multiple input files given. Defaults to STDOUT.")
    parser.add_argument(
        "--cwd",
        default=False,
        action='store_true',
        help="If creating multiple output files from multiple inputs, "
             "create output files in current directory, not in directory "
             "with input files. (By default, a suffix is appended to the "
             "full path of the input file)")
    parser.add_argument("input_files", nargs="*",
                        default=[sys.stdin, ],
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
    if cwd set to False (for multipleinpus) create output files in
        same dir as inputs. Otherwise, create files in current dir
        with input names plus suffix

    yields pairs of input and output streams as 2-element tuples
    """
    outfile_name = arguments.output_file
    infiles = [InputFile(f) for f in arguments.input_files]
    if len(infiles) == 1:
        inhandle = infiles[0]
        infile_name = inhandle.name
        if outfile_name is None:
            logger.info("IO: %s -> STDOUT" % infile_name)
            yield (inhandle, sys.stdout)
        else:
            outhandle = open(outfile_name, 'w')
            logger.info("IO: %s -> %s" % (infile_name, outfile_name))
            yield(inhandle, outhandle)
            outhandle.close()
        inhandle.close()
    else:
        for inhandle in infiles:
            infile_name = inhandle.name
            if outfile_name is None:
                logger.info("IO: %s -> STDOUT" % (infile_name))
                yield (inhandle, sys.stdout)
            else:
                # use outfile_name as suffix
                if arguments.cwd:
                    # strip path info first
                    (infile_path, infile_name) = os.path.split(infile_name)
                    infile_name = "./" + infile_name
                outhandle = open("%s%s" % (infile_name, outfile_name), 'w')
                logger.info("IO: %s -> %s%s" %
                            (infile_name, infile_name, outfile_name))
                yield (inhandle, outhandle)
                outhandle.close()
            inhandle.close()


def add_universal_arguments(parser, addQuiet=True):
    parser.add_argument('-V', '--version',
                        action="version", version=VERSION)
    parser.add_argument("-v", "--verbose",
                        action="count", dest="verbose", default=1,
                        help="Print log messages. Use twice for debugging")
    if addQuiet:
        parser.add_argument(
            "-q",
            '--quiet',
            dest='verbose',
            action="store_const",
            const=0,
            help="Suppress warnings. Only print fatal messages")


DEFAULT_LOGGER_FORMAT =\
    ':%(asctime)s::%(levelname)s:%(name)s:%(funcName)s:\n%(message)s'


def setup_logging(parsed_args, stream=sys.stderr,
                  format=DEFAULT_LOGGER_FORMAT):
    """
    Do some basic setup common to all scripts.

    Given:
        an parsed_arguments object with:
            an integer "verbose" value between 0(silent) and 3(debug)
    Set up a logger. Accepts stream and format key word arguments
    """
    verbose = parsed_args.verbose
    if verbose == 0:
        loglevel = logging.ERROR
    elif verbose == 1:
        loglevel = logging.WARN
    elif verbose == 2:
        loglevel = logging.INFO
    elif verbose >= 3:
        loglevel = logging.DEBUG
    logging.basicConfig(stream=stream, level=loglevel, format=format)
    logging.info("Log level set to %r(%d)" % (loglevel, verbose))


###############
# Methods for parsing tables to list (e.g. for screening)
##############
def add_screen_arguments(parser, defaults={}, accs=False):
    parser.add_argument(
        "-l",
        "--listFile",
        dest="listFile",
        metavar="LISTFILE",
        help="List of names to filter with"),
    parser.add_argument("-k", "--keep",
                        action="store_true", dest="keep", default=False,
                        help="Keep listed reads instead of removing")
    parser.add_argument(
        "-D",
        "--listDelim",
        dest="listDelim",
        default=None,
        help="list delimiter. If listColumn set, default is any "
             "whitespace, otherwise, the whole line (stripped of whietspace "
             "at the ends) is used. '\\t' will split on tab characters.",
        metavar="DELIM")
    parser.add_argument(
        "-C",
        "--listColumn",
        dest="listColumn",
        default=None,
        type=int,
        help="Column in listFile to get names from. Defaults to 0 if "
             "a delimiter is set.")
    parser.add_argument(
        "-G",
        "--galaxy",
        default=False,
        action="store_true",
        help="Column indices should start with 1")
    if accs:
        parser.add_argument("-a", "--accs",
                            action="store_true", dest="accs", default=False,
                            help="parse accession from read name in fasta")


dotRE = re.compile(r'(\.\d+)$')


def _get_key_type_translator(translate=False, regexp=dotRE):
    if translate:
        # hack the keyType to do our regex subs
        logging.debug("Stripping versions")
        return lambda r: dotRE.sub('', r)
    else:
        return None


def get_screen_list(arguments, accs=False):
    """
    Use the listFile, listDelim, and listColumn fields in arguments to
    parse a text file into a dict of names all mapped to True

    If accs==True, remove version (s/\\.\\d+$//) from accessions.

    Uses parse_list_to_set() to do the actual parsing
    """

    translator = _get_key_type_translator(accs, dotRE)

    if arguments.listDelim is not None or arguments.listColumn is not None:
        if arguments.listColumn is None:
            arguments.listColumn = 0
        if arguments.galaxy:
            arguments.listColumn -= 1

    return parse_list_to_set(
        arguments.listFile,
        col=arguments.listColumn,
        delim=arguments.listDelim,
        keyType=translator)


def parseAcc(read):
    m = accessionRE.search(read)
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
    for kid in sorted(node.children, **kwargs):
        if kid is not node:
            for n in treeGenerator(kid, kidsFirst=kidsFirst, **kwargs):
                yield n
    if kidsFirst:
        yield node


def returnSelf(x):
    return x


def pairwise(items, sortKey=returnSelf):
    """
    iterate over unique pairs of items
    """
    itemList = sorted(items, key=sortKey)
    for i in range(len(itemList)):
        for j in range(i + 1, len(itemList)):
            yield itemList[i], itemList[j]


def ascii_histogram(
        histogram,
        log=False,
        width=60,
        label='length',
        maxLabelWidth=10):
    (values, edges) = histogram[:2]

    maxValue = max(values)

    centers = [int(float(sum(edges[i:i + 2])) / 2.)
               for i in range(len(values))]
    largestLabel = max(max([len(str(c)) for c in centers]), len(label))
    if largestLabel < 6:
        largestLabel = 6
    elif largestLabel > maxLabelWidth:
        largestLabel = maxLabelWidth

    plotWidth = width - largestLabel + 1

    midPoint = npexp((nplog(maxValue) - nplog(.5)) / 2) \
        if log else maxValue / 2
    midPointStr = str(int(midPoint))
    padding_a_width = int(plotWidth / 2) - len(midPointStr) - len("count")
    padding_a = "".join([" " for i in range(padding_a_width)])
    maxValueStr = str(int(maxValue))
    padding_b_width = int(ceil(plotWidth / 2.)) - 1 - len(maxValueStr)
    padding_b = "".join([" " for i in range(padding_b_width)])
    output = "%s|count%s%s|%s%s|\n" % (rightPad(label, largestLabel),
                                       padding_a,
                                       midPointStr,
                                       padding_b,
                                       maxValueStr,
                                       )
    # output+="%s|%s\n" % ("".join(["_" for i in range(largestLabel)]),
    #                     "".join(["_" for i in range(plotWidth)]),
    #                     )
    for i, v in enumerate(values):
        output += "%s|%s\n" % (rightPad(str(centers[i]), largestLabel),
                               getBarString(v, maxValue, plotWidth, log))
    return output


logChars = ['-', '~', '=', '#']


def getBarString(value, maxValue, maxWidth, log):
    """
    return string of various signs (-,~,=,#) based on value and scale
    """
    if log:
        value = nplog(value) - nplog(.5) if value > 0 else 0
        maxValue = nplog(maxValue) - nplog(.5)
    width = maxWidth * (value / float(maxValue))
    if width < 1:
        return ''
    char = logChars[0]
    s = char
    while len(s) < width:
        if log:
            # print "s: %s, mw: %s, w: %s" % (s, maxWidth, width)
            char = logChars[
                int(ceil(len(logChars) * len(s) / float(maxWidth)) - 1)]
        s += char
    return s


def rightPad(name, width):
    if width < 6:
        width = 6
        logger.warn("Can't force names to be fewer than 6 characters")
    if len(name) > width:
        # remove middle and insert elipsis to get to -width- characters
        return name[:width - 4] + '***' + name[-1:]
    while len(name) < width:
        # pad with trailing space to get to 13 characters
        name += ' '
    return name

######
# manipulating large collections


def head(iterable, N=10):
    """
    yields the first N(=10) items of the given iterator or collection
    """

    iterator = iter(iterable)
    while N > 0:
        N -= 1
        try:
            yield(next(iterator))
        except StopIteration:
            break


def indexed_sample_generator(records, N, P=None):
    """
    Generate a random sample of N records from iterator "records".

    If the total numer of records in terator (P) is not specified,
    then reservoir sampling (storing entire sample in RAM) will be used!
    """
    if P > 0:
        # use random.sample to get list of indices ahead of time
        indexes_to_return = set(random.choice(range(P), N, replace=False))
        logger.debug("returning {} indices:\n{}"
                     .format(N, repr(indexes_to_return)))

        # return just records with tose indices
        record_count = 0
        yield_count = 0
        for i, record in enumerate(records):
            record_count += 1
            if i in indexes_to_return:
                yield_count += 1
                yield record
        logger.debug("Returned {} of {} records".format(yield_count,
                                                        record_count))

    else:
        # We don't know how many records to expect, so use reservoir
        for record in reservoir_sample(records, N=N, return_count=False):
            yield record


def reservoir_sample(iterator, N=100, return_count=False):
    """
    Randomly sample N items from an iterator of unknown size. Returns
    a list of the
    selected items.

    If the iterator returns fewer than N items, all will be in the sample,
    but the sample will not have N elements.

    IF return_count set to true, return value is tuple:
        (sample, totalItemsFromIterator)
    """

    sample = []

    numerator = float(N)
    for i, item in enumerate(iterator):
        n = i + 1
        val = random.rand()
        if n <= N:
            # Fill up sample with first N elements
            sample.append(item)
        elif val < numerator / n:
            # Replace random item in sample
            # with next element with probabliity N/i
            sample[random.random_integers(0, N - 1)] = item

    if return_count:
        return (sample, n)
    else:
        return sample


class ReservoirSamplingList(list):
    """
    Behaves mostly as a generic python container with one primary difference:
       it has a maximum size and
       only keeps a random subsample of added items once that size is reached.

    Other properties:

     * Insert() at an index is not supported.
     * Deleting items is not supported.
     * Only append() is allowed to change the contents
     * The total_added veriable tracks how many iterms were added
     * The order of items is not guaranteed to match the order in which
       they were added.
     * The order is nevertheless stable and elements can be retrieved by
       index or slice
    """

    def __init__(self, sample_size=100, iterable=[], preserve_order=False):
        """
        Create a new reservoirSampleingList with the given sample size
         (defaults to 100).

        If an iterable object is given as a second object,
         the reservoirSamplingList will be initialized with the items
         returned by iterating ofer the iterable.
        """
        self.N = int(sample_size)
        self.total_added = 0
        if preserve_order:
            self.append = self._append_preserve_order
        else:
            self.append = self._append
        for item in iterable:
            self.append(item)

    def __delitem__(self, key):
        raise Exception(
            "This container can only changed by appending items."
            "Deleting is not supported")

    def _append(self, item):
        """
        Add an item to the reservoirSamplingList.
        If the current size of the list is less than the sample_size, this
        item will simply be added.
        As the total added surpasses the sample_size, each new item
        will have a decreasing probability
         of replacing an existing item.
        This will simulate a subsample of the total set of added items.
        """
        self.total_added += 1
        if self.total_added <= self.N:
            # Fill up container with first N elements
            super(ReservoirSamplingList, self).append(item)
            return

        val = random.rand()
        if val < self.N / float(self.total_added):
            # Replace random item in sample
            # with next element with probabliity N/i
            self[random.random_integers(0, self.N - 1)] = item

    def _append_preserve_order(self, item):
        """
        Add an item to the reservoirSamplingList.
        If the current size of the list is less than the sample_size,
        this item will simply be added.
        As the total added surpasses the sample_size, each new item
        will have a decreasing probability
         of replacing an existing item.
        This will simulate a subsample of the total set of added items.
        """
        self.total_added += 1
        if self.total_added <= self.N:
            # Fill up container with first N elements
            super(ReservoirSamplingList, self).append(item)
            return

        val = random.rand()
        # Replace random item in sample
        # with next element with probabliity N/i
        if val < self.N / float(self.total_added):
            # pick random element to replace
            drop_index = random.random_integers(0, self.N - 1)
            # shift following items over
            for i in xrange(drop_index, len(self) - 1):
                self[i] = self[i + 1]
            # add new one at the end
            self[i + 1] = item

    def extend(self, iterable):
        for item in iterable:
            self.add(item)
