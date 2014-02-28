#!/usr/bin/python
#$ -cwd
#$ -V
#$ -p -10
#$ -S /usr/bin/python
#$ -l greedy=24

import sys, threading, logging, tempfile, subprocess, shutil
from edl.batch import fileTypeMap, getSizePerChunk, fragmentInputBySize, formatCommand, getFragmentPath
from edl.util import addUniversalOptions, setupLogging

binDir='/usr/local/bin'
#binDir='/common/bin'
lastBin='%s/lastal' % binDir
tantanBin='%s/tantan' % binDir
#scriptDir='/Users/jmeppley/work/delong/projects/scripts'
scriptDir='/minilims/home/jmeppley/work/delong/projects/scripts'
#scriptDir='/common/bin/scripts'
tmpDirRoot='/minilims/tmp'
#tmpDirRoot='/tmp'

protDefaults={'-b':'1','-x':'15','-y':'7','-z':'25'}
nuclDefaults={}

def main():
    ## set up CLI
    description = """
    Usage: lastWrapper.py [OPTIONS] -o OUTPUT_FILE LASTDB INPUT_FILE

    The last argument is taken as the input file and is fragmented into chunks on a local disk. Each chunk is run through lastal with all the given options (except output is writeen to the local disk). The results of each are comcateneted into the requested output file (the -o argument).

    Breaks input file into fragments to run last in pseudo multithreaded state. All lastal options (except -n) are accepted. Run lastal -h to see them. Additional options modify the batch behavior and post-porocessing. By default, output is converted to blast m8 like (aka 'gene') format and grouped by read.

The input file may be fasta or fastq. Fasta files are fragmented using the ">" character. FastqFiles are assumed to have four lines per record.

    The recommeded lastal options for reproducing BLASTX results are -b 1 -x 15 -y 7 -z 25, and these are invoked by default if the -F flag is used. To use different values, set them explicitly. This script will also mask the input fasta or fastq using tantan and pass '-u 2' to lastal. If the reads are already masked or to disable, supply a value for -u to this script. You MUST speicfy a frameshift penalty with -F if the database is protein. 15 is a good value.

    Batch Behavior:
    -C CHUNK_SIZE               Set the number of reads per chunk (defers to -N)
    -N NUM_CHUNKS               Set the number of threads (defaults to 4)
    For detailed options for fragmenting fasta, run fragmentRecords.py -h

    Post Processing:
    -f FORMAT                   'gene' (the default), 'blast', or 'liz' for blast-like m8.
                                '0' or '1' for lastal formats
    -O                          Original order. By default, the tabular formats
                                (ie '0','gene','blast','liz') are grouped by read and
                                sorted by score within reads.
    -n HITS_PER_READ            Maximim number of hits per read to keep
                                Defaults to 10. Set to -1 to turn off.

    Help/Info:
    -A, --about, -h, --help     This message

    """

    (options,args)=parseArgs()

    setupLogging(options, description, stream=sys.stdout)

    # Some basic checks
    if len(args)<2:
        raise Exception("Command does not seem long enough for lastal!")

    # Get last argument as input file name
    infile=args.pop(-1)
    logging.info("Reading sequences from: " + infile)
    outfile=options.outfile
    logging.info("Writing output to: %s " % outfile)

    #if options.verbose>1:
    #    args.insert(-1,'-v')
    if options.format=='1':
        if options.sort or options.maxHits > 0:
            logging.warn("Cannot sort or filter raw last output. Leaving untouched.")
    elif options.maxHits > 0:
        if not options.sort and options.format == '0':
            sys.exit("Cannot limit hits unless sorting or converting to M8 ('blast', 'gene', or 'liz')")

    # Apply any defaults not set by user
    if "-F" in options.userFlags:
        # this was a protein search
        defaultDict = protDefaults
    else:
        # this is a nucleotide search
        defaultDict = nuclDefaults
    for key in defaultDict:
        if key not in options.userFlags:
            args.insert(-1,key)
            args.insert(-1, defaultDict[key])

    ##
    # Fragment input file to temporary local dir
    if options.fastq:
        # fastq
        fileType=fileTypeMap['fastq']
    else:
        fileType=fileTypeMap['fasta']
    fragPref='fragment'
    insuff='.in'
    outsuff='.out'

    # create local tmp dir
    localdir = tempfile.mkdtemp(suffix='lastWrapper',dir=tmpDirRoot)

    # make sure we know how big to make chunks
    if options.chunk is None:
        if options.splits is None:
            logging.info("Defaulting to 4 chunks")
            options.splits=4
        options.chunk = getSizePerChunk(infile,options.splits,fileType,splitOnSize=options.splitOnSize)

    # mask with tantan unless mask is already set
    if options.mask is None:
        # add -u 2 to the command to use tantan results
        args.insert(-1,'-u')
        args.insert(-1,'2')
        # setup tantan command to pipe through fragmentInput
        command=[tantanBin,infile]
        logging.info("Masking with tantan")
        logging.debug(command)
        p=subprocess.Popen(command,stdout=subprocess.PIPE)
        instream=p.stdout
    else:
        # no masking (user has taken care of it)
        p=None
        instream=infile

    # fragment
    num=fragmentInputBySize(instream, localdir, options.chunk, fileType, fragPref, splitOnSize=options.splitOnSize, suffix=insuff)
    logging.info("Created %d fragments in %s" % (num, localdir))

    # check masking exit code if used
    if p is not None:
        ttCode = p.wait()
        if ttCode != 0:
            sys.exit("Tantan exited with code %d" % (ttCode))

    ## Run jobs
    # setup threads
    threads=[]
    for i in range(num):
        inFrag=getFragmentPath(localdir, fragPref, i+1, insuff)
        outFrag=getFragmentPath(localdir, fragPref, i+1, outsuff)
        # clone argument list and create command for this fragment
        cmd=list(args)
        # if post processing is needed, change command to string and pipe
        logging.debug("Sort check: %r %r" % (options.sort,options.format))
        if options.format=='0' and (options.sort or options.maxHits>0):
            cmd.append(inFrag)
            # sort (and possibly filter) last-formatted hit table
            cmd = "%s | %s" % (getCommandString(cmd),
                               getSortCommand(options.sort,options.maxHits))
            useShell=True
        elif options.format in ('blast','gene','liz'):
            cmd.append(inFrag)
            # convert to m8 (and sort)
            cmd = "%s | %s" % (getCommandString(cmd), getConvertCommand(options.format,
                                                                        options.sort,
                                                                        options.maxHits))
            useShell=True
        else:
            cmd.insert(-1,'-o')
            cmd.insert(-1,outFrag)
            cmd.append(inFrag)
            useShell=False

        # create thread
        threads.append(CommandThread(cmd,outFrag,shell=useShell))

    # start jobs
    for thread in threads:
        thread.start()

    # wait and collect output
    exitcode=0
    output=None
    for thread in threads:
        thread.join()
        if output is None:
            if outfile is None:
                output=sys.stdout
            else:
                output=open(outfile,'w')
        if thread.exitcode != 0:
            if thread.shell:
                logging.error("Command '%s' returned %s" % (thread.cmd, thread.exitcode))
            else:
                logging.error("Command '%s' returned %s" % (formatCommand(thread.cmd), thread.exitcode))
            exitcode=thread.exitcode
        else:
            logging.info("Thread %s completed!" % (str(thread)))

            threadstream=open(thread.outfile)
            # just copy
            for line in threadstream:
                output.write(line)
            threadstream.close()
    output.close()

    if options.verbose<=1:
        shutil.rmtree(localdir)

    sys.exit(exitcode)

##############
# Classes
##############
class CommandThread( threading.Thread ):
    def __init__(self, cmd, outfile, shell=False):
        self.exitcode=None
        self.cmd=cmd
        self.outfile=outfile
        self.shell=shell
        threading.Thread.__init__ ( self )

    def run(self):
        logging.debug("Launching thread: %s" % (self.cmd))
        if self.shell:
            outstream=open(self.outfile,'w')
        else:
            outstream=None

        self.exitcode=subprocess.call(self.cmd,shell=self.shell,stdout=outstream)
        if self.shell:
            outstream.close()

class Options():
    def __init__(self):
        self.verbose=1
        self.about=False
        self.chunk=None
        self.splits=None
        self.outfile=None
        self.splitOnSize=False
        self.help=False
        self.mask=None
        self.sort=True
        self.format='gene'
        self.userFlags=[]
        self.maxHits=10
        self.fastq=False

##############
# Methods
##############
def getOutfile(args):
    """
    Find the '-o' in the command string and get the output file name. Remove and return name.
    """
    for i,arg in enumerate(args):
        if arg.strip()=='-o':
            args.pop(i)
            return args.pop(i)
    else:
        return None

def parseArgs():
    options=Options()
    args=sys.argv
    args[0]=lastBin

    index=1
    while index < len(args):
        arg=args[index].strip()
        if arg=='-v':
            options.verbose+=1
        elif arg=='-C':
            args.pop(index)
            options.chunk=int(args.pop(index))
            continue
        elif arg=='-N':
            args.pop(index)
            options.splits=int(args.pop(index))
            continue
        elif arg=='-A' or arg=='--about':
            options.about=True
            break
        elif arg=='-o':
            args.pop(index)
            options.outfile=args.pop(index)
            continue
        elif arg=='-u':
            index+=1
            options.mask=args[index]
        elif arg=='-S':
            args.pop(index)
            options.splitOnSize=True
            continue
        elif arg=='-f':
            index+=1
            format=args[index]
            options.format=format
            if format in ('blast','gene','liz'):
                args[index]="1"
        elif arg=='-O':
            args.pop(index)
            options.sort=False
            continue
        elif arg=='-n':
            args.pop(index)
            options.maxHits=int(args.pop(index))
            continue
        elif arg=='-h' or arg=='--help':
            #options.help=True
            #printHelp([lastBin,'-h'])
            options.about=True
            break
        elif len(arg)>1 and arg[0]=='-':
            options.userFlags.append(arg[:2])
            if len(arg)==2:
                # if it's just a flag (eg '-x'),
                #  make sure we skip the value that follows
                index+=1
                # (Otherise, value is part of this string, eg. '-x34')
                val=args[index]
            else:
                val=arg[2:]
            if arg[:2]=='-Q' and val!='0':
                options.fastq=True
        index+=1

    return (options,args)

def printHelp(cmd):
    ## TODO
    pass

def convertOutput(raw,out,format,sort=True,maxHits=-1):
    """
    pipe output through liz's parser and then through sort
    """
    command=getConvertCommand(format,sort,maxHits)
    logging.debug("processing output with: %s" % (command))
    tmpstr=tempfile.TemporaryFile()
    try:
        subprocess.check_call(command, shell=True, stdin=raw, stdout=out, stderr=tmpstr)
    except:
        tmpstr.seek(0)
        logging.error("Conversion to M8 failed:\n%s" % "".join(tmpstr.readlines()))
        tmpstr.close()
        raise
    tmpstr.seek(0)
    logging.info("Conversion log:\n" + "".join(tmpstr.readlines()))
    tmpstr.close()

def getConvertCommand(format, sort, maxHits):
    """
    Return a command string to convert lastal output to (filtered and/or sorted) m8
    """

    command = "%s/LAST_output_converter.pl -q -pipe -f %s" % (scriptDir, format)
    if sort:
        if format=='liz':
            scoreCol=10
        elif format=='blast':
            scoreCol=9
        else:
            scoreCol=11

        command += " | sort -T %s -t $'\t' -k 1,1 -k %drn,%d" % (tmpDirRoot, scoreCol,scoreCol)

        if maxHits >= 0:
            command += " | %s/filter_blast_m8.py -H %s -f %s" % (scriptDir, maxHits, format)

    return command

def getSortCommand(sort,maxHits):
    """
    return command string to (filter and) sort lastal tabular output
    """
    if maxHits <= 0:
        maxHits=""
        if not sort:
            raise Exception("Code shouldn't get here if sort is False and maxHits<0")
    else:
        maxHits="| %s/filter_blast_m8.py -f last -H %s" % (scriptDir,maxHits)

    if sort:
        return "grep -v '^#' | sort -T %s -t $'\t' -k 7,7 -k 1rn,1 %s" % (tmpDirRoot,maxHits)
    else:
        return maxhits[2:]

def sortOutput(raw,out,sort=True,maxHits=-1):
    """
    Run output through grep and sort
    """
    subprocess.check_call(getSortCommand(sort,maxHits), shell=True, stdin=raw, stdout=out)

def getCommandString(cmd):
    """
    Given a list of strings (suitable for launching a subprocess), return a single long
    string for launching the same process as a shell command.
    E.G. ["ls","-la"] becomes "ls -la"

    Attempts to put quotes around elements with spaces in them
    """
    return ' '.join([cleanToken(t) for t in cmd])

def cleanToken(token):
    """
    Helper for getCommandString. Adds quotes and escapes as necessary for the shell to interpret as a single token.
    """
    if ' ' not in token:
        return token
    if '"' not in token:
        return '"%s"' % token
    if "'" not in token:
        return "'%s'" % token
    raise Exception("Sorry, cannot handle command tokens with both kinds of quotes!")

##############
if __name__ == '__main__':
    print sys.argv
    main()
