#!/usr/bin/env python
"""
Given a list of runs and an optional regular expression for sample names
run the samples through the Metagenomic pipeline in Galaxy
"""
import sys, logging, re, os
import shutil, requests
import edl.galaxy
from edl.util import addUniversalOptions,setupLogging
from optparse import OptionParser

def main():
    usage = "usage: %prog [OPTIONS] EXPR [EXPR ...]"
    description="""
Given a list of regular expressions pulls datasets from each matching history. Any datasets matching the dataset_regex (-d) or one of the dataset_indices (-D) will be retrieved. If this is run on a machine that can access the galaxy files directory, it will create symlinks to the originals unles the "-c" option is given.
"""
    parser = OptionParser(usage, description=description)
    parser.add_option('-u', '--api_url', 
                      default="https://localhost/api",
                      help="URL of Galaxy API")
    parser.add_option('-k', '--api_key', default=None,
                      help="Galaxy API key for connecting. REQUIRED!")
    parser.add_option('-d', '--dataset_index', dest='dataset_indices',action='append', default=[], type='int',
                      help="which dataset to get in each history, can be invoked multiple times for multiple datasets")
    parser.add_option('-D', '--dataset_regex', default=None,
                      help="Expression for matching dataset names. ")
    parser.add_option("-S","--sameDir", default=False, action='store_true',
            help="Save all datasets to same directory, don't created subdirs for each history")
    parser.add_option("-O",'--saveDir', default=".",
            help="Direcotry in which to save files. A subdirectory will be created for each matchin history. Defaults to the current directory")
    parser.add_option("-o", "--outfileName", default=None,
                     help="If set, give this name to output files, otherwise pull name from galaxy")
    parser.add_option("-C", "--chunk_size", default=1024, type='int', 
            help="Chunk size for file downloads. Bigger should speed up the download of large files, but slow down lots of small files. Defaults to 1024")
    parser.add_option("-l", "--force_link", default=False, action='store_true',
            help="Symlink data files even if URL is not local")
    parser.add_option("-c", "--force_copy", default=False, action='store_true',
            help="Copy data files even if symlinks are possible")

    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options, description)

    if options.force_copy and options.force_link:
        parser.error("Please use only one of -l (--force_link) and -c (--force_copy)") 

    if options.api_key is None:
        key = edl.galaxy.getApiKey()
        if key is None:
            parser.error("You must speicfy an API key with the -k flag!")
        else:
            options.api_key = key

    # set up dataset search values
    if options.dataset_regex is None:
        if len(options.dataset_indices)==0:
            parser.error('Please supply a dataset number or dataset string with -d or -D!')

    logging.debug("Looking for dataset: %r/%r" % (options.dataset_regex, 
                                                  options.dataset_indices))

    # Are we on the same machine as galaxy?
    filesAreLocal = re.search(r'://(localhost|127\.0\.0\.1)',options.api_url) is not None or options.force_link
    logging.debug("URL seems local, will attempt to link")

    for (history,dataset) in edl.galaxy.findDatasets(options.api_key, args,
                                                     dsName=options.dataset_regex,
                                                     dsNums=options.dataset_indices,
                                                     apiURL=options.api_url):

        # create dir for history if it doesn't already exist
        hdir = _get_output_dir(options, history)

        # generate name for downloaded/linked file
        out_file_name = _get_output_file(options,
                hdir, dataset)

        # Create symlink if we can find the file locally
        if filesAreLocal and 'file_path' in dataset and not options.force_copy:
            logging.debug("Linking dataset")
            originalFile=dataset['file_path']
            os.symlink(originalFile, out_file_name)
        else:
            # Copy to local file from download URL
            url = dataset['download_url']
            response = requests.get(url, stream=True, verify=False)
            logging.debug(repr(dataset))
            logging.debug("Copying data from:\n%s" % (url))
            logging.info("Copyting to: %s" % (out_file_name))
            #with open(out_file_name, 'wb') as out_file:
            #    shutil.copyfileobj(response.raw, out_file)
            with open(out_file_name, 'wb') as out_file:
                for chunk in response.iter_content(options.chunk_size):
                    out_file.write(chunk)
            del response
                                   
def _get_output_dir(options, history):
    """
    Helper method to create an output directory name for the hirstory
    """

    # Create history dir, add id if necessary
    if options.sameDir:
        hdir = options.saveDir
    else:
        hname=_sanitize(history[u'name'])
        hname = hname + "_" + history[u'id']
        hdir=os.path.sep.join((options.saveDir,hname))
    if not os.path.exists(hdir):
        os.makedirs(hdir)
    return hdir

def _get_output_file(options, hdir, dataset):
    """
    Helper method to create an output file name
    """

    # get file name, add index if necessary
    if options.outfileName is None:
        dsname=_sanitize(dataset[u'name'])
    else:
        dsname=options.outfileName
    out_file_name=os.path.sep.join((hdir, dsname))
    index=0
    while os.path.exists(out_file_name):
        index+=1
        out_file_name=os.path.sep.join((hdir, "%s_%d" % (dsname,index)))

    return out_file_name

def _sanitize(string):
    return re.sub(r'[^\w.-_+:]','_',string)

if __name__ == '__main__':
    main()
