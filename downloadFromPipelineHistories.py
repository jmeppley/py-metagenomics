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
Given a list of regular expressions pulls the Nth (default is 24th) dataset from each matching history. If this is run on a machine that can access the galaxy files directory, it will create symlinks to the originals unles the "-c" option is given.
"""
    parser = OptionParser(usage, description=description)
    parser.add_option('-u', '--api_url', 
                      default="https://localhost/api",
                      help="URL of Galaxy API")
    parser.add_option('-k', '--api_key', default=None,
                      help="Galaxy API key for connecting. REQUIRED!")
    parser.add_option('-d', '--dataset_index', default=24, type='int',
                      help="which dataset to get in each history")
    parser.add_option('-D', '--dataset_regex', default=None,
                      help="Expression for matching dataset names. If set, overrides the dataset_index (-d)")
    parser.add_option("-O",'--saveDir', default=".",
            help="Direcotry in which to save files. A subdirectory will be created for each matchin history. Defaults to the current directory")
    parser.add_option("-o", "--outfileName", default=None,
                     help="If set, give this name to output files, otherwise pull name from galaxy")
    parser.add_option("-C", "--chunk_size", default=1024, type='int', 
            help="Chunk size for file downloads. Bigger should speed up the download of large files, but slow down lots of small files. Defaults to 1024")
    parser.add_option("-c", "--force_copy", default=False, action='store_true',
            help="Copy data files even if symlinks are possible")

    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options, description)

    if options.api_key is None:
        parser.error("You must speicfy an API key with the -k flag!")

    # set up dataset search values
    if options.dataset_regex is None:
        dsNum=options.dataset_index
        dsName=None
    else:
        dsNum=None
        dsName=re.compile(options.dataset_regex)

    # Are we on the same machine as galaxy?
    filesAreLocal = re.search(r'://localhost',options.api_url) is not None

    for (history,dataset) in edl.galaxy.findDatasets(options.api_key, args,
                                                     dsName=dsName,
                                                     dsNum=dsNum,
                                                     apiURL=options.api_url):

        # create dir for history if it doesn't already exist
        hdir = _get_output_dir(options, history)

        # generate name for downloaded/linked file
        out_file_name = _get_output_file(options,
                hdir, dataset)

        # Create symlink if we can find the file locally
        if filesAreLocal and 'file_name' in dataset:
            originalFile=dataset['file_name']
            os.symlink(originalFile, out_file_name)
        else:
            # Copy to local file from download URL
            url = dataset['download_url']
            response = requests.get(url, stream=True, verify=False)
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
