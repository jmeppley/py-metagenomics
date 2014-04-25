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
Given a list of regular expressions pulls the Nth (default is 24th) dataset from each matching history.
"""
    parser = OptionParser(usage, description=description)
    parser.add_option('-u', '--api_url', 
                      default="http://edminilims.mit.edu/api",
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

    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options, description)

    # set up dataset search values
    if options.dataset_regex is None:
        dsNum=options.dataset_index
        dsName=None
    else:
        dsNum=None
        dsName=re.compile(options.dataset_regex)

    histories=set()
    for pattern in args:
        # loop over matching histories
        logging.debug("Looking for histories that match '%s'" % pattern)
        for history in edl.galaxy.getHistories(options.api_key,
                                               options.api_url,
                                               re.compile(pattern),
                                               returnDict=True):
            
            historyId=history[u'id']
            if history[u'id'] in histories:
                logging.warn("Skipping history '%s(%s)', it was processed in a previous regex" % (history[u'name'],history[u'id']))
                continue
            histories.add(history[u'id'])

            logging.debug("Looking at history %s" % history)
            hdir = _get_output_dir(options, history)

            # loop over matching datasets
            for dataset in edl.galaxy.getDatasetData(options.api_key,
                                                     options.api_url,
                                                     historyId=historyId,
                                                     datasetNumber=dsNum,
                                                     datasetName=dsName):
                
                # Copy to local file from download URL
                url = dataset['download_url']
                out_file_name = _get_output_file(options,
                        hdir, dataset)
                response = requests.get(url, stream=True)
                logging.debug("Copying data from:\n%s" % (url))
                logging.info("Copyting to: %s" % (out_file_name))
                with open(out_file_name, 'wb') as out_file:
                    shutil.copyfileobj(response.raw, out_file)
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
