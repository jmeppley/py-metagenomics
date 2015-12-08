#!/usr/bin/env python
"""
Given a list of runs and an optional regular expression for sample names
run the samples through the Metagenomic pipeline in Galaxy
"""
import sys, logging, re, os
import edl.galaxy
from edl.util import addUniversalOptions,setupLogging
from optparse import OptionParser
logger=logging.getLogger(__name__)

def main():
    usage = "usage: %prog [OPTIONS] RUN_NAME [RUN_NAME ...]"
    description="""
Given a list of runs and an optional regular expression for sample names
run the samples through the Metagenomic pipeline in Galaxy
"""
    parser = OptionParser(usage, description=description)
    parser.add_option('-p','--pipelineVersion', default="",
                      help="VErsion of the 'MG Pipeline' to run")
    parser.add_option('-P', '--pipelineName', default=None,
                      help="Name of pipeline (wihout version string)")
    parser.add_option('-u', '--api_url', 
                      default="https://localhost/api",
                      help="URL of Galaxy API")
    parser.add_option('-k', '--api_key', default=None,
                      help="Galaxy API key for connecting. REQUIRED!")
    parser.add_option('-s', '--sample_regex', default=None,
                      help="Expression for matching sample names. If empty, all samples in given runs are processed")
    parser.add_option('-c', '--chemistry', dest='chem', default=None,
                      help="Override chemistry from SampleSheet with this value. One of 'truseq','scriptseq',or 'nextera'",
                      choices=['truseq','scriptseq','nextera'])

    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options, description)

    if options.sample_regex is None:
        sampleRE=None
    else:
        sampleRE=re.compile(options.sample_regex)

    if options.pipelineName is None:
        parser.error("Please supply a workflow name!")
    wfName=options.pipelineName + options.pipelineVersion
    hpref=re.sub(r'[^A-Z0-9]','',options.pipelineName)
    if options.pipelineVersion != "":
        hpref += '.' + re.sub(r'[^A-Z0-9]','',options.pipelineVersion)

    logger.debug("SampleRE:\t%r\nworkflow:\t%s\nhistPref:\t%s" % (sampleRE, wfName, hpref))

    # attempt to stifle security warnings
    try:
        logging.captureWarnings(True)
    except:
        pass

    if options.api_key is None:
        key = edl.galaxy.getApiKey()
        if key is None:
            parser.error("You must speicfy an API key with the -k flag!")
        else:
            options.api_key = key

    total=0
    for runName in args:
        wfcount=0
        for r in edl.galaxy.launchWorkflowOnSamples(options.api_key, runName,
                                                    workflowName=wfName, 
                                                    sampleRE=sampleRE,
                                                    apiURL=options.api_url,
                                                    historyPrefix=hpref,
                                                    chemistry=options.chem):
            logger.debug(repr(r))
            if 'error' in r:
                logger.warn(r['error'])
            else:
                wfcount+=1

        logger.info("Launched workflow on %d samples in %s" % (wfcount,runName))
        total+=wfcount

    logger.info("Launched workflow %s on %d samples from %s runs" % (wfName, total, len(args)))

if __name__ == '__main__':
    main()
