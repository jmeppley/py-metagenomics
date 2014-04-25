#!/usr/bin/env python
"""
Given a list of runs and an optional regular expression for sample names
run the samples through the Metagenomic pipeline in Galaxy
"""
import sys, logging, re
import edl.galaxy
from edl.util import addUniversalOptions,setupLogging
from optparse import OptionParser

def main():
    usage = "usage: %prog [OPTIONS] RUN_NAME [RUN_NAME ...]"
    description="""
Given a list of runs and an optional regular expression for sample names
run the samples through the Metagenomic pipeline in Galaxy
"""
    parser = OptionParser(usage, description=description)
    parser.add_option('-p','--pipelineVersion', default="Beta 011",
                      help="VErsion of the 'MG Pipeline' to run")
    parser.add_option('-P', '--pipelineName', default="MG Pipeline ",
                      help="Name of pipeline (wihout version string)")
    parser.add_option('-u', '--api_url', 
                      default="http://edminilims.mit.edu/api",
                      help="URL of Galaxy API")
    parser.add_option('-k', '--api_key', default=None,
                      help="Galaxy API key for connecting. REQUIRED!")
    parser.add_option('-s', '--sample_regex', default=None,
                      help="Expression for matching sample names. If empty, all samples in given runs are processed")

    addUniversalOptions(parser)

    (options, args) = parser.parse_args()

    setupLogging(options, description)

    if options.sample_regex is None:
        sampleRE=None
    else:
        sampleRE=re.compile(options.sample_regex)
    wfName=options.pipelineName + options.pipelineVersion
    hpref=re.sub(r'[^A-Z0-9]','',options.pipelineName) + '.' + \
          re.sub(r'[^A-Z0-9]','',options.pipelineVersion)

    logging.debug("SampleRE:\t%r\nworkflow:\t%s\nhistPref:\t%s" % (sampleRE, wfName, hpref))

    total=0
    for runName in args:
        wfcount=0
        for r in edl.galaxy.launchWorkflowOnSamples(options.api_key, runName,
                                                    workflowName=wfName, 
                                                    sampleRE=sampleRE,
                                                    apiURL=options.api_url,
                                                    historyPrefix=hpref):
            wfcount+=1
            logging.debug(repr(r))

        logging.info("Launched workflow on %d samples in %s" % (wfcount,runName))
        total+=wfcount

    logging.info("Launched workflow %s on %d samples from %s runs" % (wfName, total, len(args)))

if __name__ == '__main__':
    main()
