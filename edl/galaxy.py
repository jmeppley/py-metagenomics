import logging
logger=logging.getLogger(__name__)

import simplejson as json
import urllib2, re, os, sys
from bioblend.galaxy import GalaxyInstance

# For retrieving data
def getDatasetFile(apiKey, apiURL, historyName, datasetNumber):
    """
    Given the URL and KEY for a Galaxy instance's API and:
        A history name
        A history item index
    return a file handle to that file

    Under the hood, this is an HTTP connection using urllib2. If there are multiple histories with the same name, it will only look in the first one returned by the API.
    """
    hurl=apiURL+"/histories?key="+apiKey
    logger.debug("Searching for histories at url: %s" % hurl)
    for history in json.loads(urllib2.urlopen(hurl).read()):
        if history[u'name']==historyName:
            break
    else:
        raise Exception("Cannot find history: %s" % historyName)

    hurl=apiURL+"/histories/"+history[u'id']+"/contents?key="+apiKey
    logger.debug("Searching for datasets in url: %s" % hurl)
    for dataset in json.loads(urllib2.urlopen(hurl).read()):
        logger.debug("Getting details with url: %s" % hurl)
        hurl = apiURL+"/histories/" + history[u'id']+"/contents/" + dataset[u'id'] + "?key=" + apiKey
        details = json.loads(urllib2.urlopen(hurl).read())
        if details[u'hid']==datasetNumber:
            break
    else:
        raise Exception("Cannot find dataset %d in history %s" % (datasetNumber, historyName))
    if 'download_url' in details:
        durl = details['download_url'][4:]
        if "?" in durl:
            qstringsep="&"
        else:
            qstringsep="?"
        hurl = apiURL + durl + qstringsep + "key=" + apiKey
        logger.debug("Download with URL: %s" % hurl)
        return urllib2.urlopen(hurl)
    else:
        raise Exception("Dataset does not have a download url!")


# Running workflows on MiSeq data 
#  These are pretty specific to our setup

def locateDatasets(runName, galaxyInstance, libraryNameTemplate = "MiSeq Run: %s", sampleRE=None):
    """
    Parses a MiSeq run in the Galaxy shared library area and returns a dictionary keyed on the sample number (as a string, eg: '1'). Each entry is a dictionary of the form: {'name': __, 'R1': __, 'R2': __}, where the last two blanks are encoded Galaxy IDs for use in the API.
    """
    files={}
    libID=galaxyInstance.libraries.get_libraries(name=libraryNameTemplate % (runName))[0][u'id']
    for item in galaxyInstance.libraries.show_library(libID,contents=True):
        # skip folders and other things that are not files
        if item[u'type']!=u'file':
            continue
        
        # Parse name into components
        m = re.search(r'([^\/]+)_S(\d+)_L\d+_(R[12])_.+\.fastq', item[u'name'])
        if m:
            (name, number, direction) = m.groups()
            if sampleRE is not None:
                if sampleRE.search(name) is None:
                    continue
            fileId=item[u'id']
            if number in files:
                files[number][direction]=fileId
            else:
                files[number]={'name':name,direction:fileId}
    return files

def parseSampleSheet(runName,**kwargs):
    """
    Parses a SampleSheet from the file system for the given run. The run name should correspond to a subdirectory of the 'dataDir' which defaults to /minilims/data/incoming/miseq. 

    Returns a two-element tuple:
        chemistry: either 'nextera' or 'truseq'
        barcodes: dictionary from sample number (as integer) to two-element list of barcodes. Second element will be an empty string for truseq.
    """
    dataDir=kwargs.get('dataDir','/minilims/data/incoming/miseq')
    sampleSheet=os.path.sep.join([dataDir,runName,'SampleSheet.csv'])
    chemistry='truseq'
    barcodes={}
    with open(sampleSheet) as f:
        
        # Skip ahead to sample table
        line = ""
        while re.match(r'\[Data\]',line) is None:
            line = f.next()
            
        # parse first line as headers
        headers = f.next().split(',')
        indexIndex = headers.index('index')
        if 'index2' in headers:
            # overly simplistic logic: two barcodes => Nextera
            chemistry = 'nextera'
            index2Index = headers.index('index2')
            
        # get the barcode for each sample
        for i,line in enumerate(f):
            cells=line.split(',')
            if chemistry=='nextera':
                barcodes[i+1]=[cells[indexIndex],cells[index2Index]]
            else:
                # put an empty string for unused barcode
                barcodes[i+1]=[cells[indexIndex],""]
    return (chemistry, barcodes)


def launchWorkflowOnSamples(apiKey, runName, workflowID=None, workflowName=None, historyPrefix='MGP.b011', apiURL=u'http://localhost:8443/api', **kwargs):
    """
    Given:
        apiKey: the connection code for the Galaxy API
        runName: the folder name of a MiSeq run
        either the workflowName or encoded workflowID
        optionally: 
            apiURL or historyPrefix for naming the created histories
            sampleRE (compiled regex obj) for filtering samples

    Requiremnts:
        run must be loaded into galaxy
        workflow must have two inputs (for forward and reverse reads)
        workflow must have the workflow parameter 'name'
        workflow must have a 'CreatePrimerFile' step

    Runs the given workflow on each sample in the MiSeq run. The results are put into a new history with the name: "historyPrefix: runName: sampleName".
    """
    historyTemplate='%s: %s: %s'
    galaxyInstance = GalaxyInstance(apiURL, apiKey)
    logger.debug("workflowId: %s" % workflowID)
    if workflowID is None:
        if workflowName is None:
            raise Exception('Please supply a workflow ID or Name!')
        logger.debug("Looking up workflow: %s" % workflowName)
        workflow = getWorkflow(galaxyInstance, workflowName)
        workflowID = workflow[u'id']
        logger.info("Found workflow: %s" % workflowID)
    workflowInputs = getWorkflowInputs(galaxyInstance, workflowID)
    files = locateDatasets(runName, galaxyInstance,**kwargs)
    (chemistry,barcodes) = parseSampleSheet(runName,**kwargs)
    responses=[]
    for (sample, sampleData) in files.iteritems():
        response={'sample':sample}
        try:
            sample = int(sample)
            if sample not in barcodes:
                #barcode = "NNNNNN"
                raise Exception("No barcode for sample %s" % (sample))
            else:
                (barcode1,barcode2) = barcodes[sample]
            sampleName=sampleData[u'name']

            # check if history exists
            historyName = historyTemplate % (historyPrefix, runName, sampleName)
            for history in galaxyInstance.histories.get_histories(name=historyName):
                raise Exception("History already exists: %s" % historyName)
            
            dsMap={}
            for (direction, inputid) in workflowInputs.iteritems():
                if direction not in sampleData:
                    raise Exception("Missing direction in input: %s" % (direction))
                dsMap[inputid] = {'src':'ld','id':sampleData[direction]}
    
            data={'workflow_id':workflowID,
                  'history':historyName,
                  'ds_map':dsMap,
                  'parameters':{'wf_parm':{'name':sampleName},
                                u'CreatePrimerFile':{'chemistry':chemistry,'barcode1':barcode1,'barcode2':barcode2}}
            }
            response['data']=data
            response['response']=submit(apiKey, apiURL+"/workflows",data)
        except Exception as e:
            response['error']=e
        yield response

def getWorkflowInputs(galaxyInstance, workflowId):
    workflowInputs = {}
    inputs = galaxyInstance.workflows.show_workflow(workflowId)[u'inputs']
    for (i,inputId) in enumerate(sorted(inputs.keys())):
        direction = u'R%d' % (i+1)
        workflowInputs[direction] = inputId
    return workflowInputs

def getWorkflow(galaxyInstance, name):
    """ Return the first workflow with this name """
    for workflow in galaxyInstance.workflows.get_workflows():
        if workflow[u'name']==name:
            break
    else:
        raise Excepton("Cannot find workflow: %s" % (name))
    return workflow 

####
# Copied from galaxy api common.py
def make_url( api_key, url, args=None ):
    # Adds the API Key to the URL if it's not already there.
    if args is None:
        args = []
    argsep = '&'
    if '?' not in url:
        argsep = '?'
    if '?key=' not in url and '&key=' not in url:
        args.insert( 0, ( 'key', api_key ) )
    return url + argsep + '&'.join( [ '='.join( t ) for t in args ] )

def post( api_key, url, data ):
    # Do the actual POST.
    url = make_url( api_key, url )
    req = urllib2.Request( url, headers = { 'Content-Type': 'application/json' }, data = json.dumps( data ) )
    return json.loads( urllib2.urlopen( req ).read() )

def submit( api_key, url, data, return_formatted=True ):
    # Sends an API POST request and acts as a generic formatter for the JSON response.
    # 'data' will become the JSON payload read by Galaxy.
    try:
        r = post( api_key, url, data )
    except urllib2.HTTPError, e:
        if return_formatted:
            print e
            print e.read( 1024 )
            sys.exit( 1 )
        else:
            return 'Error. '+ str( e.read( 1024 ) )
    if not return_formatted:
        return r
    print 'Response'
    print '--------'
    if type( r ) == list:
        # Currently the only implemented responses are lists of dicts, because
        # submission creates some number of collection elements.
        for i in r:
            if type( i ) == dict:
                if 'url' in i:
                    print i.pop( 'url' )
                else:
                    print '----'
                if 'name' in i:
                    print '  name: %s' % i.pop( 'name' )
                for k, v in i.items():
                    print '  %s: %s' % ( k, v )
            else:
                print i
    else:
        print r


