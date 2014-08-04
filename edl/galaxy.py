import logging
logger=logging.getLogger(__name__)

import simplejson as json
import urllib2, re, os, sys
from bioblend.galaxy import GalaxyInstance

# For retrieving data
def findDatasets(apiKey, patterns, dsName=None, dsNum=None, apiURL='http://localhost/api'):
    histories=set()
    for pattern in patterns:
        # loop over matching histories
        logger.debug("Looking for histories that match '%s'" % pattern)
        for history in getHistories(apiKey,
                                    apiURL,
                                    re.compile(pattern),
                                    returnDict=True):
            
            historyId=history[u'id']
            if history[u'id'] in histories:
                logger.warn("Skipping history '%s(%s)', it was processed in a previous regex" % (history[u'name'],history[u'id']))
                continue
            histories.add(history[u'id'])

            logger.debug("Looking at history %s" % history)

            # loop over matching datasets
            for dataset in getDatasetData(apiKey,
                                          apiURL,
                                          historyId=historyId,
                                          datasetNumber=dsNum,
                                          datasetName=dsName):
                yield (history, dataset)

def getDatasetFile(apiKey, apiURL, historyName, datasetNumber, returnURL=False, returnDict=False):
    """
    Given the URL and KEY for a Galaxy instance's API and:
        A history name
        A history item index
    return a file handle to the first matched file (or just the URL)

   Uses getDatasetData() to find a given file. Only takes hitoryName and datasetNumber. Can also return just the URL (instead of an open file-like-objet)
    """
    for dataset in getDatasetData(apiKey, apiURL, historyName=historyName, datasetNumber=datasetNumber):
        durl=dataset['download_url']
        if returnURL:
            return durl
        else:
            return urllib2.urlopen(durl)

    else:
        raise Exception("No match found for item %d in history '%s'" % (datasetNumber, historyName))

def getDatasetData(apiKey, apiURL, historyName=None, historyId=None, datasetNumber=None, datasetName=None):
    """
    Given the URL and KEY for a Galaxy instance's API and:
        A history (name or ID)
        A history item (index integer, name, or regex to match name) 
    return the data dicts (as an iterator/generator) for any matching datasets

    Under the hood, this is an HTTP connection using urllib2. 

    HistoryName can be a string or compiled re object. It will be ignored if historyId is not None. The same is true for datasetNumber and datasetName.
    """

    # Collate history IDs
    historyIds=[]
    if historyId is None:
        for history in getHistories(apiKey, apiURL, historyName, 
                returnDict=True):
            historyIds.append(history[u'id'])
        if len(historyIds)==0:
            raise Exception("Cannot find history: %s" % historyName)
    else:
        historyIds.append(historyId)

    # Loop over found histories
    for historyId in historyIds:
        hurl=apiURL+"/histories/"+historyId+"/contents?key="+apiKey
        logger.debug("Searching for datasets in url: %s" % hurl)
        count=0
        for dataset in json.loads(urllib2.urlopen(hurl).read()):
            logger.debug("Getting details with url: %s" % hurl)
            hurl = apiURL+"/histories/" + historyId + "/contents/" + dataset[u'id'] + "?key=" + apiKey
            details = json.loads(urllib2.urlopen(hurl).read())
            if _dataset_match(details, datasetNumber, datasetName):
                durl = details['download_url'][4:]
                if "?" in durl:
                    qstringsep="&"
                else:
                    qstringsep="?"
                details['download_url'] = apiURL + durl + qstringsep + "key=" + apiKey
                count+=1
                yield details

        if count==0:
            logger.warn("No matching dataset in history: %s" % historyId)

def _dataset_match(dataset, datasetNumber, datasetName):
    """
    helper method to check for dataset number or name in dataset details
    """
    if datasetNumber is not None:
        return dataset[u'hid']==datasetNumber
    elif isinstance(datasetName,str):
        return dataset[u'name']==datasetName
    else:
        return datasetName.search(dataset[u'name'])!=None

def getHistories(apiKey, apiURL, nameRE=None, returnDict=False):
    """
    Return an iterator over histories.
    
    Optionally filter list with a regular expression (nameRE). If nameRE is a string, it must match exactly.
    """
    hurl=apiURL+"/histories?key="+apiKey
    logger.debug("Searching for histories at url: %s" % hurl)
    for history in json.loads(urllib2.urlopen(hurl).read()):
        if nameRE is not None:
            if isinstance(nameRE, str):
                # just do exact match
                if nameRE != history[u'name']:
                    continue
            else:
                m=nameRE.search(history[u'name'])
                if m is None:
                    continue

        if u'deleted' in history and history[u'deleted']:
            # Skip delted histories
            continue

        if returnDict:
            yield history
        else:
            yield history[u'name']

# Running workflows on MiSeq data 
#  These are pretty specific to our setup

def locateDatasets(runName, galaxyInstance, libraryNameTemplate = "MiSeq Run: %s", sampleRE=None):
    """
    Parses a MiSeq run in the Galaxy shared library area and returns a dictionary keyed on the sample number (as a string, eg: '1'). Each entry is a dictionary of the form: {'name': __, 'R1': __, 'R2': __}, where the last two blanks are encoded Galaxy IDs for use in the API.
    """
    files={}
    libName=libraryNameTemplate % (runName)
    libs=galaxyInstance.libraries.get_libraries(name=libName)
    if len(libs)==0:
        raise Exception("No library named %s" % (libName))
    libID=libs[0][u'id']
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
    chemistry='scriptseq'
    barcodes={}
    with open(sampleSheet) as f:

        # find Assay line
        line = ""
        while re.match(r'Assay,',line) is None:
            line = f.next()

        # not the fastest approach, but should be clear
        if re.search(r'[Nn]extera',line) is not None:
            print (line)
            chemistry='nextera'
        elif re.search(r'[Tt]rue?[Ss]eq',line) is not None:
            print (line)
            chemistry='truseq'
        elif re.search(r'[Ss]cript?[Ss]eq',line) is not None:
            print (line)
            chemistry='scriptseq'
        
        # Skip ahead to sample table
        line = ""
        while re.match(r'\[Data\]',line) is None:
            line = f.next()
            
        # parse first line as headers
        headers = f.next().split(',')
        indexIndex = headers.index('index')
        if 'index2' in headers:
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

def launchWorkflowOnSamples(apiKey, runName, workflowID=None, workflowName=None, historyPrefix='MGP.b011', apiURL=u'http://edminilims.mit.edu/api', chemistry=None, **kwargs):
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
    primerToolID=getPrimerToolId(galaxyInstance, workflowID)
    workflowInputs = getWorkflowInputs(galaxyInstance, workflowID)
    files = locateDatasets(runName, galaxyInstance,**kwargs)
    (ssChemistry,barcodes) = parseSampleSheet(runName,**kwargs)
    if chemistry is None:
        chemistry=ssChemistry
    else:
        logger.info("Using user specified chemistry: %s" % chemistry)
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
                logger.warn("History already exists: %s" % historyName)
                raise Exception("History already exists: %s" % historyName)
            
            dsMap={}
            for (direction, inputid) in workflowInputs.iteritems():
                if direction not in sampleData:
                    raise Exception("Missing direction in input: %s" % (direction))
                dsMap[inputid] = {'src':'ld','id':sampleData[direction]}
    
            data={'workflow_id':workflowID,
                  'history':historyName,
                  'ds_map':dsMap,
                  'replacement_params':{'name':sampleName},
                  'parameters':{primerToolID:{'chemistry':chemistry,
                                              'barcode1':barcode1,
                                              'barcode2':barcode2}}
            }
            response['data']=data
            response['response']=submit(apiKey, apiURL+"/workflows",data)
        except Exception as e:
            response['error']=e
        yield response

def getPrimerToolId(gi, workflowID):
    for stepid, step in gi.workflows.show_workflow(workflowID)[u'steps'].iteritems():
        if u'tool_id' not in step or step[u'tool_id'] is None:
            continue
        if re.search(r'CreatePrimerFile',step[u'tool_id']):
            return step[u'tool_id']
    else:
        raise Exception("Could not find CreatePrimerFile tool in workflow!")

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


