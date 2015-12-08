import logging
logger=logging.getLogger(__name__)

import simplejson as json
import urllib2, re, os, sys
from bioblend.galaxy import GalaxyInstance
from httplib import IncompleteRead

def getApiKey(apiKeyFile = '.galaxy_api_key', apiEnvVar = 'GALAXY_API_KEY'):
    """
    Check the filesystem for a .galaxy_api_key file (in . and ~)
    Then check the GALAXY_API_KEY env variable
    """
    #  look for .galaxy_api_key file
    api_key=None
    if os.path.exists(apiKeyFile):
        with open(apiKeyFile) as f:
            api_key = f.next().strip()
    elif os.path.exists(os.environ.get('HOME','.') + os.path.sep + apiKeyFile):
        with open(os.environ.get('HOME','.') + os.path.sep + apiKeyFile) as f:
            api_key = f.next().strip()
    elif apiEnvVar in os.environ:
        api_key = os.environ[apiEnvVar].strip()

    logging.debug("Found API key: %s" % api_key)
    return api_key

# For retrieving data
def findDatasets(apiKey, patterns, 
                 dsName=None, 
                 dsNum=None, 
                 apiURL='http://localhost/api',
                 skipDeleted=True):
    """
    Return an iterator over (history,dataset) tuples where the history name matches one of the given patterns and the dataset is identified by name or number.
    """
    histories=set()
    for pattern in patterns:
        # loop over matching histories
        logger.debug("Looking for histories that match '%s'" % pattern)
        for history in getHistories(apiKey,
                                    apiURL,
                                    re.compile(pattern),
                                    returnDict=True,
                                    skipDeleted=skipDeleted):
            
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
                                          datasetName=dsName,
                                          skipDeleted=skipDeleted):
                yield (history, dataset)

def getDatasetFile(apiKey, apiURL, historyName, datasetNumber, 
                   skipDeleted=True, returnURL=False):
    """
    Given the URL and KEY for a Galaxy instance's API and:
        A history name
        A history item index
    return a file handle to the first matched file (or just the URL)

   Uses getDatasetData() to find a given file. Only takes hitoryName and datasetNumber. Can also return just the URL (instead of an open file-like-objet)
    """
    for dataset in getDatasetData(apiKey, apiURL, 
                                  historyName=historyName, 
                                  datasetNumber=datasetNumber,
                                  skipDeleted=skipDeleted):
        durl=dataset['download_url']
        if returnURL:
            return durl
        else:
            return urllib2.urlopen(durl)

    else:
        raise Exception("No match found for item %d in history '%s'" % (datasetNumber, historyName))

def getDatasetData(apiKey, apiURL, 
                   historyName=None, 
                   historyId=None, 
                   datasetNumber=None, 
                   datasetName=None, 
                   skipDeleted=True):
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
            if (details['deleted'] or details['state']=='error') and skipDeleted:
                continue
            if _dataset_match(details, datasetNumber, datasetName):
                if 'download_url' not in details:
                    logger.warn("Can't find the download URL in %r\nHistory: %r" % (details,historyId))
                    continue
                details['download_url']=fixDownloadUrl(details['download_url'],
                                                       apiURL,
                                                       apiKey)
                count+=1
                yield details

        if count==0:
            logger.warn("No matching dataset in history: %s" % historyId)

def fixDownloadUrl(downloadUrl, apiUrl, apiKey):
    durl = downloadUrl[4:]
    if "?" in durl:
        qstringsep="&"
    else:
        qstringsep="?"
    return apiUrl + durl + qstringsep + "key=" + apiKey

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

def getHistories(apiKey, apiURL, 
                 nameRE=None, 
                 returnDict=False, 
                 skipDeleted=True):
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

        if u'deleted' in history and history[u'deleted'] and skipDeleted:
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
    Parses a MiSeq or NextSeq run in the Galaxy shared library area and returns a dictionary keyed on the sample number (as a string, eg: '1'). Each entry is a dictionary of the form: 
    {'name': __, 'lanes':{
        '001':{'R1': __, 'R2': __},
        ...}
    }
    where the last two blanks are encoded Galaxy IDs for use in the API. MiSeq runs will only have the 001 lane, NextSeq runs should have up to 004.
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
        m = re.search(r'([^\/]+)_S(\d+)_(?:L(\d+)_)?(R[12])_.+\.fastq', 
                      item[u'name'])
        if m:
            (name, number, lane, direction) = m.groups()

            logger.debug("Matched %s:%s" % (number, name))
            if sampleRE is not None:
                if sampleRE.search(name) is None:
                    logger.debug("%s not found in %s" % (sampleRE.pattern, 
                                                         name))
                    continue
            fileId=item[u'id']
            if number in files:
                if lane in files[number]['lanes']:
                    files[number]['lanes'][lane][direction]=fileId
                else:
                    files[number]['lanes'][lane]={direction:fileId}
            else:
                files[number]={'name':name,'lanes':{lane:{direction:fileId}}}

        # look for sample sheet, too
        m = re.search(r'SampleSheet.csv$',item[u'name'])
        if m:
            files['sampleSheet'] = item[u'id']

        logger.debug("Found %d files" % (len(files)))
    
    return files

def loadSampleSheetFromGalaxy(sampleSheetId, galaxyInstance):
    """
    The old bioblend GalaxyInstance can't do this without help.

    I borrowed this code from bioblend.galaxy.objects
    """
    logging.debug("Attempting to download SampleSheet %s from galaxy" % (sampleSheetId))
    base_url=galaxyInstance._make_url(galaxyInstance.libraries) + \
            "/datasets/download/uncompressed"
    kwargs = {'stream': True,
              'params': {'ld_ids%5B%5D': sampleSheetId},
              }
    r = galaxyInstance.make_get_request(base_url, **kwargs)
    if r.status_code == 500:
        # compatibility with older Galaxy releases
        kwargs['params'] = {'ldda_ids%5B%5D': sampleSheetId}
        r = galaxyInstance.make_get_request(base_url, **kwargs)
    r.raise_for_status()
    return (r.iter_lines(), r.close)

def loadSampleSheetFromFileSystem(runName,**kwargs):
    dataDir=kwargs.get('dataDir','/minilims/data/incoming/miseq')
    sampleSheet=os.path.sep.join([dataDir,runName,'SampleSheet.csv'])
    ssin = open(sampleSheet)
    return ssin, ssin.close

def parseSampleSheet(runName, sampleSheetId, galaxyInstance, **kwargs):
    """
    Parses a SampleSheet from the file system for the given run. The run name should correspond to a subdirectory of the 'dataDir' which defaults to /minilims/data/incoming/miseq. 

    Returns a two-element tuple:
        chemistry: either 'nextera' or 'truseq'
        barcodes: dictionary from sample number (as integer) to two-element list of barcodes. Second element will be an empty string for truseq.
    """
    chemistry='scriptseq'
    barcodes={}

    if sampleSheetId is None:
        f,closeStream = loadSampleSheetFromFileSystem(runName, **kwargs)
    else:
        f,closeStream = loadSampleSheetFromGalaxy(sampleSheetId, galaxyInstance)

    try:
        # find Assay line
        line = ""
        while re.match(r'\[Reads\]',line) is None:
            line = f.next()

            # not the fastest approach, but should be clear
            if re.match(r'Assay,', line) is not None:
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
        closeStream()
    except StopIteration:
        closeStream()
        raise Exception("Unexpected end of SampleSheet!")

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
    logger.info("Found %d files in %s" % (len(files), runName))
    sampleSheetId = files.pop('sampleSheet',None)
    if sampleSheetId is None:
        logger.warn("No sample sheet!")
    (ssChemistry,barcodes) = parseSampleSheet(runName,
                                              sampleSheetId,
                                              galaxyInstance,
                                              **kwargs)
    logger.debug("Chemistry is %s for %d samples" % (ssChemistry,
                                                     len(barcodes)))
    if chemistry is None:
        chemistry=ssChemistry
    else:
        logger.info("Using user specified chemistry: %s" % chemistry)
    responses=[]
    for (sample, sampleData) in files.iteritems():
        response={'sample':sample}
        logger.debug("Checking sample %s" % (sample))
        try:
        #if True:
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
                #logger.warn("History already exists: %s" % historyName)
                raise Exception("History already exists: %s" % historyName)
            
            dsMap={}
            lanes=sampleData['lanes']
            if len(lanes)==1:
                # MiSeq data has only one lane, so we can just launch
                #  the workflow directly
                # (Also, some NextSeq runs may already have lanes merged)
                laneData = lanes.values()[0]
                for (direction, inputid) in workflowInputs.iteritems():
                    if direction not in laneData:
                        raise Exception("Missing direction in input: %s" % (direction))
                    dsMap[inputid] = {'src':'ld','id':laneData[direction]}
            else:
                # NextSeq data has 4 lanes that need to be merged first
                historyId,datasets = launchNextSeqLaneMerge(lanes,
                                                         historyName,
                                                         galaxyInstance,
                                                         sampleName=sampleName)

                # create dataset map
                for (direction, inputid) in workflowInputs.iteritems():
                    dsMap[inputid] = {'src':'hda','id':datasets[direction]}

                # This is how we tell galaxy to use an existing hitory
                historyName="hist_id=%s" % (historyId)
                logger.debug("Using same history: %s" % (historyName))

            data={'workflow_id':workflowID,
                  'history':historyName,
                  'ds_map':dsMap,
                  'replacement_params':{'name':sampleName},
                  'parameters':{primerToolID:{'chemistry':chemistry,
                                              'barcode1':barcode1,
                                              'barcode2':barcode2}}
            }
            response['data']=data
            logger.debug(data)
            response['response']=submit(apiKey, apiURL+"/workflows",data)
        except Exception as e:
            response['error']=e
        yield response

def launchNextSeqLaneMerge(laneDatasets, historyName,
                           galaxyInstance, 
                           mergeWorkflowName="NextSeq Lane Merge",
                           mergeWorkflowID=None,
                           sampleName="Merge"):
    """
    Runs the "NextSeq Lane Merge" workflow on the 4 lanes of data.
    Returns the id of the created history and Ca datasetMap to the 
    merged read files
    """

    # Locate workflow
    if mergeWorkflowID is None:
        logger.debug("Looking up workflow: %s" % mergeWorkflowName)
        workflow = getWorkflow(galaxyInstance, mergeWorkflowName)
        mergeWorkflowID = workflow[u'id']
        logger.info("Found workflow: %s" % mergeWorkflowID)

    # make datasetMap for merge workflow
    dsMap={}
    workflowInputs = getNextSeqMergeInputs(galaxyInstance, mergeWorkflowID)
    for (direction, laneInputs) in workflowInputs.iteritems():
        for lane, inputid in laneInputs.iteritems():
            dsMap[inputid] = {'src':'ld','id':laneDatasets[lane][direction]}

    # run merge workflow
    data={'workflow_id':mergeWorkflowID,
          'history':historyName,
          'ds_map':dsMap,
          'replacement_params':{'name':sampleName}
    }
    logger.debug("Launching merge workflow")
    response=submit(galaxyInstance.key, 
                    galaxyInstance.base_url+"/workflows",
                    data,
                    return_formatted=False)
    logger.debug(response)
    historyID = response[u'history']

    # find merge outputs 
    mergedDatasets={}
    mergedRE = re.compile(r'(R[12])\s+merged')
    for dataset in galaxyInstance.histories.show_history(historyID, 
                                                         contents=True):
        m = mergedRE.search(dataset[u'name'])
        if m:
            mergedDatasets[m.group(1)]=dataset[u'id']

    logger.debug("Merging in process: %s, %r" % (historyID, mergedDatasets))
    return (historyID,mergedDatasets)

inputRE = re.compile(r'Lane\s+(\d{3})\s+(R\d)')
def getNextSeqMergeInputs(galaxyInstance, workflowId):
    inputs={'R1':{},'R2':{}}
    wfInputs = galaxyInstance.workflows.show_workflow(workflowId)['inputs']
    for inputId, inputInfo in wfInputs.iteritems():
        label=inputInfo[u'label']
        lane,direction = inputRE.search(label).groups()
        inputs[direction][lane]=inputId
    return inputs

def getPrimerToolId(galaxyInstance, workflowID):
    for stepid, step in galaxyInstance.workflows.show_workflow(workflowID)[u'steps'].iteritems():
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
        raise Exception("Cannot find workflow: %s" % (name))
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
        logging.debug("SUBMIT: %s\n%r" % (url,data))
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

######
# Methods for manipulating libraries
######
def findLibrariesWithName(name, galaxyInstance, libs_cache=None, createMissing=True):
    """
    Look for a library with the given name (on the server or in the given cached list).
    Return list of libraries with matching names.
    If nothing found, create the library (turn off with createMissing=False)
    """
    # list for returned data
    libraries=[]

    # Library list from cache or from server API
    if libs_cache is None:
        libs=galaxyInstance.libraries.get_libraries()
    else:
        libs=libs_cache

    # brute force search for matching name(s)
    for lib in libs:
        if lib[u'name']==name:
            lib[u'contents']=galaxyInstance.libraries.show_library(lib[u'id'],contents=True)
            libraries.append(lib)

    # return new library if not found
    if len(libraries)==0 and createMissing:
        # Create new library
        lib=galaxyInstance.libraries.create_library(name)
        if isinstance(lib,list):
            raise Exception("create_library returned a list, this usually means your request was redirected and the library wasn't created. Try switching the root api url from http to https or vice versa")
        # Add class and contents to object to match what you get from get_libraries()
        lib[u'model_class']=u'Library'
        lib[u'contents']=galaxyInstance.libraries.show_library(lib[u'id'],contents=True)
        # Add to return list
        libraries.append(lib)

        if libs_cache is not None:
            # Add to cached list, if given
            libs.append(lib)

    return libraries

def createOrFindFolderInLibrary(folder, library, galaxyInstance):
    """
    Given a folder name, look for it in the given library. 
    First check the 'contents' entry, if missing, create folder in gi and add to contents
    Returns tuple with folder dict and boolean for whether or not it had to be created.
    The folder dict has 'library' entry added that points back to the enclosing library.
    """
    path=u'/'+folder

    # loop over items
    for item in library[u'contents']:
        if u'type' in item and item[u'type']==u'folder':
            if item[u'name']==path:
                item[u'library']=library
                return (item,False)
            if item[u'name']==u'/':
                # save root folder in case we need to create the folder
                rootFolder=item
    else:
        # Folder not found, create it
        try:
            libid=library[u'id']
            data=galaxyInstance.libraries.create_folder(libid, folder, base_folder_id=rootFolder[u'id'])
        except IncompleteRead as inst:
            # This error usually happens after the folder is
            # created. Let's see if it exists...
            for item in galaxyInstance.libraries.show_library(library[u'id'], contents=True):
                if u'type' in item and item[u'type']==u'folder' and item[u'name']==path:
                    # found it!
                    library[u'contents'].append(item)
                    item[u'library']=library
                    return (item,True)
            else:
                # nope, it was an error
                raise inst
        item=data[0]
        item[u'name']=path
        library[u'contents'].append(item)
        item[u'library']=library
        return (item,True)

def uploadFileToFolder(galaxyInstance, folder, fileName, metadata=None, file_type=None):
    """
    Check to see if this file has already been uploaded.
    Upload if it's new.
    Return data dictionary and boolean indicating whether it is new
    """
    library=folder[u'library']
    libID=library[u'id']
    runFolderID=folder[u'id']

    # search library contents for file path
    fullPath="%s/%s" % (folder[u'name'],os.path.split(fileName)[1])
    fileData=getFileDataFromLibrary(fullPath,library[u'contents'])
    if fileData is not None:
        logging.debug("Skipping file %s" % (fileName))
        return (fileData, False)

    # upload file
    logging.debug("Uploading file %s" % (fileName))
    try:
        if metadata is not None:
            upload=galaxyInstance.libraries.upload_file_from_local_path(libID,    fileName, folder_id=runFolderID, file_type=file_type, extended_metadata=metadata)
        else:
            # no metadata
            upload=galaxyInstance.libraries.upload_file_from_local_path(libID,    fileName, folder_id=runFolderID, file_type=file_type)
        fileData=upload[0]
    except IncompleteRead as inst:
        # This error usually happens after the dataset is
        # created. Let's see if it exists:
        fileData=getFileDataFromLibrary(fullPath,galaxyInstance.libraries.        show_library(library[u'id']))
        if fileData is None:
            # Nope.
            raise inst
        else:
            logging.debug("Ignored IncompleteRead in upload of %s" %
                          basecall)
    library[u'contents'].append(fileData)
    return (fileData,True)

def getFileDataFromLibrary(fullPath, libraryContents):
    """
    Look for file (and file without gz suffix) in library
    """
    paths=[fullPath]
    if fullPath[-3:]=='.gz':
        paths.append(fullPath[:-3])
    for item in libraryContents:
        if u'type' in item and item[u'type']==u'file' and item[u'name'] in paths:
            return item
    return None


