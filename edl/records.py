from edl.util import parse_list_to_set

def recordIterator(stream, separatorRE, idRE=None):
    """
    Given:

   na file-like object (any iterator over strings)
        1 or 2 regular expressions that define record boundaries and identifiers
    Return:
        an iterator over records that returns a tuple of (id, [recordLines])
    If only a separator given, it is assumed to match the record id
    """
    recordId=None
    recordLines=[]
    for line in stream:
        m=separatorRE.search(line)
        if m:
            # is there a previous record?
            if recordId is not None:
                yield (recordId, recordLines)
                recordId=None
            recordLines=[line,]

            if idRE is None:
                recordId=m.group(1)

            continue

        recordLines.append(line)

        if idRE is not None:
            m=idRE.search(line)
            if m:
                recordId=m.group(1)

    if recordId is not None:
        yield (recordId, recordLines)

def screenRecords(stream, separatorRE, idRE=None, keep=False, screen_set=None, screenFile=None):
    """
    uses recordIterator(strean, separatorRE, idRE) to parse input into records
    uses screen_set (can be read from screenFile) to identify records 
    identified records are kept or skipped based on the value of keep
    """

    if screen_set is None:
        if screenFile is None:
            raise Exception("Please supply a hash(Python map) or file of record keys")
        else:
            screen_set=parse_list_to_set(screenFile)

    for (recordId, recordLines) in recordIterator(stream,separatorRE,idRE=idRE):
        screened = recordId in screen_set
        if screened == keep:
            for line in recordLines:
                yield line

