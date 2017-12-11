import re
from edl.kegg import koRE

fastaRE = re.compile(r'^>(\S+)(?:\s+(\S.*))?$')
accessionRE = re.compile(
    r"""(?:^|\b             # start at first chr or a word-break
        [A-Za-z]{2,3}       # followed by one of the seq codes
        \|{1,2})            # and one or two vertical bars
        ([-\.0-9a-zA-Z_]+)  # the accession is word chars plus . and - and _
        (?<!\.\d\d)(?<!\.\d\d\d)(?<!\.\d)
                            # dont include version suffix (eg .1 or .11)
        (?:\.\d+)?          # allow for version suffix outside of match
        (?:\||$)            # must end with bar (|) or nothing ($)
        (?=(?:\s|$))          # either at end or before whitespace""",
    re.X)
nrOrgRE = re.compile(
    r'\[([^|\[\]]*(?:\[[^|\[\]]+\])?[^|\[\]]*)\]?\s*(?:(?=[a-z]{2,3}\|)|$)')
giRE = re.compile('^gi\|(\d+)\|.*$')
pfamRE = re.compile(r'\b((?:COG|PF|TIGR)\d+)\.\d+\b')
