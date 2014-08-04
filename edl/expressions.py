import re
from edl.kegg import koRE

fastaRE=re.compile(r'^>(\S+)(?:\s+(\S.*))?$')
accessionRE = re.compile(r'\b(?:edl|dbj|emb|gb|pdb|pir|prf|ref|sp|tpd|tpe|tpg|CDD|lcl)\|{1,2}([-\.0-9a-zA-Z_]+)(?<!\.\d)(?:\.\d+)?(?:\||$)')
nrOrgRE = re.compile(r'\[([^|\[\]]*(?:\[[^|\[\]]+\])?[^|\[\]]*)\]?\s*(?:(?=[a-z]{2,3}\|)|$)')
giRE=re.compile('^gi\|(\d+)\|.*$')
