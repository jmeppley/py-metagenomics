"""
Load the NCBI taxdump and refseq taxid map and write a pruned version of the taxdump that will be quicker to load

Also, check for missing taxids (from mismatched taxdump and refseq metadata)
and write corrected map.
"""
from Bio import Entrez
import re

def main(input, output, params):
    from edl.taxon import readTaxonomy
    from edl.util import parseMapFile
    from edl.silva import writeDumpFiles
    import os

    refseq_tax_table = str(input.taxids)
    refseq_tax_dir = os.path.dirname(str(input.names))

    taxonomy = readTaxonomy(refseq_tax_dir)
    taxid_map = parseMapFile(refseq_tax_table, valueType=int)

    # find missing taxids and get corected version
    orig_taxid_set = set(taxid_map.values())
    missing_taxids = orig_taxid_set.difference(taxonomy.idMap)
    new_taxids = {t:get_best_taxid(t, taxonomy)
				  for t in missing_taxids}
    with open(str(output.taxids), 'wt') as TAXOUT:
        for acc, taxid in taxid_map.items():
            fixed_id = new_taxids.get(taxid, taxid)
            TAXOUT.write(f'{acc}\t{fixed_id}\n')

    # prune tree and save new version
    taxid_set = \
        orig_taxid_set \
            .difference(new_taxids.keys()) \
            .union(new_taxids.values())
    prune(taxonomy.root, taxid_set)
    with open(str(output.names), 'wt') as NAMES, open(str(output.nodes), 'wt') as NODES:
        writeDumpFiles(taxonomy.root, NODES, NAMES)

    
def get_best_taxid(taxid, taxonomy):
    # first see if there is an AKA taxid in our taxonomy
    handle = Entrez.efetch(db='taxonomy', id=taxid)
    r = handle.read()
    m = re.search(r"AkaTax.+\>(\d+)\<\/Item", r)
    if m:
        aka_taxid = int(m.group(1))
        if aka_taxid != 0 and aka_taxid in taxonomy.idMap:
            return aka_taxid
    
    # next, go through the lineage and return the first one in our taxonomy
    handle = Entrez.efetch(db='taxonomy', id=taxid)
    r = handle.read()
    for candid in re.findall(r'<TaxId>(\d+)</TaxId>', r):
        candid = int(candid)
        if candid != taxid:
            if candid in taxonomy.idMap:
                best_alt = candid
    return best_alt

def prune(node, keep_set):
    keep = []
    for child in node.children:
        if child.id == node.id:
            continue
        if prune(child, keep_set):
            keep.append(child)
    if len(keep) > 0 or node.id in keep_set:
        node.children = keep
        return True
    return False

if __name__ == "__main__":
    main(snakemake.input, snakemake.output, snakemake.params)

