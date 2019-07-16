"""
Tools for parsing GTDB headers into a taxonomy model
"""

import os
import re
import sys
import logging
from edl.taxon import TaxNode, Taxonomy
from edl.silva import writeDumpFiles
from edl.util import treeGenerator

logger = logging.getLogger(__name__)

def generate_taxdump(fasta=None, table=None, dump="."):
    """ convert a GTDB faa file to ncbi style taxdumps """
    if fasta is not None:
        tax_file = fasta
        fmt = 'fasta'
    elif table is not None:
        tax_file = table
        fmt = 'table'
    else:
        raise Exception("Please supply 'fasta' or 'table' file")

    taxonomy = parse_gtdb_lineages(tax_file, fmt)
    dump_taxonomy(taxonomy, dump)

def generate_lineages(tax_file, fmt='fasta'):
    """ return acc,lineage tuple from file """
    if fmt == 'fasta':
        with open(tax_file) as fasta_h:
            for line in fasta_h:
                if line.startswith(">"):
                    # in GTDB headers, lineage is second chunk
                    acc, lineage = line[1:].split(None, 2)[:2]
                    yield (acc, lineage)
    elif fmt == 'table':
        with open(tax_file) as table_h:
            # skip header
            try:
                next(table_h)
            except StopIteration:
                raise Exception("Table is empty!\n" + tax_file)
            for line in table_h:
                org, _species, lineage = \
                        [x.strip() \
                         for x in line.split('\t', 4)[:3]]
                yield (org, lineage)
    else:
        raise Exception("Uknown format: " + fmt)

def parse_gtdb_lineages(tax_file, fmt='fasta'):
    """ returns taxonomy object """

    id_map = {}
    root = TaxNode('root', None, None)
    tree = {'root': root}

    logger.debug("Parsing %s", tax_file)

    # generate taxonomy tree
    for acc, lineage in generate_lineages(tax_file, fmt):
        # create TaxNode
        node = add_lineage_to_tree(lineage, tree)
        id_map[acc] = node

    logger.debug("Adding id numbers to %d nodes", len(tree))

    # assign numeric IDs
    i = 0
    for node in treeGenerator(root):
        i += 1
        node.id = i

    logger.debug("Added %d id numbers", i)

    return Taxonomy(id_map, None, None, tax_file, root)

RANK_DICT = {'d': 'domain', 'p': 'phylum', 'c': 'class',
             'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species'}

def add_lineage_to_tree(lineage, tree):
    """ parse given lineage
        create new TaxNode objects as needed
        return leaf node """
    last_node = tree['root']
    sub_lineage = []
    for taxon_string in lineage.split(';'):
        rank_char, taxon_name = taxon_string.split('__')
        rank_char = re.sub(r'^_', '', rank_char)
        sub_lineage.append(taxon_string)
        taxon = ';'.join(sub_lineage)
        try:
            last_node = tree[taxon]
        except KeyError:
            try:
                rank = RANK_DICT[rank_char]
            except KeyError:
                print(lineage)
                print(rank_char)
                exit(-1)
            new_node = TaxNode(taxon, last_node.id, rank)
            new_node.name = taxon_name
            new_node.setParent(last_node)
            tree[taxon] = new_node
            last_node = new_node
    return last_node

def dump_taxonomy(taxonomy, dump_path, map_file_name='gtdb.acc.to.taxid'):
    """ generate nodes.dmp and names.dmp """

    # Write dump files
    if not os.path.exists(dump_path):
        os.makedirs(dump_path)
    with open(os.path.sep.join((dump_path, 'nodes.dmp')), 'w') as nodes_h:
        with open(os.path.sep.join((dump_path, 'names.dmp')), 'w') as names_h:
            writeDumpFiles(taxonomy.root, nodes_h, names_h)

    # Write hit->tax mapping file
    with open(os.path.sep.join((dump_path, map_file_name)),
              'w') as acc_map_h:
        for (hitid, tax_node) in taxonomy.idMap.items():
            acc_map_h.write("%s\t%d\n" % (hitid, tax_node.id))

if __name__ == '__main__':
    """ convert a GTDB faa file to ncbi style taxdumps

        kw arguments to generate_taxdump passed as args like:
            python edl/grdb.py fasta=/path/to/x.faa dump=/path/to/dump

        for reference:
          generate_taxdump(fasta=None, table=None, dump="."):
    """

    kwargs = dict(w.split("=", 1) for w in sys.argv[1:])

    logging.basicConfig(level=logging.DEBUG)
    logger.debug("args are: %r from:\n%s", kwargs, sys.argv)

    # do the work:
    generate_taxdump(**kwargs)
