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

GTDB='gtdb'
GTDBTAB='gtdb_table'
PHYLODB='phylodb'

def generate_taxdump(fasta=None, table=None, dump=".", **kwargs):
    """ convert a GTDB faa file to ncbi style taxdumps """
    if fasta is not None:
        tax_file = fasta
        fmt = 'fasta'
    elif table is not None:
        tax_file = table
        fmt = 'table'
    else:
        raise Exception("Please supply 'fasta' or 'table' file")

    tax_args = {k:v for k,v in kwargs.items() if k in ['style']}
    taxonomy = parse_lineages(tax_file, fmt, **tax_args)

    dump_args = {k:v for k,v in kwargs.items() if k in ['map_file_name']}
    dump_taxonomy(taxonomy, dump, **dump_args)

def generate_gtdb_lineages_from_table(tax_file):
    """ return acc,lineage tuple from file """
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

def generate_gtdb_lineages(fasta_file):
    """ return acc,lineage tuple from file """
    with open(fasta_file) as fasta_h:
        for line in fasta_h:
            if line.startswith(">"):
                # in GTDB headers, lineage is second chunk
                acc, lineage = line[1:].split(None, 2)[:2]
                yield (acc, lineage)

def generate_phylodb_lineages(fasta_file):
    """ return acc,lineage tuple from file """
    with open(fasta_file) as fasta_h:
        for line in fasta_h:
            if line.startswith(">"):
                # in GTDB headers, lineage is second chunk
                acc, lineage = line[1:].split("\t", 2)[:2]
                yield (acc, lineage)

def parse_lineages(tax_file, fmt='fasta', style=GTDB):
    """ returns taxonomy object """

    id_map = {}
    root = TaxNode('root', None, None)
    tree = {'root': root}

    logger.debug("Parsing %s", tax_file)

    if style == GTDB:
        add_lineage_to_tree = add_gtdb_lineage_to_tree
        generate_lineages = generate_gtdb_lineages
    else:
        add_lineage_to_tree = add_phylodb_lineage_to_tree
        generate_lineages = generate_phylodb_lineages

    # generate taxonomy tree
    for acc, lineage in generate_lineages(tax_file):
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

RANK_LIST = ['domain', 'phylum', 'class',
             'order', 'family', 'genus', 'species']
def add_phylodb_lineage_to_tree(lineage, tree):
    """ parse given lineage
        create new TaxNode objects as needed
        assumes that there are 7 elements in lineage, one for each rank
        return leaf node """
    last_node = tree['root']
    sub_lineage = []
    if lineage.startswith('Euk'):
        # There is an extra entr in the PhyloDB Euk lineages
        ranks = [RANK_LIST[0], None,] + RANK_LIST[1:]
    else:
        ranks = RANK_LIST
    for rank, taxon_string in zip(ranks, lineage.split(';')):
        sub_lineage.append(taxon_string)
        taxon = ';'.join(sub_lineage)
        try:
            last_node = tree[taxon]
        except KeyError:
            new_node = TaxNode(taxon, last_node.id, rank)
            new_node.name = taxon_string
            new_node.setParent(last_node)
            tree[taxon] = new_node
            last_node = new_node
    return last_node

RANK_DICT = {'d': 'domain', 'p': 'phylum', 'c': 'class',
             'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species'}

def add_gtdb_lineage_to_tree(lineage, tree):
    """ parse given lineage
        create new TaxNode objects as needed
        assumes lineage names atart with x__ where x is a rank abbreviation
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
    if map_file_name is None:
        return
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
