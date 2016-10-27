#!/usr/bin/env python
"""
Takes in a set of gff files in two lists. Features (like tRNA and rRNA)
and then features that we don't want to overlap with them (like CDS
from prodigal)

Outputs 3 files: a master gff, a matching fna of fastq seqs, and a faa of
amino acid seqs of just the CDSs that don't overlap RNA.
"""

import sys
import re
import os
import glob
import argparse
import logging
from edl.blastm8 import GFF, generate_hits
from edl.util import add_universal_arguments, setup_logging
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("contigs_fasta",
                        help="The contigs as a fasta files")
    parser.add_argument("output_file_prefix",
                        help="Three files are created with this prefix \
                              and one of hte suffixes: '.gff','.fna','.faa'")
    parser.add_argument("-r", "--rna-gff", action='append', default=[],
                        help="GFF file containing RNA annotations. \
                              May be specified multiple times for \
                              multiple files")
    parser.add_argument("-c", "--cds-gff", action='append', default=[],
                        help="GFF file containing RNA annotations. \
                              May be specified multiple times for \
                              multiple files")

    add_universal_arguments(parser)
    args = parser.parse_args()
    setup_logging(args)

    logging.info("Generating annotations for {contigs_fasta}".
                 format(**vars(args)))

    merge_gffs(args.rna_gff,
               args.cds_gff,
               args.contigs_fasta,
               args.output_file_prefix)


def merge_gffs(rna_gff_files, cds_gff_files, contigs_fasta_file,
               output_file_prefix):

    # parse RNA GFF files
    rna_hits = {f: get_gff_hits(f) for f in rna_gff_files}

    # parse CDS files
    cds_hits = {}
    for cds_gff_file in cds_gff_files:
        for contig, hits in generate_hits(cds_gff_file, format=GFF):
            # get regions with rRNAs for this contig
            rna_regions = get_rna_regions(rna_hits, contig)

            # collect CDSs that don't overlap
            cds_hits.setdefault(contig, []).extend(
                [h for h in hits if h.checkForOverlap(rna_regions)[1] is None])

    # the source data
    hit_list_dicts = list(rna_hits.values())
    hit_list_dicts.append(cds_hits)

    # output files
    with open(output_file_prefix + ".gff", 'w') as GFFOUT:
        with open(output_file_prefix + ".fna", 'w') as FNAOUT:
            with open(output_file_prefix + ".faa", 'w') as FAAOUT:
                write_annotations_to_files(hit_list_dicts,
                                           contigs_fasta_file,
                                           GFFOUT,
                                           FNAOUT,
                                           FAAOUT)


def write_annotations_to_files(hit_list_dicts, contigs_fasta_file,
                               gff_output_handle,
                               fna_output_handle,
                               faa_output_handle):
    # #### Loop over contigs
    # Merge and sort hits by lowest position index

    # loop over contigs
    for contig in SeqIO.parse(contigs_fasta_file, 'fasta'):

        # collect hits for this contig
        contig_hits = []
        for hit_list_dict in hit_list_dicts:
            hits = hit_list_dict.get(contig.id, [])
            for hit in hits:
                contig_hits.append(hit)

        # sort by start position
        contig_hits.sort(key=lambda h: min(h.qstart, h.qend))

        # loop over hits
        for i, hit in enumerate(contig_hits):
            n = i + 1
            # output GFF as it came in
            gff_output_handle.write(hit.line)

            # get naming informatoion from gff line and contig descritpion
            length, coverage = re.search(r'length_(\d+)_cov_([0-9.]+)',
                                         contig.description).groups()
            source, feature_type, start, end, score, strand = \
                hit.line.split('\t')[1:7]
            # name gene with contig name and index
            gene_name = contig.id + "_{n}".format(n=n)
            # put everything else in the description
            gene_desc =\
                ("source={source};type={feature_type};score={score};" +
                 "strand={strand};start={start};end={end};") \
                .format(source=source, start=start, end=end,
                        feature_type=feature_type, score=score,
                        strand=strand) + \
                hit.hitDesc + \
                "contig_length={length};contig_cov={coverage}" \
                .format(length=length, coverage=coverage)

            if hit.qstart <= hit.qend:
                subsequence = contig.seq[hit.qstart - 1:hit.qend]
            else:
                subsequence = contig.seq[
                    hit.qend - 1:hit.qstart].reverse_complement()
            if strand == '-':
                subsequence = subsequence.reverse_complement()
            fna_output_handle.write(">{seqid}\t{desc}\n{seq}\n"
                                    .format(seqid=gene_name,
                                            desc=gene_desc,
                                            seq=str(subsequence)
                                            ))

            if feature_type == 'CDS':
                aa_string = str(subsequence.translate(table=11))
                faa_output_handle.write(">{seqid}\t{desc}\n{seq}\n"
                                        .format(seqid=gene_name,
                                                desc=gene_desc,
                                                seq=aa_string,
                                                ))


def get_rna_regions(rna_hit_dicts, contig):
    """
    merges all the rna hits for a contig into a list of (start,end) tuples
    """
    rna_tuples = []
    for rna_hit_dict in rna_hit_dicts.values():
        rna_regions = get_regions(rna_hit_dict.get(contig, []))
        rna_tuples = merge_regions(rna_tuples, rna_regions)
    return rna_tuples


def get_gff_hits(hit_table_gff, **filter_args):
    return {c: list(h) for c, h in generate_hits(hit_table_gff,
                                                 format=GFF,
                                                 **filter_args)}


def get_regions(hits):
    """return list of sorted start end tuples"""
    regions = []
    for hit in hits:
        # for our purposes, start is the lowest overlapped position, ignoring
        # driection/strand
        start = min(hit.qstart, hit.qend)
        end = max(hit.qstart, hit.qend)
        # add to list
        regions.append((start, end))
    return regions


def merge_regions(regions, new_regions):
    """
    Simply combine lists for now. We way want to be more intelligent about
    this in the future
    """
    regions.extend(new_regions)
    return regions

if __name__ == "__main__":
    main()
