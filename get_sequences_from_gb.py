#! /usr/bin/env python
"""
Reads one or more sequence files from stdin or the argument list and outputs
a list of sequences.
"""

from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import argparse
import sys
import traceback


def main():
    description = __doc__

    # command line options
    parser = argparse.ArgumentParser(description)
    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        metavar="OUTFILE",
        help="Write output sequences to OUTFILE.")
    parser.add_argument(
        "-f",
        "--formatIn",
        dest="formatIn",
        default="genbank",
        help="Input sequence format (see biopython docs)",
        metavar="FORMAT")
    parser.add_argument(
        "-F",
        "--formatOut",
        dest="formatOut",
        default="fasta",
        help="Output sequence format (see biopython docs)",
        metavar="FORMAT")
    parser.add_argument(
        "-c",
        "--codingSeq",
        dest="cds",
        default=False,
        action='store_true',
        help="Extract features of type CDS")
    parser.add_argument("-t", "--translate",
                        default=False, action="store_true",
                        help="Translate output sequence to amino acids")
    parser.add_argument("-r", "--refseq", default=False, action='store_true',
                        help="If input is GBK and output is FASTA, make "
                             "record id look like a proper RefSeq "
                             "entry: 'gi|XXX|ref|XXX'")
    parser.add_argument("-v", "--verbose",
                        action="store_true", dest="verbose", default=False,
                        help="Print status messages to stderr")
    parser.add_argument("-A", "--about",
                        action="store_true", dest="about", default=False,
                        help="Print description")
    parser.add_argument("args", nargs="*", metavar='FILENAME',
                        help="Files to process")

    arguments = parser.parse_args()

    if arguments.about:
        print(description)
        exit(0)

    if arguments.verbose:
        global verbose
        verbose = True

    # output
    if arguments.outfile is None:
        log("Writting %s sequences to STDOUT" % (arguments.formatOut))
        outstream = sys.stdout
    else:
        log("Writting %s sequences to %s" %
            (arguments.formatOut, arguments.outfile))
        outstream = open(arguments.outfile, 'w')

    if len(arguments.args) == 0:
        log("reading sequences from STDIN")
        instream = sys.stdin
        translateStream(
            instream,
            arguments.formatIn,
            outstream,
            arguments.formatOut,
            arguments.cds,
            arguments.translate,
            arguments.refseq)
    else:
        for name in arguments.args:
            log("reading %s sequences from %s" % (arguments.formatIn, name))
            instream = open(name, 'rU')
            try:
                translateStream(
                    instream,
                    arguments.formatIn,
                    outstream,
                    arguments.formatOut,
                    arguments.cds,
                    arguments.translate,
                    arguments.refseq)
            except Exception:
                warn("Exception parsing %s:\n-----\n" % (name))
                traceback.print_exc(file=sys.stderr)
            instream.close()


#############
# Functions #
#############
verbose = False


def log(msg):
    if verbose:
        sys.stderr.write(msg)
        sys.stderr.write("\n")


def die(msg):
    sys.stderr.write("%s\n" % msg)
    sys.exit()


def warn(msg):
    sys.stderr.write("WARNING: %s\n" % (msg))


def translateStream(
        instream,
        inf,
        outstream,
        outf,
        cds,
        translate,
        makeRefSeq):
    log("translating records from %s to %s (%s,%s)" %
        (inf, outf, cds, translate))
    records = SeqIO.parse(instream, inf)
    for record in records:
        if cds:
            # get coding sequences if requested
            translations = getCodingSequences(record, makeRefSeq)
        else:
            if makeRefSeq and 'gi' in record.annotations:
                record.id = "gi|%s|ref|%s|" % \
                    (record.annotations['gi'], record.id)
            translations = (record,)

        if translations is None or len(translations) == 0:
            warn("record %s has no features!" % (record.id))
            continue

        # change alphabet
        if translate:
            translations = [translateRecord(t) for t in translations]

        # write in new format
        for t in translations:
            log("writing %s" % (str(t)))
            SeqIO.write([t], outstream, outf)


def getCodingSequences(record, makeRefSeq):
    try:
        org = " [%s]" % (record.annotations['organism'])
    except KeyError:
        org = None

    seqs = []
    gene_count = 0
    for f in record.features:
        if f.type == 'CDS':
            gene_count += 1
            if makeRefSeq:
                if 'protein_id' in f.qualifiers:
                    acc = f.qualifiers['protein_id'][0]
                else:
                    continue
                if 'db_xref' in f.qualifiers:
                    for ref in f.qualifiers['db_xref']:
                        if ref[0:2] == 'GI':
                            gi = ref[3:]
                            break
                    else:
                        continue
                if 'translation' in f.qualifiers:
                    translation = f.qualifiers['translation'][0]
                else:
                    continue
                seq = Seq(translation, IUPAC.protein)
                r = SeqRecord.SeqRecord(seq,
                                        id="gi|%s|acc|%s|" % (gi, acc),
                                        name=acc)
            else:
                seq = f.extract(record.seq)
                r = SeqRecord.SeqRecord(seq)
                foundName = True
                if 'protein_id' in f.qualifiers:
                    r.id = f.qualifiers['protein_id'][0]
                    r.name = r.id
                elif 'db_xref' in f.qualifiers:
                    for ref in f.qualifiers['db_xref']:
                        if ref[0:2] == 'GI':
                            r.id = ref
                            r.name = ref
                            break
                    else:
                        foundName = False
                else:
                    foundName = False

                if not foundName:
                    for q in ('locus_tag', 'name', 'id'):
                        if q in f.qualifiers:
                            r.id = f.qualifiers[q][0]
                            r.name = r.id
                            break
                    else:
                        r.name = '%s_GENE_%s' % (record.id, gene_count)
                        r.id = r.name

            desc = r.name
            for q in ('product', 'gene', 'note'):
                if q in f.qualifiers:
                    desc = f.qualifiers[q][0]
                    break

            if org is not None:
                desc += org

            r.description = desc

            if 'db_xrefs' in f.qualifiers:
                r.dbxrefs = f.qualifiers['db_xrefs']

            log("created:\n%s\n from:\n%s" % (repr(r), repr(f)))
            seqs.append(r)

    return seqs


def translateRecord(record):
    record.seq = record.seq.translate()
    return record


if __name__ == '__main__':
    main()
