from edl.expressions import *

def test_accession_re():
    with open('test/data/sample.1.blastx.b50.m8') as F:
        try:
            for line in F:
                acc = accessionRE.search(line).group(1)
        except AttributeError:
            # There should have been a match in every line of this file
            assert False


    test_data = {
            'ref|YP_002498923.1|': 'YP_002498923',
            'gi|109900248|ref|YP_663503.1|': 'YP_663503',
            }
    for data, acc in test_data.items():
        new_acc = accessionRE.search(data).group(1)
        assert acc == new_acc


def test_fasta_re():
    file_data = {
            'test/data/test.gbk.faa': 3941,
            'test/data/test.gbk.fna': 3941,
            'test/data/createPrimerNextera.fasta': 2,
            'test/data/createPrimerTruseq.fasta': 4,
            'test/data/HOT_100_reads.8.fasta': 8,
            'test/data/HOT_100_reads.fasta': 100,
    }
    count=0
    for file_name, expected_count in file_data.items():
        new_count = _count_re_hits(file_name, fastaRE)
        assert new_count == expected_count
        count+=1
    assert count == len(file_data)
    assert count == 6


def _count_re_hits(file_name, regex):
    count=0
    with open(file_name) as INF:
        for line in INF:
            if regex.search(line):
                count+=1
    return count


