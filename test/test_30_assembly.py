from edl import assembly
import pandas
from math import floor

fasta_file = 'test/data/HOT_100_reads.fasta'

def test_get_stats_from_contigs():
    stats = assembly.get_stats_from_contigs(fasta_file)

    # its a data frame
    assert isinstance(stats, pandas.core.frame.DataFrame)

    # 100 x 2
    assert stats.shape == (100,2)

    # test a couple arbitrary values
    assert stats.loc['000005_1741_3371','GC'] == 32.8125
    assert stats.loc['000483_1123_3166','length'] == 269


def test_get_column_stats():
    stats = assembly.get_stats_from_contigs(fasta_file)
    data = stats['GC']
    stats = assembly.get_column_stats(data)

    assert floor(stats['mean'])==39
    assert floor(stats['median'])==37


def test_contig_length_stats():
    assert isinstance(assembly.calc_stats(fasta_file,return_type='report'),str)
    data = assembly.calc_stats(fasta_file,return_type='data')
    assert data['N75'] == 142
    assert data['count'] == 100
    assert data['mean'] == 157.97

