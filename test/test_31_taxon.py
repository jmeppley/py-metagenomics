from edl import taxon


def test_compare_ranks():
    """compareRanks returns cmp-style -1/0/1 from the recorded comparisons."""
    taxon.comparisons.clear()
    taxon.comparisons[('genus', 'species')] = 3
    try:
        assert taxon.compareRanks('genus', 'species') == 1
        assert taxon.compareRanks('species', 'genus') == -1
        assert taxon.compareRanks('order', 'class') == 0
    finally:
        taxon.comparisons.clear()
