import edl.batch

def test_chunk_size_calcs():
    assert edl.batch.calculateChunkSize(100,100,10) == 10
    assert edl.batch.calculateChunkSize(1000,100,10) == 100
    assert edl.batch.calculateChunkSize(100,100,9) == 12 
    assert edl.batch.calculateChunkSize(1000,1000,9) == 112
    


def test_fasta_record_sizer():
    """fastaRecordSizer counts sequence characters, skipping the header."""
    record_lines = [">seq1 description\n", "ACGT\n", "  GG  \n"]
    assert edl.batch.fastaRecordSizer(record_lines) == 6
