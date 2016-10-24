#!/usr/bin/env bats

@test "Testing get_seqs_from_gb: faa" {
    run ./get_sequences_from_gb.py test/data/test.gbk -c -o test/data/.tst.test.gbk.faa -f genbank -t
    [ "$status" = 0 ]
    run diff test/data/test.gbk.faa test/data/.tst.test.gbk.faa
    [ "${#lines[@]}" = 0 ]
}

@test "Testing get_seqs_from_gb: fna" {
    run ./get_sequences_from_gb.py test/data/test.gbk -c -o test/data/.tst.test.gbk.fna -f genbank
    [ "$status" = 0 ]
    run diff test/data/test.gbk.fna test/data/.tst.test.gbk.fna
    [ "$status" = 0 ]
}
