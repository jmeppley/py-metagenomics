#!/usr/bin/env bats

@test "Testing screen_list" {
    run ./screen_list.py test/data/HOT_100_reads.fasta -l test/data/HOT_100_reads.8.names -k -o test/data/.tst.HOT_100_reads.8.fasta
    [ "$status" = 0 ]
    run diff test/data/HOT_100_reads.8.fasta test/data/.tst.HOT_100_reads.8.fasta
    [ "$status" = 0 ]
}
