#!/usr/bin/env bats

@test "Testing get_seqs_from_m8" {
    run ./get_sequences_from_m8.py -i test/data/HOT_100_reads.fasta test/data/HOT_100_reads.fasta.v.SAGs.f0.lastx -o test/data/.tst.HOT_100_reads.fasta.v.SAGs.f0.lastx.faa -F 0 -f last -t 
    [ "$status" = 0 ]
    run diff test/data/HOT_100_reads.fasta.v.SAGs.f0.lastx.faa test/data/.tst.HOT_100_reads.fasta.v.SAGs.f0.lastx.faa
    [ "$status" = 0 ]
}
