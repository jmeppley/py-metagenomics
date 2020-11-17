#!/usr/bin/env bats

@test "Testing basic compilation" {
    run ./fragment_records.py -h
    [ "$status" = 0 ]
}

@test "C10 and fasta parsing" {
    run ./fragment_records.py -i test/data/HOT_100_reads.fasta -T fasta -C 10 -o test/data/.tst.frag.C10.fasta -Z 2
    [ "$status" = 0 ]
    [ -e test/data/.tst.frag.C10.10.fasta ]
    [ -e test/data/.tst.frag.C10.01.fasta ]
    run grep -c "^>" test/data/.tst.frag.C10.10.fasta
    [ "${lines[0]}" = "10" ]
    run grep -c "^>" test/data/.tst.frag.C10.01.fasta
    [ "${lines[0]}" = "10" ]
}

@test "N10 and fasta extension recognition" {
    run ./fragment_records.py -i test/data/HOT_100_reads.fasta -N 10 -o test/data/.tst.frag.N10.fasta
    [ "$status" = 0 ]
    [ -e test/data/.tst.frag.N10.10.fasta ]
    [ -e test/data/.tst.frag.N10.01.fasta ]
    run grep -c "^>" test/data/.tst.frag.N10.10.fasta
    [ "${lines[0]}" = "10" ]
    run grep -c "^>" test/data/.tst.frag.N10.01.fasta
    [ "${lines[0]}" = "10" ]
}

@test "C5 + N5" {
    run rm test/data/.tst.frag.N5.C5.?.fasta
    run ./fragment_records.py -i test/data/HOT_100_reads.fasta -C 5 -N 5 -o test/data/.tst.frag.N5.C5.fasta
    [ "$status" = 0 ]
    [ -e test/data/.tst.frag.N5.C5.5.fasta ]
    [ -e test/data/.tst.frag.N5.C5.1.fasta ]
    [ ! -e test/data/.tst.frag.N5.C5.6.fasta ]
    run grep -c "^>" test/data/.tst.frag.N5.C5.1.fasta
    [ "$status" = 0 ]
    [ "${lines[0]}" = "5" ]
    run grep -c "^>" test/data/.tst.frag.N5.C5.5.fasta
    [ "$status" = 0 ]
    [ "${lines[0]}" = "5" ]
}

@test "evening out fragments" {
    run ./fragment_records.py -i test/data/HOT_100_reads.fasta -C 40 -o test/data/.tst.frag.C40E.fasta -v -v -E -Z 1
    [ "$status" = 0 ]
    [ -e test/data/.tst.frag.C40E.1.fasta ]
    run grep -c "^>" test/data/.tst.frag.C40E.1.fasta
    [ "${lines[0]}" = "34" ]
}

