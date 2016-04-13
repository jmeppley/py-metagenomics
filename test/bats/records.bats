#!/usr/bin/env bats

@test "Testing sample_records.py -h" {
    run python sample_records.py -h
    [ "$status" = 0 ]
}

@test "Testing sample_records.py -s 10 with count" {
    run bash -c "python sample_records.py -s 10 test/data/HOT_100_reads.fasta | grep -c '>'"
    [ "$status" = 0 ]
    [ "$output" = "10" ]
}

@test "Testing sample_records.py -s 10 with uniq" {
    run bash -c "python sample_records.py -s 10 test/data/HOT_100_reads.fasta | sort | uniq | grep -c '>'"
    [ "$status" = 0 ]
    [ "$output" = "10" ]
}

