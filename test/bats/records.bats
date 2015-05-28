#!/usr/bin/env bats

@test "Testing sample_records.py -h" {
    run python sample_records.py -h
    [ "$status" = 0 ]
}

@test "Testing sample_records.py -s 10" {
    run bash -c "python sample_records.py -s 10 test/data/HOT_100_reads.fasta | grep -c '>'"
    [ "$status" = 0 ]
    [ "$output" = "10" ]
}

