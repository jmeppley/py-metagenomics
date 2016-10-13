#!/usr/bin/env bats

@test "Testing count_hits.py compilation" {
    run ./count_hits.py -h
    [ "$status" = 0 ]
}

@test "Testing count_hits.py" {
    run ./count_hits.py -i test/data/read.assignments.sample -o test/data/.tst.read.assignments.sample.counts -H 1 -a first -v
    [ "$status" = 0 ]
}
@test "Checking output from count_hits.py" {
    run diff test/data/read.assignments.sample.counts test/data/.tst.read.assignments.sample.counts
    [ "$status" = 0 ]
}
