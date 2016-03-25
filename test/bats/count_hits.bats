#!/usr/bin/env bats

@test "Testing count_hits.py" {
    run ./countHits.py -i test/data/read.assignments.sample -o test/data/.tst.read.assignments.sample.counts -H 1 -a first -v
    [ "$status" = 0 ]
    run diff test/data/read.assignments.sample.counts test/data/.tst/read.assignments.sample.counts
    [ "$status" = 0 ]
}
