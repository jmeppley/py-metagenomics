#!/usr/bin/env bats

@test "Testing generating sunburst JSON" {
    run ./sunburstFromMeganCSV.py -n test/data test/data/taxcounts.csv --figsize=6,6 -c Bacteria=r,Proteobacteria=g,Cyanobacteria=c -C 0.01 -o test/data/.tst.taxcounts.sunburst.json -J
    [ "$status" = 0 ]
    run diff -w <(head test/data/taxcounts.sunburst.json) <(head test/data/.tst.taxcounts.sunburst.json)
    [ "$status" = 0 ]
}

@test "Testing generating sunburst PDF" {
    run ./sunburstFromMeganCSV.py -n test/data test/data/taxcounts.csv --figsize=6,6 -c Bacteria=r,Proteobacteria=g,Cyanobacteria=c -C 0.01 -o test/data/.tst.taxcounts.sunburst.pdf
    [ "$status" = 0 ]
}
