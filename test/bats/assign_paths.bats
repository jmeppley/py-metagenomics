#!/usr/bin/env bats

@test "Testing assignPaths squash" {
    run ./assignPaths.py -H test/data/kobrite/ko00001.keg -l ko -l 3  -m test/data/ko.map.partial -v -s test/data/contig.CDSs.faa.vs.KEGG.lastal -o test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash
    [ "$status" = 0 ]
    run diff <(head test/data/contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash) <(head test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash)
    [ "$status" = 0 ]
    run diff <(tail test/data/contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash) <(tail test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash)
    [ "$status" = 0 ]
}

@test "Testing assignPaths split" {
    run ./assignPaths.py -H test/data/kobrite/ko00001.keg -l ko -l 3  -m test/data/ko.map.partial -v test/data/contig.CDSs.faa.vs.KEGG.lastal -o test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.split
    [ "$status" = 0 ]
    run diff <(head test/data/contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.split) <(head test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.split)
    [ "$status" = 0 ]
    run diff <(tail test/data/contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.split) <(tail test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.split)
    [ "$status" = 0 ]
}
