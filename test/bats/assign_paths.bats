#!/usr/bin/env bats

@test "Testing assign_paths compilation" {
    run ./assign_paths.py -h
    [ "$status" = 0 ]
}

@test "Testing assign_paths squash" {
    run ./assign_paths.py -H test/data/kobrite/ko00001.keg -l ko -l 3  -m test/data/ko.map.partial -v -s test/data/contig.CDSs.faa.vs.KEGG.lastal -o test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash
    [ "$status" = 0 ]
}
@test "Checking assign_paths squash output" {
    run diff test/data/contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash
    [ "$status" = 0 ]
}

@test "Testing assign_paths split" {
    run ./assign_paths.py -H test/data/kobrite/ko00001.keg -l ko -l 3  -m test/data/ko.map.partial -v test/data/contig.CDSs.faa.vs.KEGG.lastal -o test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.split
    [ "$status" = 0 ]
}
@test "Checking assign_paths split output" {
    run diff test/data/contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.split test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.split
    [ "$status" = 0 ]
}
