#!/usr/bin/env bats

@test "Testing add_pathway compilation" {
    run ./add_kegg_pathways.py -h
    [ "$status" = 0 ]
}

@test "Testing add_kegg_pathways" {
    run ./add_kegg_pathways.py -k test/data/kobrite/ko00001.keg -l PATHWAY test/data/contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash -c 2 -o test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash.add_path
    [ "$status" = 0 ]
}
@test "Checking add_kegg_pathways output" {
    run diff test/data/contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash.add_path test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash.add_path
    [ "$status" = 0 ]
}
@test "Testing add_kegg_pathways level 1" {
    run ./add_kegg_pathways.py -k test/data/kobrite/ko00001.keg -l 1 test/data/contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash -c 2 -C -1 -o test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash.add_k1
    [ "$status" = 0 ]
}
@test "Checking add_kegg_pathways level 1 output" {
    run diff test/data/contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash.add_k1 test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash.add_k1
    [ "$status" = 0 ]
}
@test "Testing add_kegg_pathways" {
    run ./add_kegg_pathways.py -k test/data/kobrite/ko00001.keg -L test/data/contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash -c 2 -C -1 -o test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash.add_path_long
    [ "$status" = 0 ]
}
@test "Checking add_kegg_pathways output" {
    run diff test/data/contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash.add_path_long test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash.add_path_long
    [ "$status" = 0 ]
}
