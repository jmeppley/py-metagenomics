#!/usr/bin/env bats

@test "Testing assign_paths compilation" {
    run ./assign_paths.py -h
    [ "$status" = 0 ]
}

@test "Testing assign_paths squash" {
    run ./assign_paths.py -H test/data/kobrite/ko00001.keg -l ko -l 3  -m test/data/ko.map.partial -v -S test/data/contig.CDSs.faa.vs.KEGG.lastal -o test/data/.tst.contig.CDSs.faa.vs.KEGG.lastal.assignments.ko.3.squash
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

@test "Testing assign_paths consensus" {
	run ./assign_paths.py -f hmmsearchdom --sort score -C consensus -F 5 -r -p hitid test/data/reads.non-rRNA.vs.COG.tbl.reads.head50 -o test/data/.tst.reads.non-rRNA.vs.COG.tbl.reads.head50.COGs.consensus
	[ "$status" = 0 ]
}
@test "Checking assign_paths consensus" {
    run diff <(sort test/data/reads.non-rRNA.vs.COG.tbl.reads.head50.COGs.consensus) <(sort test/data/.tst.reads.non-rRNA.vs.COG.tbl.reads.head50.COGs.consensus)
    [ "$status" = 0 ]
}

@test "Testing assign_paths most" {
	run ./assign_paths.py -f hmmsearchdom --sort score -C most -F 5 -r -p hitid test/data/reads.non-rRNA.vs.COG.tbl.reads.head50 -o test/data/.tst.reads.non-rRNA.vs.COG.tbl.reads.head50.COGs.most
	[ "$status" = 0 ]
}
@test "Checking assign_paths most" {
    run diff <(sort test/data/reads.non-rRNA.vs.COG.tbl.reads.head50.COGs.most) <(sort test/data/.tst.reads.non-rRNA.vs.COG.tbl.reads.head50.COGs.most)
    [ "$status" = 0 ]
}

@test "Checking assign_paths with sort by evalue" {
    run ./assign_paths.py -f hmmsearchdom -s evalue -C first test/data/genes.v.tigr.tbl -o test/data/.tst.genes.v.tigr.tbl.paths
    [ "$status" = 0 ]
}
@test "Checking assign_paths output from sort by evalue" {
    run diff test/data/genes.v.tigr.tbl.paths test/data/.tst.genes.v.tigr.tbl.paths
    [ "$status" = 0 ]
}
