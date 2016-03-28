#!/usr/bin/env bats

@test "Testing translateColumn" {
    run ./translateColumn.py test/data/sample.?.kegg.blastx.b50.m8 -m test/data/ko.map.partial -f NA -c 3 -o .addko.tst.
    [ "$status" = 0 ]
}
@test "Checking translateColumn output 1" {
    run diff test/data/sample.1.kegg.blastx.b50.m8.addko test/data/sample.1.kegg.blastx.b50.m8.addko.tst.
    [ "$status" = 0 ]
}
@test "Checking translateColumn output 2" {
    run diff test/data/sample.2.kegg.blastx.b50.m8.addko test/data/sample.2.kegg.blastx.b50.m8.addko.tst.
    [ "$status" = 0 ]
}
@test "Checking translateColumn output 3" {
    run diff test/data/sample.3.kegg.blastx.b50.m8.addko test/data/sample.3.kegg.blastx.b50.m8.addko.tst.
    [ "$status" = 0 ]
}

