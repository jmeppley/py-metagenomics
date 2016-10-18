#!/usr/bin/env bats

@test "Testing translate_column compilation" {
    run ./translate_column.py -h
    [ "$status" = 0 ]
}

@test "Testing translate_column" {
    run ./translate_column.py test/data/sample.?.kegg.blastx.b50.m8 -m test/data/ko.map.partial -f NA -c 3 -o .addko.tst.
    [ "$status" = 0 ]
}
@test "Checking translate_column output 1" {
    run diff test/data/sample.1.kegg.blastx.b50.m8.addko test/data/sample.1.kegg.blastx.b50.m8.addko.tst.
    [ "$status" = 0 ]
}
@test "Checking translate_column output 2" {
    run diff test/data/sample.2.kegg.blastx.b50.m8.addko test/data/sample.2.kegg.blastx.b50.m8.addko.tst.
    [ "$status" = 0 ]
}
@test "Checking translate_column output 3" {
    run diff test/data/sample.3.kegg.blastx.b50.m8.addko test/data/sample.3.kegg.blastx.b50.m8.addko.tst.
    [ "$status" = 0 ]
}

