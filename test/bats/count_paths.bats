#!/usr/bin/env bats

@test "Testing count_paths with acc map" {
    run python ./count_paths.py test/data/sample.?.blastx.b50.m8 -o test/data/.tst.keggCounts.F0.acc.c0.tophit -c 0 -m test/data/acc.to.ko.protein.plus.filtered -C tophit -r -F 0 -p accs
    [ "$status" = 0 ]
}
@test "Checking output from count_paths with acc map" {
    run diff test/data/keggCounts.F0.acc.c0.tophit test/data/.tst.keggCounts.F0.acc.c0.tophit
    [ "$status" = 0 ]
}

@test "Testing count_paths with kegg map and multiple levels and a cutoff" {
    run python ./count_paths.py test/data/sample.?.kegg.blastx.b50.m8 -H test/data/kobrite/ko00001.keg -m test/data/ko.map.partial -c 0.15 -l ko -o test/data/.tst.small_files.kegg.c.15.out -l 2 -l 3
    [ "$status" = 0 ]
}
@test "Checking level 3 output from count_paths with kegg map" {
    run diff <(sort test/data/small_files.kegg.c.15.out.3) <(sort test/data/.tst.small_files.kegg.c.15.out.3)
    [ "$status" = 0 ]
}
@test "Checking level 2 output from count_paths with kegg map" {
    run diff test/data/small_files.kegg.c.15.out.2 test/data/.tst.small_files.kegg.c.15.out.2
    [ "$status" = 0 ]
}
@test "Checking level ko output from count_paths with kegg map" {
    run diff test/data/small_files.kegg.c.15.out.ko test/data/.tst.small_files.kegg.c.15.out.ko
    [ "$status" = 0 ]
}

@test "Testing count_paths with kegg map and multiple levels with no cutoff" {
    run python count_paths.py test/data/sample.?.kegg.blastx.b50.m8     -H test/data/kobrite/ko00001.keg -m test/data/ko.map.partial -c 0. -l ko -o test/data/.tst.small_files.kegg.out -l 2 -l 3
    [ "$status" = 0 ]
}
@test "Checking level 3 output from count_paths with kegg map" {
    run diff <(sort test/data/small_files.kegg.out.3) <(sort test/data/.tst.small_files.kegg.out.3)
    [ "$status" = 0 ]
}
@test "Checking level 2 output from count_paths with kegg map" {
    run diff test/data/small_files.kegg.out.2 test/data/.tst.small_files.kegg.out.2
    [ "$status" = 0 ]
}
@test "Checking level 2 output from count_paths with kegg map" {
    run diff test/data/small_files.kegg.out.ko test/data/.tst.small_files.kegg.out.ko
    [ "$status" = 0 ]
}

