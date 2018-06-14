#!/usr/bin/env bats

@test "Testing count_paths with acc map" {
    run ./count_paths.py test/data/sample.?.blastx.b50.m8 -o test/data/.tst.keggCounts.F0.acc.c0.tophit -c 0 -m test/data/acc.to.ko.protein.plus.filtered -C tophit -r -F 0 -p accs
    [ "$status" = 0 ]
}
@test "Checking output from count_paths with acc map" {
    run diff test/data/keggCounts.F0.acc.c0.tophit test/data/.tst.keggCounts.F0.acc.c0.tophit
    [ "$status" = 0 ]
}

@test "Testing sample naming" {
    run python ./count_paths.py S1=test/data/sample.1.blastx.b50.m8 S2=test/data/sample.2.blastx.b50.m8 S3=test/data/sample.3.blastx.b50.m8 -o test/data/.tst.keggCounts.F0.acc.c0.tophit.named -c 0 -m test/data/acc.to.ko.protein.plus.filtered -C tophit -r -F 0 -p accs
    [ "$status" = 0 ]
}
@test "Checking output from count_paths with sample naming" {
    run diff test/data/keggCounts.F0.acc.c0.tophit.named test/data/.tst.keggCounts.F0.acc.c0.tophit.named
    [ "$status" = 0 ]
}

@test "Testing count_paths with kegg map and multiple levels and a cutoff" {
    run python ./count_paths.py test/data/sample.?.kegg.blastx.b50.m8 -H test/data/kobrite/ko00001.keg -m test/data/ko.map.partial -c 0.15 -l ko -o test/data/.tst.small_files.kegg.c.15.out -l 2 -l 3
    [ "$status" = 0 ]   # run
    run diff <(sort test/data/small_files.kegg.c.15.out.3) <(sort test/data/.tst.small_files.kegg.c.15.out.3)
    [ "$status" = 0 ]   # level 3
    run diff <(test/data/small_files.kegg.c.15.out.2) <(test/data/.tst.small_files.kegg.c.15.out.2)
    [ "$status" = 0 ]   # level 2
    run diff test/data/small_files.kegg.c.15.out.ko test/data/.tst.small_files.kegg.c.15.out.ko
    [ "$status" = 0 ]   # level ko
}

@test "Testing count_paths with kegg map and multiple levels with no cutoff" {
    run python count_paths.py test/data/sample.?.kegg.blastx.b50.m8     -H test/data/kobrite/ko00001.keg -m test/data/ko.map.partial -c 0. -l ko -o test/data/.tst.small_files.kegg.out -l 2 -l 3
    [ "$status" = 0 ]   # run
    run diff <(sort test/data/small_files.kegg.out.3) <(sort test/data/.tst.small_files.kegg.out.3)
    [ "$status" = 0 ]   # level 3
    run diff <(sort test/data/small_files.kegg.out.2) <(sort test/data/.tst.small_files.kegg.out.2)
    [ "$status" = 0 ]   # level 2
    run diff test/data/small_files.kegg.out.ko test/data/.tst.small_files.kegg.out.ko
    [ "$status" = 0 ]   # level ko
}

@test "Testing count_paths -C first with kegg map and multiple levels with no cutoff" {
    run python count_paths.py test/data/sample.?.kegg.blastx.b50.m8 -C first -H test/data/kobrite/ko00001.keg -m test/data/ko.map.partial -c 0. -l ko -o test/data/.tst.small_files.kegg.Cfirst.out -l 2 -l 3
    [ "$status" = 0 ]   # run
    run diff <(sort test/data/small_files.kegg.Cfirst.out.3) <(sort test/data/.tst.small_files.kegg.Cfirst.out.3)
    [ "$status" = 0 ]   # level 3
    run diff <(test/data/small_files.kegg.Cfirst.out.2) <(test/data/.tst.small_files.kegg.Cfirst.out.2)
    [ "$status" = 0 ]   # level 2
    run diff test/data/small_files.kegg.Cfirst.out.ko test/data/.tst.small_files.kegg.Cfirst.out.ko
    [ "$status" = 0 ]   # level ko
}

