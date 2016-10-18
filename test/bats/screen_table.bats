#!/usr/bin/env bats
@test "Testing screen_table compilation" {
    run ./screen_table.py -h
    [ "$status" = 0 ]
}

@test "Testing screen_table: intersection" {
    run python ./screen_table.py test/data/table.1 -k -l test/data/table.2 -o test/data/.tst.table.in.both -D "\t" -v
    [ "$status" = 0 ]
    run diff test/data/table.in.both test/data/.tst.table.in.both
    [ "$status" = 0 ]
}

@test "Testing screen_table: unique" {
    run python ./screen_table.py test/data/table.1 -l test/data/table.2 -o test/data/.tst.table.1.uniq -D "\t" -v
    [ "$status" = 0 ]
    run diff test/data/table.1.uniq test/data/.tst.table.1.uniq
    [ "$status" = 0 ]
}
