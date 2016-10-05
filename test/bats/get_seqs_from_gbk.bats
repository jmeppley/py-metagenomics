#!/usr/bin/env bats

@test "Testing getSeqsFromGBK: faa" {
    run ./getSequencesFromGbk.py test/data/test.gbk -c -o test/data/.tst.test.gbk.faa -f genbank -t
    [ "$status" = 0 ]
    run diff test/data/test.gbk.faa test/data/.tst.test.gbk.faa
    [ "${#lines[@]}" = 0 ]
}

@test "Testing getSeqsFromGBK: fna" {
    run ./getSequencesFromGbk.py test/data/test.gbk -c -o test/data/.tst.test.gbk.fna -f genbank
    [ "$status" = 0 ]
    run diff test/data/test.gbk.fna test/data/.tst.test.gbk.fna
    [ "$status" = 0 ]
}
