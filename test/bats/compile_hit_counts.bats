#!/usr/bin/env bats

@test "Testing compile_hit_counts.py compilation" {
    run ./compile_hit_counts.py -h
    [ "$status" = 0 ]
}

@test "Testing compile_hit_counts to full table" {
    run ./compile_hit_counts.py -S -1 test/data/reads.annotations.taxon_rank.order.tsv -2 test/data/reads.annotations.gene_family.KEGG.tsv -o test/data/.tst.reads.annotations.order.KEGG.tab -v
    [ "$status" = 0 ]
}
@test "Checking output from compile_hits to full table" {
    run diff test/data/reads.annotations.order.KEGG.tab test/data/.tst.reads.annotations.order.KEGG.tab
    [ "$status" = 0 ]
}

@test "Testing compile_hit_counts to two-index table for sparse data" {
    run ./compile_hit_counts.py -L -1 test/data/reads.annotations.taxon_rank.order.tsv -2 test/data/reads.annotations.gene_family.KEGG.tsv -o test/data/.tst.reads.annotations.order.KEGG.tsv -S -v
    [ "$status" = 0 ]
    run diff test/data/reads.annotations.order.KEGG.tsv test/data/.tst.reads.annotations.order.KEGG.tsv
    [ "$status" = 0 ]
}

@test "Testing compile_hit_counts to two-index table for sparse data with known read count" {
    run ./compile_hit_counts.py -L -1 test/data/reads.annotations.taxon_rank.order.tsv -2 test/data/reads.annotations.gene_family.KEGG.tsv -o test/data/.tst.reads.annotations.order.KEGG.wUnk.tsv -S -T 900
    [ "$status" = 0 ]
    run diff test/data/reads.annotations.order.KEGG.wUnk.tsv test/data/.tst.reads.annotations.order.KEGG.wUnk.tsv
    [ "$status" = 0 ]
}

