#!/usr/bin/env bats
@test "Testing count_taxa compilation" {
    run ./count_taxa.py -h
    [ "$status" = 0 ]
}

@test "Testing countTaxa: simplehit id counts" {
    run ./countTaxaFromBlasts.py S1=test/data/sample.1.kegg.blastx.b50.m8 S2=sample.2.kegg.blastx.b50.m8 S3=sample.3.kegg.blastx.b50.m8 -p hitid -F 0 -o test/data/.tst.sample.kegg.hitid.counts
    [ "$status" = 0 ]
    run diff test/data/sample.kegg.hitid.counts test/data/.tst.sample.kegg.hitid.counts
    [ "$status" = 0 ]
}

@test "Testing countTaxa: org lookup and lineages" {
    run ./countTaxaFromBlasts.py S1=test/data/sample.1.blastx.b50.m8       S2=sample.2.blastx.b50.m8 S3=sample.3.blastx.b50.m8 -o test/data/sample.org.lineage.counts -p orgs -F 0  -n test/data -C LCA -c 0.025 -r organism -r phylum -r genus -r species -R domain -R phylum
    [ "$status" = 0 ]
    run diff test/data/sample.org.lineage.counts.organism test/data/.tst.sample.org.lineage.counts.organism
    [ "$status" = 0 ]
    run diff test/data/sample.org.lineage.counts.genus test/data/.tst.sample.org.lineage.counts.genus
    [ "$status" = 0 ]
    run diff test/data/sample.org.lineage.counts.phylum test/data/.tst.sample.org.lineage.counts.phylum
    [ "$status" = 0 ]
    run diff test/data/sample.org.lineage.counts.species test/data/.tst.sample.org.lineage.counts.species
    [ "$status" = 0 ]
}

@test "Testing countTaxa: acc lookup" {
    run ./countTaxaFromBlasts.py S1=test/data/sample.1.blastx.b50.m8       S2=sample.2.blastx.b50.m8 S3=sample.3.blastx.b50.m8 -o test/data/sample.org.acc.counts -p accs -c 0.0 -C first -n test/data -m test/data/acc.to.taxid.proetin.plus.filtered
    [ "$status" = 0 ]
    run diff test/data/sample.org.acc.counts test/data/.tst.sample.org.acc.counts
    [ "$status" = 0 ]
}

@test "Testing countTaxa: top hit assignment" {
    run ./countTaxaFromBlasts.py S1=test/data/sample.1.blastx.b50.m8       S2=sample.2.blastx.b50.m8 S3=sample.3.blastx.b50.m8 -o test/data/sample.topacc.counts -p accs -c 0.0 -C tophit -R Nonr
    [ "$status" = 0 ]
    run diff test/data/sample.topacc.counts test/data/.tst.sample.topacc.counts
    [ "$status" = 0 ]
}

