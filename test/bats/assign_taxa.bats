#!/usr/bin/env bats

@test "Testing assign_taxa compilation" {
    run ./assign_taxa.py -h
    [ "$status" = 0 ]
}

@test "Testing assign_taxa accs" {
    run ./assign_taxa.py test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos -F 0 -f gene -p accs -o test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.accs
    [ "$status" = 0 ]
}
@test "Checking assign_taxa accs output" {
    run diff test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.accs test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.accs
    [ "$status" = 0 ]
}

@test "Testing assign_taxa org" {
    run ./assign_taxa.py -F 0 -f gene test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc -p accs -m test/data/acc.to.taxid.proetin.plus.filtered -n test/data -o test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc.orgs
    [ "$status" = 0 ]
}
@test "Checking assign_taxa org output" {
    run diff test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc.orgs test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc.orgs
    [ "$status" = 0 ]
}
@test "Testing assign_taxa phylum" {
    run ./assign_taxa.py -F 0 -f gene test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc -p accs -m test/data/acc.to.taxid.proetin.plus.filtered -n test/data -o test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc.phyla -r phylum
    [ "$status" = 0 ]
}
@test "Checking assign_taxa phylum output" {
    run diff test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc.phyla test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc.phyla
    [ "$status" = 0 ]
}
