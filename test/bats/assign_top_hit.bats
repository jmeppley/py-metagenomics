#!/usr/bin/env bats

@test "Testing assign_top_hit by top acc abundance" {
    run ./assign_top_hit.py test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos -C toporg -o test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc
    [ "$status" = 0 ]
}
@test "Checking output from assign_top_hit by top acc abundance" {
    run diff <(sort test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc) <(sort test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc)
    [ "$status" = 0 ]
}

@test "Testing assign_top_hit by top org abundance" {
    run ./assign_top_hit.py test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos -n test/data -m test/data/acc.to.taxid.proetin.plus.filtered -C toporg -o test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrg
    [ "$status" = 0 ]
}
@test "Checking output from assign_top_hit by top org abundance" {
    run diff <(sort test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrg) <(sort test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrg)
    [ "$status" = 0 ]
}

@test "Testing assign_top_hit by proportional org abundance" {
    run ./assign_top_hit.py test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos -n test/data -m test/data/acc.to.taxid.proetin.plus.filtered -C toporg -P -o test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrgP
    [ "$status" = 0 ]
}
@test "Checking output from assign_top_hit by proportional org abundance" {
    run diff <(sort test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrgP) <(sort test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrgP)
    [ "$status" = 0 ]
}

@test "Testing assign_top_hit on multiple inputs" {
    run ./assign_top_hit.py test/data/sample.?.blastx.b50.m8 -C tophit -o .tophit.tst.
    [ "$status" = 0 ]
}
@test "Checking output 1" {
    run diff <(sort test/data/sample.1.blastx.b50.m8.tophit) <(sort test/data/sample.1.blastx.b50.m8.tophit.tst.)
    [ "$status" = 0 ]
}
@test "Checking output 2" {
    run diff <(sort test/data/sample.2.blastx.b50.m8.tophit) <(sort test/data/sample.2.blastx.b50.m8.tophit.tst.)
    [ "$status" = 0 ]
}
@test "Checking output 3" {
    run diff <(sort test/data/sample.3.blastx.b50.m8.tophit) <(sort test/data/sample.3.blastx.b50.m8.tophit.tst.)
    [ "$status" = 0 ]
}



