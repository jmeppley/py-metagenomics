#!/usr/bin/env bats

@test "Testing assignTopHit by top acc abundance" {
    run ./assignTopHit.py test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos -C toporg -o test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc
    [ "$status" = 0 ]
}
@test "Checking output from assignTopHit by top acc abundance" {
    run diff <(sort test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc) <(sort test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc)
    [ "$status" = 0 ]
}

@test "Testing assignTopHit by top org abundance" {
    run ./assignTopHit.py test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos -n test/data -m test/data/acc.to.taxid.proetin.plus.filtered -C toporg -o test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrg
    [ "$status" = 0 ]
}
@test "Checking output from assignTopHit by top org abundance" {
    run diff <(sort test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrg) <(sort test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrg)
    [ "$status" = 0 ]
}

@test "Testing assignTopHit by proportional org abundance" {
    run ./assignTopHit.py test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos -n test/data -m test/data/acc.to.taxid.proetin.plus.filtered -C toporg -P -o test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrgP
    [ "$status" = 0 ]
}
@test "Checking output from assignTopHit by proportional org abundance" {
    run diff <(sort test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrgP) <(sort test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrgP)
    [ "$status" = 0 ]
}

@test "Testing assignTopHit on multiple inputs" {
    run ./assignTopHit.py test/data/sample.?.blastx.b50.m8 -C tophit -o .tophit.tst.
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



