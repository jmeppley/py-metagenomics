#!/usr/bin/env bats

@test "Testing assignTopHit by top acc abundance" {
    run ./assignTopHit.py test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos -C toporg | sort > test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc
    [ "$status" = 0 ]
    run diff <(head test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc) <(head test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc)
    [ "$status" = 0 ]
    run diff <(tail test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc) <(tail test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc)
    [ "$status" = 0 ]
}

@test "Testing assignTopHit by top org abundance" {
    run ./assignTopHit.py test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos -n test/data -m test/data/acc.to.taxid.proetin.plus.filtered -C toporg | sort > test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrg
    [ "$status" = 0 ]
    run diff <(head test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrg) <(head test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrg)
    [ "$status" = 0 ]
    run diff <(tail test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrg) <(tail test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrg)
    [ "$status" = 0 ]
}

@test "Testing assignTopHit by proportional org abundance" {
    run ./assignTopHit.py test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos -n test/data -m test/data/acc.to.taxid.proetin.plus.filtered -C toporg -P | sort > test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrgP
    [ "$status" = 0 ]
    run diff <(head test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrgP) <(head test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrgP)
    [ "$status" = 0 ]
    run diff <(tail test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrgP) <(tail test/data/.tst.HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByOrgP)
    [ "$status" = 0 ]
}

@test "Testing assignTopHit on multiple inputs" {
    run ./assignTopHit.py test/data/sample.?.blastx.b50.m8 -C tophit -o .tophit
    [ "$status" = 0 ]
    run for F in sample.?.blastx.b50.m8.tophit; do mv $F $F.tmp; sort $F.tmp > $F; rm $F.tmp; done
    [ "$status" = 0 ]
    run diff test/data/sample.1.blastx.b50.m8 test/data/.tst.sample.1.blastx.b50.m8
    [ "$status" = 0 ]
    run diff test/data/sample.2.blastx.b50.m8 test/data/.tst.sample.2.blastx.b50.m8
    [ "$status" = 0 ]
    run diff test/data/sample.3.blastx.b50.m8 test/data/.tst.sample.3.blastx.b50.m8
    [ "$status" = 0 ]
}



