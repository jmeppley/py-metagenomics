#!/usr/bin/env bats

@test "Testing createPrimer.py" {
    run bash -c"./createPrimerFile.py primerTemplates/truseq.1x.primers.fasta CAT > test/data/.tst.createPrimerTruseq.fasta"
    [ "$status" = 0 ]
    run diff test/data/createPrimerTruseq.fasta test/data/.tst.createPrimerTruseq.fasta
    [ "$status" = 0 ]
}

@test "Testing createPrimer.py" {
    run bash -c "./createPrimerFile.py primerTemplates/nextera.2x.primers.fasta CAT TAG > test/data/.tst.createPrimerNextera.fasta"
    [ "$status" = 0 ]
    run diff test/data/createPrimerNextera.fasta test/data/.tst.createPrimerNextera.fasta
    [ "$status" = 0 ]
}

