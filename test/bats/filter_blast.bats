#!/usr/bin/env bats

@test "Testing filter_blast_m8.py -h" {
    run ./filter_blast_m8.py -h
    [ "$status" = 0 ]
}

@test "Testing filter_blast_m8.py with no filter" {
    OF=test/data/.tst.contigs.lastn
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn > $OF"
    [ "$status" = 0 ]
    [ "$output" = "" ]
    [ -e $OF ]
    [ -s $OF ]
    OFCOMP=test/data/contigs.lastn
    run diff <(sort $OFCOMP) <(sort $OF)
    [ "$status" = 0 ]
}

@test "Testing filter_blast with non overlapping and sort" {
    OF='test/data/.tst.contigs.lastn.L100.U.H1'
    OFCOMP='test/data/contigs.lastn.L100.U.H1'
    run ./filter_blast_m8.py -f blast test/data/contigs.lastn --nonoverlapping --sort pctid -L 100 -H 1 -o $OF
    [ "$status" = 0 ]
    run diff <(sort $OFCOMP) <(sort $OF)
    [ "$status" = 0 ]
}

@test "Testing filter_blast with non overlapping with buffer" {
    OF='test/data/.tst.contigs.lastn.L100.U200.H1'
    OFCOMP='test/data/contigs.lastn.L100.U200.H1'
    run ./filter_blast_m8.py -f blast test/data/contigs.lastn --nonoverlapping 200 --sort pctid -L 100 -H 1 -o $OF
    [ "$status" = 0 ]
    run diff <(sort $OFCOMP) <(sort $OF)
    [ "$status" = 0 ]
}

@test "Testing sam pctid filtering" {
    SAMF=test/data/ALOHA_XVII_500_SSU_hits.sam
    OF=test/data/.tst.ALOHA_XVII_500_SSU_hits.filtered
    run ./filter_blast_m8.py -f sam $SAMF -o $OF -L 70
    [ "$status" = 0 ]
    run grep -c . $OF
    [ "$output" = 355 ]
    run ./filter_blast_m8.py -f sam $SAMF -o $OF -L 70 -I 97
    [ "$status" = 0 ]
    run grep -c . $OF
    [ "$output" = 265 ]
}

@test "Testing filter_blast auto output name" {
    OF=test/data/contigs.lastn.I80.L100.B50.P10.H1
    rm -f $OF
    run ./filter_blast_m8.py test/data/contigs.lastn -f blast -H 1 -P 10 -I 80 -B 50 -L 100 -O
    [ "$status" = 0 ]
    [ -e $OF ]
    [ -s $OF ]
    run grep -c . $OF
    [ "$output" = "6" ]
    rm $OF
}

@test "Checking a whole slew of filter combos" {
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 0 -B 50 -L 50 | grep -c ."
    [ "$output" = "356" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 0 -B 50 -L 100 | grep -c ."
    [ "$output" = "260" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 0 -B 50 -L 500 | grep -c ."
    [ "$output" = "98" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 0 -B 100 -L 50 | grep -c ."
    [ "$output" = "157" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 0 -B 100 -L 100 | grep -c ."
    [ "$output" = "155" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 0 -B 100 -L 500 | grep -c ."
    [ "$output" = "98" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 0 -B 500 -L 50 | grep -c ."
    [ "$output" = "30" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 0 -B 500 -L 100 | grep -c ."
    [ "$output" = "30" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 0 -B 500 -L 500 | grep -c ."
    [ "$output" = "30" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 0 -B 1000 -L 50 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 0 -B 1000 -L 100 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 0 -B 1000 -L 500 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 1 -B 50 -L 50 | grep -c ."
    [ "$output" = "170" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 1 -B 50 -L 100 | grep -c ."
    [ "$output" = "124" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 1 -B 50 -L 500 | grep -c ."
    [ "$output" = "48" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 1 -B 100 -L 50 | grep -c ."
    [ "$output" = "75" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 1 -B 100 -L 100 | grep -c ."
    [ "$output" = "74" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 1 -B 100 -L 500 | grep -c ."
    [ "$output" = "48" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 1 -B 500 -L 50 | grep -c ."
    [ "$output" = "15" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 1 -B 500 -L 100 | grep -c ."
    [ "$output" = "15" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 1 -B 500 -L 500 | grep -c ."
    [ "$output" = "15" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 1 -B 1000 -L 50 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 1 -B 1000 -L 100 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 0 -H 1 -B 1000 -L 500 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 0 -B 50 -L 50 | grep -c ."
    [ "$output" = "8" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 0 -B 50 -L 100 | grep -c ."
    [ "$output" = "8" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 0 -B 50 -L 500 | grep -c ."
    [ "$output" = "7" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 0 -B 100 -L 50 | grep -c ."
    [ "$output" = "7" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 0 -B 100 -L 100 | grep -c ."
    [ "$output" = "7" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 0 -B 100 -L 500 | grep -c ."
    [ "$output" = "7" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 0 -B 500 -L 50 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 0 -B 500 -L 100 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 0 -B 500 -L 500 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 0 -B 1000 -L 50 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 0 -B 1000 -L 100 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 0 -B 1000 -L 500 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 1 -B 50 -L 50 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 1 -B 50 -L 100 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 1 -B 50 -L 500 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 1 -B 100 -L 50 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 1 -B 100 -L 100 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 1 -B 100 -L 500 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 1 -B 500 -L 50 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 1 -B 500 -L 100 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 1 -B 500 -L 500 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 1 -B 1000 -L 50 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 1 -B 1000 -L 100 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 1 -H 1 -B 1000 -L 500 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 0 -B 50 -L 50 | grep -c ."
    [ "$output" = "63" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 0 -B 50 -L 100 | grep -c ."
    [ "$output" = "60" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 0 -B 50 -L 500 | grep -c ."
    [ "$output" = "38" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 0 -B 100 -L 50 | grep -c ."
    [ "$output" = "50" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 0 -B 100 -L 100 | grep -c ."
    [ "$output" = "50" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 0 -B 100 -L 500 | grep -c ."
    [ "$output" = "38" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 0 -B 500 -L 50 | grep -c ."
    [ "$output" = "20" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 0 -B 500 -L 100 | grep -c ."
    [ "$output" = "20" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 0 -B 500 -L 500 | grep -c ."
    [ "$output" = "20" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 0 -B 1000 -L 50 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 0 -B 1000 -L 100 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 0 -B 1000 -L 500 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 1 -B 50 -L 50 | grep -c ."
    [ "$output" = "26" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 1 -B 50 -L 100 | grep -c ."
    [ "$output" = "26" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 1 -B 50 -L 500 | grep -c ."
    [ "$output" = "18" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 1 -B 100 -L 50 | grep -c ."
    [ "$output" = "23" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 1 -B 100 -L 100 | grep -c ."
    [ "$output" = "23" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 1 -B 100 -L 500 | grep -c ."
    [ "$output" = "18" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 1 -B 500 -L 50 | grep -c ."
    [ "$output" = "10" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 1 -B 500 -L 100 | grep -c ."
    [ "$output" = "10" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 1 -B 500 -L 500 | grep -c ."
    [ "$output" = "10" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 1 -B 1000 -L 50 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 1 -B 1000 -L 100 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 10 -P 10 -H 1 -B 1000 -L 500 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 0 -B 50 -L 50 | grep -c ."
    [ "$output" = "237" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 0 -B 50 -L 100 | grep -c ."
    [ "$output" = "206" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 0 -B 50 -L 500 | grep -c ."
    [ "$output" = "98" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 0 -B 100 -L 50 | grep -c ."
    [ "$output" = "157" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 0 -B 100 -L 100 | grep -c ."
    [ "$output" = "155" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 0 -B 100 -L 500 | grep -c ."
    [ "$output" = "98" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 0 -B 500 -L 50 | grep -c ."
    [ "$output" = "30" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 0 -B 500 -L 100 | grep -c ."
    [ "$output" = "30" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 0 -B 500 -L 500 | grep -c ."
    [ "$output" = "30" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 0 -B 1000 -L 50 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 0 -B 1000 -L 100 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 0 -B 1000 -L 500 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 1 -B 50 -L 50 | grep -c ."
    [ "$output" = "114" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 1 -B 50 -L 100 | grep -c ."
    [ "$output" = "98" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 1 -B 50 -L 500 | grep -c ."
    [ "$output" = "48" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 1 -B 100 -L 50 | grep -c ."
    [ "$output" = "75" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 1 -B 100 -L 100 | grep -c ."
    [ "$output" = "74" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 1 -B 100 -L 500 | grep -c ."
    [ "$output" = "48" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 1 -B 500 -L 50 | grep -c ."
    [ "$output" = "15" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 1 -B 500 -L 100 | grep -c ."
    [ "$output" = "15" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 1 -B 500 -L 500 | grep -c ."
    [ "$output" = "15" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 1 -B 1000 -L 50 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 1 -B 1000 -L 100 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 0 -H 1 -B 1000 -L 500 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 0 -B 50 -L 50 | grep -c ."
    [ "$output" = "8" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 0 -B 50 -L 100 | grep -c ."
    [ "$output" = "8" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 0 -B 50 -L 500 | grep -c ."
    [ "$output" = "7" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 0 -B 100 -L 50 | grep -c ."
    [ "$output" = "7" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 0 -B 100 -L 100 | grep -c ."
    [ "$output" = "7" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 0 -B 100 -L 500 | grep -c ."
    [ "$output" = "7" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 0 -B 500 -L 50 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 0 -B 500 -L 100 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 0 -B 500 -L 500 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 0 -B 1000 -L 50 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 0 -B 1000 -L 100 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 0 -B 1000 -L 500 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 1 -B 50 -L 50 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 1 -B 50 -L 100 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 1 -B 50 -L 500 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 1 -B 100 -L 50 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 1 -B 100 -L 100 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 1 -B 100 -L 500 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 1 -B 500 -L 50 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 1 -B 500 -L 100 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 1 -B 500 -L 500 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 1 -B 1000 -L 50 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 1 -B 1000 -L 100 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 1 -H 1 -B 1000 -L 500 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 0 -B 50 -L 50 | grep -c ."
    [ "$output" = "53" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 0 -B 50 -L 100 | grep -c ."
    [ "$output" = "53" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 0 -B 50 -L 500 | grep -c ."
    [ "$output" = "38" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 0 -B 100 -L 50 | grep -c ."
    [ "$output" = "50" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 0 -B 100 -L 100 | grep -c ."
    [ "$output" = "50" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 0 -B 100 -L 500 | grep -c ."
    [ "$output" = "38" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 0 -B 500 -L 50 | grep -c ."
    [ "$output" = "20" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 0 -B 500 -L 100 | grep -c ."
    [ "$output" = "20" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 0 -B 500 -L 500 | grep -c ."
    [ "$output" = "20" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 0 -B 1000 -L 50 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 0 -B 1000 -L 100 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 0 -B 1000 -L 500 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 1 -B 50 -L 50 | grep -c ."
    [ "$output" = "23" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 1 -B 50 -L 100 | grep -c ."
    [ "$output" = "23" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 1 -B 50 -L 500 | grep -c ."
    [ "$output" = "18" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 1 -B 100 -L 50 | grep -c ."
    [ "$output" = "23" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 1 -B 100 -L 100 | grep -c ."
    [ "$output" = "23" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 1 -B 100 -L 500 | grep -c ."
    [ "$output" = "18" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 1 -B 500 -L 50 | grep -c ."
    [ "$output" = "10" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 1 -B 500 -L 100 | grep -c ."
    [ "$output" = "10" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 1 -B 500 -L 500 | grep -c ."
    [ "$output" = "10" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 1 -B 1000 -L 50 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 1 -B 1000 -L 100 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-10 -P 10 -H 1 -B 1000 -L 500 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 0 -B 50 -L 50 | grep -c ."
    [ "$output" = "184" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 0 -B 50 -L 100 | grep -c ."
    [ "$output" = "174" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 0 -B 50 -L 500 | grep -c ."
    [ "$output" = "98" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 0 -B 100 -L 50 | grep -c ."
    [ "$output" = "157" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 0 -B 100 -L 100 | grep -c ."
    [ "$output" = "155" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 0 -B 100 -L 500 | grep -c ."
    [ "$output" = "98" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 0 -B 500 -L 50 | grep -c ."
    [ "$output" = "30" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 0 -B 500 -L 100 | grep -c ."
    [ "$output" = "30" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 0 -B 500 -L 500 | grep -c ."
    [ "$output" = "30" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 0 -B 1000 -L 50 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 0 -B 1000 -L 100 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 0 -B 1000 -L 500 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 1 -B 50 -L 50 | grep -c ."
    [ "$output" = "90" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 1 -B 50 -L 100 | grep -c ."
    [ "$output" = "85" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 1 -B 50 -L 500 | grep -c ."
    [ "$output" = "48" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 1 -B 100 -L 50 | grep -c ."
    [ "$output" = "75" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 1 -B 100 -L 100 | grep -c ."
    [ "$output" = "74" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 1 -B 100 -L 500 | grep -c ."
    [ "$output" = "48" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 1 -B 500 -L 50 | grep -c ."
    [ "$output" = "15" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 1 -B 500 -L 100 | grep -c ."
    [ "$output" = "15" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 1 -B 500 -L 500 | grep -c ."
    [ "$output" = "15" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 1 -B 1000 -L 50 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 1 -B 1000 -L 100 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 0 -H 1 -B 1000 -L 500 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 0 -B 50 -L 50 | grep -c ."
    [ "$output" = "7" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 0 -B 50 -L 100 | grep -c ."
    [ "$output" = "7" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 0 -B 50 -L 500 | grep -c ."
    [ "$output" = "7" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 0 -B 100 -L 50 | grep -c ."
    [ "$output" = "7" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 0 -B 100 -L 100 | grep -c ."
    [ "$output" = "7" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 0 -B 100 -L 500 | grep -c ."
    [ "$output" = "7" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 0 -B 500 -L 50 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 0 -B 500 -L 100 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 0 -B 500 -L 500 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 0 -B 1000 -L 50 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 0 -B 1000 -L 100 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 0 -B 1000 -L 500 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 1 -B 50 -L 50 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 1 -B 50 -L 100 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 1 -B 50 -L 500 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 1 -B 100 -L 50 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 1 -B 100 -L 100 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 1 -B 100 -L 500 | grep -c ."
    [ "$output" = "3" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 1 -B 500 -L 50 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 1 -B 500 -L 100 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 1 -B 500 -L 500 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 1 -B 1000 -L 50 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 1 -B 1000 -L 100 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 1 -H 1 -B 1000 -L 500 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 0 -B 50 -L 50 | grep -c ."
    [ "$output" = "50" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 0 -B 50 -L 100 | grep -c ."
    [ "$output" = "50" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 0 -B 50 -L 500 | grep -c ."
    [ "$output" = "38" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 0 -B 100 -L 50 | grep -c ."
    [ "$output" = "50" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 0 -B 100 -L 100 | grep -c ."
    [ "$output" = "50" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 0 -B 100 -L 500 | grep -c ."
    [ "$output" = "38" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 0 -B 500 -L 50 | grep -c ."
    [ "$output" = "20" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 0 -B 500 -L 100 | grep -c ."
    [ "$output" = "20" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 0 -B 500 -L 500 | grep -c ."
    [ "$output" = "20" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 0 -B 1000 -L 50 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 0 -B 1000 -L 100 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 0 -B 1000 -L 500 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 1 -B 50 -L 50 | grep -c ."
    [ "$output" = "23" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 1 -B 50 -L 100 | grep -c ."
    [ "$output" = "23" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 1 -B 50 -L 500 | grep -c ."
    [ "$output" = "18" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 1 -B 100 -L 50 | grep -c ."
    [ "$output" = "23" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 1 -B 100 -L 100 | grep -c ."
    [ "$output" = "23" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 1 -B 100 -L 500 | grep -c ."
    [ "$output" = "18" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 1 -B 500 -L 50 | grep -c ."
    [ "$output" = "10" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 1 -B 500 -L 100 | grep -c ."
    [ "$output" = "10" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 1 -B 500 -L 500 | grep -c ."
    [ "$output" = "10" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 1 -B 1000 -L 50 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 1 -B 1000 -L 100 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 1e-15 -P 10 -H 1 -B 1000 -L 500 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 0 -B 50 -L 50 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 0 -B 50 -L 100 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 0 -B 50 -L 500 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 0 -B 100 -L 50 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 0 -B 100 -L 100 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 0 -B 100 -L 500 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 0 -B 500 -L 50 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 0 -B 500 -L 100 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 0 -B 500 -L 500 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 0 -B 1000 -L 50 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 0 -B 1000 -L 100 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 0 -B 1000 -L 500 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 1 -B 50 -L 50 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 1 -B 50 -L 100 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 1 -B 50 -L 500 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 1 -B 100 -L 50 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 1 -B 100 -L 100 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 1 -B 100 -L 500 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 1 -B 500 -L 50 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 1 -B 500 -L 100 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 1 -B 500 -L 500 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 1 -B 1000 -L 50 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 1 -B 1000 -L 100 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 0 -H 1 -B 1000 -L 500 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 0 -B 50 -L 50 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 0 -B 50 -L 100 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 0 -B 50 -L 500 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 0 -B 100 -L 50 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 0 -B 100 -L 100 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 0 -B 100 -L 500 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 0 -B 500 -L 50 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 0 -B 500 -L 100 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 0 -B 500 -L 500 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 0 -B 1000 -L 50 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 0 -B 1000 -L 100 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 0 -B 1000 -L 500 | grep -c ."
    [ "$output" = "2" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 1 -B 50 -L 50 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 1 -B 50 -L 100 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 1 -B 50 -L 500 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 1 -B 100 -L 50 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 1 -B 100 -L 100 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 1 -B 100 -L 500 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 1 -B 500 -L 50 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 1 -B 500 -L 100 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 1 -B 500 -L 500 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 1 -B 1000 -L 50 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 1 -B 1000 -L 100 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 1 -H 1 -B 1000 -L 500 | grep -c ."
    [ "$output" = "1" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 0 -B 50 -L 50 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 0 -B 50 -L 100 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 0 -B 50 -L 500 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 0 -B 100 -L 50 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 0 -B 100 -L 100 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 0 -B 100 -L 500 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 0 -B 500 -L 50 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 0 -B 500 -L 100 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 0 -B 500 -L 500 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 0 -B 1000 -L 50 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 0 -B 1000 -L 100 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 0 -B 1000 -L 500 | grep -c ."
    [ "$output" = "12" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 1 -B 50 -L 50 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 1 -B 50 -L 100 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 1 -B 50 -L 500 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 1 -B 100 -L 50 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 1 -B 100 -L 100 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 1 -B 100 -L 500 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 1 -B 500 -L 50 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 1 -B 500 -L 100 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 1 -B 500 -L 500 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 1 -B 1000 -L 50 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 1 -B 1000 -L 100 | grep -c ."
    [ "$output" = "6" ]
    run bash -c "./filter_blast_m8.py -f blast test/data/contigs.lastn -E 0 -P 10 -H 1 -B 1000 -L 500 | grep -c ."
    [ "$output" = "6" ]
}
