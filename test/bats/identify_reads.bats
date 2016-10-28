#!/usr/bin/env bats

@test "Testing identify reads" {
    run ./identifyReads.py test/data/sample.?.blastx.b50.m8 -o .Cyanos.tst. -F 0 -a -g Cyanobacteria -n test/data -m test/data/acc.to.taxid.proetin.plus.filtered 
    [ "$status" = 0 ]
    run diff test/data/sample.1.blastx.b50.m8.Cyanos test/data/sample.1.blastx.b50.m8.Cyanos.tst.
    [ "$status" = 0 ]
    run diff test/data/sample.2.blastx.b50.m8.Cyanos test/data/sample.2.blastx.b50.m8.Cyanos.tst.
    [ "$status" = 0 ]
    run diff test/data/sample.3.blastx.b50.m8.Cyanos test/data/sample.3.blastx.b50.m8.Cyanos.tst.
    [ "$status" = 0 ]
}
