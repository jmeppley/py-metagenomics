#!/usr/bin/env bats

@test "run nosetests" {
    run nosetests
    [ "$status" = 0 ]
}
