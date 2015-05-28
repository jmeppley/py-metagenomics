#!/usr/bin/env bats

@test "PEP8 tests for CLI scripts" {
    run pep8 *.py
    [ "$status" = 0 ]
}

@test "PEP8 tests for edl module " {
    run pep8 edl
    [ "$status" = 0 ]
}
