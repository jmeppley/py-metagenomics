#!/usr/bin/env bats

@test "PEP8 tests for CLI scripts" {
    run pep8 --ignore E731,E126,E124 *.py
    [ "$status" = 0 ]
}

@test "PEP8 tests for edl module " {
    run pep8 --ignore E731,E126,E124 edl/*.py
    [ "$status" = 0 ]
}
