#!/usr/bin/env bats

@test "run pytest" {
    run python -m pytest -q test
    [ "$status" = 0 ]
}
