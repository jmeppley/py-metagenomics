"""Checks that the packaging metadata in pyproject.toml stays consistent."""
import glob
import os
import re

import pytest

import edl

try:
    import tomllib
except ImportError:  # python 3.10
    tomllib = pytest.importorskip("tomli")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_pyproject():
    """Return the parsed contents of the repository's pyproject.toml."""
    with open(os.path.join(REPO_ROOT, "pyproject.toml"), "rb") as handle:
        return tomllib.load(handle)


def test_script_files_match_top_level_scripts():
    """Every top-level *.py script is installed, and nothing else is listed."""
    listed_scripts = set(load_pyproject()["tool"]["setuptools"]["script-files"])
    top_level_scripts = {os.path.basename(path)
                         for path in glob.glob(os.path.join(REPO_ROOT, "*.py"))}
    assert listed_scripts == top_level_scripts


def test_version_is_valid():
    """edl.__version__ is a plain X.Y.Z release number."""
    assert re.fullmatch(r"\d+\.\d+\.\d+", edl.__version__)
