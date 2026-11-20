#!/bin/bash
set -euo pipefail

# build dependencies come from the host environment, not an isolated one
$PYTHON -m pip install --no-deps --no-build-isolation --ignore-installed -vv .

# database download/format workflows are not part of the python package;
# ship them under $PREFIX/share/py-metagenomics
outdir=$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM
mkdir -p $outdir
rm -f $PREFIX/share/$PKG_NAME
ln -s $PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM $PREFIX/share/$PKG_NAME
cp -R ./databases $outdir/
