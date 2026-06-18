#!/bin/sh
#
# Run with any make arguments e.g.:
#    ./capture_autotools_commands.sh -j10

set -eu

make clean
./bootstrap.sh
./configure
echo ::::::::::::::::::::
echo ::::::::::::::::::::
echo ::::::::::::::::::::
make V=1 "$@" > autotools_commands.txt
