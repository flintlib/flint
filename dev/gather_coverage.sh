#!/usr/bin/env bash

# This file was modified from GAP-system's dev/ci-gather-coverage.sh.  Visit
# their website at <https://www.gap-system.org/>.

# NOTE: This needs to be executed from the FLINT directory *after* checks has
# been made.

set -ex

find build -type f -name '*.gcno' -exec gcov -pb {} +
