#!/bin/sh
set -e

rm -rf autom4te.cache

autoreconf -f -i -v -Wall

# The following lines are from Semigroups/Semigroups, written by Max Horn.
if ! test -x config.guess -a -x config.sub ;
then
    automake -acf 2> /dev/null || :
fi

# There is a second bug in autoconf 2.69 where the generated configure
# script complains about install-sh not being there (even though it does
# not actually need it). As a workaround, we just provide an empty file
# instead. Since newer autoconf versions such as 2.71 are not affected
# by the bug, we add a test to limit when this workaround is applied
if fgrep -q ac_aux_dir/install-sh configure ;
then
    touch config/install-sh
fi

rm -rf autom4te.cache
