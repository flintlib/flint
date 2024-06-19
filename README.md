[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://flintlib.org/doc/#)
[![CI](https://github.com/flintlib/flint/actions/workflows/CI.yml/badge.svg)](https://github.com/flintlib/flint/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/flintlib/flint/graph/badge.svg?token=qKHycWTpUP)](https://codecov.io/gh/flintlib/flint)

# FLINT (Fast Library for Number Theory)

Website: https://flintlib.org

Mailing list: https://groups.google.com/g/flint-devel

## Overview

FLINT is a C library in support of computations in number theory. It's also a
research project into algorithms in number theory. FLINT consists mainly of fast
scalar and polynomial arithmetic, factorization and linear algebra over many
basic rings (integers, rationals, reals, finite fields, number fields, p-adics).
It includes some higher-level functionality for algebraic and analytic number
theory.

FLINT 2, released in 2011 was a complete rewrite of FLINT 1.x from scratch.
FLINT 3, released in 2023, incorporates the [Arb](https://arblib.org/),
[Antic](https://github.com/flintlib/antic),
[Calcium](https://fredrikj.net/calcium/) and
[Generic-Rings](https://github.com/fredrik-johansson/generic-rings) libraries,
formerly developed separately.

## Documentation

For FLINT's online documentation, see https://flintlib.org/doc/.

## Building from source

This example assumes that [GMP](https://gmplib.org/), [MPFR](https://www.mpfr.org/)
and the [GNU build system](https://www.gnu.org/software/automake/manual/html_node/GNU-Build-System.html)
are already installed. To install them on a Ubuntu system, write

    apt install libgmp-dev libmpfr-dev make autoconf libtool-bin

possibly with super-user privileges.

To download, bootstrap, configure and build everything, write

    git clone https://github.com/flintlib/flint.git && cd flint
    ./bootstrap.sh
    ./configure                        # ./configure --help for more options
    make
    make check                         # optional
    make install                       # optional
    make examples                      # optional
    cd doc && make html && cd ..       # optional: documentation

See FLINT's documentation for further instructions on how to build FLINT.

## Authors

FLINT was started in 2007 by
[David Harvey](https://web.maths.unsw.edu.au/~davidharvey/) and
[William Hart](https://www.dpmms.cam.ac.uk/person/wh369). Maintenance was later
taken over solely by William Hart who remained in charge of the project
until 2022. A large number of authors have contributed to FLINT over the years;
for a complete list, see https://flintlib.org/authors.html or the `AUTHORS` file.

The current maintainers are:

* [Fredrik Johansson](https://fredrikj.net/) (fredrik.johansson@gmail.com) (project leader since 2022)
* [Albin Ahlbäck](https://albinahlback.gitlab.io/) (albin.ahlback@gmail.com)

## License

FLINT is distributed under LGPL (GNU Lesser General Public License) version 3 or
later. See the `COPYING.LESSER` and `COPYING` files.
