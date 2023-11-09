.. _introduction:

**Introduction**
===============================================================================

What is Flint?
-------------------------------------------------------------------------------

FLINT is a C library of functions for doing basic arithmetic in support of
computational number theory and other areas of computer algebra. It is highly
optimised and can be compiled on numerous platforms.

FLINT provides highly optimised implementations of basic rings, such as the
integers, rationals, `p`-adics, finite fields, etc., and linear algebra and
univariate and multivariate polynomials over most of these rings.

FLINT also has some multithreading capabilities. To this end, the library is
threadsafe, with few exceptions noted in the appropriate place, and a number of
key functions have multithreaded implementations.

Maintainers and Authors
-------------------------------------------------------------------------------

FLINT is currently maintained by Fredrik Johansson of INRIA Bordeaux.

FLINT was originally designed by William Hart and David Harvey. Since then
FLINT was rewritten as FLINT 2 by William Hart, Fredrik Johansson and
Sebastian Pancratz. Many other substantial contributions have been made
by other authors, e.g. Tom Bachmann, Mike Hansen, Daniel Schultz and Andy
Novocin. There have been a great number of other contributors, listed on
the main Flint website and the contributors section of this documentation.

Requirements
-------------------------------------------------------------------------------

FLINT and following should compile on any machine with GCC and a standard
GNU toolchain, though GCC 4.8 and following are recommended.

Flint is specially optimised for x86 (32 and 64 bit) machines. There is also
limited optimisation for ARM machines.

As of version 3.0, FLINT requires GMP 6.2.1 or later, and MPFR 4.1.0 or later.
Note that earlier, MPIR, a fork of GMP, was supported. However, as of FLINT 3.0,
this support has been dropped.

It is also required that the platform provide a ``uint64_t`` type if a
native 64 bit type is not available. Full C99 compliance is not required.

Structure of Flint
-----------------------------------------------------------------------------

FLINT is supplied as a set of modules, ``fmpz``, ``fmpz_poly``, etc.,
each of which can be linked to a C program making use of their functionality.

All of the functions in FLINT have a corresponding test function provided
in an appropriately named test file.  For example, the function
``fmpz_poly_add`` located in ``src/fmpz_poly/add.c`` has test code in the
file ``src/fmpz_poly/test/t-add.c``.

Some modules have a ``profile`` directory in which profile programs can be
found.

Documentation exists in the ``doc/source`` directory in a series of ``.rst``
files.

License
-----------------------------------------------------------------------------

FLINT is distributed under the LGPL License, version 2.1+. There is a copy
of the license included in the repository and distribution tarballs.

