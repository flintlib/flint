---
title: Introduction
summary: What FLINT is, what it contains, how the library and its documentation are organised, what it requires, and its license.
order: 1
---

FLINT (Fast Library for Number Theory) is a C library for computations in
number theory and computer algebra. It provides arithmetic in basic
rings and fields, among them integers,
rationals, integers modulo $n$, finite fields, $p$-adic numbers and
number fields, together with univariate and multivariate polynomials,
power series and matrices, with factorisation, linear
algebra and number theoretic functions. It also provides arbitrary
precision real and complex arithmetic with rigorous error bounds (ball
arithmetic), special functions, exact real and complex numbers, and a
generic interface to all of these rings. FLINT is highly optimised, and
it is portable: it is written in C11, with assembly and SIMD code for
x86-64 and ARM64.

FLINT is also a research project in the algorithms of computer algebra.
The algorithm pages of this documentation describe the methods that are
implemented, with references to the literature.

## Contents of the library {#contents}

The library is divided into modules, one per header file. The home page
lists them by topic. In outline:

- *Integers and rationals*: `fmpz` (integers of any size, stored in one
  word when small), `fmpq`, and integer factorisation and primality
  proving.
- *Word-size arithmetic*: `ulong_extras` (number theory on words) and
  `nmod` (arithmetic modulo a word-size $n$).
- *Polynomials and matrices* over $\Z$, $\Q$, $\Z/n\Z$, finite fields,
  $p$-adic numbers and number fields: for example `fmpz_poly`,
  `nmod_poly`, `fq_poly`, `fmpz_mat`, `nmod_mat` and the multivariate
  `fmpz_mpoly` and `nmod_mpoly`.
- *Real and complex numbers*: the ball arithmetic of `arb` and `acb`
  with polynomials, matrices and special functions (the former Arb
  library), floating-point types (`arf`, `nfloat`, `mp_real`), and exact
  real and complex numbers in `ca` and `qqbar` (the former Calcium
  library).
- *Generic rings*: the `gr` interface, through which the same code works
  over any of the rings above, with generic polynomials, matrices and
  series in `gr_poly`, `gr_mat` and `gr_series`; the guide
  [](gr/getting-started) introduces it.

Many operations automatically use several threads when the user allows
it; see [](guide/conventions/threading).

## Organisation {#organisation}

A module `foo` consists of the header `src/foo.h` and the directory
`src/foo`, which holds one source file per function or small group of
functions, a `test` directory with one test program `t-*.c` per function
or group, and in some modules `profile` and `tune` directories with
benchmark and tuning programs. For example, `fmpz_poly_add` is in
`src/fmpz_poly/add.c` and is tested by `src/fmpz_poly/test/t-add.c`.
Programs using FLINT include the headers as `<flint/foo.h>` and link
against the single library `libflint`.

This documentation has four parts: the User guide (this part), the
reference documentation of each module (every public function with its
prototype; the module pages also contain guides to using the module),
the algorithm pages, and the examples and tutorials. Conventions
shared by all modules (types, naming, memory management, aliasing) are
described in [](guide/conventions/types-and-naming) and the following
pages.

## Requirements {#requirements}

FLINT requires

- GMP, version 6.2.1 or later (https://gmplib.org/);
- MPFR, version 4.1.0 or later (https://mpfr.org/);
- a C compiler supporting C11, such as GCC or Clang, and on Windows also
  MSVC;
- for building, GNU Make (with the Autotools when building from a git
  checkout), or CMake on Windows.

MPIR, a fork of GMP that FLINT supported before version 3.0, is no longer
supported. Optional dependencies are a BLAS library,
the Boehm–Demers–Weiser garbage collector, and
NTL (only for testing the conversion functions of the C++ header
`NTL-interface.h`). FLINT is optimised in
particular for x86-64 and ARM64 processors, and also runs on other 64-bit
and 32-bit platforms. See [](guide/getting-started/building).

## History and authors {#history}

FLINT was started by William Hart and David Harvey in 2007, and William
Hart led the project until 2022. FLINT 2, released in 2011, was a
complete rewrite by William Hart, Sebastian Pancratz and Fredrik
Johansson; the 2.x series added, among much else, finite fields (Mike
Hansen) and multivariate polynomials (Daniel Schultz). FLINT 3, released
in 2023, merged the formerly separate libraries Arb (ball arithmetic) and
Calcium (exact real and complex numbers), started by Fredrik Johansson in
2012 and 2020, Antic (algebraic number theory), started by William Hart
in 2013, and Generic-Rings. FLINT is currently maintained by Fredrik
Johansson (project leader) and Albin Ahlbäck. Many other people have
contributed; see [](guide/development/contributors) and the contributors
listed for each release in [](history).

## License {#license}

FLINT is free software, distributed under the GNU Lesser General Public
License (LGPL), version 3 or later. A copy of the license is included in
the source distribution (files `COPYING` and `COPYING.LESSER`). FLINT
versions 2.6 to 3.0 were distributed under the LGPL version 2.1 or later,
and earlier versions under the GPL version 2 or later.

## Citing FLINT {#citing}

If FLINT contributes to published research, please cite it as "The
FLINT team. FLINT: Fast Library for Number Theory, Version 3.x.y,
https://flintlib.org", with the version that was used. Specific algorithms often have their own
papers, given in the references of the corresponding pages.
