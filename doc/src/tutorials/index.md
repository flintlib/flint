---
title: Examples & tutorials
summary: Worked examples that build complete programs step by step, and a guide to the example programs shipped with FLINT. Every program shown is compiled and run, and its output checked.
---

The examples on function pages are deliberately minimal: they show one
function at work. The tutorials in this chapter are longer. Each one
solves a small but realistic problem, explains the mathematics behind it
and the choices made in the code, and points to the reference
documentation and the algorithm pages for details. The programs are
complete: they can be compiled and linked as shown in
[](tutorials/getting-started/first-program). The documentation checks
compile and run every program and compare its output with the output
shown.

FLINT's source tree also contains a directory `examples` of larger
programs: command-line tools, demonstrations and benchmarks, some of
which were used for record computations. They are described in
[](tutorials/getting-started/example-programs), and several tutorials
explain the techniques they use with shorter programs.

The tutorials are grouped by topic:

- [](tutorials/getting-started): a first program with integers and
  polynomials, compiling and linking, a program that combines several
  modules, and the example programs shipped with FLINT. If FLINT is new
  to you, start with [](tutorials/getting-started/first-program).
- [](tutorials/number-theory): integers, primes, modular arithmetic and
  elliptic curves.
- [](tutorials/algebra): programs built on polynomial arithmetic, such as
  error-correcting codes over finite fields.
- [](tutorials/linear-algebra): matrix computations over exact rings.
- [](tutorials/numerics): rigorous numerical computation with ball
  arithmetic: constants to many digits, special functions, roots of
  functions and polynomials, integrals and the zeros of the Riemann zeta
  function.
- [](tutorials/exact): exact computation with algebraic numbers and exact
  real and complex numbers.

The tutorials assume some familiarity with C and with the mathematics of
their topic, but not with FLINT. For the conventions that all modules
share (types, memory management, aliasing, threading) see
[](guide/conventions); for the meaning of the results of ball arithmetic
see [](guide/types/ball-arithmetic).
