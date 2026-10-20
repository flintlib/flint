---
title: Ball arithmetic
short: Ball arithmetic
summary: Rigorous real and complex arithmetic with midpoint-radius balls: representation and error propagation, dot products, polynomial and matrix multiplication, power series, roots of polynomials and real functions, linear systems, eigenvalues and integration.
order: 11
---

This section describes the algorithms behind FLINT's ball arithmetic: the
modules [`mag`](mag), [`arf`](arf), [`arb`](arb) and [`acb`](acb), and
the polynomial, matrix and calculus modules built on them (`arb_poly`,
`acb_poly`, `arb_mat`, `acb_mat`, `arb_calc`, `acb_calc` and
`arb_fmpz_poly`). Each page gives the algorithm, the
argument that its result is an enclosure, and the functions that
implement it. Start with the overview, which explains the design and
lists the pages by topic. Special functions have their own group,
[](algorithms/special-functions), and the other floating-point types of
FLINT are compared in [](algorithms/floating-point/overview).

::: summary algorithms/ball-arithmetic/overview :::
