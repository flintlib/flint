---
title: Integer arithmetic
short: Integers
summary: Multiplication, division, GCDs, roots and radix conversion of multiprecision integers, and binary splitting.
order: 2
---

This section describes how FLINT computes with multiprecision integers:
the algorithms behind the `fmpz` functions and the limb-level routines of
`mpn_extras` on which they are built, from schoolbook methods to fast
Fourier transforms and Newton iteration. Multiplication comes first,
since the cost of every other operation is measured in multiplications:
division and square roots by Newton iteration cost a constant factor
more, and GCDs, radix conversion and binary splitting a logarithmic
factor more. FLINT uses GMP for part of this work (Toom–Cook
multiplication in a middle range of sizes, GCDs, $n$-th roots, most
modular exponentiation, radix conversion in a middle range of sizes and
in bases other than 10) and its
own code elsewhere; each page says which.

::: summary algorithms/integers/multiplication :::
