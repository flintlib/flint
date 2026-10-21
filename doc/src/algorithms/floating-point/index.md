---
title: Floating-point arithmetic
short: Floating point
summary: Fixed-precision multiprecision floating-point and decimal arithmetic, double-word and interval arithmetic with machine doubles.
order: 13
---

This section describes how FLINT implements floating-point arithmetic
outside the `arb` ball arithmetic (which has its own group,
[](algorithms/ball-arithmetic/overview)): the fixed-precision `nfloat`
type and its kernels, the limb-level `mp_real` balls and their
elementary functions, constants and hypergeometric series, double-word
to quad-double expansions and balls in `dfloat`, and interval
arithmetic with doubles. Start with the overview, which compares all
floating-point types in FLINT and says which to use.
