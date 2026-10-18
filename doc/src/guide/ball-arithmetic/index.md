---
title: Ball arithmetic
short: Ball arithmetic
summary: How to use FLINT's rigorous real and complex arithmetic and how to interpret its results.
order: 3
---

FLINT computes with real and complex numbers in *ball arithmetic*: every
number is a floating-point midpoint together with an error bound, and
every operation returns a ball that is guaranteed to contain the exact
result. The modules `mag`, `arf`, `arb`, `acb` and those built on them
(`arb_poly`, `arb_hypgeom`, `acb_theta`, ...) provide this arithmetic and
the functions computed with it. (They were developed as the separate
library Arb, which was merged into FLINT in version 3.0, and keep its
module names and interfaces.)

The pages below explain what ball arithmetic is, how to write programs
with it, how the precision relates to the accuracy of the results, and
the conventions specific to ball arithmetic. The general rules of
[](guide/conventions) and [](guide/correctness) apply as well. The module pages
[arb](arb) and [acb](acb) and the guide [](arb/getting-started) give
the programming interface, and the algorithms are described in
[](algorithms/ball-arithmetic) and [](algorithms/special-functions).
