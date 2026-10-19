---
title: Special functions
short: Special functions
summary: Rigorous evaluation of elementary functions, constants, the gamma function, hypergeometric functions and their special cases, zeta and L-functions, and elliptic, modular and theta functions in ball arithmetic.
order: 12
---

The sections below describe how FLINT evaluates special functions with
rigorous error bounds, for real and complex balls and for power
series: the elementary functions and constants of [`arb`](arb) and
[`acb`](acb), the gamma and hypergeometric functions of
[`arb_hypgeom`](arb_hypgeom) and [`acb_hypgeom`](acb_hypgeom), the zeta
and $L$-functions of [`acb_dirichlet`](acb_dirichlet), and the
elliptic, modular and theta functions of [`acb_elliptic`](acb_elliptic),
[`acb_modular`](acb_modular) and [`acb_theta`](acb_theta). Each section
gives the formulas, the error bounds with proofs or references, the way
the method and the working precision are chosen, and the functions that
implement them. The overview explains the principles they share —
truncation bounds, cancellation-aware precision, propagated error
bounds, rectangular and binary splitting — and maps every function
family to its methods.

::: summary algorithms/special-functions/overview :::

The sections on the gamma function and on hypergeometric functions
describe the two main families of algorithms; those after them
treat the special cases that need more than the general machinery.

::: summary algorithms/special-functions/gamma :::

::: summary algorithms/special-functions/hypergeometric :::
