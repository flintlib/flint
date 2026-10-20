---
title: Univariate polynomials
short: Polynomials
summary: Dense univariate polynomials over generic rings, $\Z/n\Z$, $\Z$ and finite fields: multiplication, Newton iteration, GCDs and resultants, evaluation, interpolation and composition, Hensel lifting, squarefree factorisation, factoring over finite fields and over $\Z$, and roots.
order: 6
---

The dense univariate polynomial modules are the generic
[`gr_poly`](gr_poly) and type-specific modules such as `nmod_poly`,
`fmpz_mod_poly`, `fmpz_poly`, `fmpq_poly` and `fq_poly`, with their
factorisation modules. Division, power series, GCDs, evaluation,
interpolation and composition are reduced to multiplication so that they
inherit its quasi-linear complexity. Over $\Z/n\Z$ the basic operations
rest on the vector arithmetic of [](algorithms/linear-algebra) and the
modular arithmetic of [](algorithms/arithmetic).

Start with the overview, which describes how the generic algorithms are
organised and where the type-specific modules take over.

::: summary algorithms/polynomials/overview :::

::: summary algorithms/polynomials/multiplication :::

::: summary algorithms/polynomials/multiplication-integers :::

::: summary algorithms/polynomials/gcd :::
