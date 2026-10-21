---
title: Modular arithmetic
short: Modular arithmetic
summary: Arithmetic modulo word-size and multiword integers: division with precomputed inverses, Shoup, Montgomery, Barrett and Lemire reduction, divisibility tests, simultaneous inversion, multiword modular multiplication and delayed reduction, and how to choose among them.
order: 1
---

Almost every computation over $\Z/n\Z$, in `nmod`, `fmpz_mod`, `mpn_mod`,
the vector, polynomial and matrix modules built on them and the `gr`
contexts for them, spends its time multiplying and reducing modulo $n$.
A hardware division is too slow to do this directly, so FLINT replaces
it with multiplications by precomputed data. Most of the sections below
treat a word-size modulus; those on multiword moduli describe the
analogous methods for moduli of several limbs: a precomputed inverse of
the modulus, and sums of products accumulated exactly and reduced once.
The two overviews compare the methods and say which to use; the other
sections state each method and prove it correct.

::: summary algorithms/arithmetic/choosing :::

::: summary algorithms/arithmetic/multiword-overview :::
