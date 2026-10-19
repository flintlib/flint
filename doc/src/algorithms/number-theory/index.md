---
title: Primes and number theory
short: Number theory
summary: Primality testing and proving, integer factoring, square roots and discrete logarithms modulo primes, Dirichlet characters, rational reconstruction and continued fractions, number fields, Gaussian integers, binary quadratic forms, and Bernoulli, partition, Stirling and Bell numbers.
order: 3
---

The sections below describe FLINT's algorithms for the problems of
computational number theory. For word-size integers, most methods are
made deterministic with precomputed tables, and the choice between them
is driven by constant factors rather than asymptotic complexity; the
sections on multiprecision integers cover primality proofs by APRCL and
ECPP and factoring by the elliptic curve method and the quadratic sieve.
The group also covers computations modulo primes (square roots,
discrete logarithms, Dirichlet characters), rational reconstruction and
continued fractions, arithmetic in [number fields](algorithms/number-theory/number-fields)
and in the [Gaussian integers](algorithms/number-theory/gaussian-integers),
binary quadratic forms, and the combinatorial numbers: Bernoulli
numbers, partitions, Stirling and Bell numbers. The overviews below
compare the methods for primality testing and factoring and say which
functions use them.

::: summary algorithms/number-theory/primality-word :::

::: summary algorithms/number-theory/factoring-word :::

::: summary algorithms/number-theory/integer-factoring :::
