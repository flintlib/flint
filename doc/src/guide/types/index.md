---
title: Types and mathematical domains
short: Types
summary: Which FLINT type to use for integers, residues, finite fields, polynomials, matrices, p-adic, algebraic, real and complex numbers; ball arithmetic and exact real and complex numbers.
order: 2.7
---

FLINT often offers several types for the same mathematical object: four
representations of $\Z/n\Z$, three of finite fields, and ball,
floating-point and exact representations of real numbers, besides the
generic rings of [`gr`](gr). The sections below compare the types by
domain and explain how to compute with rigorous ball arithmetic. The
naming and calling conventions common to all types are in
[](guide/conventions/types-and-naming), and the ways in which one
algorithm serves many types in [](guide/generic).
