---
title: Types and generic programming
short: Types
summary: Which FLINT type to use for integers, residues, finite fields, polynomials, matrices, p-adic, algebraic, real and complex numbers, and how FLINT implements the same algorithm for many types.
order: 2.7
---

FLINT often offers several types for the same mathematical object: four
representations of $\Z/n\Z$, three of finite fields, ball, floating-point
and exact representations of real numbers, and a generic version of
almost everything. The first page below compares them and says which one
to choose; the second explains the ways in which FLINT shares one
algorithm between many types, from parallel copies of a module to the
generic rings of [`gr`](gr), and what this means for code that uses
them. The naming and calling conventions common to all types are in
[](guide/conventions/types-and-naming).
