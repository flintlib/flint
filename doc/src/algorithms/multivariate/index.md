---
title: Multivariate polynomials
short: Multivariate
summary: Sparse multivariate polynomials: the packed representation, heap and dense multiplication and division, powering and square roots, modular GCD algorithms with dense and sparse interpolation, Hensel lifting, and factorisation over the integers and over finite fields.
order: 7
---

The sections below describe the algorithms behind the sparse
multivariate polynomial modules `fmpz_mpoly`, `fmpq_mpoly`, `nmod_mpoly`,
`fmpz_mod_mpoly`, `fq_nmod_mpoly`, `fq_zech_mpoly` and `gr_mpoly`, their
factorisation modules and the rational function modules. The overview
compares the methods and says which functions use them; the others give
the algorithms with proofs, starting from the representation of
polynomials, through arithmetic, to GCDs and factorisation. Those on GCDs
over small finite fields, factoring over finite fields and bivariate
resultants treat what is specific to finite fields.
