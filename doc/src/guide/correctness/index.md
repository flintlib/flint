---
title: Correctness
short: Correctness
summary: What FLINT guarantees about its results, how it is tested, which inputs are valid, size limits, and why FLINT must not be used for cryptography or on untrusted data.
order: 2.5
---

FLINT is used to obtain mathematical results, and a user needs to know
how far these results can be trusted. The sections below explain what FLINT
promises and what it does not:

- the algorithms are provably correct, except where a function
  explicitly says otherwise (probable primes, Monte Carlo methods,
  floating-point approximations), and some answers are allowed to be
  "unknown" rather than wrong;
- the library is extensively tested, but not formally verified, so
  results should be checked independently when they matter;
- most functions assume that their input is valid and check nothing,
  and sizes such as lengths, shifts and precisions are not checked
  against overflow;
- FLINT is meant for research, education and experimentation; it is
  not suitable for implementing cryptographic protocols or for
  processing untrusted input outside a sandbox.

The conventions that define valid input (aliasing, normalisation, the
meaning of a leading underscore) are described in
[](guide/conventions/types-and-naming), and the semantics of ball
arithmetic, in which every result is a rigorous enclosure, in
[](guide/types/ball-arithmetic).
