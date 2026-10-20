---
title: Generic algorithms
short: Generic rings
summary: How the generic rings framework is designed: contexts and method tables, the three levels of operations and specialisation, inexact elements and truth values, domain and unable errors, and randomised testing of ring axioms.
order: 4
---

FLINT's generic rings framework [`gr`](gr) lets one implementation of an
algorithm run over every ring, from $\Z/n\Z$ to power series over
complex balls, while the type-specific modules keep their speed by
overriding the critical operations. The sections below describe the design
decisions behind the framework and the guarantees that generic
algorithms rely on. The design is described in detail
by Johansson [@Joh2025gr]; the sections below follow the current source.

The algorithms themselves (polynomial multiplication, GCDs, linear
algebra, series) are described in their own groups, which start from
the generic version: see [](algorithms/polynomials/overview) and
[](algorithms/linear-algebra/overview).
