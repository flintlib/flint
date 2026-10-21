---
title: Exact real and complex numbers
short: Exact numbers
summary: Algebraic numbers and exact real and complex numbers: representations, arithmetic and decision procedures.
order: 14
---

FLINT computes exactly with algebraic numbers, represented by minimal
polynomials and isolating balls ([](qqbar)), and with real and complex
numbers built from rational numbers by arithmetic, roots, exponentials,
logarithms and other functions, represented as elements of formal fields
([](ca)). The sections below describe how these representations
are constructed and certified, how arithmetic is carried out, how
equality is tested, and how
polynomial roots and linear algebra are computed on top of them. The
[overview](algorithms/exact/overview) compares the representations and
says which algorithm each operation uses. The design is described in
[@Joh2021calcium].
