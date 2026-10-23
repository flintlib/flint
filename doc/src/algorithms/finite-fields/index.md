---
title: Finite fields
short: Finite fields
summary: Arithmetic in finite fields of non-prime order: the polynomial and Zech logarithm representations, reduction modulo sparse defining polynomials, the Frobenius map, norm and trace, embeddings between fields, and the choice of defining polynomials.
order: 9
---

This section describes how the modules `fq_default`, `fq`, `fq_nmod` and
`fq_zech` and their vector, matrix, polynomial and embedding modules
compute in finite fields $\mathbb{F}_{p^d}$. The overview compares the
representations and lists the methods used for each task; the other
pages state the algorithms and prove them correct. Polynomials and
matrices over finite fields mostly use the generic algorithms of the
[polynomial](algorithms/polynomials/overview) and
[linear algebra](algorithms/linear-algebra/overview) groups, and their
factorisation is described in
[](algorithms/polynomials/factoring-finite-fields).
