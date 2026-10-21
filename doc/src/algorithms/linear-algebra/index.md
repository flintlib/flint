---
title: Vectors and matrices
short: Linear algebra
summary: Dot products and matrix multiplication modulo word-size integers, over $\Z$ and in floating point; Gaussian elimination, fraction-free and multimodular algorithms, characteristic polynomials, lattice reduction, Hermite and Smith normal forms and matrices of polynomials, over generic rings and their specialisations.
order: 5
---

This section describes the algorithms behind the matrix modules: the
generic [`gr_mat`](gr_mat) and the type-specific modules such as
`nmod_mat`, `fmpz_mod_mat`, `fmpz_mat`, `fmpq_mat` and the modules of
polynomial matrices. Over $\Z/n\Z$ for word-size $n$, most of the time
is spent in dot products and matrix products, and the algorithms keep
modular reductions out of the inner loops by accumulating exact sums in
one, two or three words, or by computing in floating point when the
result is guaranteed to be exact. Over $\Z$ and $\Q$, results are
computed modulo many primes or by $p$-adic lifting and then
reconstructed, or by fraction-free elimination; lattice reduction and
the Hermite and Smith normal forms complete the group.

Start with the overview, which describes the generic algorithms and
where the type-specific modules take over. The topic pages on matrix
multiplication say which method is used for which sizes, modulo $n$ and
over $\Z$.

::: summary algorithms/linear-algebra/overview :::

::: summary algorithms/linear-algebra/matrix-multiplication :::

::: summary algorithms/linear-algebra/integer-matrix-multiplication :::
