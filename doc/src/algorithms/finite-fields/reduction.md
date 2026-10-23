---
title: Reduction modulo the defining polynomial
summary: Elements of $\mathbb{F}_p[x]/(f)$ are reduced after each product either by eliminating the high coefficients with the few nonzero terms of a sparse $f$, or by Newton division with a precomputed inverse of the reversal of $f$; sums of products are reduced only once.
complexity: $(\ell - d)(t - 1)$ multiplications in $\mathbb{F}_p$ for a polynomial of length $\ell$ and a defining polynomial with $t$ nonzero terms; $2M(d) + O(d)$ for Newton division.
order: 10
---

In the polynomial representations of `fq` and `fq_nmod`, an element of
$\mathbb{F}_q = \mathbb{F}_p[x]/(f)$, $\deg f = d$, is stored as its
representative of degree less than $d$. A product of two elements is
computed as a product of polynomials, of length at most $2d - 1$, and
then reduced modulo $f$. For the small degrees that are common in
practice the reduction costs about as much as the product, so FLINT
chooses defining polynomials with few nonzero terms when it can
([](algorithms/finite-fields/defining-polynomials)) and exploits them.

Throughout, let $\tilde f = f / \operatorname{lc}(f) = x^d + \sum_{k}
a_k x^{j_k}$, with $t$ nonzero terms in total, so that the sum runs over
$t - 1$ exponents $0 \le j_k < d$. The context stores the coefficients
$a_k$ and the exponents $j_k$ of the monic polynomial $\tilde f$, which
generates the same ideal as $f$.

## Sparse elimination {#sparse}

::: algorithm Sparse reduction {#alg-sparse}
**Input:** $R = \sum_{i < \ell} r_i x^i$.

**Output:** $R \bmod f$, in the coefficients $r_0, \ldots, r_{d-1}$.

1. For $i = \ell - 1, \ell - 2, \ldots, d$:
   for every $k$, $r_{j_k + i - d} \gets r_{j_k + i - d} - r_i a_k$;
   then $r_i \gets 0$.
2. Reduce $r_0, \ldots, r_{d-1}$ modulo $p$.
:::

::: proposition {#prop-sparse}
[](#alg-sparse) computes $R \bmod f$ with $(\ell - d)(t - 1)$
multiplications.
:::

::: proof
Step 1 replaces $R$ by $R - r_i x^{i-d} \tilde f$, which does not change
its class modulo $f$ and clears the coefficient of $x^i$, since
$x^{i-d} \tilde f$ has leading term $x^i$; the coefficients modified
have indices $j_k + i - d < i$. After the iteration for $i$, all
coefficients of index at least $i$ are zero, so at the end $\deg R < d$.
Each iteration costs $t - 1$ multiplications.
:::

For a product, $\ell = 2d - 1$, so a trinomial costs $2(d - 1)$
multiplications, much less than a product of polynomials of degree
$d - 1$. `_fq_sparse_reduce` works with unreduced integer coefficients
(`fmpz_submul`) and reduces modulo $p$ only at the end, so the
coefficients grow by about $\log_2 p$ bits per elimination that touches
them; this is cheaper than reducing after every step for the small $t$
for which the method is used. `_fq_nmod_sparse_reduce` reduces after each
step and has four variants: for $p = 2$ the update is an exclusive or of
bits; if all $a_k$ are 1 (as for $x^d + x^j + 1$ and other polynomials
with unit coefficients), only subtractions are needed; for
$p < 2^{\bits/2}$ the update $r + r_i (p - a_k)$ fits in one word and is
reduced with a Barrett reduction (`n_mod_barrett`); otherwise
`nmod_addmul` is used.

## Newton division {#newton}

For a dense defining polynomial, the reduction is a division with
remainder by $f$ using a precomputed power series inverse. Let
$\operatorname{rev}_n(g) = x^{n-1} g(1/x)$ denote the reversal of a
polynomial of length $n$. When the context is created, the series
$v = \operatorname{rev}_{d+1}(f)^{-1} \bmod x^{d+1}$ is computed and
stored (`ctx->inv`).

::: proposition {#prop-newton}
Let $\deg R = \ell - 1$ with $d \le \ell - 1 \le 2d - 1$, and let
$R = Qf + S$ with $\deg S < d$. Then
$\operatorname{rev}_{\ell - d}(Q) = \operatorname{rev}_\ell(R) \cdot v
\bmod x^{\ell - d}$.
:::

::: proof
Reversing $R = Qf + S$ gives $\operatorname{rev}_\ell(R) =
\operatorname{rev}_{\ell - d}(Q) \operatorname{rev}_{d+1}(f) +
x^{\ell - d} \operatorname{rev}_d(S)$. Modulo $x^{\ell - d}$ the last
term vanishes, and multiplying by $v$, the inverse of
$\operatorname{rev}_{d+1}(f)$ modulo $x^{d+1}$ and hence modulo
$x^{\ell - d}$ because $\ell - d \le d + 1$, gives the claim.
:::

The quotient thus costs one truncated product of length $\ell - d \le d$
and the remainder $S = R - Qf$ another product, about $2M(d)$ operations
in all, independently of the number of terms of $f$
[@GG2003, §9.1]. `_fq_dense_reduce` and `_fq_nmod_dense_reduce` call
`_fmpz_mod_poly_divrem_newton_n_preinv` and
`_nmod_poly_divrem_newton_n_preinv`, which require
$\ell \le 2d$; this is the precondition of `fq_reduce` and
`fq_nmod_reduce` for dense moduli. The sparse method has no such limit,
and `fq_set_fmpz_poly` and `fq_nmod_set_nmod_poly` use a general division
for longer inputs.

`_fq_reduce` and `_fq_nmod_reduce` use [](#alg-sparse) if $f$ has at most
five nonzero terms (the flag `sparse_modulus`, set when the context is
created) and Newton division otherwise. With $t \le 5$ and $\ell = 2d - 1$,
the sparse method costs at most $4(d - 1)$ multiplications.

## Delayed reduction {#delayed}

Since reduction is linear, a sum $\sum_i a_i b_i$ of products of
elements can be reduced once: the products of the representatives are
accumulated as an unreduced polynomial of length $2d - 1$, over $\Z$ for
`fq` or modulo $p$ for `fq_nmod`, and only the sum is reduced. The dot
products of `_fq_vec_dot` and `_fq_nmod_vec_dot` (implemented in the
`gr` contexts), and hence the classical matrix and polynomial
multiplications, which are built on them, save all but one reduction per
output entry. Kronecker substitution ([](fq_poly/mul),
`fq_mat_mul_KS`) goes further: the whole product is computed over $\Z$
and each output coefficient is reduced once when it is unpacked.

## Use in FLINT {#use}

- `_fq_reduce`, `_fq_sparse_reduce`, `_fq_dense_reduce` and `fq_reduce`
  in `fq`; `_fq_nmod_reduce`, `_fq_nmod_sparse_reduce`,
  `_fq_nmod_dense_reduce` and `fq_nmod_reduce` in `fq_nmod`.
- `fq_mul`, `fq_sqr`, `_fq_pow`, `fq_bit_unpack` and the corresponding
  `fq_nmod` functions reduce every product.
- The `fq_zech` multiplication algorithm for polynomials reduces the
  blocks of the substituted product with its own copy of
  [](#alg-sparse), accumulating several updates in one word when
  $(d - 1)(p - 1)^2 + p < 2^{\bits}$.
