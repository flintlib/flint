---
title: fq_nmod_embed.h
header: fq_nmod_embed.h
template-headers: fq_embed_templates
description: Embeddings between finite fields in the word-size characteristic.
category: Finite fields
order: 16
---

This module computes embeddings $\mathbb{F}_{p^m} \to \mathbb{F}_{p^n}$, $m \mid n$, between fields in the representation of [](fq_nmod) (a prime $p < 2^{\bits}$, elements stored as `nmod_poly`), as matrices over $\mathbb{F}_p$ (`nmod_mat_t`), and the change-of-basis matrices they rest on. It is generated from the same template as [](fq_embed) and has the same functions, with `fq_nmod` in place of `fq` and `nmod_poly` and `nmod_mat` in place of `fmpz_mod_poly` and `fmpz_mod_mat`.

::: summary algorithms/finite-fields/embeddings :::

## Embeddings {#embeddings}

::: function
void fq_nmod_embed_gens(fq_nmod_t gen_sub, fq_nmod_t gen_sup, nmod_poly_t minpoly, const fq_nmod_ctx_t sub_ctx, const fq_nmod_ctx_t sup_ctx)

Given contexts for fields $k$ = `sub_ctx` of degree $m$ and $K$ =
`sup_ctx` of degree $n$ with $m \mid n$, computes a generator `gen_sub`
of $k$ over $\mathbb{F}_p$, its minimal polynomial `minpoly`, and a root
`gen_sup` of `minpoly` in $K$. These data define the embedding
$k \to K$ that maps `gen_sub` to `gen_sup`.

::: details
Generated from the same template as `fq_embed_gens`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_embed_gens_naive(fq_nmod_t gen_sub, fq_nmod_t gen_sup, nmod_poly_t minpoly, const fq_nmod_ctx_t sub_ctx, const fq_nmod_ctx_t sup_ctx)

Computes the data of `fq_nmod_embed_gens` as follows: `gen_sub` is the
generator $x$ of `sub_ctx`, `minpoly` is the defining polynomial of
`sub_ctx`, and `gen_sup` is a root of it in `sup_ctx`, found with
`fq_nmod_poly_factor_split_single`, that is, by root finding for a
polynomial of degree $m$ over $K$. Requires $m > 1$ to divide $n$.

::: details
Generated from the same template as `_fq_embed_gens_naive`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_embed_matrices(nmod_mat_t embed, nmod_mat_t project, const fq_nmod_t gen_sub, const fq_nmod_ctx_t sub_ctx, const fq_nmod_t gen_sup, const fq_nmod_ctx_t sup_ctx, const nmod_poly_t gen_minpoly)

Given the output of `fq_nmod_embed_gens` for fields $k$ = `sub_ctx` and $K$ =
`sup_ctx` of degrees $m \mid n$, computes the $n \times m$ matrix
`embed` of the embedding $\varphi: k \to K$ that maps `gen_sub` to
`gen_sup`, and an $m \times n$ matrix `project` of a $k$-linear
projection $\pi: K \to k$ with $\pi \circ \varphi = \mathrm{id}$, so that
$\mathit{project} \cdot \mathit{embed}$ is the $m \times m$ identity.
The matrices act on column vectors of coordinates in the monomial bases.

::: details
Generated from the same template as `fq_embed_matrices`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_embed_trace_matrix(nmod_mat_t res, const nmod_mat_t basis, const fq_nmod_ctx_t sub_ctx, const fq_nmod_ctx_t sup_ctx)

Given contexts for $k$ = `sub_ctx` and $K$ = `sup_ctx` of degrees
$m \mid n$ and the $n \times m$ matrix `basis` of an embedding
$\varphi: k \to K$, sets `res` to the $m \times n$ matrix of the map
$y \mapsto \varphi^{-1}(\operatorname{Tr}_{K/\varphi(k)}(y))$. It is
computed as $D_k \cdot \mathit{basis}^T \cdot M_K$, with $D_k$ =
`fq_nmod_embed_dual_to_mono_matrix` for $k$ and $M_K$ =
`fq_nmod_embed_mono_to_dual_matrix` for $K$. For $m = n$ the result is the
inverse of `basis`.

::: details
Generated from the same template as `fq_embed_trace_matrix`, whose documentation has the details.
:::
:::

## Matrices of elements {#matrices}

::: function
void fq_nmod_embed_composition_matrix(nmod_mat_t matrix, const fq_nmod_t gen, const fq_nmod_ctx_t ctx)
void fq_nmod_embed_composition_matrix_sub(nmod_mat_t matrix, const fq_nmod_t gen, const fq_nmod_ctx_t ctx, slong trunc)

Set `matrix` to the composition matrix of the element $a$ = `gen` of a
field of degree $d$: the $d \times d$ matrix whose columns are the
coordinates of $a^0, a^1, \ldots, a^{d-1}$, which maps the coordinates of
a polynomial $g$ of degree less than $d$ to those of $g(a)$. The `_sub`
version computes only the first `trunc` $\le d$ columns, into a
$d \times \mathit{trunc}$ matrix.

::: details
Generated from the same template as `fq_embed_composition_matrix`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_embed_mul_matrix(nmod_mat_t matrix, const fq_nmod_t gen, const fq_nmod_ctx_t ctx)

Sets the $d \times d$ matrix `matrix` to the matrix of multiplication by
$a$ = `gen`, whose columns are the coordinates of
$a, ax, \ldots, ax^{d-1}$. The defining polynomial need not be monic.

::: details
Generated from the same template as `fq_embed_mul_matrix`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_embed_mono_to_dual_matrix(nmod_mat_t res, const fq_nmod_ctx_t ctx)

Sets the $d \times d$ matrix `res` to the change of basis from the
monomial basis $1, x, \ldots, x^{d-1}$ to its dual basis with respect to
the trace form $(u, v) \mapsto \operatorname{Tr}(uv)$: the Hankel matrix
with entries $\operatorname{Tr}(x^{i+j})$, computed from the power series
expansion of $\operatorname{rev}(f')/\operatorname{rev}(f)$ (Newton's
identities). Mapping a coordinate vector of $y$ by this matrix gives the
traces $\operatorname{Tr}(y x^i)$.

::: details
Generated from the same template as `fq_embed_mono_to_dual_matrix`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_embed_dual_to_mono_matrix(nmod_mat_t res, const fq_nmod_ctx_t ctx)

Sets the $d \times d$ matrix `res` to the inverse of
`fq_nmod_embed_mono_to_dual_matrix`. By Euler's formula the dual basis of
$1, x, \ldots, x^{d-1}$ is $b_0 / f'(x), \ldots, b_{d-1} / f'(x)$, where
$f(X)/(X - x) = \sum b_j X^j$, so the matrix is the product of the
multiplication matrix of $1/f'(x)$ and the Hankel matrix of the
coefficients of $f$.

::: details
Generated from the same template as `fq_embed_dual_to_mono_matrix`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_modulus_pow_series_inv(nmod_poly_t res, const fq_nmod_ctx_t ctx, slong trunc)

Sets `res` to the inverse of the reversal of the defining polynomial, as
a power series modulo $x^{\mathit{trunc}}$.

::: details
Generated from the same template as `fq_modulus_pow_series_inv`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_modulus_derivative_inv(fq_nmod_t m_prime, fq_nmod_t m_prime_inv, const fq_nmod_ctx_t ctx)

Sets `m_prime` to $f'(x)$, the derivative of the defining polynomial
evaluated at the generator, and `m_prime_inv` to its inverse. $f'(x) \ne
0$ because $f$ is irreducible, hence separable.

::: details
Generated from the same template as `fq_modulus_derivative_inv`, whose documentation has the details.
:::
:::

## Example {#example}
